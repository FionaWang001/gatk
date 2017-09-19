package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVLocation;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD;

/**
 * This deals with the special case where a contig's multiple (> 2) alignments has head and tail mapped to the same chr.
 * For the case where the head and tail mapped to different chromosome, we could decide to emit all BND records, but
 * that could be dealt with later.
 */
final class CpxVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> localAssemblyContigs, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final SAMSequenceDictionary refSequenceDictionary,
                                   final Logger toolLogger){
        // preprocess AI list

        // after the preprocessing, the configuration might be different from the input configuration,
        // so again divert the contigs into different units

        // extract reference ordered jumping locations on reference

        // segment affect reference regions by jumping locations

        // make sense of event, i.e. provide interpretation, and extract corresponding alt haplotype

        // output VCF
    }

    /**
     * Preprocess provided {@code originalConfiguration} of a particular contig and return a configuration after operation.
     * Note that "original" is meant to be possibly different from the returned configuration, 
     * but DOES NOT mean the alignments of the contig as given by the aligner, i.e. the configuration should be 
     * one of the best given by {@link FilterLongReadAlignmentsSAMSpark#pickBestConfigurations(AlignedContig, Set)}.
     * 
     * Note that the configuration after this preprocessing may not even have the same number of alignment as in input configuration.
     */
    private static List<AlignmentInterval> preprocessAlignments(final List<AlignmentInterval> originalConfiguration,
                                                                final int mapQThresholdInclusive,
                                                                final int uniqReadLenInclusive) {

        // two pass, each focusing on removing the alignments of a contig that offers low uniqueness

        // first pass is for removing alignments with low ref uniqueness, using low mapping quality as the criteria
        final List<AlignmentInterval> firstPassSurvivors =
                originalConfiguration.stream().filter(ai -> ai.mapQual >= mapQThresholdInclusive).collect(Collectors.toList());

        // second pass, the slower one, is to remove alignments offering low read uniqueness,
        // i.e. with only a very short part of the read being explained by this particular alignment;
        // the steps are:
        //      search bi-directionally until cannot find overlap any more, subtract from it all overlaps.
        //      This gives unique read region it explains. If this unique read region is "short"
        //      (e.g. less than half of its aligner-assigned read consumption length, or shorter than 30 bp),
        //      drop it.
        final List<Integer> idxToRemove = new ArrayList<>(firstPassSurvivors.size());
        for (int i = 0; i < firstPassSurvivors.size(); ++i) { // implementation is less-than efficient because we are computing the overlap twice, but prototyping for now
            final AlignmentInterval cur = firstPassSurvivors.get(i);
            int maxOverlap = -1;
            for (int j = 0; j != i && j < firstPassSurvivors.size(); ++j) { // bi-directional
                final int overlap = AlignmentInterval.overlapOnContig(cur, firstPassSurvivors.get(j));
                if (overlap > 0)
                    maxOverlap = Math.max(maxOverlap, overlap);
                else // following ones, as guaranteed by the ordering of alignments in the contig, cannot overlap
                    break;
            }
            // TODO: 9/21/17 clip a copy of cur here and check if it is too short
            if (cur.endInAssembledContig - cur.startInAssembledContig + 1 < uniqReadLenInclusive)
                idxToRemove.add(i);
        }

        if ( idxToRemove.isEmpty() )
            return firstPassSurvivors;

        // removing in reverse order so that iterators are not invalidated if we were to remove from start
        final ListIterator<Integer> rit = idxToRemove.listIterator(idxToRemove.size());
        while (rit.hasPrevious()) {
            firstPassSurvivors.remove( rit.previous().intValue() );
        }

        return firstPassSurvivors;
    }

    /**
     * Each pair of neighboring reference locations are meant to be used closed, i.e. [a, b].
     */
    private static List<SVLocation> extractReferenceOrdereredJumpLocations(final List<AlignmentInterval> alignmentConfiguration) {

        // A configuration has a series jumps on the reference as indicated by the chimeric alignments.
        // A jump has a starting and landing ref location.
        // A jump can be: 1) gapped, meaning a part of read is uncovered by neighboring AI's;
        //                2) connected, meaning neighboring--but not overlapping on the read--AI's leave no base on the read uncovered;
        //                3) retracting, meaning neighboring AI's overlap on the read, pointing to homology between their ref span
        // Among them, retracting jumps are the most difficult to deal with, mainly due to how to have a consistent
        // homology-yielding scheme

        return null;
    }

    private static List<SVInterval> segmentReference(final List<SVLocation> jumpingLocations ) {
        return null;
    }

    /**
     * Splits input alignments into ones that are intended to be used for chimeric alignments and
     * ones that are not reliable for that purpose,
     * based on provided MQ and unique ref/read span length threshold.
     */
    private static Tuple2<List<AlignmentInterval>, List<AlignmentInterval>> classifyAlignments(final Iterator<AlignmentInterval> iterator,
                                                                                               final int uniqueRefSpanThreshold,
                                                                                               final int uniqueReadSpanThreshold,
                                                                                               final boolean filterWhollyContainedAlignments) {

        final List<AlignmentInterval> good = new ArrayList<>(10); // 10 is a blunt guess
        final List<AlignmentInterval> bad  = new ArrayList<>(10);

        AlignmentInterval current = iterator.next();
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            final AlnPairUniqueLength alnPairUniqueLength = new AlnPairUniqueLength(current, next);

            if (alignmentIsNonInformative(current.mapQual, CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD,
                    alnPairUniqueLength.oneUniqRefLen, uniqueRefSpanThreshold, alnPairUniqueLength.oneUniqReadLen, uniqueReadSpanThreshold)){
                bad.add(current);
                current = next;
            } else if (alignmentIsNonInformative(next.mapQual, CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD,
                    alnPairUniqueLength.twoUniqRefLen, uniqueRefSpanThreshold, alnPairUniqueLength.twoUniqReadLen, uniqueReadSpanThreshold)) {
                bad.add(next);
            } else {
                good.add(current);
                current = next;
            }
        }

        return new Tuple2<>(good, bad);
    }

    private static boolean alignmentIsNonInformative(final int mapQ, final int mapqThresholdInclusive,
                                                     final int uniqRefSpanLen, final int uniqueRefSpanThreshold,
                                                     final int uniqReadSpanLen, final int uniqueReadSpanThreshold) {
        return mapQ < mapqThresholdInclusive
                || uniqRefSpanLen < uniqueRefSpanThreshold
                || uniqReadSpanLen < uniqueReadSpanThreshold;

    }

    private static final class AlnPairUniqueLength {
        final int oneUniqRefLen;
        final int oneUniqReadLen;
        final int twoUniqRefLen;
        final int twoUniqReadLen;

        AlnPairUniqueLength(final AlignmentInterval one, final AlignmentInterval two) {
            Utils.validateArg(one.startInAssembledContig <= two.startInAssembledContig,
                    "assumption that input alignments are order along read is violated");

            final int overlapOnRefSpan = AlignmentInterval.overlapOnRefSpan(one, two);
            final int overlapOnRead = AlignmentInterval.overlapOnContig(one, two);

            if (overlapOnRead == 0) {
                oneUniqRefLen = one.referenceSpan.size() - overlapOnRefSpan;
                twoUniqRefLen = two.referenceSpan.size() - overlapOnRefSpan;
                oneUniqReadLen = one.endInAssembledContig - one.startInAssembledContig + 1;
                twoUniqReadLen = two.endInAssembledContig - two.startInAssembledContig + 1;
            } else {
                // TODO: 10/16/17 hardclip offset
                final int i = SvCigarUtils.computeAssociatedDistOnRef(one.cigarAlong5to3DirectionOfContig, two.startInAssembledContig, overlapOnRead);
                final int j = SvCigarUtils.computeAssociatedDistOnRef(two.cigarAlong5to3DirectionOfContig, two.startInAssembledContig, overlapOnRead);
                oneUniqRefLen = one.referenceSpan.size() - Math.max(i, overlapOnRefSpan);
                twoUniqRefLen = two.referenceSpan.size() - Math.max(j, overlapOnRefSpan);
                oneUniqReadLen = one.endInAssembledContig - one.startInAssembledContig + 1 - overlapOnRead;
                twoUniqReadLen = two.endInAssembledContig - two.startInAssembledContig + 1 - overlapOnRead;
            }
        }
    }
}
