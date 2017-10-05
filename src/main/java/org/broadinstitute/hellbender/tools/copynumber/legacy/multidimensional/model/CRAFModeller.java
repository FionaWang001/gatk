package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.model;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model.AlleleFractionPrior;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model.AlleleFractionSegmentedData;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.model.CopyRatioModeller;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.model.CopyRatioSegmentedData;
import org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation.CRAFSegmentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;
import org.broadinstitute.hellbender.utils.mcmc.ParameterWriter;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents a segmented model for copy ratio and allele fraction.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CRAFModeller {
    private static final Logger logger = LogManager.getLogger(CRAFModeller.class);

    private static final String DOUBLE_FORMAT = "%6.6f";

    //use 95% HPD interval to construct {@link PosteriorSummary} for segment means and minor allele fractions
    private static final double CREDIBLE_INTERVAL_ALPHA = 0.05;

    private final String sampleName;
    private final CopyRatioCollection denoisedCopyRatios;
    private final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector;
    private final AllelicCountCollection allelicCounts;
    private final OverlapDetector<AllelicCount> allelicCountOverlapDetector;
    private final AlleleFractionPrior alleleFractionPrior;

    private CopyRatioModeller copyRatioModeller;
    private AlleleFractionModeller alleleFractionModeller;

    private List<SimpleInterval> currentSegments;
    private final List<ModeledSegment> modeledSegments = new ArrayList<>();

    //similar-segment merging may leave model in a state where it is not completely fit (i.e., posterior modes will be improperly specified)
    private boolean isModelFit;

    private final int numSamplesCopyRatio;
    private final int numBurnInCopyRatio;
    private final int numSamplesAlleleFraction;
    private final int numBurnInAlleleFraction;
    private final JavaSparkContext ctx;

    /**
     * Constructs a copy-ratio and allele-fraction modeller, specifying number of total samples
     * and number of burn-in samples for Markov-Chain Monte Carlo model fitting.
     * An initial model fit is performed.
     * @param ctx   JavaSparkContext, used for kernel density estimation in {@link PosteriorSummary}
     */
    public CRAFModeller(final CRAFSegmentCollection crafSegments,
                        final CopyRatioCollection denoisedCopyRatios,
                        final AllelicCountCollection allelicCounts,
                        final AlleleFractionPrior alleleFractionPrior,
                        final int numSamplesCopyRatio,
                        final int numBurnInCopyRatio,
                        final int numSamplesAlleleFraction,
                        final int numBurnInAlleleFraction,
                        final JavaSparkContext ctx) {
        Utils.validateArg(Stream.of(
                Utils.nonNull(crafSegments).getSampleName(),
                Utils.nonNull(denoisedCopyRatios).getSampleName(),
                Utils.nonNull(allelicCounts).getSampleName()).distinct().count() == 1,
                "Sample names from all inputs must match.");
        ParamUtils.isPositive(crafSegments.size(), "Number of segments must be positive.");
        sampleName = crafSegments.getSampleName();
        currentSegments = crafSegments.getIntervals();
        this.denoisedCopyRatios = denoisedCopyRatios;
        copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
        this.allelicCounts = allelicCounts;
        allelicCountOverlapDetector = allelicCounts.getOverlapDetector();
        this.alleleFractionPrior = Utils.nonNull(alleleFractionPrior);
        this.numSamplesCopyRatio = numSamplesCopyRatio;
        this.numBurnInCopyRatio = numBurnInCopyRatio;
        this.numSamplesAlleleFraction = numSamplesAlleleFraction;
        this.numBurnInAlleleFraction = numBurnInAlleleFraction;
        this.ctx = ctx;
        logger.info("Fitting initial model...");
        fitModel();
    }

    public ModeledSegmentCollection getModeledSegments() {
        return new ModeledSegmentCollection(sampleName, modeledSegments);
    }

    /**
     * Performs Markov-Chain Monte Carlo model fitting using the
     * number of total samples and number of burn-in samples specified at construction.
     */
    private void fitModel() {
        //perform MCMC to generate posterior samples
        logger.info("Fitting copy-ratio model...");
        final CopyRatioSegmentedData copyRatioSegmentedData =
                new CopyRatioSegmentedData(denoisedCopyRatios, currentSegments);
        copyRatioModeller = new CopyRatioModeller(copyRatioSegmentedData);
        final AlleleFractionSegmentedData alleleFractionSegmentedData =
                new AlleleFractionSegmentedData(allelicCounts, currentSegments);
        copyRatioModeller.fitMCMC(numSamplesCopyRatio, numBurnInCopyRatio);
        logger.info("Fitting allele-fraction model...");
        alleleFractionModeller = new AlleleFractionModeller(alleleFractionSegmentedData, alleleFractionPrior);
        alleleFractionModeller.fitMCMC(numSamplesAlleleFraction, numBurnInAlleleFraction);

        //update list of ModeledSegment with new PosteriorSummaries
        modeledSegments.clear();
        final List<PosteriorSummary> segmentMeansPosteriorSummaries =
                copyRatioModeller.getSegmentMeansPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final List<PosteriorSummary> minorAlleleFractionsPosteriorSummaries =
                alleleFractionModeller.getMinorAlleleFractionsPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        for (int segmentIndex = 0; segmentIndex < currentSegments.size(); segmentIndex++) {
            final SimpleInterval segment = currentSegments.get(segmentIndex);
            final int numPointsCopyRatio = copyRatioMidpointOverlapDetector.getOverlaps(segment).size();
            final int numPointsAlleleFraction = allelicCountOverlapDetector.getOverlaps(segment).size();
            final PosteriorSummary segmentMeansPosteriorSummary = segmentMeansPosteriorSummaries.get(segmentIndex);
            final PosteriorSummary minorAlleleFractionPosteriorSummary = minorAlleleFractionsPosteriorSummaries.get(segmentIndex);
            modeledSegments.add(new ModeledSegment(
                    segment, numPointsCopyRatio, numPointsAlleleFraction, segmentMeansPosteriorSummary, minorAlleleFractionPosteriorSummary));
        }
        isModelFit = true;
    }

    /**
     * @param numSmoothingIterationsPerFit  if this is zero, no refitting will be performed between smoothing iterations
     */
    public void smoothSegments(final int maxNumSmoothingIterations,
                               final int numSmoothingIterationsPerFit,
                               final double smoothingCredibleIntervalThresholdCopyRatio,
                               final double smoothingCredibleIntervalThresholdAlleleFraction) {
        ParamUtils.isPositiveOrZero(maxNumSmoothingIterations,
                "The maximum number of smoothing iterations must be non-negative.");
        ParamUtils.isPositiveOrZero(smoothingCredibleIntervalThresholdCopyRatio,
                "The number of smoothing iterations per fit must be non-negative.");
        ParamUtils.isPositiveOrZero(smoothingCredibleIntervalThresholdAlleleFraction,
                "The allele-fraction credible-interval threshold for segmentation smoothing must be non-negative.");
        logger.info(String.format("Initial number of segments before smoothing: %d", modeledSegments.size()));
        //perform iterations of similar-segment merging until all similar segments are merged
        for (int numIterations = 1; numIterations <= maxNumSmoothingIterations; numIterations++) {
            logger.info(String.format("Smoothing iteration: %d", numIterations));
            final int prevNumSegments = modeledSegments.size();
            if (numSmoothingIterationsPerFit > 0 && numIterations % numSmoothingIterationsPerFit == 0) {
                //refit model after this merge iteration
                performSmoothingIteration(smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction, true);
            } else {
                //do not refit model after this merge iteration (posterior modes will be identical to posterior medians)
                performSmoothingIteration(smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction, false);
            }
            if (modeledSegments.size() == prevNumSegments) {
                break;
            }
        }
        if (!isModelFit) {
            //make sure final model is completely fit (i.e., posterior modes are specified)
            fitModel();
        }
        logger.info(String.format("Final number of segments after smoothing: %d", modeledSegments.size()));
    }

    /**
     * Performs one iteration of similar-segment merging on the list of {@link ModeledSegment} held internally.
     * Markov-Chain Monte Carlo model fitting is optionally performed after each iteration using the
     * number of total samples and number of burn-in samples specified at construction.
     * @param intervalThresholdSegmentMean         threshold number of credible intervals for segment-mean similarity
     * @param intervalThresholdMinorAlleleFraction threshold number of credible intervals for minor-allele-fraction similarity
     * @param doModelFit                           if true, refit MCMC model after merging
     */
    private void performSmoothingIteration(final double intervalThresholdSegmentMean,
                                           final double intervalThresholdMinorAlleleFraction,
                                           final boolean doModelFit) {
        logger.info("Number of segments before smoothing iteration: " + modeledSegments.size());
        final List<ModeledSegment> mergedSegments = SimilarSegments.mergeSimilarSegments(modeledSegments, intervalThresholdSegmentMean, intervalThresholdMinorAlleleFraction);
        logger.info("Number of segments after smoothing iteration: " + mergedSegments.size());
        currentSegments = mergedSegments.stream().map(ModeledSegment::getInterval).collect(Collectors.toList());
        if (doModelFit) {
            fitModel();
        } else {
            modeledSegments.clear();
            modeledSegments.addAll(mergedSegments);
            isModelFit = false;
        }
    }

    /**
     * Writes posterior summaries for the global model parameters to a file.
     */
    public void writeModelParameterFiles(final File copyRatioParameterFile,
                                         final File alleleFractionParameterFile) {
        Utils.nonNull(copyRatioParameterFile);
        Utils.nonNull(alleleFractionParameterFile);
        ensureModelIsFit();
        logger.info("Writing posterior summaries for copy-ratio global parameters to " + copyRatioParameterFile);
        writeModelParameterFile(copyRatioModeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx), copyRatioParameterFile);
        logger.info("Writing posterior summaries for allele-fraction global parameters to " + alleleFractionParameterFile);
        writeModelParameterFile(alleleFractionModeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx), alleleFractionParameterFile);
    }

    private void ensureModelIsFit() {
        if (!isModelFit) {
            logger.warn("Attempted to write ACNV results to file when model was not completely fit. Performing model fit now.");
            fitModel();
        }
    }

    private <T extends Enum<T> & ParameterEnum> void writeModelParameterFile(final Map<T, PosteriorSummary> parameterPosteriorSummaries,
                                                                             final File outFile) {
        try (final ParameterWriter<T> writer = new ParameterWriter<>(outFile, DOUBLE_FORMAT)) {
            writer.writeAllRecords(parameterPosteriorSummaries.entrySet());
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }

    /**
     * Contains private methods for similar-segment merging.
     */
    private static final class SimilarSegments {
        /**
         * Returns a new, modifiable list of segments with similar segments (i.e., adjacent segments with both
         * segment-mean and minor-allele-fractions posteriors similar; posteriors are similar if the difference between
         * posterior central tendencies is less than intervalThreshold times the posterior credible interval of either summary)
         * merged.  The list of segments is traversed once from beginning to end, and each segment is checked for similarity
         * with the segment to the right and merged until it is no longer similar.
         * @param intervalThresholdSegmentMean         threshold number of credible intervals for segment-mean similarity
         * @param intervalThresholdMinorAlleleFraction threshold number of credible intervals for minor-allele-fraction similarity
         */
        private static List<ModeledSegment> mergeSimilarSegments(final List<ModeledSegment> segments,
                                                                 final double intervalThresholdSegmentMean,
                                                                 final double intervalThresholdMinorAlleleFraction) {
            final List<ModeledSegment> mergedSegments = new ArrayList<>(segments);
            int index = 0;
            while (index < mergedSegments.size() - 1) {
                final ModeledSegment segment1 = mergedSegments.get(index);
                final ModeledSegment segment2 = mergedSegments.get(index + 1);
                if (segment1.getContig().equals(segment2.getContig()) &&
                        SimilarSegments.areSimilar(segment1, segment2,
                                intervalThresholdSegmentMean, intervalThresholdMinorAlleleFraction)) {
                    mergedSegments.set(index, SimilarSegments.merge(segment1, segment2));
                    mergedSegments.remove(index + 1);
                    index--; //if merge performed, stay on current segment during next iteration
                }
                index++; //if no merge performed, go to next segment during next iteration
            }
            return mergedSegments;
        }

        //checks similarity of posterior summaries to within a credible-interval threshold;
        //posterior summaries are similar if the difference between posterior central tendencies is less than
        //intervalThreshold times the credible-interval width of either summary
        private static boolean areSimilar(final ModeledSegment.SimplePosteriorSummary summary1,
                                          final ModeledSegment.SimplePosteriorSummary summary2,
                                          final double intervalThreshold) {
            if (Double.isNaN(summary1.getDecile50()) || Double.isNaN(summary2.getDecile50())) {
                return true;
            }
            final double absoluteDifference = Math.abs(summary1.getDecile50() - summary2.getDecile50());
            return absoluteDifference < intervalThreshold * (summary1.getDecile90() - summary1.getDecile10()) ||
                    absoluteDifference < intervalThreshold * (summary2.getDecile90() - summary2.getDecile10());
        }

        //checks similarity of modeled segments to within credible-interval thresholds for segment mean and minor allele fraction
        private static boolean areSimilar(final ModeledSegment segment1,
                                          final ModeledSegment segment2,
                                          final double intervalThresholdSegmentMean,
                                          final double intervalThresholdMinorAlleleFraction) {
            return areSimilar(segment1.getLog2CopyRatioSimplePosteriorSummary(), segment2.getLog2CopyRatioSimplePosteriorSummary(), intervalThresholdSegmentMean) &&
                    areSimilar(segment1.getMinorAlleleFractionSimplePosteriorSummary(), segment2.getMinorAlleleFractionSimplePosteriorSummary(), intervalThresholdMinorAlleleFraction);
        }

        //merges posterior summaries naively by approximating posteriors as normal
        private static ModeledSegment.SimplePosteriorSummary merge(final ModeledSegment.SimplePosteriorSummary summary1,
                                                                   final ModeledSegment.SimplePosteriorSummary summary2) {
            if (Double.isNaN(summary1.getDecile50()) && !Double.isNaN(summary2.getDecile50())) {
                return summary2;
            }
            if ((!Double.isNaN(summary1.getDecile50()) && Double.isNaN(summary2.getDecile50())) ||
                    (Double.isNaN(summary1.getDecile50()) && Double.isNaN(summary2.getDecile50()))) {
                return summary1;
            }
            //use credible half-interval as standard deviation
            final double standardDeviation1 = (summary1.getDecile90() - summary1.getDecile10()) / 2.;
            final double standardDeviation2 = (summary2.getDecile90() - summary2.getDecile10()) / 2.;
            final double variance = 1. / (1. / Math.pow(standardDeviation1, 2.) + 1. / Math.pow(standardDeviation2, 2.));
            final double mean =
                    (summary1.getDecile50() / Math.pow(standardDeviation1, 2.) + summary2.getDecile50() / Math.pow(standardDeviation2, 2.))
                            * variance;
            final double standardDeviation = Math.sqrt(variance);
            //we simply use the naive mean as the mode
            return new ModeledSegment.SimplePosteriorSummary(mean, mean, mean - standardDeviation, mean + standardDeviation);
        }

        private static ModeledSegment merge(final ModeledSegment segment1,
                                            final ModeledSegment segment2) {
            return new ModeledSegment(mergeSegments(segment1.getInterval(), segment2.getInterval()),
                    segment1.getNumPointsCopyRatio() + segment2.getNumPointsCopyRatio(),
                    segment1.getNumPointsAlleleFraction() + segment2.getNumPointsAlleleFraction(),
                    merge(segment1.getLog2CopyRatioSimplePosteriorSummary(), segment2.getLog2CopyRatioSimplePosteriorSummary()),
                    merge(segment1.getMinorAlleleFractionSimplePosteriorSummary(), segment2.getMinorAlleleFractionSimplePosteriorSummary()));
        }

        private static SimpleInterval mergeSegments(final SimpleInterval segment1,
                                                    final SimpleInterval segment2) {
            Utils.validateArg(segment1.getContig().equals(segment2.getContig()),
                    String.format("Cannot join segments %s and %s on different chromosomes.", segment1.toString(), segment2.toString()));
            final int start = Math.min(segment1.getStart(), segment2.getStart());
            final int end = Math.max(segment1.getEnd(), segment2.getEnd());
            return new SimpleInterval(segment1.getContig(), start, end);
        }
    }
}
