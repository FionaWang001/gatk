package org.broadinstitute.hellbender.engine;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class AssemblyRegionIteratorUnitTest extends BaseTest {

    @DataProvider
    public Object[][] testCorrectRegionsHaveCorrectReadsAndSizeData() {
        return new Object[][] {
                { NA12878_20_21_WGS_bam, b37_reference_20_21, new SimpleInterval("20", 10000000, 10100000), 50, 300, 100 }
        };
    }

    /*
     * This test checks that over a 100,000 base interval, all assembly regions created by the AssemblyRegionIterator
     * have the correct reads, that the reads are stored in the correct order, and that region padding and other
     * settings are respected.
     *
     * We determine this by going by to the original BAM and doing a fresh query for each region, to ensure that
     * the query results match the reads actually in the region.
     */
    @Test(dataProvider = "testCorrectRegionsHaveCorrectReadsAndSizeData")
    public void testRegionsHaveCorrectReadsAndSize( final String reads, final String reference, final SimpleInterval shardInterval, final int minRegionSize, final int maxRegionSize, final int assemblyRegionPadding ) throws IOException {
        try ( final ReadsDataSource readsSource = new ReadsDataSource(IOUtils.getPath(reads));
              final ReferenceDataSource refSource = ReferenceDataSource.of(new File(reference)) ) {
            final SAMSequenceDictionary readsDictionary = readsSource.getSequenceDictionary();
            final LocalReadShard readShard = new LocalReadShard(shardInterval, shardInterval.expandWithinContig(assemblyRegionPadding, readsDictionary), readsSource);
            final AssemblyRegionEvaluator evaluator = new HaplotypeCallerEngine(new HaplotypeCallerArgumentCollection(), false, false, readsSource.getHeader(), new CachingIndexedFastaSequenceFile(new File(b37_reference_20_21)));
            final ReadCoordinateComparator readComparator = new ReadCoordinateComparator(readsSource.getHeader());

            final List<ReadFilter> readFilters = new ArrayList<>(2);
            readFilters.add(new WellformedReadFilter());
            readFilters.add(new ReadFilterLibrary.MappedReadFilter());
            final CountingReadFilter combinedReadFilter = CountingReadFilter.fromList(readFilters, readsSource.getHeader());
            readShard.setReadFilter(combinedReadFilter);

            final AssemblyRegionIterator iter = new AssemblyRegionIterator(readShard, readsSource.getHeader(), refSource, null, evaluator, minRegionSize, maxRegionSize, assemblyRegionPadding, 0.002, 50);

            AssemblyRegion previousRegion = null;
            while ( iter.hasNext() ) {
                final AssemblyRegion region = iter.next();

                Assert.assertTrue(region.getSpan().size() <= maxRegionSize, "region size " + region.getSpan().size() + " exceeds the configured maximum: " + maxRegionSize);

                final int regionContigLength = readsDictionary.getSequence(region.getSpan().getContig()).getSequenceLength();
                final int expectedLeftRegionPadding = region.getSpan().getStart() - assemblyRegionPadding > 0 ? assemblyRegionPadding : region.getSpan().getStart() - 1;
                final int expectedRightRegionPadding = region.getSpan().getEnd() + assemblyRegionPadding <= regionContigLength ? assemblyRegionPadding : regionContigLength - region.getSpan().getEnd();
                Assert.assertEquals(region.getSpan().getStart() - region.getExtendedSpan().getStart(), expectedLeftRegionPadding, "Wrong amount of padding on the left side of the region");
                Assert.assertEquals(region.getExtendedSpan().getEnd() - region.getSpan().getEnd(), expectedRightRegionPadding, "Wrong amount of padding on the right side of the region");
                final SimpleInterval regionInterval = region.getExtendedSpan();
                final List<GATKRead> regionActualReads = region.getReads();

                if ( previousRegion != null ) {
                    Assert.assertTrue(IntervalUtils.isBefore(previousRegion.getSpan(), region.getSpan(), readsDictionary), "Previous assembly region's span is not before the current assembly region's span");
                    Assert.assertEquals(previousRegion.getSpan().getEnd(), region.getSpan().getStart() - 1, "previous and current regions are not contiguous");
                }

                GATKRead previousRead = null;
                for ( final GATKRead currentRead : regionActualReads ) {
                    if ( previousRead != null ) {
                        Assert.assertTrue(readComparator.compare(previousRead, currentRead) <= 0, "Reads are out of order within the assembly region");
                    }
                    previousRead = currentRead;
                }

                try ( final ReadsDataSource innerReadsSource = new ReadsDataSource(IOUtils.getPath(reads)) ) {
                    final List<GATKRead> regionExpectedReads = Lists.newArrayList(innerReadsSource.query(regionInterval)).stream().filter(combinedReadFilter).collect(Collectors.toList());

                    final List<GATKRead> actualNotInExpected = new ArrayList<>();
                    final List<GATKRead> expectedNotInActual = new ArrayList<>();
                    for ( final GATKRead expectedRead : regionExpectedReads ) {
                        if ( ! regionActualReads.contains(expectedRead) ) {
                            expectedNotInActual.add(expectedRead);
                        }
                    }

                    for ( final GATKRead actualRead : regionActualReads ) {
                        if ( ! regionExpectedReads.contains(actualRead) ) {
                            actualNotInExpected.add(actualRead);
                        }
                    }

                    Assert.assertEquals(regionActualReads.size(), regionExpectedReads.size(), "Wrong number of reads in region " + region + " for extended interval " + regionInterval +
                            ". Expected reads not in actual reads: " + expectedNotInActual + ". Actual reads not in expected reads: " + actualNotInExpected);

                    Assert.assertEquals(regionActualReads, regionExpectedReads, "Wrong reads in region " + region + " for extended interval " + regionInterval +
                            ". Expected reads not in actual reads: " + expectedNotInActual + ". Actual reads not in expected reads: " + actualNotInExpected);
                }
            }
        }
    }
}
