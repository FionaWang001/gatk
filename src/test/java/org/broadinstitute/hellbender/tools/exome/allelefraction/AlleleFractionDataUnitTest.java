package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * @author David Benjamin
 */
public final class AlleleFractionDataUnitTest {

    @Test
    public void testData() {
        final List<AllelicCount> ac = new ArrayList<>();
        final List<SimpleInterval> segments= new ArrayList<>();

        // segment 0: hets 0-2
        segments.add(new SimpleInterval("chr", 1, 5));
        ac.add(new AllelicCount(new SimpleInterval("chr", 1, 1), 0, 5));
        ac.add(new AllelicCount(new SimpleInterval("chr", 2, 2), 5, 0));
        ac.add(new AllelicCount(new SimpleInterval("chr", 3, 3), 5, 5));

        // segment 1: hets 3-4
        segments.add(new SimpleInterval("chr", 10, 15));
        ac.add(new AllelicCount(new SimpleInterval("chr", 10, 10), 1, 1));
        ac.add(new AllelicCount(new SimpleInterval("chr", 11, 11), 2, 2));

        final Genome genome = new Genome(AlleleFractionSimulatedData.TRIVIAL_TARGETS, ac);

        final SegmentedGenome segmentedGenome = new SegmentedGenome(segments, genome);

        final AlleleFractionData dc = new AlleleFractionData(segmentedGenome);
        Assert.assertEquals(dc.getNumSegments(), 2);

        Assert.assertEquals(dc.getCountsInSegment(0).get(1).getRefReadCount(), 5);
        Assert.assertEquals(dc.getCountsInSegment(0).get(1).getAltReadCount(), 0);

        Assert.assertEquals(dc.getNumPoints(), 5);
    }
}