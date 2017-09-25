package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.OptionalDouble;

public class CRAFSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPointsCopyRatio;
    private final int numPointsAlleleFraction;
    private final double meanLog2CopyRatio;
    private final double meanMinorAlleleFraction;

    public CRAFSegment(final SimpleInterval interval,
                       final int numPointsCopyRatio,
                       final int numPointsAlleleFraction,
                       final double meanLog2CopyRatio,
                       final double meanMinorAlleleFraction) {
        Utils.nonNull(interval);
        Utils.validateArg(numPointsCopyRatio > 0 || numPointsAlleleFraction > 0,
                "Number of copy-ratio points or number of allele-fraction points must be positive.");
        this.interval = interval;
        this.numPointsCopyRatio = numPointsCopyRatio;
        this.numPointsAlleleFraction = numPointsAlleleFraction;
        this.meanLog2CopyRatio = meanLog2CopyRatio;
        this.meanMinorAlleleFraction = meanMinorAlleleFraction;
    }

    public CRAFSegment(final SimpleInterval interval,
                       final List<Double> denoisedCopyRatios,
                       final List<Double> alternateAlleleFractions) {
        Utils.nonNull(interval);
        Utils.nonNull(denoisedCopyRatios);
        Utils.nonNull(alternateAlleleFractions);
        this.interval = interval;
        numPointsCopyRatio = denoisedCopyRatios.size();
        numPointsAlleleFraction = alternateAlleleFractions.size();
        final OptionalDouble optionalMeanLog2CopyRatio = denoisedCopyRatios.stream().mapToDouble(Double::doubleValue).average();
        meanLog2CopyRatio = optionalMeanLog2CopyRatio.isPresent() ? optionalMeanLog2CopyRatio.getAsDouble() : Double.NaN;
        final OptionalDouble optionalMeanMinorAlleleFraction = alternateAlleleFractions.stream().mapToDouble(aaf -> Math.min(aaf, 1. - aaf)).average();
        meanMinorAlleleFraction = optionalMeanMinorAlleleFraction.isPresent() ? optionalMeanMinorAlleleFraction.getAsDouble() : Double.NaN;
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public int getNumPointsCopyRatio() {
        return numPointsCopyRatio;
    }

    public int getNumPointsAlleleFraction() {
        return numPointsAlleleFraction;
    }

    public double getMeanLog2CopyRatio() {
        return meanLog2CopyRatio;
    }

    public double getMeanMinorAlleleFraction() {
        return meanMinorAlleleFraction;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final CRAFSegment that = (CRAFSegment) o;

        return numPointsCopyRatio == that.numPointsCopyRatio &&
                numPointsAlleleFraction == that.numPointsAlleleFraction &&
                Double.compare(that.meanLog2CopyRatio, meanLog2CopyRatio) == 0 &&
                Double.compare(that.meanMinorAlleleFraction, meanMinorAlleleFraction) == 0 &&
                interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + numPointsCopyRatio;
        result = 31 * result + numPointsAlleleFraction;
        temp = Double.doubleToLongBits(meanLog2CopyRatio);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(meanMinorAlleleFraction);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
