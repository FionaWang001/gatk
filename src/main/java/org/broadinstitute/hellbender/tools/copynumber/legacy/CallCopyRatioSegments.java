package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller.ReCapSegCaller;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.LegacyCopyNumberArgument;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Calls segments as amplified, deleted or copy number neutral given files containing denoised copy ratios
 * and a list of segments.
 *
 * @author David Benjamin
 *
 * <h3>Examples</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CallCopyRatioSegments \
 *   --denoisedCopyRatios tumor.denoisedCR.tsv \
 *   --segments tumor.seg \
 *   --output tumor.called
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Call copy-ratio segments as amplified, deleted, or copy number neutral.",
        oneLineSummary = "Call copy-ratio segments as amplified, deleted, or copy number neutral.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class CallCopyRatioSegments extends CommandLineProgram {

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_FULL_NAME,
            shortName = LegacyCopyNumberArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME
    )
    private File inputDenoisedCopyRatiosFile;

    @Argument(
            doc = "Input file containing copy-ratio segments (output of ModelSegments).",
            fullName = LegacyCopyNumberArgument.SEGMENTS_FILE_FULL_NAME,
            shortName = LegacyCopyNumberArgument.SEGMENTS_FILE_SHORT_NAME
    )
    protected File segmentsFile;

    @Argument(
            doc = "Output file for called segments.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected File outFile;

    @Override
    protected Object doWork() {
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
        final CopyRatioSegmentCollection segments = new CopyRatioSegmentCollection(segmentsFile);

        final CalledCopyRatioSegmentCollection calledSegments = new ReCapSegCaller(denoisedCopyRatios, segments).makeCalls();
        calledSegments.write(outFile, denoisedCopyRatios.getSampleName());
        return "SUCCESS";
    }
}
