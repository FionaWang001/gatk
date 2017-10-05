package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.VersatileAnnotatedRegionParser;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class UnionSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/copynumber/utils/");
    public static final String SEG1 = TEST_RESOURCE_DIR.getAbsolutePath() + "/seg1.tsv";
    public static final String SEG1_DIFFERENT_HEADERS = TEST_RESOURCE_DIR.getAbsolutePath() + "/seg1_different_headers.tsv";

    @Test
    public void testRunWithExactSegments() throws IOException {
        // Segment intervals are the same in the input files.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("ob_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG1);
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG1_DIFFERENT_HEADERS);
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("MEAN_LOG2_COPY_RATIO");
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("CALL");
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Mean");
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Call");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final VersatileAnnotatedRegionParser parser = new VersatileAnnotatedRegionParser();
        final List<SimpleAnnotatedGenomicRegion> regions = parser.readAnnotatedRegions(outputFile, Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", "Segment_Mean", "Segment_Call"));
        Assert.assertEquals(regions.size(), 4);

        //TODO: More tests here
    }

    @Test
    public void testRunWithExactSameFiles() throws IOException {
        // Input files are exactly the same.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("ob_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG1);
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG1);
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("MEAN_LOG2_COPY_RATIO");
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("CALL");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final VersatileAnnotatedRegionParser parser = new VersatileAnnotatedRegionParser();
        final List<SimpleAnnotatedGenomicRegion> regions = parser.readAnnotatedRegions(outputFile, Sets.newHashSet("MEAN_LOG2_COPY_RATIO_1", "CALL_1", "MEAN_LOG2_COPY_RATIO_2", "CALL_2"));
        Assert.assertEquals(regions.size(), 4);

        // TODO: More tests here (mimic the ones above but with different column names
    }

    // TODO: Add tests for annotations that do not match.
}
