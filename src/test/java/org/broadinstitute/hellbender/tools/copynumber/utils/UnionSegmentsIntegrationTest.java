package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Sets;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.VersatileAnnotatedRegionParser;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class UnionSegmentsIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_RESOURCE_DIR = new File("src/test/resources/org/broadinstitute/hellbender/tools/copynumber/utils/");
    public static final String SEG1 = TEST_RESOURCE_DIR.getAbsolutePath() + "/seg1.tsv";
    public static final String SEG_GT = TEST_RESOURCE_DIR.getAbsolutePath() + "/seg_fake_gt.tsv";
    public static final String SEG1_DIFFERENT_HEADERS = TEST_RESOURCE_DIR.getAbsolutePath() + "/seg1_different_headers.tsv";

    @Test
    public void testRunWithExactSegments() throws IOException {
        // Segment intervals are the same in the input files.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("unionseg_", ".tsv");
        final Set<String> columnSet = Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", "Segment_Mean", "Segment_Call");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG1);
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG1_DIFFERENT_HEADERS);

        Utils.stream(columnSet.iterator()).forEach(s -> {
            arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
            arguments.add(s);
        });

        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final VersatileAnnotatedRegionParser parser = new VersatileAnnotatedRegionParser();
        final List<SimpleAnnotatedGenomicRegion> regions = parser.readAnnotatedRegions(outputFile, Sets.newHashSet("MEAN_LOG2_COPY_RATIO", "CALL", "Segment_Mean", "Segment_Call"));
        Assert.assertEquals(regions.size(), 4);
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().size() == columnSet.size()));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().keySet().containsAll(columnSet)));
        Assert.assertTrue(regions.stream().noneMatch(r -> r.getAnnotations().values().contains("")));
    }

    @Test
    public void testRunWithExactSameFiles() throws IOException {
        // Input files are exactly the same.  Therefore, the union should only generate more columns.
        final File outputFile = File.createTempFile("unionseg_", ".tsv");
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

        final Set<String> gtColumnSet = Sets.newHashSet("MEAN_LOG2_COPY_RATIO_1", "CALL_1", "MEAN_LOG2_COPY_RATIO_2", "CALL_2");
        Assert.assertEquals(regions.size(), 4);
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().size() == gtColumnSet.size()));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().keySet().containsAll(gtColumnSet)));
        Assert.assertTrue(regions.stream().noneMatch(r -> r.getAnnotations().values().contains("")));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().get("MEAN_LOG2_COPY_RATIO_1").equals(r.getAnnotations().get("MEAN_LOG2_COPY_RATIO_2"))));
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().get("CALL_1").equals(r.getAnnotations().get("CALL_2"))));
    }
    @Test
    public void testRunWithNotExactOverlaps() throws IOException {
        // This test is a bit more like the real world
        final File outputFile = File.createTempFile("unionseg_", ".tsv");
        final List<String> arguments = new ArrayList<>();
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG1_DIFFERENT_HEADERS);
        arguments.add("-" + ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME);
        arguments.add(SEG_GT);
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Mean");
        arguments.add("-" + UnionSegments.COLUMNS_OF_INTEREST_SHORT_NAME);
        arguments.add("Segment_Call");
        arguments.add("-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME);
        arguments.add(outputFile.getAbsolutePath());
        runCommandLine(arguments);

        Assert.assertTrue(outputFile.exists());

        final VersatileAnnotatedRegionParser parser = new VersatileAnnotatedRegionParser();
        final List<SimpleAnnotatedGenomicRegion> regions = parser.readAnnotatedRegions(outputFile, Sets.newHashSet("Segment_Mean_1", "Segment_Call_1", "Segment_Mean_2", "Segment_Call_2"));
        Assert.assertEquals(regions.size(), 13);
        Assert.assertTrue(regions.stream().allMatch(r -> r.getAnnotations().size() == 4));

        // Painstakingly made by hand
        Assert.assertEquals(regions.get(0), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 5000, 10000),
                ImmutableMap.of("Segment_Call_1", "", "Segment_Call_2", "0",	"Segment_Mean_1", "", "Segment_Mean_2", "-0.04")));
        Assert.assertEquals(regions.get(1), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 10001, 10500),
                ImmutableMap.of("Segment_Call_1", "0", "Segment_Call_2", "0",	"Segment_Mean_1", "-0.03", "Segment_Mean_2", "-0.04")));
        Assert.assertEquals(regions.get(2), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 10501, 52499),
                ImmutableMap.of("Segment_Call_1", "", "Segment_Call_2", "0",	"Segment_Mean_1", "", "Segment_Mean_2", "-0.04")));
        Assert.assertEquals(regions.get(3), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 52500, 60000),
                ImmutableMap.of("Segment_Call_1", "-", "Segment_Call_2", "0",	"Segment_Mean_1", "-0.76", "Segment_Mean_2", "-0.04")));
        Assert.assertEquals(regions.get(4), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 60001, 69999),
                ImmutableMap.of("Segment_Call_1", "-", "Segment_Call_2", "",	"Segment_Mean_1", "-0.76", "Segment_Mean_2", "")));
        Assert.assertEquals(regions.get(5), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 70000, 100000),
                ImmutableMap.of("Segment_Call_1", "-", "Segment_Call_2", "-",	"Segment_Mean_1", "-0.76", "Segment_Mean_2", "-0.8")));
        Assert.assertEquals(regions.get(6), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 100001, 109750),
                ImmutableMap.of("Segment_Call_1", "-", "Segment_Call_2", "",	"Segment_Mean_1", "-0.76", "Segment_Mean_2", "")));
        Assert.assertEquals(regions.get(7), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 109751, 119999),
                ImmutableMap.of("Segment_Call_1", "0", "Segment_Call_2", "",	"Segment_Mean_1", "-0.10", "Segment_Mean_2", "")));
        Assert.assertEquals(regions.get(8), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 120000, 220000),
                ImmutableMap.of("Segment_Call_1", "0", "Segment_Call_2", "0",	"Segment_Mean_1", "-0.10", "Segment_Mean_2", "-0.1")));
        Assert.assertEquals(regions.get(9), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 220001, 229999),
                ImmutableMap.of("Segment_Call_1", "0", "Segment_Call_2", "",	"Segment_Mean_1", "-0.10", "Segment_Mean_2", "")));
        Assert.assertEquals(regions.get(10), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 230000, 230500),
                ImmutableMap.of("Segment_Call_1", "0", "Segment_Call_2", "-",	"Segment_Mean_1", "-0.10", "Segment_Mean_2", "-0.8")));
        Assert.assertEquals(regions.get(11), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 230501, 258500),
                ImmutableMap.of("Segment_Call_1", "-", "Segment_Call_2", "-",	"Segment_Mean_1", "-0.60", "Segment_Mean_2", "-0.8")));
        Assert.assertEquals(regions.get(12), new SimpleAnnotatedGenomicRegion(new SimpleInterval("1", 258501, 300000),
                ImmutableMap.of("Segment_Call_1", "", "Segment_Call_2", "-",	"Segment_Mean_1", "", "Segment_Mean_2", "-0.8")));
    }
}