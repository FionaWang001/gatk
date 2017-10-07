# Workflow for running the GATK CNV pipeline on a single normal or tumor BAM. Supports both WGS and WES.
#
# Notes:
#
# - The target file (targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_length (default 10000) will be used.
#
# - The sites file (common_sites) should be a Picard or GATK-style interval list.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_somatic_bam_workflow.wdl myParameters.json
#   See cnv_somatic_bam_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks

workflow CNVSomaticBAMWorkflow {
    File? targets
    File common_sites
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File read_count_pon
    String gatk_jar

    # If no target file is input, then do WGS workflow
    Boolean is_wgs = select_first([targets, ""]) == ""

    String gatk_docker

    if (!is_wgs) {
        call CNVTasks.PadTargets {
            input:
                targets = select_first([targets, ""]),
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker
        }
    }

    call CNVTasks.CollectReadCounts {
        input:
            padded_targets = PadTargets.padded_targets,
            bam = bam,
            bam_idx = bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }
    
    call CNVTasks.CollectAllelicCounts {
        input:
            common_sites = common_sites,
            bam = bam,
            bam_idx = bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call DenoiseReadCounts {
        input:
            entity_id = CollectReadCounts.entity_id,
            read_counts = if is_wgs then CollectReadCounts.read_counts_hdf5 else CollectReadCounts.read_counts,
            read_count_pon = read_count_pon,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call ModelSegments {
        input:
            entity_id = CollectReadCounts.entity_id,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            allelic_counts = CollectAllelicCounts.allelic_counts,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call CallCopyRatioSegments {
        input:
            entity_id = CollectReadCounts.entity_id,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            copy_ratio_segments = ModelSegments.copy_ratio_segments,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call PlotDenoisedCopyRatios {
        input:
            entity_id = CollectReadCounts.entity_id,
            standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratios,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    call PlotModeledSegments {
        input:
            entity_id = CollectReadCounts.entity_id,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            het_allelic_counts = ModelSegments.het_allelic_counts,
            modeled_segments = ModelSegments.modeled_segments,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    output {
        String entity_id = CollectReadCounts.entity_id
        File read_counts = CollectReadCounts.read_counts
        File allelic_counts = CollectAllelicCounts.allelic_counts
        File standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratios
        File denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios
        File het_allelic_counts = ModelSegments.het_allelic_counts
        File copy_ratio_segments = ModelSegments.copy_ratio_segments
        File allele_fraction_segments = ModelSegments.allele_fraction_segments
        File modeled_segments = ModelSegments.modeled_segments
        File copy_ratio_parameters = ModelSegments.copy_ratio_parameters
        File allele_fraction_parameters = ModelSegments.allele_fraction_parameters
        File called_copy_ratio_segments = CallCopyRatioSegments.called_copy_ratio_segments
        File denoised_copy_ratios_plot = PlotDenoisedCopyRatios.denoised_copy_ratios_plot
        File denoised_copy_ratios_lim_4_plot = PlotDenoisedCopyRatios.denoised_copy_ratios_lim_4_plot
        File modeled_segments_plot = PlotModeledSegments.modeled_segments_plot
    }
}

task DenoiseReadCounts {
    String entity_id
    File read_counts
    File read_count_pon
    Int? number_of_eigensamples #use all eigensamples in panel by default
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default="4" mem}g -jar ${gatk_jar} DenoiseReadCounts \
            --input ${read_counts} \
            --readCountPanelOfNormals ${read_count_pon} \
            ${"--numberOfEigensamples " + number_of_eigensamples} \
            --standardizedCopyRatios ${entity_id}.standardizedCR.tsv \
            --denoisedCopyRatios ${entity_id}.denoisedCR.tsv
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(read_count_pon, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File standardized_copy_ratios = "${entity_id}.standardizedCR.tsv"
        File denoised_copy_ratios = "${entity_id}.denoisedCR.tsv"
    }
}

task ModelSegments {
    String entity_id
    File denoised_copy_ratios
    File allelic_counts
    Int? max_num_segments_per_chromosome
    Int? min_total_allele_count
    Float? genotyping_p_value_threshold
    Float? genotyping_base_error_rate
    Float? kernel_variance_copy_ratio
    Float? kernel_variance_allele_fraction
    Int? kernel_approximation_dimension
    Array[Int]? window_sizes = [8, 16, 32, 64, 128, 256]
    Float? num_changepoints_penalty_factor_copy_ratio
    Float? num_changepoints_penalty_factor_allele_fraction
    Float? minor_allele_fraction_prior_alpha
    Int? num_samples_copy_ratio
    Int? num_burn_in_copy_ratio
    Int? num_samples_allele_fraction
    Int? num_burn_in_allele_fraction
    Float? smoothing_threshold_copy_ratio
    Float? smoothing_threshold_allele_fraction
    Int? max_num_smoothing_iterations
    Int? num_smoothing_iterations_per_fit
    String? output_dir
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # If optional output_dir not specified, use "."
    String output_dir_ = select_first([output_dir, "."])

    command {
        java -Xmx${default="4" mem}g -jar ${gatk_jar} ModelSegments \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            --allelicCounts ${allelic_counts} \
            --maxNumSegmentsPerChromosome ${default="500" max_num_segments_per_chromosome} \
            --minTotalAlleleCount ${default="30" min_total_allele_count} \
            --genotypingPValueThreshold ${default="0.001" genotyping_p_value_threshold} \
            --genotypingBaseErrorRate ${default="0.01" genotyping_base_error_rate} \
            --kernelVarianceCopyRatio ${default="0.0" kernel_variance_copy_ratio} \
            --kernelVarianceAlleleFraction ${default="0.01" kernel_variance_allele_fraction} \
            --kernelApproximationDimension ${default="100" kernel_approximation_dimension} \
            --windowSizes ${sep= " --windowSizes " window_sizes} \
            --numChangepointsPenaltyFactorCopyRatio ${default="1.0" num_changepoints_penalty_factor_copy_ratio} \
            --numChangepointsPenaltyFactorAlleleFraction ${default="1.0" num_changepoints_penalty_factor_allele_fraction} \
            --minorAlleleFractionPriorAlpha ${default="25.0" minor_allele_fraction_prior_alpha} \
            --numSamplesCopyRatio ${default=100 num_samples_copy_ratio} \
            --numBurnInCopyRatio ${default=50 num_burn_in_copy_ratio} \
            --numSamplesAlleleFraction ${default=100 num_samples_allele_fraction} \
            --numBurnInAlleleFraction ${default=50 num_burn_in_allele_fraction} \
            --smoothingThresholdCopyRatio ${default="4.0" smoothing_threshold_copy_ratio} \
            --smoothingThresholdAlleleFraction ${default="2.0" smoothing_threshold_allele_fraction} \
            --maxNumSmoothingIterations ${default=10 max_num_smoothing_iterations} \
            --numSmoothingIterationsPerFit ${default=0 num_smoothing_iterations_per_fit} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File het_allelic_counts = "${output_dir_}/${entity_id}.hets.tsv"
        File copy_ratio_segments = "${output_dir_}/${entity_id}.cr.seg"
        File allele_fraction_segments = "${output_dir_}/${entity_id}.af.seg"
        File combined_segments = "${output_dir_}/${entity_id}.craf.seg"
        File modeled_segments = "${output_dir_}/${entity_id}.modelFinal.seg"
        File copy_ratio_parameters = "${output_dir_}/${entity_id}.modelFinal.cr.param"
        File allele_fraction_parameters = "${output_dir_}/${entity_id}.modelFinal.af.param"
    }
}

task CallCopyRatioSegments {
    String entity_id
    File denoised_copy_ratios
    File copy_ratio_segments
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default="4" mem}g -jar ${gatk_jar} CallCopyRatioSegments \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            --segments ${copy_ratio_segments} \
            --output ${entity_id}.called.seg
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File called_copy_ratio_segments = "${entity_id}.called.seg"
    }
}

task PlotDenoisedCopyRatios {
    String entity_id
    File standardized_copy_ratios
    File denoised_copy_ratios
    File ref_fasta_dict
    Int? minimum_contig_length
    String? output_dir
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # If optional output_dir not specified, use "."
    String output_dir_ = select_first([output_dir, "."])

    command {
        mkdir -p ${output_dir_}; \
        java -Xmx${default="4" mem}g -jar ${gatk_jar} PlotDenoisedCopyRatios \
            --standardizedCopyRatios ${standardized_copy_ratios} \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            -SD ${ref_fasta_dict} \
            --minimumContigLength ${default="1000000" minimum_contig_length} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File denoised_copy_ratios_plot = "${output_dir_}/${entity_id}.denoised.png"
        File denoised_copy_ratios_lim_4_plot = "${output_dir_}/${entity_id}.denoisedLimit4.png"
    }
}

task PlotModeledSegments {
    String entity_id
    File denoised_copy_ratios
    File het_allelic_counts
    File modeled_segments
    File ref_fasta_dict
    Int? minimum_contig_length
    String? output_dir
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # If optional output_dir not specified, use "."
    String output_dir_ = select_first([output_dir, "."])

    command {
        mkdir -p ${output_dir_}; \
        java -Xmx${default="4" mem}g -jar ${gatk_jar} PlotModeledSegments \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            --allelicCounts ${het_allelic_counts} \
            --segments ${modeled_segments} \
            -SD ${ref_fasta_dict} \
            --minimumContigLength ${default="1000000" minimum_contig_length} \
            --output ${output_dir_} \
            --outputPrefix ${entity_id}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File modeled_segments_plot = "${output_dir_}/${entity_id}.modeled.png"
    }
}