# Workflow for running the GATK CNV pipeline on a single normal or tumor BAM. Supports both WGS and WES.
#
# Notes:
#
# - The target file (targets) is required for the WES workflow and should be a TSV file with the column headers:
#    contig    start    stop    name
#   These targets will be padded on both sides by the amount specified by PadTargets.padding (default 250).
#
# - If a target file is not provided, then the WGS workflow will be run instead and the specified value of
#   wgs_bin_size (default 10000) will be used.
#
# - The sites file (common_sites) is must be a Picard interval list.
#
# - Example invocation:
#    java -jar cromwell.jar run cnv_somatic_bam_workflow.wdl myParameters.json
#   See cnv_somatic_bam_workflow_template.json for a template json file to modify with your own parameters (please save
#   your modified version with a different filename and do not commit to the gatk repository).
#
#############

import "cnv_common_tasks.wdl" as CNVTasks
import "cnv_somatic_oncotate.wdl" as Oncotate

workflow CNVSomaticBAMWorkflow {
    # Workflow input files
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

    Boolean is_run_oncotator = false

    # docker images
    String gatk_docker
    String oncotator_docker = "broadinstitute/oncotator:1.9.3.0-eval-gatk-protected"

    if (!is_wgs) {
        call CNVTasks.PadTargets {
            input:
                # The task will fail if targets is not defined when it gets here, but that should not be allowed to happen.
                targets = select_first([targets, ""]),
                gatk_jar = gatk_jar,
                gatk_docker = gatk_docker
        }
    }

    call CNVTasks.CollectReadCounts {
        input:
            padded_targets = padded_targets,
            bam = bam,
            bam_idx = bam_idx,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker,
            is_wgs = is_wgs
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
            read_counts = CollectReadCounts.read_counts,
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

    call PlotDenoisedCopyRatios  {
        input:
            entity_id = CollectReadCounts.entity_id,
            standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratios,
            denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios,
            ref_fasta_dict = ref_fasta_dict,
            gatk_jar = gatk_jar,
            gatk_docker = gatk_docker
    }

    if (is_run_oncotator) {
        call Oncotate.CNVOncotateCalledSegments as OncotateCalledCNVWorkflow {
            input:
                 called_file = CallCopyRatioSegments.called_copy_ratio_segments,
                 oncotator_docker = oncotator_docker
        }
    }

    output {
        String entity_id = CollectReadCounts.entity_id
        File read_counts = CollectReadCounts.read_counts
        File allelic_counts = CollectAllelicCounts.allelic_counts
        File standardized_copy_ratios = DenoiseReadCounts.standardized_copy_ratio
        File denoised_copy_ratios = DenoiseReadCounts.denoised_copy_ratios
        File copy_ratio_segments = ModelSegments.copy_ratio_segments
        File allele_fraction_segments = ModelSegments.allele_fraction_segments
        #File unioned_segments = ModelSegments.unioned_segments
        File called_copy_ratio_segments = CallCopyRatioSegments.called_copy_ratio_segments
        File denoised_copy_ratios_plot = PlotDenoisedCopyRatios.denoised_copy_ratios_plot
        File denoised_copy_ratios_lim_4_plot = PlotDenoisedCopyRatios.denoised_copy_ratios_lim_4_plot

        File? oncotated_called_file = OncotateCalledCNVWorkflow.oncotated_called_file
    }
}

# Denoise the coverage
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

# Segment the denoised copy-ratio profile
task ModelSegments {
    String entity_id
    File denoised_copy_ratios
    File allelic_counts
    Int? max_num_segments_per_chromosome
    Float? kernel_variance_copy_ratio
    Float? kernel_variance_allele_fraction
    Int? kernel_approximation_dimension
    Array[Int]? window_sizes = [8, 16, 32, 64, 128, 256]
    Float? num_changepoints_penalty_factor_copy_ratio
    Float? num_changepoints_penalty_factor_allele_fraction
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default="4" mem}g -jar ${gatk_jar} ModelSegments \
            --denoisedCopyRatios ${denoised_copy_ratios} \
            --allelicCounts ${allelic_counts} \
            --maxNumSegmentsPerChromosome ${default="50" max_num_segments_per_chromosome} \
            --kernelVarianceCopyRatio ${default="0.0" kernel_variance_copy_ratio} \
            --kernelVarianceAlleleFraction ${default="1.0" kernel_variance_allele_fraction} \
            --kernelApproximationDimension ${default="100" kernel_approximation_dimension} \
            --windowSizes ${sep= " --windowSizes " window_sizes} \
            --numChangepointsPenaltyFactorCopyRatio ${default="1.0" num_changepoints_penalty_factor_copy_ratio} \
            --numChangepointsPenaltyFactorAlleleFraction ${default="1.0" num_changepoints_penalty_factor_allele_fraction} \
            --outputPrefix ${entity_id}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File copy_ratio_segments = "${entity_id}.cr.seg"
        File allele_fraction_segments = "${entity_id}.af.seg"
    }
}

# Make calls (amplified, neutral, or deleted) on each segment
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

# Create plots of denoised copy ratios
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
        File denoised_copy_ratios_plot = "${output_dir_}/${entity_id}_Before_After.png"
        File denoised_copy_ratios_lim_4_plot = "${output_dir_}/${entity_id}_Before_After_CR_Lim_4.png"
    }
}