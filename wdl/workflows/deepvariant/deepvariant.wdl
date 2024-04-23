version 1.0

import "../../structs.wdl"

workflow deepvariant {
  meta {
    description: "Call variants from aligned HiFi reads using DeepVariant"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex"
    }
    aligned_bams: {
      name: "Aligned BAM"
    }
    aligned_bam_indices: {
      name: "Aligned BAI"
    }
    regions_bed: {
      name: "Regions BED"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FAI"
    }
    ref_name: {
      name: "Reference name"
    }
    ref_par_regions_bed: {
      name: "Reference PAR regions BED"
    }
    ref_allosomes: {
      name: "Comma delimited list of allosome names"
    }
    deepvariant_version: {
      name: "DeepVariant Version"
    }
    custom_deepvariant_model_tar: {
      name: "Custom DeepVariant Model tar"
    }
    gpu: {
      name: "Use GPU for DeepVariant call_variants"
    }
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
    vcf: {
      name: "VCF"
    }
    vcf_index: {
      name: "VCF index"
    }
    gvcf: {
      name: "GVCF"
    }
    gvcf_index: {
      name: "GVCF index"
    }
  }

  input {
    String sample_id
    String? sex
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    File? regions_bed

    File ref_fasta
    File ref_index
    String ref_name
    File? ref_par_regions_bed
    String? ref_allosomes

    String deepvariant_version
    File? custom_deepvariant_model_tar

    Boolean gpu

    RuntimeAttributes default_runtime_attributes
  }

  Int total_deepvariant_tasks = 64
  Int num_shards              = 8
  Int tasks_per_shard         = total_deepvariant_tasks / num_shards

  Boolean adjust_allosome_ploidy = defined(ref_par_regions_bed) && defined(ref_allosomes) && defined(sex) && (sex == "MALE")

  String docker_image = if (default_runtime_attributes.backend == "AWS-HealthOmics") then default_runtime_attributes.container_registry else "google" + "/deepvariant:~{deepvariant_version}"

  scatter (shard_index in range(num_shards)) {
    Int task_start_index = shard_index * tasks_per_shard

    call deepvariant_make_examples {
      input:
        sample_id               = sample_id,
        aligned_bams            = aligned_bams,
        aligned_bam_indices     = aligned_bam_indices,
        regions_bed             = regions_bed,
        ref_fasta               = ref_fasta,
        ref_index               = ref_index,
        task_start_index        = task_start_index,
        tasks_per_shard         = tasks_per_shard,
        total_deepvariant_tasks = total_deepvariant_tasks,
        docker_image            = docker_image,
        runtime_attributes      = default_runtime_attributes
    }
  }

  call deepvariant_call_variants {
    input:
      sample_id                    = sample_id,
      ref_name                     = ref_name,
      example_tfrecord_tars        = deepvariant_make_examples.example_tfrecord_tar,
      custom_deepvariant_model_tar = custom_deepvariant_model_tar,
      total_deepvariant_tasks      = total_deepvariant_tasks,
      docker_image                 = docker_image + if gpu then "-gpu" else "",
      gpu                          = gpu,
      runtime_attributes           = default_runtime_attributes
  }

  call deepvariant_postprocess_variants {
    input:
      sample_id                     = sample_id,
      tfrecords_tar                 = deepvariant_call_variants.tfrecords_tar,
      nonvariant_site_tfrecord_tars = deepvariant_make_examples.nonvariant_site_tfrecord_tar,
      ref_fasta                     = ref_fasta,
      ref_index                     = ref_index,
      ref_name                      = ref_name,
      ref_par_regions_bed           = ref_par_regions_bed,
      ref_allosomes                 = ref_allosomes,
      total_deepvariant_tasks       = total_deepvariant_tasks,
      docker_image                  = docker_image,
      adjust_allosome_ploidy        = adjust_allosome_ploidy,
      runtime_attributes            = default_runtime_attributes
  }

  output {
    File vcf        = deepvariant_postprocess_variants.vcf
    File vcf_index  = deepvariant_postprocess_variants.vcf_index
    File gvcf       = deepvariant_postprocess_variants.gvcf
    File gvcf_index = deepvariant_postprocess_variants.gvcf_index
  }
}

task deepvariant_make_examples {
  meta {
    description: "Run DeepVariant make_examples step"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    aligned_bams: {
      name: "Aligned BAM"
    }
    aligned_bam_indices: {
      name: "Aligned BAM index"
    }
    regions_bed: {
      name: "Regions BED"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    task_start_index: {
      name: "Task start index"
    }
    tasks_per_shard: {
      name: "Tasks per shard"
    }
    total_deepvariant_tasks: {
      name: "Total DeepVariant tasks"
    }
    docker_image: {
      name: "Docker image URL"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    example_tfrecord_tar: {
      name: "Example TFRecord tar"
    }
    nonvariant_site_tfrecord_tar: {
      name: "Nonvariant Site TFRecord tar"
    }
  }

  input {
    String sample_id
    Array[File] aligned_bams
    Array[File] aligned_bam_indices
    File? regions_bed

    File ref_fasta
    File ref_index

    Int task_start_index
    Int tasks_per_shard

    Int total_deepvariant_tasks
    String docker_image

    RuntimeAttributes runtime_attributes
  }

  Int task_end_index = task_start_index + tasks_per_shard - 1
  Int disk_size      = ceil(size(aligned_bams, "GB") * 2 + size(ref_fasta, "GB") + 20)
  Int mem_gb         = tasks_per_shard * 4

  command <<<
    set -euo pipefail

    mkdir example_tfrecords nonvariant_site_tfrecords

    echo "DeepVariant version: $VERSION"

    seq ~{task_start_index} ~{task_end_index} \
    | parallel \
      --jobs ~{tasks_per_shard} \
      --halt 2 \
      /opt/deepvariant/bin/make_examples \
        --norealign_reads \
        --vsc_min_fraction_indels 0.12 \
        --pileup_image_width 199 \
        --track_ref_reads \
        --phase_reads \
        --partition_size=25000 \
        --max_reads_per_partition=600 \
        --alt_aligned_pileup=diff_channels \
        --add_hp_channel \
        --sort_by_haplotypes \
        --parse_sam_aux_fields \
        --min_mapping_quality=1 \
        --mode calling \
        --ref ~{ref_fasta} \
        ~{if defined(regions_bed) then "--regions " + regions_bed else ""} \
        --reads ~{sep="," aligned_bams} \
        --examples example_tfrecords/~{sample_id}.examples.tfrecord@~{total_deepvariant_tasks}.gz \
        --gvcf nonvariant_site_tfrecords/~{sample_id}.gvcf.tfrecord@~{total_deepvariant_tasks}.gz \
        --task {}

    tar -zcvf ~{sample_id}.~{task_start_index}.example_tfrecords.tar.gz example_tfrecords \
      && rm -rf example_tfrecords
    tar -zcvf ~{sample_id}.~{task_start_index}.nonvariant_site_tfrecords.tar.gz nonvariant_site_tfrecords \
      && rm -rf nonvariant_site_tfrecords
  >>>

  output {
    File example_tfrecord_tar         = "~{sample_id}.~{task_start_index}.example_tfrecords.tar.gz"
    File nonvariant_site_tfrecord_tar = "~{sample_id}.~{task_start_index}.nonvariant_site_tfrecords.tar.gz"
  }

  runtime {
    docker: docker_image
    cpu: tasks_per_shard
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}

task deepvariant_call_variants {
  meta {
    description: "Run DeepVariant call_variants step"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    ref_name: {
      name: "Reference name"
    }
    example_tfrecord_tars: {
      name: "Example TFRecord tars"
    }
    total_deepvariant_tasks: {
      name: "Total DeepVariant tasks"
    }
    docker_image: {
      name: "Docker image URL"
    }
    custom_deepvariant_model_tar: {
      name: "Custom DeepVariant Model tar"
    }
    gpu: {
      name: "Use GPU for DeepVariant call_variants"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    tfrecords_tar: {
      name: "TFRecords tar"
    }
  }

  input {
    String sample_id
    String ref_name
    Array[File] example_tfrecord_tars

    File? custom_deepvariant_model_tar
    Int total_deepvariant_tasks
    String docker_image

    Boolean gpu = false

    RuntimeAttributes runtime_attributes
  }

  Int threads        = if gpu then 8 else total_deepvariant_tasks
  Int writer_threads = if gpu then 4 else 8
  Int mem_gb         = if gpu then 32 else (total_deepvariant_tasks * 4)
  Int disk_size      = ceil(size(example_tfrecord_tars, "GB") * 2 + 100)

  command <<<
    set -euo pipefail

    while read -r tfrecord_tar || [[ -n "${tfrecord_tar}" ]]; do
      tar -zxvf "${tfrecord_tar}"
    done < ~{write_lines(example_tfrecord_tars)}

    if ~{defined(custom_deepvariant_model_tar)}; then
      mkdir -p ./custom_deepvariant_model
      tar --no-same-owner -zxvf ~{custom_deepvariant_model_tar} -C ./custom_deepvariant_model
      DEEPVARIANT_MODEL="./custom_deepvariant_model"
    else
      DEEPVARIANT_MODEL="/opt/models/pacbio"
    fi

    echo "DeepVariant version: $VERSION"
    echo "DeepVariant model: $DEEPVARIANT_MODEL"

    /opt/deepvariant/bin/call_variants \
      --writer_threads ~{writer_threads} \
      --outfile ~{sample_id}.~{ref_name}.call_variants_output.tfrecord.gz \
      --examples "example_tfrecords/~{sample_id}.examples.tfrecord@~{total_deepvariant_tasks}.gz" \
      --checkpoint "${DEEPVARIANT_MODEL}"

    tar -zcvf ~{sample_id}.~{ref_name}.call_variants_output.tar.gz ~{sample_id}.~{ref_name}.call_variants_output*.tfrecord.gz \
      && rm ~{sample_id}.~{ref_name}.call_variants_output*.tfrecord.gz \
      && rm -rf example_tfrecords \
      && rm -rf ./custom_deepvariant_model
  >>>

  output {
    File tfrecords_tar = "~{sample_id}.~{ref_name}.call_variants_output.tar.gz"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    gpu: gpu
    gpuCount: if gpu then 1 else 0
    gpuType: if gpu then runtime_attributes.gpuType else ""
    acceleratorCount: if gpu then 1 else 0  # !UnknownRuntimeKey
    acceleratorType: if gpu then runtime_attributes.gpuType else ""  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task deepvariant_postprocess_variants {
  meta {
    description: "Run DeepVariant postprocess_variants step"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    tfrecords_tar: {
      name: "TFRecords tar"
    }
    nonvariant_site_tfrecord_tars: {
      name: "Nonvariant Site TFRecord tars"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    ref_name: {
      name: "Reference name"
    }
    ref_par_regions_bed: {
      name: "Reference PAR regions BED"
    }
    ref_allosomes: {
      name: "Comma delimited list of allosome names"
    }
    total_deepvariant_tasks: {
      name: "Total DeepVariant tasks"
    }
    docker_image: {
      name: "Docker image URL"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    vcf : {
      name: "VCF"
    }
    vcf_index : {
      name: "VCF index"
    }
    gvcf : {
      name: "gVCF"
    }
    gvcf_index : {
      name: "gVCF index"
    }
  }

  input {
    String sample_id
    File tfrecords_tar
    Array[File] nonvariant_site_tfrecord_tars

    File ref_fasta
    File ref_index
    String ref_name
    File? ref_par_regions_bed
    String? ref_allosomes

    Int total_deepvariant_tasks
    String docker_image

    Boolean adjust_allosome_ploidy = false
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 6
  Int mem_gb    = 64
  Int disk_size = ceil((size(tfrecords_tar, "GB") + size(ref_fasta, "GB") + size(nonvariant_site_tfrecord_tars, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    tar -zxvf "~{tfrecords_tar}"

    while read -r nonvariant_site_tfrecord_tar || [[ -n "${nonvariant_site_tfrecord_tar}" ]]; do
      tar -zxvf "${nonvariant_site_tfrecord_tar}"
    done < ~{write_lines(nonvariant_site_tfrecord_tars)}

    echo "DeepVariant version: $VERSION"

    /opt/deepvariant/bin/postprocess_variants \
      --cpus ~{threads} \
      --vcf_stats_report=false \
      --ref ~{ref_fasta} \
      --infile ~{sample_id}.~{ref_name}.call_variants_output.tfrecord.gz \
      --outfile ~{sample_id}.~{ref_name}.small_variants.vcf.gz \
      --nonvariant_site_tfrecord_path "nonvariant_site_tfrecords/~{sample_id}.gvcf.tfrecord@~{total_deepvariant_tasks}.gz" \
      --gvcf_outfile ~{sample_id}.~{ref_name}.small_variants.g.vcf.gz

      # ~{if adjust_allosome_ploidy then "--par_regions_bed " + ref_par_regions_bed else ""} \
      # ~{if adjust_allosome_ploidy then "--haploid_contigs " + ref_allosomes else ""} \

    rm ~{sample_id}.~{ref_name}.call_variants_output*.tfrecord.gz \
      && rm -rf nonvariant_site_tfrecords
  >>>

  output {
    File vcf        = "~{sample_id}.~{ref_name}.small_variants.vcf.gz"
    File vcf_index  = "~{sample_id}.~{ref_name}.small_variants.vcf.gz.tbi"
    File gvcf       = "~{sample_id}.~{ref_name}.small_variants.g.vcf.gz"
    File gvcf_index = "~{sample_id}.~{ref_name}.small_variants.g.vcf.gz.tbi"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}