version 1.0

import "../structs.wdl"

task methbat {
  meta {
    description: "Profile methylation regions with methbat."
  }

  parameter_meta {
    sample_prefix: {
      name: "Sample prefix"
    }
    methylation_pileup_beds: {
      name: "Methylation pileup BED files"
    }
    region_tsv: {
      name: "Regions TSV file"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
  }

  input {
    String sample_prefix

    Array[File] methylation_pileup_beds

    File region_tsv

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil((size(methylation_pileup_beds, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    methbat --version

    methbat profile \
      --input-prefix ~{sample_prefix} \
      --input-regions ~{region_tsv} \
      --output-region-profile ~{out_prefix}.methbat.profile.tsv
  >>>

  output {
    File profile = "~{out_prefix}.methbat.profile.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/methbat@sha256:7cf8f17ced418a3b8e1bec7522acc2d84cc73fb29776c5aaf2a02784ede7c52b"
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform
  }
}
