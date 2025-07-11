version 1.0

import "../structs.wdl"

task pbmm2_align_wgs {
  meta {
    description: "Align HiFi reads to a reference genome with pbmm2."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    prealigned_bam: {
      name: "Name of Output directory where the aligned bam is"
    }
    prealigned_bam_bai: {
      name: "Name of Output directory where the aligned bam is"
    }
    bam: {
      name: "HiFi reads (BAM)"
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
    strip_kinetics: {
      name: "Strip kinetics tags"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    bam_stats: {
      name: "BAM stats"
    }
  }

  input {
    String sample_id
    File bam
    File prealigned_bam
    File prealigned_bam_bai
    File ref_fasta
    File ref_index
    String ref_name

    Boolean strip_kinetics = true

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 1
  Int mem_gb    = ceil(threads * 4)
  Int disk_size = ceil(size(bam, "GB") * 3 + size(ref_fasta, "GB") + 70)
  String movie = basename(bam, ".bam")

  command <<<
    set -eo pipefail
    cat << EOF > extract_read_length_and_qual.py
    import math, pysam
    MAX_QV = 60
    save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
    bamin = pysam.AlignmentFile('~{bam}', check_sq=False)
    pysam.set_verbosity(save)  # restore warnings
    for b in bamin:
      errorrate = 1.0 - b.get_tag('rq')
      if math.isnan(b.get_tag('rq')):
        print(f'Warning: read {b.query_name} has tag rq:f:nan.', file=sys.stderr)
        continue
      readqv = MAX_QV if errorrate == 0 else math.floor(-10 * math.log10(errorrate))
      print(f"{b.query_name.split('/')[0]}\t{b.query_name}\t{len(b.query_sequence)}\t{readqv}")
    bamin.close()
    EOF

    python3 ./extract_read_length_and_qual.py \
    | gzip --stdout > ~{sample_id}.~{movie}.read_length_and_quality.tsv.gz &
    BAM_STATS_PID=$!

    cat << EOF > detect_bam_tags.py
    import json, pysam
    def check_bam_file(bam_file_path, n_records):
      output = dict()
      save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
      with pysam.AlignmentFile(bam_file_path, 'rb', check_sq=False) as bam_file:
        pysam.set_verbosity(save)  # restore warnings
        aligned = bool(bam_file.nreferences)
        unique_tags = set()
        for i, record in enumerate(bam_file):
          if i >= n_records: break
          unique_tags.update(tag[0] for tag in record.tags)
      output['kinetics'] = bool(unique_tags & {'fi', 'ri', 'fp', 'rp', 'ip', 'pw'})
      output['base_modification'] = bool(unique_tags & {'MM', 'ML', 'Mm', 'Ml'})
      output['aligned'] = aligned
      output['haplotagged'] = bool(unique_tags & {'HP', 'PS', 'PC'})
      return output
    print(json.dumps(check_bam_file('~{bam}', 10000)))
    EOF

    read -r kinetics base_modification aligned haplotagged <<< "$(python3 ./detect_bam_tags.py | jq -r '. | [.kinetics, .base_modification, .aligned, .haplotagged] | @tsv')"

    if [ "$aligned" = true ]; then
      echo "Input ~{basename(bam)} is already aligned.  Alignments and and haplotype tags will be stripped."
    fi

    if [ "$base_modification" = false ]; then
      echo "Input ~{basename(bam)} does not contain base modification tags.  5mCpG pileups will not be generated."
    fi

    if [ "$kinetics" = true ]; then
      echo "Input ~{basename(bam)} contains consensus kinetics tags."
      if [ "~{strip_kinetics}" = true ]; then
        echo "Kinetics will be stripped from the output."
      fi
    fi

    # cat << EOF > pbmm2Slurm.sh
    # module load smrtlink-sequel2
    #echo ~(readline)
    #job_Id=`head -n 1 slurm_singularity.log.txt | cut -d' ' -f 3`
    #SLURM_TMPDIR="/localscratch/felixant.${job_Id}.0"

    # pbmm2 --version
    # export TMPDIR=/mnt/tmp/  #Using local tmpdir files will greatly accelerate the sorting process

    # echo /mnt/tmp/
    # cp ~{ref_fasta} /mnt/tmp/
    # cp ~{ref_index} /mnt/tmp/
    # cp ~{bam} /mnt/tmp/
    # bam_file=`echo ~{bam} | rev | cut -d/ -f 1 | rev` #name of the file without the path
    # ref_fasta_file=`echo ~{ref_fasta} | rev | cut -d/ -f 1 | rev`

    # pbmm2 align \
    #         --num-threads 24 \
    #         --sort-memory 6G \
    #         --preset HIFI \
    #         --sample GM232447_SPRQ \
    #         --log-level INFO \
    #         --log-level DEBUG \
    #         --sort \
    #         --unmapped \
    #         /mnt/tmp/${ref_fasta_file} \
    #         /mnt/tmp/${bam_file} \
    #         /mnt/tmp/aligned.bam

    # cp /mnt/tmp//aligned.bam .
    #EOF


    #srun --time=12:00:00 --account=def-rallard --cpus-per-task=24 --mem=90G /home/felixant/scratch/MINIWDL/pbmm2Slurm.sh ~{ref_fasta} ~{bam}
    # pbmm2 --version

    # pbmm2 align \
    #   --num-threads ~{threads} \
    #   --sort-memory 4G \
    #   --preset HIFI \
    #   --sample ~{sample_id} \
    #   --log-level INFO \
    #   --sort \
    #   ~{true='--strip' false='' strip_kinetics} \
    #   --unmapped \
    #   ~{ref_fasta} \
    #   ~{bam} \
    #   aligned.bam

    # if [ "$haplotagged" = true ]; then
    #   # remove haplotype tags
    #   samtools view \
    #     ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
    #     --bam --no-PG \
    #     --remove-tag HP,PS,PC \
    #     -o ~{sample_id}.~{movie}.~{ref_name}.aligned.bam \
    #     aligned.bam \
    #   && rm --verbose aligned.bam aligned.bam.bai
    #   samtools index \
    #     ~{if threads > 1 then "-@ " + (threads - 1) else ""} \
    #     ~{sample_id}.~{movie}.~{ref_name}.aligned.bam
    # else
    #   mv --verbose aligned.bam ~{sample_id}.~{movie}.~{ref_name}.aligned.bam
    #   mv --verbose aligned.bam.bai ~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai
    # fi

    cp --verbose ~{prealigned_bam} ./aligned.bam
    cp --verbose ~{prealigned_bam_bai} ./aligned.bam.bai
    mv --verbose aligned.bam ~{sample_id}.~{movie}.~{ref_name}.aligned.bam
    mv --verbose aligned.bam.bai ~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai

    wait ${BAM_STATS_PID}
  >>>

  output {
    File aligned_bam       = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
    File aligned_bam_index = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai"
    File bam_stats         = "~{sample_id}.~{movie}.read_length_and_quality.tsv.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:24218cb5cbc68d1fd64db14a9dc38263d3d931c74aca872c998d12ef43020ef0"
    cpu: threads
    memory: mem_gb + " GB"
    time_minutes: "30"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: 0
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}