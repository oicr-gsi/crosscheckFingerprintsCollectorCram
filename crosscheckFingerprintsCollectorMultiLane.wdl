version 1.0

# ============================================================
#  crosscheckFingerprintsCollectorMultiLane
#
#  Aligned-input, multi-lane variant of
#  crosscheckFingerprintsCollector.
#
#  The general-purpose workflow accepts fastq, bam or cram and
#  therefore carries six optional File inputs plus inputType /
#  aligner discriminators, and pulls in the bwaMem and star
#  subworkflows. This version takes an already-aligned merged-
#  lanes file - either a cram or a bam - so the input signature
#  is flat and no alignment imports are needed.
#
#  The input is always assumed to hold more than one read group:
#  it is filtered to the fingerprint intervals, split by read
#  group, and every lane is then processed independently
#  (optional duplicate marking, then fingerprint + metrics).
#  There is no lane-level input path; use the general-purpose
#  workflow for an already-split file.
# ============================================================

struct OutputGroup {
   String limsId
   File outputVcf
   File outputTbi
   File json
   File samstats
}

struct GenomeResources {
    String refFasta
    String refHapMap
    String intervalBed
    String intervalsToParallelizeByString
    String filterBamModules
    String splitLanesModules
    String markDuplicatesModules
    String mergeBamsModules
    String alignmentMetricsModules
    String extractFingerprintModules
}

workflow crosscheckFingerprintsCollectorMultiLane {
   input {
        File? cram
        File? cramIndex
        File? bam
        File? bamIndex
        Boolean markDups
        Boolean filterBam
        String outputFileNamePrefix
        String reference
        Int maxReads = 0
        String sampleId
   }
   parameter_meta {
        cram: "merged-lanes cram file; supply this with cramIndex, or bam with bamIndex"
        cramIndex: "index (.crai) for the cram file"
        bam: "merged-lanes bam file; supply this with bamIndex, or cram with cramIndex"
        bamIndex: "index (.bai) for the bam file"
        markDups: "should the alignment be duplicate marked?, generally yes"
        filterBam: "should filterBam prefilter the input to the fingerprint intervals before splitting? Generally true"
        outputFileNamePrefix: "Optional output prefix for the output"
        reference: "the reference genome for input sample"
        maxReads: "recorded in the metrics json only; no downsampling is done for aligned input"
        sampleId: "value that will be used as the sample identifier in the vcf fingerprint"
   }

Map[String,GenomeResources] resources = {
  "hg38": {
    "refFasta" : "$HG38_ROOT/hg38_random.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_intervals.bed",
    "intervalsToParallelizeByString" : "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM",
    "filterBamModules" : "crosscheckfingerprints-haplotype-map/20230324 samtools/1.14 hg38/p12",
    "splitLanesModules" : "samtools/1.15 hg38/p12",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.15 hg38/p12",
    "mergeBamsModules" : "gatk/4.2.0.0",
    "alignmentMetricsModules" : "samtools/1.15 hg38/p12",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 hg38/p12 crosscheckfingerprints-haplotype-map/20230324"
  },
  "hg19": {
    "refFasta" : "$HG19_ROOT/hg19_random.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_intervals.bed",
    "intervalsToParallelizeByString" : "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM",
    "filterBamModules" : "crosscheckfingerprints-haplotype-map/20230324 samtools/1.14 hg19/p13",
    "splitLanesModules" : "samtools/1.15 hg19/p13",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.15 hg19/p13",
    "mergeBamsModules" : "gatk/4.2.0.0",
    "alignmentMetricsModules" : "samtools/1.15 hg19/p13",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 hg19/p13 crosscheckfingerprints-haplotype-map/20230324"
  }}

   # -------------------------------------------------------
   # Stage 0: Resolve the alignment input.
   #
   # Exactly one of (cram + cramIndex) or (bam + bamIndex) is expected; cram
   # wins if both are given. Each index is resolved inside the same conditional
   # as its alignment file, so a cram can never end up paired with a .bai. The
   # single-element select_first fails loudly when the matching index is missing.
   # -------------------------------------------------------
   if (defined(cram)) {
     File cramInput      = select_first([cram])
     File cramInputIndex = select_first([cramIndex])
   }
   if (!defined(cram)) {
     File bamInput      = select_first([bam])
     File bamInputIndex = select_first([bamIndex])
   }
   File alignFile  = select_first([cramInput,      bamInput])
   File alignIndex = select_first([cramInputIndex, bamInputIndex])

   # -------------------------------------------------------
   # Stage 1: Reduce the merged-lanes input to per-lane files.
   #
   # Filtering to the fingerprint intervals happens before the split so
   # splitLanes works on a much smaller file. splitLanes always emits bam,
   # whichever format came in.
   # -------------------------------------------------------
   if (filterBam) {
     call filterBam as filterBamPreSplit {
       input:
         inputBam = alignFile,
         inputBai = alignIndex,
         intervalBed = resources[reference].intervalBed,
         refFasta = resources[reference].refFasta,
         outputFileNamePrefix = outputFileNamePrefix,
         modules = resources[reference].filterBamModules
     }
   }

   call splitLanes {
     input:
       inputBam = select_first([filterBamPreSplit.bam,      alignFile]),
       inputBai = select_first([filterBamPreSplit.bamIndex, alignIndex]),
       refFasta = resources[reference].refFasta,
       outputFileNamePrefix = outputFileNamePrefix,
       modules = resources[reference].splitLanesModules
   }

   # One element per read group found in the input.
   Array[File] bamsToProcess = splitLanes.laneBams
   Array[File] baisToProcess = splitLanes.laneBamIndexes

   # -------------------------------------------------------
   # Stage 2: Per-lane processing (scattered in parallel)
   #   filterBam -> markDuplicates (chr-scatter) -> merge chr bams
   # -------------------------------------------------------
   call splitStringToArray {
     input:
       str = resources[reference].intervalsToParallelizeByString
   }
   Array[Array[String]] intervalsToParallelizeBy = splitStringToArray.out

   scatter (idx in range(length(bamsToProcess))) {
     # Derive a clean per-lane prefix: strip the .bam extension left by splitLanes, then
     # replace any remaining dots with underscores. Output-provisioning derives file
     # identity from the name and mishandles base names containing multiple dots, which
     # stalls provision-out for that lane.
     String laneBase   = sub(basename(bamsToProcess[idx]), "\\.bam$", "")
     String lanePrefix = sub(laneBase, "\\.", "_")

     # Redundant when filterBamPreSplit already ran (it always does when filterBam is
     # true), but kept so the per-lane knobs stay available and the filtered lane file
     # is what feeds duplicate marking.
     if (filterBam) {
       call filterBam as filterBamLane {
         input:
           inputBam = bamsToProcess[idx],
           inputBai = baisToProcess[idx],
           intervalBed = resources[reference].intervalBed,
           refFasta = resources[reference].refFasta,
           outputFileNamePrefix = lanePrefix,
           modules = resources[reference].filterBamModules
       }
     }

     if (markDups) {
       scatter (intervals in intervalsToParallelizeBy) {
         call markDuplicates {
           input:
             inputBam = select_first([filterBamLane.bam, bamsToProcess[idx]]),
             inputBai = select_first([filterBamLane.bamIndex, baisToProcess[idx]]),
             outputFileNamePrefix = lanePrefix,
             intervals = intervals,
             refFasta = resources[reference].refFasta,
             modules = resources[reference].markDuplicatesModules
         }
       }
       call mergeBams as mergeIntervalBams {
         input:
           bams = markDuplicates.bam,
           outputFileName = lanePrefix,
           suffix = "",
           modules = resources[reference].mergeBamsModules
       }
     }

     # Best available output for this lane, in priority order:
     #   markDups merged > filterBam filtered > original lane file
     File laneResult      = select_first([mergeIntervalBams.mergedBam,      filterBamLane.bam,      bamsToProcess[idx]])
     File laneResultIndex = select_first([mergeIntervalBams.mergedBamIndex, filterBamLane.bamIndex, baisToProcess[idx]])

     # -------------------------------------------------------
     # Stage 3: Metrics and fingerprint - one per lane
     # -------------------------------------------------------
     call alignmentMetrics {
       input:
          inputBam = laneResult,
          inputBai = laneResultIndex,
          outputFileNamePrefix = lanePrefix,
          markDups = markDups,
          maxReads = maxReads,
          refFasta = resources[reference].refFasta,
          modules = resources[reference].alignmentMetricsModules
     }

     call extractFingerprint {
       input:
          inputBam = laneResult,
          inputBai = laneResultIndex,
          haplotypeMap = resources[reference].refHapMap,
          refFasta = resources[reference].refFasta,
          outputFileNamePrefix = lanePrefix,
          sampleId = sampleId,
          modules = resources[reference].extractFingerprintModules
      }

     # Extract the read group ID (RGID) from this lane's file; used as limsId.
     call fingerprintReadgroupInfo {
       input:
         inputBam = laneResult,
         refFasta = resources[reference].refFasta,
         modules  = resources[reference].alignmentMetricsModules
     }

     OutputGroup laneOutput = {
       "limsId":    fingerprintReadgroupInfo.readgroupId,
       "outputVcf": extractFingerprint.vgz,
       "outputTbi": extractFingerprint.tbi,
       "json":      alignmentMetrics.json,
       "samstats":  alignmentMetrics.samstats
     }
   }

   output {
      Array[OutputGroup] outputFingerprints = laneOutput
   }

    meta {
     author: "Lawrence Heisler, Gavin Peng"
     email: "lawrence.heisler@oicr.on.ca, gpeng@oicr.on.ca"
     description: "Multi-lane crosscheckFingerprintsCollector for aligned input. Takes a merged-lanes cram or bam, splits it by read group, and generates a genotype fingerprint per lane using gatk ExtractFingerprint. Outputs are vcf files that can be processed through gatk CrosscheckFingerprints\n##"
     dependencies: [
      {
        name: "gatk/4.2.0.0",
        url: "https://gatk.broadinstitute.org"
      },
      {
        name: "tabix/0.2.6",
        url: "http://www.htslib.org"
      },
      {
        name: "samtools/1.14",
        url: "http://www.htslib.org/"
      },
      {
        name: "samtools/1.15",
        url: "http://www.htslib.org/"
      },
      { name: "gsi crosscheckfingerprints-haplotype-map module : crosscheckfingerprints-haplotype-map/20230324",
        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
      },
      {
        name: "gsi hg38 modules : hg38/p12",
        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
      },
      {
        name: "gsi hg19 modules : hg19/p13",
        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
      }
     ]
     output_meta: {
     outputFingerprints: {
         description: "per-lane output groups; each carries the lane read group ID (limsId), the crosscheck fingerprint vcf.gz and its .tbi index, the alignment metrics json, and the samstats summary"
     }
     }
  }
}


# ==========================================
#  Split a merged cram/bam into per-lane bam
#  files using samtools split (by read group)
# ==========================================

task splitLanes {
  input {
    File inputBam
    File inputBai
    String refFasta
    String outputFileNamePrefix
    String modules
    Int jobMemory = 16
    Int timeout = 24
  }
  parameter_meta {
    inputBam: "input .cram or .bam file to split by read group"
    inputBai: "index for the input file (.crai or .bai)"
    refFasta: "path to reference FASTA (required for CRAM input)"
    outputFileNamePrefix: "prefix for output lane bam files"
    modules: "Names and versions of modules"
    jobMemory: "memory allocated for Job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -euo pipefail
    EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
    if [ "$EXT" = "cram" ]; then
      # samtools split does not support -T; convert CRAM to BAM first
      ln -s ~{inputBam} input.cram
      ln -s ~{inputBai} input.cram.crai
      samtools view -b -T ~{refFasta} -o input_converted.bam input.cram
      samtools index input_converted.bam
      samtools split -f "~{outputFileNamePrefix}_%!.bam" input_converted.bam
    else
      ln -s ~{inputBam} input.bam
      ln -s ~{inputBai} input.bam.bai
      samtools split -f "~{outputFileNamePrefix}_%!.bam" input.bam
    fi

    # samtools split silently writes nothing when the input carries no @RG
    # headers. The rest of the workflow scatters over these files, so an empty
    # split has to fail here rather than yield zero fingerprints.
    shopt -s nullglob
    lanes=(~{outputFileNamePrefix}_*.bam)
    if [[ ${#lanes[@]} -eq 0 ]]; then
      echo "ERROR: samtools split produced no per-lane bam. Does the input have @RG headers?" >&2
      exit 1
    fi
    echo "split into ${#lanes[@]} lane(s): ${lanes[*]}" >&2

    for f in "${lanes[@]}"; do samtools index "$f"; done
  >>>

  output {
    Array[File] laneBams       = glob("~{outputFileNamePrefix}_*.bam")
    Array[File] laneBamIndexes = glob("~{outputFileNamePrefix}_*.bam.bai")
  }

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
}


# ==========================================
#  Filter Cram/Bam to Intervals
# ==========================================

task filterBam {
 input{
  File inputBam
  File inputBai
  String intervalBed
  String refFasta
  String modules
  String outputFileNamePrefix
  Int jobMemory = 16
  Int overhead = 6
  Int timeout = 24
 }
 parameter_meta {
  inputBam: "input .cram or .bam file"
  inputBai: "index of the input file"
  intervalBed: "bed file of the fingerprint intervals to keep"
  outputFileNamePrefix: "prefix for making names for output files"
  refFasta: "path to reference FASTA (required for CRAM input, harmless for BAM)"
  jobMemory: "memory allocated for Job"
  overhead: "memory allocated to overhead of the job other than used in the filter command"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }

command <<<
  set -euo pipefail
  EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
  ln -s ~{inputBam} input.$EXT
  if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
  else                          ln -s ~{inputBai} input.bam.bai
  fi
  samtools view -b -T ~{refFasta} -L ~{intervalBed} input.$EXT > ~{outputFileNamePrefix}.filtered.bam
  samtools index ~{outputFileNamePrefix}.filtered.bam
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File bam      = "~{outputFileNamePrefix}.filtered.bam"
  File bamIndex = "~{outputFileNamePrefix}.filtered.bam.bai"
 }
}


# ==========================================
#  Split a string to array
# ==========================================
task splitStringToArray {
  input {
    String str
    String lineSeparator = ","
    String recordSeparator = "+"

    Int jobMemory = 1
    Int threads = 1
    Int timeout = 1
    String modules = ""
  }

  command <<<
    set -euo pipefail

    echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
  >>>

  output {
    Array[Array[String]] out = read_tsv(stdout())
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{threads}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    str: "Interval string to split (e.g. chr1,chr2,chr3+chr4)."
    lineSeparator: "Interval group separator - these are the intervals to split by."
    recordSeparator: "Interval interval group separator - this can be used to combine multiple intervals into one group."
    jobMemory: "Memory allocated to job (in GB)."
    threads: "The number of threads to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}


# ==========================================
#  Duplicate Marking
# ==========================================

task markDuplicates {
 input{
  File inputBam
  File inputBai
  String refFasta
  String modules
  String outputFileNamePrefix
  Array[String] intervals
  Int jobMemory = 16
  Int overhead = 6
  Int timeout = 24
 }
 parameter_meta {
  inputBam: "input .cram or .bam file"
  inputBai: "index of the input file"
  refFasta: "path to reference FASTA (required for CRAM input, harmless for BAM)"
  outputFileNamePrefix: "prefix for making names for output files"
  intervals: "intervals to restrict this shard of duplicate marking to"
  jobMemory: "memory allocated for Job"
  overhead: "memory allocated to overhead of the job other than used in markDuplicates command"
  modules: "Names and versions of modules"
  timeout: "Timeout in hours, needed to override imposed limits"
 }

command <<<
  set -euo pipefail
  EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
  ln -s ~{inputBam} input.$EXT
  if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
  else                          ln -s ~{inputBai} input.bam.bai
  fi
  samtools view -b -T ~{refFasta} input.$EXT \
        ~{sep=" " intervals} > intervalBam.bam
  samtools index intervalBam.bam intervalBam.bam.bai

  $GATK_ROOT/bin/gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
                      -I intervalBam.bam \
                      --METRICS_FILE ~{outputFileNamePrefix}.dupmetrics \
                      --VALIDATION_STRINGENCY SILENT \
                      --CREATE_INDEX true \
                      -O ~{outputFileNamePrefix}.dupmarked.bam
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File bam = "~{outputFileNamePrefix}.dupmarked.bam"
 }
}


# ==========================================
#  Merge Bam files
# ==========================================

task mergeBams {
  input {
    Array[File] bams
    String outputFileName
    String suffix = ".merge"
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int threads = 1
    Int timeout = 6
    String modules = "gatk/4.2.0.0"
  }

  command <<<
    set -euo pipefail

    $GATK_ROOT/bin/gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
    --INPUT ~{sep=" --INPUT " bams} \
    --OUTPUT "~{outputFileName}~{suffix}.bam" \
    --CREATE_INDEX true \
    --SORT_ORDER coordinate \
    --ASSUME_SORTED false \
    --USE_THREADING true \
    --VALIDATION_STRINGENCY SILENT \
    ~{additionalParams}
  >>>

  output {
    File mergedBam      = "~{outputFileName}~{suffix}.bam"
    File mergedBamIndex = "~{outputFileName}~{suffix}.bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{threads}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    bams: "Array of bam files to merge together."
    outputFileName: "Output files will be prefixed with this."
    suffix: "Suffix appended to outputFileName before the .bam extension."
    additionalParams: "Additional parameters to pass to GATK MergeSamFiles."
    jobMemory: "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    threads: "The number of threads to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}


# ==========================================
#  coverage metrics from the file used for fingerprint analysis
# ==========================================

 task alignmentMetrics{
   input{
    File inputBam
    File inputBai
    String refFasta
    String modules
    String outputFileNamePrefix
    Boolean markDups
    Int jobMemory = 8
    Int timeout = 24
    Int maxReads
   }
   parameter_meta {
    inputBam: "input .cram or .bam file"
    inputBai: "index of the input file"
    refFasta: "path to reference FASTA (required for CRAM input, harmless for BAM)"
    outputFileNamePrefix: "prefix for making names for output files"
    markDups: "whether duplicate marking was run; recorded in the json"
    jobMemory: "memory allocated for Job"
    modules: "Names and versions of modules"
    maxReads: "the maximum number of reads used; recorded in the json"
    timeout: "Timeout in hours, needed to override imposed limits"
   }

command <<<
  set -euo pipefail

  EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
  ln -s ~{inputBam} input.$EXT
  if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
  else                          ln -s ~{inputBai} input.bam.bai
  fi

  ### samtools stats
  $SAMTOOLS_ROOT/bin/samtools stats --reference ~{refFasta} input.$EXT > ~{outputFileNamePrefix}.samstats.txt
  reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "raw total sequences:" | cut -f3`
  mapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads mapped:" | cut -f 3`
  unmapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads unmapped:" | cut -f 3`
  mapped_bases=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "bases mapped:" | cut -f 3`
  reads_duplicated=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads duplicated:" | cut -f 3`

  ### samtools coverage, with duplicates
  $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL --reference ~{refFasta} input.$EXT > ~{outputFileNamePrefix}.coverage.txt
  mean_cvg=`cat ~{outputFileNamePrefix}.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }'`

  ### samtools coverage, deduplicated
  $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL,DUP --reference ~{refFasta} input.$EXT > ~{outputFileNamePrefix}.dedup.coverage.txt
  mean_dedup_cvg=`cat ~{outputFileNamePrefix}.dedup.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }'`

  ### json file
  echo \{\"reads\":$reads,\"mapped_reads\":$mapped_reads,\"unmapped_reads\":$unmapped_reads,\"mapped_bases\":$mapped_bases,\"reads_duplicated\":$reads_duplicated,\"mean_raw_cvg\":$mean_cvg,\"mean_dedup_cvg\":$mean_dedup_cvg\,\"markDups\":~{markDups},\"maxReads\":~{maxReads}} > ~{outputFileNamePrefix}.json
>>>

  runtime {
   memory:  "~{jobMemory} GB"
   modules: "~{modules}"
   timeout: "~{timeout}"
  }

  output {
    File json     = "~{outputFileNamePrefix}.json"
    File samstats = "~{outputFileNamePrefix}.samstats.txt"
  }
}


# ==========================================
#  configure and run extractFingerprintsCollector
# ==========================================

task extractFingerprint {
input {
 File inputBam
 File inputBai
 String modules
 String refFasta
 String outputFileNamePrefix
 String haplotypeMap
 String sampleId
 Int jobMemory = 8
 Int timeout = 24
}
parameter_meta {
 inputBam: "input .cram or .bam file"
 inputBai: "index of the input file"
 refFasta: "Path to reference FASTA file"
 outputFileNamePrefix: "prefix for making names for output files"
 haplotypeMap: "Hotspot SNPs are the locations of variants used for genotyping"
 sampleId : "value used as the sample identifier in the vcf fingerprint"
 jobMemory: "memory allocated for Job"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
  set -euo pipefail

 # GATK needs the index beside the alignment file; Cromwell may localize
 # them to different directories, so link both into the task dir.
 EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
 ln -s ~{inputBam} input.$EXT
 if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
 else                          ln -s ~{inputBai} input.bam.bai
 fi

 $GATK_ROOT/bin/gatk ExtractFingerprint \
                    -R ~{refFasta} \
                    -H ~{haplotypeMap} \
                    -I input.$EXT \
                    -O ~{outputFileNamePrefix}.vcf \
                    --SAMPLE_ALIAS ~{sampleId}

 $TABIX_ROOT/bin/bgzip -c ~{outputFileNamePrefix}.vcf > ~{outputFileNamePrefix}.vcf.gz
 $TABIX_ROOT/bin/tabix -p vcf ~{outputFileNamePrefix}.vcf.gz
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File vcf = "~{outputFileNamePrefix}.vcf"
  File vgz = "~{outputFileNamePrefix}.vcf.gz"
  File tbi = "~{outputFileNamePrefix}.vcf.gz.tbi"
 }
}


# ==========================================
#  Extract the read group ID (RGID) from the
#  @RG header line of a lane's cram/bam file
# ==========================================

task fingerprintReadgroupInfo {
  input {
    File inputBam
    String refFasta
    String modules
    Int jobMemory = 8
    Int timeout = 24
  }
  parameter_meta {
    inputBam:  "lane-level cram/bam file used for fingerprint extraction"
    refFasta:  "path to reference FASTA (required for CRAM decoding)"
    modules:   "Names and versions of modules"
    jobMemory: "memory allocated for job"
    timeout:   "timeout in hours"
  }

  command <<<
    set -euo pipefail
    samtools view -H -T ~{refFasta} ~{inputBam} \
      | awk -F'\t' '!found && /^@RG/ { for (i=1; i<=NF; i++) if ($i ~ /^ID:/) { sub(/^ID:/, "", $i); print $i; found=1 } }'
  >>>

  output {
    String readgroupId = read_string(stdout())
  }

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
}
