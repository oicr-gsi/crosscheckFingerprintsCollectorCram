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
#  The input is always assumed to hold more than one read group.
#  Rather than splitting it into per-lane files on disk, the read
#  groups are read out of the header and each lane is then
#  streamed from the merged file - one read group, restricted to
#  the fingerprint intervals - directly into the tool that
#  consumes it. No cram-to-bam conversion and no per-lane bam is
#  written; the one exception is duplicate marking, which reads
#  its input twice and so cannot take a stream.
#
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
    String readGroupIdsModules
    String markDuplicatesModules
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
        filterBam: "should the per-lane read streams be restricted to the fingerprint intervals? Generally true"
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
    "readGroupIdsModules" : "samtools/1.15 hg38/p12",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.15 hg38/p12 crosscheckfingerprints-haplotype-map/20230324",
    "alignmentMetricsModules" : "samtools/1.15 hg38/p12 crosscheckfingerprints-haplotype-map/20230324",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 samtools/1.15 hg38/p12 crosscheckfingerprints-haplotype-map/20230324"
  },
  "hg19": {
    "refFasta" : "$HG19_ROOT/hg19_random.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_intervals.bed",
    "readGroupIdsModules" : "samtools/1.15 hg19/p13",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.15 hg19/p13 crosscheckfingerprints-haplotype-map/20230324",
    "alignmentMetricsModules" : "samtools/1.15 hg19/p13 crosscheckfingerprints-haplotype-map/20230324",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 samtools/1.15 hg19/p13 crosscheckfingerprints-haplotype-map/20230324"
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
   # Stage 1: Enumerate the lanes.
   #
   # Only the header is read - the lane set is the set of @RG IDs. Nothing is
   # split out here: every task below re-reads the merged input through
   # samtools and keeps one read group with -r.
   # -------------------------------------------------------
   call readGroupIds {
     input:
       inputBam = alignFile,
       refFasta = resources[reference].refFasta,
       modules  = resources[reference].readGroupIdsModules
   }

   # -------------------------------------------------------
   # Stage 2: Per-lane processing (scattered in parallel)
   # -------------------------------------------------------
   scatter (readGroup in readGroupIds.ids) {
     # Read group IDs routinely carry dots and other punctuation. Output
     # provisioning derives file identity from the name and mishandles base
     # names containing multiple dots, which stalls provision-out for that lane,
     # so the name-forming copy of the ID is sanitized. limsId below keeps the
     # ID verbatim.
     String lanePrefix = outputFileNamePrefix + "_" + sub(readGroup, "[^A-Za-z0-9_-]", "_")

     # The only step that cannot consume a stream: MarkDuplicates passes over
     # its input twice. The lane bam it needs is written inside the task and
     # never leaves it, and holds just this read group over the intervals.
     if (markDups) {
       call markDuplicates {
         input:
           inputBam = alignFile,
           inputBai = alignIndex,
           readGroup = readGroup,
           filterToIntervals = filterBam,
           intervalBed = resources[reference].intervalBed,
           refFasta = resources[reference].refFasta,
           outputFileNamePrefix = lanePrefix,
           modules = resources[reference].markDuplicatesModules
       }
     }

     # What this lane's reads are streamed out of. When duplicates were marked
     # that is the dup-marked lane bam - already one read group over the
     # intervals, so the -r / -L filtering below simply passes it through.
     # Otherwise the stream comes straight off the merged input.
     File laneSource      = select_first([markDuplicates.bam,      alignFile])
     File laneSourceIndex = select_first([markDuplicates.bamIndex, alignIndex])

     # -------------------------------------------------------
     # Stage 3: Metrics and fingerprint - one per lane
     # -------------------------------------------------------
     call alignmentMetrics {
       input:
          inputBam = laneSource,
          inputBai = laneSourceIndex,
          readGroup = readGroup,
          filterToIntervals = filterBam,
          intervalBed = resources[reference].intervalBed,
          outputFileNamePrefix = lanePrefix,
          markDups = markDups,
          maxReads = maxReads,
          refFasta = resources[reference].refFasta,
          modules = resources[reference].alignmentMetricsModules
     }

     call extractFingerprint {
       input:
          inputBam = laneSource,
          inputBai = laneSourceIndex,
          readGroup = readGroup,
          filterToIntervals = filterBam,
          intervalBed = resources[reference].intervalBed,
          haplotypeMap = resources[reference].refHapMap,
          refFasta = resources[reference].refFasta,
          outputFileNamePrefix = lanePrefix,
          sampleId = sampleId,
          modules = resources[reference].extractFingerprintModules
      }

     OutputGroup laneOutput = {
       "limsId":    readGroup,
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
     description: "Multi-lane crosscheckFingerprintsCollector for aligned input. Takes a merged-lanes cram or bam, reads the lane set from its @RG headers, and for each lane streams the reads of that read group that overlap the fingerprint intervals straight into gatk ExtractFingerprint - no cram-to-bam conversion and no per-lane bam on disk. Outputs are vcf files that can be processed through gatk CrosscheckFingerprints\n##"
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
#  List the read group IDs of the merged
#  input; header only, no data is decoded
# ==========================================

task readGroupIds {
  input {
    File inputBam
    String refFasta
    String modules
    Int jobMemory = 4
    Int timeout = 1
  }
  parameter_meta {
    inputBam: "merged-lanes .cram or .bam file; only the header is read"
    refFasta: "path to reference FASTA (required for CRAM)"
    modules: "Names and versions of modules"
    jobMemory: "memory allocated for Job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -euo pipefail

    samtools view -H -T "~{refFasta}" "~{inputBam}" \
      | awk -F'\t' '/^@RG/ { for (i = 1; i <= NF; i++) if ($i ~ /^ID:/) { sub(/^ID:/, "", $i); print $i } }' \
      | sort -u > readGroups.txt

    # Every lane downstream is keyed by read group, so a header without @RG
    # lines would silently yield zero fingerprints. Fail here instead.
    if [[ ! -s readGroups.txt ]]; then
      echo "ERROR: no @RG ID found in the input header. Does the input have read groups?" >&2
      exit 1
    fi
    echo "found $(wc -l < readGroups.txt) read group(s)" >&2
  >>>

  output {
    Array[String] ids = read_lines("readGroups.txt")
  }

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }
}


# ==========================================
#  Duplicate Marking, one lane at a time
# ==========================================

task markDuplicates {
 input{
  File inputBam
  File inputBai
  String readGroup
  Boolean filterToIntervals
  String intervalBed
  String refFasta
  String modules
  String outputFileNamePrefix
  Int jobMemory = 16
  Int overhead = 6
  Int timeout = 24
 }
 parameter_meta {
  inputBam: "merged-lanes .cram or .bam file"
  inputBai: "index of the input file"
  readGroup: "@RG ID of the lane to mark"
  filterToIntervals: "restrict the lane to the fingerprint intervals before marking"
  intervalBed: "bed file of the fingerprint intervals"
  refFasta: "path to reference FASTA (required for CRAM input, harmless for BAM)"
  outputFileNamePrefix: "prefix for making names for output files"
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

  # -M -L seeks to the intervals through the index; a bare -L would decode the
  # whole file just to filter it.
  regions=()
  if [[ "~{filterToIntervals}" = "true" ]]; then
    regions+=(-M -L "~{intervalBed}")
  fi

  # MarkDuplicates reads its input twice, so this is the one step that needs a
  # file rather than a stream. The lane bam is task-local and never provisioned.
  samtools view -b -T "~{refFasta}" -r "~{readGroup}" ${regions[@]+"${regions[@]}"} input.$EXT > lane.bam

  $GATK_ROOT/bin/gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
                      -I lane.bam \
                      --METRICS_FILE ~{outputFileNamePrefix}.dupmetrics \
                      --VALIDATION_STRINGENCY SILENT \
                      --CREATE_INDEX true \
                      -O ~{outputFileNamePrefix}.dupmarked.bam

  # --CREATE_INDEX writes <prefix>.bai; the streams downstream look for the
  # index beside the bam it belongs to.
  mv ~{outputFileNamePrefix}.dupmarked.bai ~{outputFileNamePrefix}.dupmarked.bam.bai
>>>

 runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
 }

 output {
  File bam      = "~{outputFileNamePrefix}.dupmarked.bam"
  File bamIndex = "~{outputFileNamePrefix}.dupmarked.bam.bai"
 }
}


# ==========================================
#  coverage metrics for the reads used for
#  fingerprint analysis
# ==========================================

 task alignmentMetrics{
   input{
    File inputBam
    File inputBai
    String readGroup
    Boolean filterToIntervals
    String intervalBed
    String refFasta
    String modules
    String outputFileNamePrefix
    Boolean markDups
    Int jobMemory = 8
    Int timeout = 24
    Int maxReads
   }
   parameter_meta {
    inputBam: "merged-lanes .cram or .bam file, or the dup-marked lane bam"
    inputBai: "index of the input file"
    readGroup: "@RG ID of the lane to report on"
    filterToIntervals: "restrict the lane stream to the fingerprint intervals"
    intervalBed: "bed file of the fingerprint intervals"
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

  regions=()
  if [[ "~{filterToIntervals}" = "true" ]]; then
    regions+=(-M -L "~{intervalBed}")
  fi

  # Each metric re-streams the lane instead of sharing a temporary bam. With
  # -M -L that is an index seek over the fingerprint intervals, not a full pass,
  # and nothing is written to disk.
  stream_lane () {
    $SAMTOOLS_ROOT/bin/samtools view -u -T "~{refFasta}" -r "~{readGroup}" ${regions[@]+"${regions[@]}"} input.$EXT
  }

  ### samtools stats
  stream_lane | $SAMTOOLS_ROOT/bin/samtools stats --reference ~{refFasta} - > ~{outputFileNamePrefix}.samstats.txt
  reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "raw total sequences:" | cut -f3`
  mapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads mapped:" | cut -f 3`
  unmapped_reads=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads unmapped:" | cut -f 3`
  mapped_bases=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "bases mapped:" | cut -f 3`
  reads_duplicated=`cat ~{outputFileNamePrefix}.samstats.txt | grep ^SN | grep "reads duplicated:" | cut -f 3`

  ### samtools coverage, with duplicates
  stream_lane | $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL --reference ~{refFasta} - > ~{outputFileNamePrefix}.coverage.txt
  mean_cvg=`cat ~{outputFileNamePrefix}.coverage.txt | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }'`

  ### samtools coverage, deduplicated
  stream_lane | $SAMTOOLS_ROOT/bin/samtools coverage --ff UNMAP,SECONDARY,QCFAIL,DUP --reference ~{refFasta} - > ~{outputFileNamePrefix}.dedup.coverage.txt
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
 String readGroup
 Boolean filterToIntervals
 String intervalBed
 String modules
 String refFasta
 String outputFileNamePrefix
 String haplotypeMap
 String sampleId
 Int jobMemory = 8
 Int timeout = 24
}
parameter_meta {
 inputBam: "merged-lanes .cram or .bam file, or the dup-marked lane bam"
 inputBai: "index of the input file"
 readGroup: "@RG ID of the lane to fingerprint"
 filterToIntervals: "restrict the lane stream to the fingerprint intervals"
 intervalBed: "bed file of the fingerprint intervals"
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

 # samtools needs the index beside the alignment file; Cromwell may localize
 # them to different directories, so link both into the task dir.
 EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
 ln -s ~{inputBam} input.$EXT
 if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
 else                          ln -s ~{inputBai} input.bam.bai
 fi

 regions=()
 if [[ "~{filterToIntervals}" = "true" ]]; then
   regions+=(-M -L "~{intervalBed}")
 fi

 # The lane's reads go straight into the fingerprinter: no per-lane bam is
 # written, and the input is decoded once, over the fingerprint intervals only.
 # Streaming uncompressed bam rather than sam keeps the bgzf EOF block, so a
 # stream cut short is an error instead of a silently under-covered fingerprint.
 samtools view -u -T "~{refFasta}" -r "~{readGroup}" ${regions[@]+"${regions[@]}"} input.$EXT \
   | $GATK_ROOT/bin/gatk ExtractFingerprint \
                    -R ~{refFasta} \
                    -H ~{haplotypeMap} \
                    -I /dev/stdin \
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
