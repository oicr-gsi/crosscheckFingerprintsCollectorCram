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
#  groups are read out of the header and each lane is then read
#  back out of the merged file - one read group, seeking to the
#  fingerprint intervals through the index. There is no cram-to-bam
#  conversion of the whole input and no per-lane bam is provisioned.
#  The samtools metrics consume that as a stream; the two Picard
#  tools cannot (both read their input twice), so they materialize
#  the lane inside their own task, at fingerprint-interval size.
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
    String prepareLanesModules
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
  "hg38_noAlt": {
    "refFasta" : "$HG38_NOALT_ROOT/hg38_noAlt.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_intervals.bed",
    "prepareLanesModules" : "samtools/1.16.1 hg38-noalt/p12 crosscheckfingerprints-haplotype-map/20230324",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.16.1 hg38-noalt/p12",
    "alignmentMetricsModules" : "samtools/1.16.1 hg38-noalt/p12",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 samtools/1.16.1 hg38-noalt/p12 crosscheckfingerprints-haplotype-map/20230324"
  },
  "hg38": {
    "refFasta" : "$HG38_ROOT/hg38_random.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg38_intervals.bed",
    "prepareLanesModules" : "samtools/1.16.1 hg38/p12 crosscheckfingerprints-haplotype-map/20230324",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.16.1 hg38/p12",
    "alignmentMetricsModules" : "samtools/1.16.1 hg38/p12",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 samtools/1.16.1 hg38/p12 crosscheckfingerprints-haplotype-map/20230324"
  },
  "hg19": {
    "refFasta" : "$HG19_ROOT/hg19_random.fa",
    "refHapMap" : "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_chr.map",
    "intervalBed": "$CROSSCHECKFINGERPRINTS_HAPLOTYPE_MAP_ROOT/oicr_hg19_intervals.bed",
    "prepareLanesModules" : "samtools/1.16.1 hg19/p13 crosscheckfingerprints-haplotype-map/20230324",
    "markDuplicatesModules" : "gatk/4.2.0.0 samtools/1.16.1 hg19/p13",
    "alignmentMetricsModules" : "samtools/1.16.1 hg19/p13",
    "extractFingerprintModules" : "gatk/4.2.0.0 tabix/0.2.6 samtools/1.16.1 hg19/p13 crosscheckfingerprints-haplotype-map/20230324"
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
   # Stage 1: Read the header once.
   #
   # Three things come out of it: the lane set (the @RG IDs), the fingerprint
   # intervals reduced to the contigs this file actually declares, and whether
   # its header lets the lane streams seek through the index. Nothing is split
   # out here - every task below re-reads the merged input through samtools and
   # keeps one read group with -r.
   # -------------------------------------------------------
   call prepareLanes {
     input:
       inputBam = alignFile,
       inputBai = alignIndex,
       intervalBed = resources[reference].intervalBed,
       filterToIntervals = filterBam,
       refFasta = resources[reference].refFasta,
       modules  = resources[reference].prepareLanesModules
   }

   # -------------------------------------------------------
   # Stage 2: Per-lane processing (scattered in parallel)
   # -------------------------------------------------------
   scatter (readGroup in prepareLanes.ids) {
     # Read group IDs routinely carry dots and other punctuation. Output
     # provisioning derives file identity from the name and mishandles base
     # names containing multiple dots, which stalls provision-out for that lane,
     # so the name-forming copy of the ID is sanitized. limsId below keeps the
     # ID verbatim.
     String lanePrefix = outputFileNamePrefix + "_" + sub(readGroup, "[^A-Za-z0-9_-]", "_")

     # MarkDuplicates passes over its input twice, so it cannot take a stream.
     # The lane bam it needs is written inside the task and never leaves it, and
     # holds just this read group over the intervals.
     if (markDups) {
       call markDuplicates {
         input:
           inputBam = alignFile,
           inputBai = alignIndex,
           readGroup = readGroup,
           filterToIntervals = filterBam,
           useIndexSeek = prepareLanes.canSeek,
           intervalBed = prepareLanes.intervals,
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
          useIndexSeek = prepareLanes.canSeek,
          intervalBed = prepareLanes.intervals,
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
          useIndexSeek = prepareLanes.canSeek,
          intervalBed = prepareLanes.intervals,
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
     description: "Multi-lane crosscheckFingerprintsCollector for aligned input. Takes a merged-lanes cram or bam, reads the lane set from its @RG headers, and for each lane reads back only that read group over the fingerprint intervals, seeking through the index - no cram-to-bam conversion of the whole input and no per-lane bam provisioned. Outputs are vcf files that can be processed through gatk CrosscheckFingerprints\n##"
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
        name: "samtools/1.16.1",
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
#  From the merged input's header: the lane
#  set, the fingerprint intervals reduced to
#  the contigs it declares, and whether the
#  lane streams may seek through the index.
# ==========================================

task prepareLanes {
  input {
    File inputBam
    File inputBai
    String intervalBed
    Boolean filterToIntervals
    String refFasta
    String modules
    Int jobMemory = 4
    Int timeout = 1
  }
  parameter_meta {
    inputBam: "merged-lanes .cram or .bam file; only the header and one probe region are read"
    inputBai: "index of the input file"
    intervalBed: "bed file of the fingerprint intervals"
    filterToIntervals: "whether the lane streams will use the intervals; when false the emitted bed is unused"
    refFasta: "path to reference FASTA (required for CRAM)"
    modules: "Names and versions of modules"
    jobMemory: "memory allocated for Job"
    timeout: "Timeout in hours, needed to override imposed limits"
  }

  command <<<
    set -euo pipefail
    EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
    ln -s ~{inputBam} input.$EXT
    if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
    else                          ln -s ~{inputBai} input.bam.bai
    fi

    samtools view -H -T "~{refFasta}" input.$EXT > header.sam

    awk -F'\t' '/^@RG/ { for (i = 1; i <= NF; i++) if ($i ~ /^ID:/) { sub(/^ID:/, "", $i); print $i } }' \
      header.sam | sort -u > readGroups.txt

    # Every lane downstream is keyed by read group, so a header without @RG
    # lines would silently yield zero fingerprints. Fail here instead.
    if [[ ! -s readGroups.txt ]]; then
      echo "ERROR: no @RG ID found in the input header. Does the input have read groups?" >&2
      exit 1
    fi
    echo "found $(wc -l < readGroups.txt) read group(s)" >&2

    # Intervals naming a contig this file does not declare can never match, so
    # drop them and say how many went.
    awk -F'\t' '/^@SQ/ { for (i = 1; i <= NF; i++) if ($i ~ /^SN:/) { sub(/^SN:/, "", $i); print $i } }' \
      header.sam | sort -u > contigs.txt

    awk -F'\t' 'NR == FNR { known[$1]; next } /^(#|track|browser)/ { next } ($1 in known)' \
      contigs.txt "~{intervalBed}" | sort -k1,1 -k2,2n > intervals.bed

    kept=$(wc -l < intervals.bed)
    total=$(grep -c -v -e '^#' -e '^track' -e '^browser' "~{intervalBed}" || true)
    echo "fingerprint intervals: kept $kept of $total (dropped: contig not in this header)" >&2

    if [[ "~{filterToIntervals}" = "true" && "$kept" -eq 0 ]]; then
      echo "ERROR: none of the $total fingerprint intervals in ~{intervalBed} name a contig" >&2
      echo "       declared by this input. Wrong reference build, or chr-prefix mismatch?" >&2
      exit 1
    fi

    # Can the lane streams seek to the intervals through the index (-M -L), or do
    # they have to filter a full sequential pass (-L)?
    #
    # Both give identical reads; -M is the one that makes per-lane streaming
    # affordable, since it touches only the fingerprint intervals instead of
    # decoding the whole file once per stream.
    #
    # It is not available on every samtools. Up to and including 1.15, sam_view.c
    # loads the index only when an index file or a region argument was given:
    #
    #   if ( settings.fn_idx_in || nregs )                          # 1.15
    #   if ( settings.fn_idx_in || nregs || settings.multi_region )  # 1.16 and up
    #
    # so with -M -L and no region argument, 1.15 hands a NULL index to the
    # multi-region iterator and dies with "Iterator could not be created.
    # Aborting.". Plain -L needs no index and is unaffected - which is why the
    # older split-based pipeline ran fine on the same samtools.
    #
    # Probe it once here rather than let every lane task rediscover it.
    echo "true" > canSeek.txt
    if [[ "$kept" -gt 0 ]]; then
      head -1 intervals.bed > probe.bed
      if ! samtools view -c -T "~{refFasta}" -M -L probe.bed input.$EXT > /dev/null 2>probe.err; then
        echo "false" > canSeek.txt
        echo "WARNING: cannot seek to the fingerprint intervals through the index:" >&2
        sed 's/^/  /' probe.err >&2
        echo "         $(samtools --version | head -1)" >&2
        echo "         Falling back to a full sequential pass per lane stream:" >&2
        echo "         every stream decodes the whole input. samtools 1.16 or" >&2
        echo "         newer restores the fast path." >&2

        # Without -M, a read overlapping more than one region of -L is emitted
        # once per region, and duplicate records make MarkDuplicates abort.
        # Stop here rather than let it fail further downstream.
        if [[ "~{filterToIntervals}" = "true" && "$kept" -gt 1 ]]; then
          echo "ERROR: interval filtering needs the -M fast path to guarantee each read" >&2
          echo "       is emitted once. Use samtools 1.16 or newer, or run with" >&2
          echo "       filterBam = false." >&2
          exit 1
        fi
      fi
    fi
    echo "index seeking: $(cat canSeek.txt)" >&2
  >>>

  output {
    Array[String] ids  = read_lines("readGroups.txt")
    File intervals     = "intervals.bed"
    Boolean canSeek    = read_boolean("canSeek.txt")
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
  Boolean useIndexSeek
  File intervalBed
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
  useIndexSeek: "seek to the intervals through the index (-M) instead of filtering a full pass"
  intervalBed: "bed file of the fingerprint intervals, reduced to the contigs of this input by prepareLanes"
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

  # -M seeks to the intervals through the index; without it -L still filters
  # correctly, but only by reading the whole file. prepareLanes probed which one
  # this input's header supports.
  regions=()
  if [[ "~{filterToIntervals}" = "true" ]]; then
    if [[ "~{useIndexSeek}" = "true" ]]; then
      regions+=(-M -L "~{intervalBed}")
    else
      regions+=(-L "~{intervalBed}")
    fi
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
    Boolean useIndexSeek
    File intervalBed
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
    useIndexSeek: "seek to the intervals through the index (-M) instead of filtering a full pass"
    intervalBed: "bed file of the fingerprint intervals, reduced to the contigs of this input by prepareLanes"
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

  # -M seeks to the intervals through the index; without it -L still filters
  # correctly, but only by reading the whole file. prepareLanes probed which one
  # this input's header supports.
  regions=()
  if [[ "~{filterToIntervals}" = "true" ]]; then
    if [[ "~{useIndexSeek}" = "true" ]]; then
      regions+=(-M -L "~{intervalBed}")
    else
      regions+=(-L "~{intervalBed}")
    fi
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
 Boolean useIndexSeek
 File intervalBed
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
 useIndexSeek: "seek to the intervals through the index (-M) instead of filtering a full pass"
 intervalBed: "bed file of the fingerprint intervals, reduced to the contigs of this input by prepareLanes"
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

 # -M seeks to the intervals through the index; without it -L still filters
 # correctly, but only by reading the whole file. prepareLanes probed which one
 # this input's header supports.
 regions=()
 if [[ "~{filterToIntervals}" = "true" ]]; then
   if [[ "~{useIndexSeek}" = "true" ]]; then
     regions+=(-M -L "~{intervalBed}")
   else
     regions+=(-L "~{intervalBed}")
   fi
 fi

 # ExtractFingerprint cannot read a pipe: Picard opens its input more than once,
 # and on a non-seekable stream the second read fails with
 #   RuntimeIOException: Read error; BinaryCodec in readmode; streamed file
 # (-I /dev/stdin does work for Picard tools that make a single pass). So the
 # lane's reads land in a task-local bam first. It holds one read group over the
 # fingerprint intervals, is never provisioned, and the merged input is still
 # decoded just once - over those intervals only.
 samtools view -b -T "~{refFasta}" -r "~{readGroup}" ${regions[@]+"${regions[@]}"} input.$EXT > lane.bam

 # No index on lane.bam: without one the fingerprinter walks it sequentially,
 # which is what we want for a file this small.
 $GATK_ROOT/bin/gatk ExtractFingerprint \
                    -R ~{refFasta} \
                    -H ~{haplotypeMap} \
                    -I lane.bam \
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
