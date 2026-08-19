# crosscheckFingerprintsCollectorMultiLane

## Overview

Multi-lane crosscheckFingerprintsCollector for aligned input. Takes a merged-lanes cram or bam, reads the lane set from its @RG headers, and for each lane reads back only that read group over the fingerprint intervals, seeking through the index - no cram-to-bam conversion of the whole input and no per-lane bam provisioned. Outputs are vcf files that can be processed through gatk CrosscheckFingerprints
##

## Dependencies

* [gatk 4.2.0.0](https://gatk.broadinstitute.org)
* [tabix 0.2.6](http://www.htslib.org)
* [samtools 1.16.1](http://www.htslib.org/)
* [gsi crosscheckfingerprints-haplotype-map module : crosscheckfingerprints-haplotype-map 20230324](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg38 modules : hg38 p12](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [gsi hg19 modules : hg19 p13](https://gitlab.oicr.on.ca/ResearchIT/modulator)


## Usage

### Cromwell
```
java -jar cromwell.jar run crosscheckFingerprintsCollectorMultiLane.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`markDups`|Boolean|should the alignment be duplicate marked?, generally yes
`filterBam`|Boolean|should the per-lane read streams be restricted to the fingerprint intervals? Generally true
`outputFileNamePrefix`|String|Optional output prefix for the output
`reference`|String|the reference genome for input sample
`sampleId`|String|value that will be used as the sample identifier in the vcf fingerprint


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`cram`|File?|None|merged-lanes cram file; supply this with cramIndex, or bam with bamIndex
`cramIndex`|File?|None|index (.crai) for the cram file
`bam`|File?|None|merged-lanes bam file; supply this with bamIndex, or cram with cramIndex
`bamIndex`|File?|None|index (.bai) for the bam file
`maxReads`|Int|0|recorded in the metrics json only; no downsampling is done for aligned input


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`prepareLanes.jobMemory`|Int|4|memory allocated for Job
`prepareLanes.timeout`|Int|1|Timeout in hours, needed to override imposed limits
`markDuplicates.jobMemory`|Int|16|memory allocated for Job
`markDuplicates.overhead`|Int|6|memory allocated to overhead of the job other than used in markDuplicates command
`markDuplicates.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`alignmentMetrics.jobMemory`|Int|8|memory allocated for Job
`alignmentMetrics.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`extractFingerprint.jobMemory`|Int|8|memory allocated for Job
`extractFingerprint.timeout`|Int|24|Timeout in hours, needed to override imposed limits


### Outputs

Output | Type | Description | Labels
---|---|---|---
`outputFingerprints`|Array[OutputGroup]|per-lane output groups; each carries the lane read group ID (limsId), the crosscheck fingerprint vcf.gz and its .tbi index, the alignment metrics json, and the samstats summary|


## Commands
This section lists command(s) run by crosscheckFingerprintsCollectorMultiLane workflow

* Running crosscheckFingerprintsCollectorMultiLane

```
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
```
```
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
```
```
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
```
```
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
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
