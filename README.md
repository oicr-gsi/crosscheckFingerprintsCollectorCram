# crosscheckFingerprintsCollectorMultiLane

Multi-lane crosscheckFingerprintsCollector for aligned input. Takes a merged-lanes cram or bam, reads the lane set from its @RG headers, and for each lane streams the reads of that read group that overlap the fingerprint intervals straight into gatk ExtractFingerprint - no cram-to-bam conversion and no per-lane bam on disk. Outputs are vcf files that can be processed through gatk CrosscheckFingerprints
##

## Overview

## Dependencies

* [gatk 4.2.0.0](https://gatk.broadinstitute.org)
* [tabix 0.2.6](http://www.htslib.org)
* [samtools 1.15](http://www.htslib.org/)
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
`readGroupIds.jobMemory`|Int|4|memory allocated for Job
`readGroupIds.timeout`|Int|1|Timeout in hours, needed to override imposed limits
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
```
```
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
```
```
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
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
