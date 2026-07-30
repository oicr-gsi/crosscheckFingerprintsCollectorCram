# crosscheckFingerprintsCollectorMultiLane

Multi-lane crosscheckFingerprintsCollector for aligned input. Takes a merged-lanes cram or bam, splits it by read group, and generates a genotype fingerprint per lane using gatk ExtractFingerprint. Outputs are vcf files that can be processed through gatk CrosscheckFingerprints
##

## Overview

## Dependencies

* [gatk 4.2.0.0](https://gatk.broadinstitute.org)
* [tabix 0.2.6](http://www.htslib.org)
* [samtools 1.14](http://www.htslib.org/)
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
`filterBam`|Boolean|should filterBam prefilter the input to the fingerprint intervals before splitting? Generally true
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
`filterBamPreSplit.jobMemory`|Int|16|memory allocated for Job
`filterBamPreSplit.overhead`|Int|6|memory allocated to overhead of the job other than used in the filter command
`filterBamPreSplit.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`splitLanes.jobMemory`|Int|16|memory allocated for Job
`splitLanes.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`splitStringToArray.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`splitStringToArray.recordSeparator`|String|"+"|Interval interval group separator - this can be used to combine multiple intervals into one group.
`splitStringToArray.jobMemory`|Int|1|Memory allocated to job (in GB).
`splitStringToArray.threads`|Int|1|The number of threads to allocate to the job.
`splitStringToArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`splitStringToArray.modules`|String|""|Environment module name and version to load (space separated) before command execution.
`filterBamLane.jobMemory`|Int|16|memory allocated for Job
`filterBamLane.overhead`|Int|6|memory allocated to overhead of the job other than used in the filter command
`filterBamLane.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`markDuplicates.jobMemory`|Int|16|memory allocated for Job
`markDuplicates.overhead`|Int|6|memory allocated to overhead of the job other than used in markDuplicates command
`markDuplicates.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`mergeIntervalBams.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeIntervalBams.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeIntervalBams.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeIntervalBams.threads`|Int|1|The number of threads to allocate to the job.
`mergeIntervalBams.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`alignmentMetrics.jobMemory`|Int|8|memory allocated for Job
`alignmentMetrics.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`extractFingerprint.jobMemory`|Int|8|memory allocated for Job
`extractFingerprint.timeout`|Int|24|Timeout in hours, needed to override imposed limits
`fingerprintReadgroupInfo.jobMemory`|Int|8|memory allocated for job
`fingerprintReadgroupInfo.timeout`|Int|24|timeout in hours


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
```
```
  set -euo pipefail
  EXT=$(basename ~{inputBam} | rev | cut -d. -f1 | rev)
  ln -s ~{inputBam} input.$EXT
  if [ "$EXT" = "cram" ]; then ln -s ~{inputBai} input.cram.crai
  else                          ln -s ~{inputBai} input.bam.bai
  fi
  samtools view -b -T ~{refFasta} -L ~{intervalBed} input.$EXT > ~{outputFileNamePrefix}.filtered.bam
  samtools index ~{outputFileNamePrefix}.filtered.bam
```
```
    set -euo pipefail

    echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
```
```
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
```
```
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
```
```
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
```
```
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
```
```
    set -euo pipefail
    samtools view -H -T ~{refFasta} ~{inputBam} \
      | awk -F'\t' '!found && /^@RG/ { for (i=1; i<=NF; i++) if ($i ~ /^ID:/) { sub(/^ID:/, "", $i); print $i; found=1 } }'
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
