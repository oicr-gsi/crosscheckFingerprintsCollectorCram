# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-07-28
### Added
- [GBS-6863](https://jira.oicr.on.ca/browse/GBS-6863) - add cram inputs, split the merged cram to lanes (by readgroup) and generate fingerprints for each

### Changed
- [GBS-7189](https://jira.oicr.on.ca/browse/GBS-7189) - crosscheckFingerprintsCollectorMultiLane no longer splits the merged input into per-lane files. The lane set is read from the @RG headers and each lane is streamed out of the merged cram/bam by read group, restricted to the fingerprint intervals, straight into gatk ExtractFingerprint. This drops the cram-to-bam conversion, the per-lane bams and their indexes. Duplicate marking, which cannot read a stream, keeps a task-local lane bam and no longer scatters by chromosome.

