# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 5.2.0 - 2024-06-25
### Added
- Add vidarr labels to outputs (changes to medata only).
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) 

## 5.1.3 -2024-02-12
### Added
- Add minMemory parameter to RAM-scaling tasks.
- [GRD-728](https://jira.oicr.on.ca/browse/GRD-728)

## 5.1.2 -2023-12-18
### Added
- Add RAM scaling by chromosome.
- [GRD-728](https://jira.oicr.on.ca/browse/GRD-728)

## 5.1.1 -2023-10-12
### Changed
- Update bamQC to reflect previous version.
- [GBS-4314] https://jira.oicr.on.ca/browse/GBS-4314

## 5.1.0 -2023-08-17
### Added
- Create two modes lane_levl and call_ready, for different inputs. The formal run bamQC steps as before woth single lane level bam, the latter filter and multiple baminputs before other bamQC steps. 

## 5.0.2 - 2022-08-31
### Changed
- Modify python module loading during validation.

## 5.0.1 - 2022-08-31
### Changed
- Update workflow name to `bamqc_call_ready_by_tumor_group` for clinical.

## 5.0.0 - 2021-06-01
### Changed
- Migration from Niassa to Vidarr.

## 4.0.5 - 2020-11-09
### Fixed
- Bugfix to fix `coverage histogram`. 

### Removed
- Remove invalid `coverage_histogram`.
- [GR-1243](https://jira.oicr.on.ca/browse/GR-1243) 

## 4.0.4 - 2020-09-09
### Fixed
- Bugfix for empty BAM input, with test.
- [GP-2469](https://jira.oicr.on.ca/browse/GP-2469)

## 4.0.3 - 2020-06-08
### Changed
- Use the mosdepth tool for fast calculation of coverage depth.
- [GP-2397](https://jira.oicr.on.ca/browse/GP-2397)

## 4.0.2 - 2020-05-07
### Fixed
- Fix another integer overflow error.
- [GP-2374](https://jira.oicr.on.ca/browse/GP-2374)

### Added
- Update README.md with `gsi-wdl-tools`.

## 4.0.1 - 2020-04-24
### Fixed
- Fix integer overflow error in Cromwell.
- [GP-2353](https://jira.oicr.on.ca/browse/GP-2353) 

## 4.0.0 - 2020-03-18
### Added
- WDL workflow for bamQC.
- [GP-2123](https://jira.oicr.on.ca/browse/GP-2123) 
- Runs Picard MarkDuplicates and bam-qc-metrics 0.2.5
- New filtering and downsampling logic; see `filter_downsample.md`.

### Removed
- Remove decider and convert previous Niassa workflow to Niassa-WDL.

## 3.0.3 - 2020-01-16
### Changed
- Update to bam-qc-metrics 0.2.3.
- [GP-2245](https://jira.oicr.on.ca/browse/GP-2245) 

## 3.0.2 - 2019-12-04
### Changed
- Update to bam-qc-metrics 0.2.2.
- [GR-908](https://jira.oicr.on.ca/browse/GR-908)

## 3.0.1 - 2019-11-27
### Changed
- Update to bam-qc-metrics 0.2.1.
- Use provided duplicate marked bam metrics input file (if provided) 
- Update to picard 2.21.2.
- [GP-2190](https://jira.oicr.on.ca/browse/GP-2190)

## 3.0 - 2019-09-13
### Changed
- Move to bam-qc-metrics version 0.2.0. 
- Complete reimplementation of QC metrics.

## 2.6 - 2019-06-11
### Changed
- Disable downsampling for smaller input files (default = 100 MB).

## 2.5 - 2015-09-08
### Added
- Json metadata must now be provided as an ini parameter: "json_metadata".

## 2.4 - 2015-03-19
### Changed
- Upgrade to SeqWare 1.1.0, common-utilities 1.6 and workflow-utilities 1.6.
