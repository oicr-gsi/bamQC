# CHANGELOG

## 4.0.3 - 2020-06-08
- [GP-2397](https://jira.oicr.on.ca/browse/GP-2397) - Use the mosdepth tool for fast calculation of coverage depth
## 4.0.2 - 2020-05-07
- [GP-2374](https://jira.oicr.on.ca/browse/GP-2374) - Fix another integer overflow error
- Update README.md with `gsi-wdl-tools`
## 4.0.1 - 2020-04-24
- [GP-2353](https://jira.oicr.on.ca/browse/GP-2353) - Fix integer overflow error in Cromwell
## 4.0.0 - 2020-03-18
- [GP-2123](https://jira.oicr.on.ca/browse/GP-2123) - WDL workflow for bamQC
- Remove decider and convert previous Niassa workflow to Niassa-WDL
- Runs Picard MarkDuplicates and bam-qc-metrics 0.2.5
- New filtering and downsampling logic; see `filter_downsample.md`
## 3.0.3 - 2020-01-16
- [GP-2245](https://jira.oicr.on.ca/browse/GP-2245) - Update to bam-qc-metrics 0.2.3
## 3.0.2 - 2019-12-04
- [GR-908](https://jira.oicr.on.ca/browse/GR-908) - Update to bam-qc-metrics 0.2.2
## 3.0.1 - 2019-11-27
- Update to bam-qc-metrics 0.2.1
- [GP-2190](https://jira.oicr.on.ca/browse/GP-2190) - Use provided duplicate marked bam metrics input file (if provided) and update to picard 2.21.2
## 3.0 - 2019-09-13
- Move to bam-qc-metrics version 0.2.0; complete reimplementation of QC metrics
## 2.6 - 2019-06-11
- Disable downsampling for smaller input files (default = 100 MB)
## 2.5 - 2015-09-08
- Json metadata must now be provided as an ini parameter: "json_metadata"
## 2.4 - 2015-03-19
- Upgrade to SeqWare 1.1.0, common-utilities 1.6 and workflow-utilities 1.6.
