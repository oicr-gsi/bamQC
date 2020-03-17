# bamQC workflow

WDL implementation of QC metrics for BAM files.

## Overview

bamQC runs Picard MarkDuplicates and GSI bam-qc-metrics on the input BAM file.

Filtering is applied to non-primary alignments, unmapped reads, and low-quality reads to exclude them from QC. For large input files, downsampling is applied separately for MarkDuplicates and bam-qc-metrics. See `filter_downsample.md` for details.

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [picard 2.21.2](https://broadinstitute.github.io/picard/command-line-overview.html)
* [python 3.6](https://www.python.org/downloads/release/python-3610/)
* [bam-qc-metrics 0.2.5](https://github.com/oicr-gsi/bam-qc-metrics.git)

## Documentation

- `bamqc_wdl.png`: Diagram of workflow structure
- `filter_downsample.md`: Filtering and downsampling logic

## Usage

### Cromwell
```
java -jar cromwell.jar run bamQC.wdl --inputs inputs.json
```
### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bamFile`|File|Input BAM file on which to compute QC metrics
`metadata`|Map[String,String]| metadata parameters

## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with wdl_doc_gen (https://github.com/oicr-gsi/wdl_doc_gen/)_
