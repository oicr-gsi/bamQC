# bamQC

QC metrics for BAM files

## Overview

bamQC runs Picard MarkDuplicates and GSI bam-qc-metrics on the input BAM file. See `bamqc_wdl.png` for workflow structure.

Filtering is applied to non-primary alignments, unmapped reads, and low-quality reads to exclude them from QC. For large input files, downsampling is applied separately for MarkDuplicates and bam-qc-metrics. See `filter_downsample.md` for details.

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [picard 2.21.2](https://broadinstitute.github.io/picard/command-line-overview.html)
* [python 3.6](https://www.python.org/downloads/)
* [bam-qc-metrics 0.2.5](https://github.com/oicr-gsi/bam-qc-metrics.git)


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
`metadata`|Map[String,String]|JSON file containing metadata


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|"bamQC"|Prefix for output files


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filter.minQuality`|Int|30|Minimum alignment quality to pass filter
`filter.modules`|String|"samtools/1.9"|required environment modules
`filter.jobMemory`|Int|16|Memory allocated for this job
`filter.threads`|Int|4|Requested CPU threads
`filter.timeout`|Int|4|hours before task timeout
`updateMetadata.modules`|String|"python/3.6"|required environment modules
`updateMetadata.jobMemory`|Int|16|Memory allocated for this job
`updateMetadata.threads`|Int|4|Requested CPU threads
`updateMetadata.timeout`|Int|4|hours before task timeout
`countInputReads.modules`|String|"samtools/1.9"|required environment modules
`countInputReads.jobMemory`|Int|16|Memory allocated for this job
`countInputReads.threads`|Int|4|Requested CPU threads
`countInputReads.timeout`|Int|4|hours before task timeout
`findDownsampleParams.targetReads`|Int|100000|Desired number of reads in downsampled output
`findDownsampleParams.minReadsAbsolute`|Int|10000|Minimum value of targetReads to allow pre-downsampling
`findDownsampleParams.minReadsRelative`|Int|2|Minimum value of (inputReads)/(targetReads) to allow pre-downsampling
`findDownsampleParams.precision`|Int|8|Number of decimal places in fraction for pre-downsampling
`findDownsampleParams.preDSMultiplier`|Float|1.5|Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative).
`findDownsampleParams.modules`|String|"python/3.6"|required environment modules
`findDownsampleParams.jobMemory`|Int|16|Memory allocated for this job
`findDownsampleParams.threads`|Int|4|Requested CPU threads
`findDownsampleParams.timeout`|Int|4|hours before task timeout
`findDownsampleParamsMarkDup.threshold`|Int|10000000|Minimum number of reads to conduct downsampling
`findDownsampleParamsMarkDup.chromosomes`|Array[String]|["chr12", "chr13", "chrXII", "chrXIII"]|Array of chromosome identifiers for downsampled subset
`findDownsampleParamsMarkDup.baseInterval`|Int|15000|Base width of interval in each chromosome, for very large BAMs
`findDownsampleParamsMarkDup.intervalStart`|Int|100000|Start of interval in each chromosome, for very large BAMs
`findDownsampleParamsMarkDup.customRegions`|String|""|Custom downsample regions; overrides chromosome and interval parameters
`findDownsampleParamsMarkDup.modules`|String|"python/3.6"|required environment modules
`findDownsampleParamsMarkDup.jobMemory`|Int|16|Memory allocated for this job
`findDownsampleParamsMarkDup.threads`|Int|4|Requested CPU threads
`findDownsampleParamsMarkDup.timeout`|Int|4|hours before task timeout
`downsample.downsampleSuffix`|String|"downsampled.bam"|Suffix for output file
`downsample.randomSeed`|Int|42|Random seed for pre-downsampling (if any)
`downsample.modules`|String|"samtools/1.9"|required environment modules
`downsample.jobMemory`|Int|16|Memory allocated for this job
`downsample.threads`|Int|4|Requested CPU threads
`downsample.timeout`|Int|4|hours before task timeout
`downsampleRegion.modules`|String|"samtools/1.9"|
`downsampleRegion.jobMemory`|Int|16|
`downsampleRegion.threads`|Int|4|
`downsampleRegion.timeout`|Int|4|
`markDuplicates.opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`markDuplicates.picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`markDuplicates.modules`|String|"picard/2.21.2"|required environment modules
`markDuplicates.jobMemory`|Int|16|Memory allocated for this job
`markDuplicates.threads`|Int|4|Requested CPU threads
`markDuplicates.timeout`|Int|4|hours before task timeout
`bamQCMetrics.refFasta`|String|None|Path to human genome FASTA reference
`bamQCMetrics.refSizesBed`|String|None|Path to human genome BED reference with chromosome sizes
`bamQCMetrics.workflowVersion`|String|None|Workflow version string
`bamQCMetrics.normalInsertMax`|Int|1500|Maximum of expected insert size range
`bamQCMetrics.modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`bamQCMetrics.jobMemory`|Int|16|Memory allocated for this job
`bamQCMetrics.threads`|Int|4|Requested CPU threads
`bamQCMetrics.timeout`|Int|4|hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`result`|File|None


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
