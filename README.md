# bamQC

QC metrics for BAM files

## Overview
bamQC runs the following tools:
- *Picard MarkDuplicates:* Identifies duplicate reads
- *mosdepth:* Fast tool to compute depth of coverage
- *bam-qc-metrics:* Package developed at GSI to compute an assortment of samtools, bedtools, and custom metrics on BAM files.

See `bamqc_wdl.png` for workflow structure.

Filtering is applied to non-primary alignments, unmapped reads, and low-quality reads to exclude them from QC. For large input
files, downsampling is applied separately for MarkDuplicates and bam-qc-metrics. See `filter_downsample.md` for details.

## Dependencies

* [samtools 1.9](https://github.com/samtools/samtools)
* [picard 2.21.2](https://broadinstitute.github.io/picard/command-line-overview.html)
* [python 3.6](https://www.python.org/downloads/)
* [bam-qc-metrics 0.2.5](https://github.com/oicr-gsi/bam-qc-metrics.git)
* [mosdepth 0.2.9](https://github.com/brentp/mosdepth)


## Usage

### Cromwell
```
java -jar cromwell.jar run bamQC.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputGroups`|Array[InputGroup]|Array of objects describing sets of bams to merge together and on which to compute QC metrics
`metadata`|Map[String,String]|JSON file containing metadata
`bamQCMetrics.refFasta`|String|Path to human genome FASTA reference
`bamQCMetrics.refSizesBed`|String|Path to human genome BED reference with chromosome sizes
`bamQCMetrics.workflowVersion`|String|Workflow version string


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|"bamQC"|Prefix for output files
`intervalsToParallelizeByString`|String|"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4).


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitStringToArray.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`splitStringToArray.recordSeparator`|String|"+"|Interval interval group separator - this can be used to combine multiple intervals into one group.
`splitStringToArray.jobMemory`|Int|1|Memory allocated to job (in GB).
`splitStringToArray.cores`|Int|1|The number of cores to allocate to the job.
`splitStringToArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`splitStringToArray.modules`|String|""|Environment module name and version to load (space separated) before command execution.
`filter.minQuality`|Int|30|Minimum alignment quality to pass filter
`filter.modules`|String|"samtools/1.9"|required environment modules
`filter.jobMemory`|Int|16|Memory allocated for this job
`filter.threads`|Int|4|Requested CPU threads
`filter.timeout`|Int|4|hours before task timeout
`mergeSplitByIntervalFiles.suffix`|String|".merge"|suffix to use for merged bam
`mergeSplitByIntervalFiles.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeSplitByIntervalFiles.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeSplitByIntervalFiles.cores`|Int|1|The number of cores to allocate to the job.
`mergeSplitByIntervalFiles.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mergeSplitByIntervalFiles.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`mergeFiles.suffix`|String|".merge"|suffix to use for merged bam
`mergeFiles.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeFiles.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeFiles.cores`|Int|1|The number of cores to allocate to the job.
`mergeFiles.timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mergeFiles.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
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
`downsampleRegion.modules`|String|"samtools/1.9"|required environment modules
`downsampleRegion.jobMemory`|Int|16|Memory allocated for this job
`downsampleRegion.threads`|Int|4|Requested CPU threads
`downsampleRegion.timeout`|Int|4|hours before task timeout
`markDuplicates.opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`markDuplicates.picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`markDuplicates.modules`|String|"picard/2.21.2"|required environment modules
`markDuplicates.jobMemory`|Int|16|Memory allocated for this job
`markDuplicates.threads`|Int|4|Requested CPU threads
`markDuplicates.timeout`|Int|4|hours before task timeout
`bamQCMetrics.normalInsertMax`|Int|1500|Maximum of expected insert size range
`bamQCMetrics.modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`bamQCMetrics.jobMemory`|Int|16|Memory allocated for this job
`bamQCMetrics.threads`|Int|4|Requested CPU threads
`bamQCMetrics.timeout`|Int|4|hours before task timeout
`runMosdepth.modules`|String|"mosdepth/0.2.9"|required environment modules
`runMosdepth.jobMemory`|Int|16|Memory allocated for this job
`runMosdepth.threads`|Int|4|Requested CPU threads
`runMosdepth.timeout`|Int|4|hours before task timeout
`cumulativeDistToHistogram.modules`|String|"python/3.6"|required environment modules
`cumulativeDistToHistogram.jobMemory`|Int|8|Memory allocated for this job
`cumulativeDistToHistogram.threads`|Int|4|Requested CPU threads
`cumulativeDistToHistogram.timeout`|Int|1|hours before task timeout
`collateResults.modules`|String|"python/3.6"|required environment modules
`collateResults.jobMemory`|Int|8|Memory allocated for this job
`collateResults.threads`|Int|4|Requested CPU threads
`collateResults.timeout`|Int|1|hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`result`|File|json file that contains metrics and meta data described in https://github.com/oicr-gsi/bam-qc-metrics/blob/master/metrics.md



 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
