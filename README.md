# bamQC

bamQC workflow collects a number of metrics which are computed using several methods (by employing third-party software tools along with some custom code) and outputs the results in JSON format. The output also contains metadata, such as the instrument and lane names. bamQC supports downsampling for faster analysis.

## Overview

## Dependencies

* [samtools 1.16.1](https://github.com/samtools/samtools)
* [samblaster 0.1.26](https://github.com/GregoryFaust/samblaster)
* [python 3.6](https://www.python.org/downloads/)
* [bam-qc-metrics 0.2.7](https://github.com/oicr-gsi/bam-qc-metrics.git)
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
`reference`|String|Reference id, we need it to pick the right reference file


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBed`|String?|None|Path to optional target bed file
`outputFileNamePrefix`|String|"bamQC"|Prefix for output files
`downsampleToReads`|Int|500000|Downsample to this many reads when running unique read count, duplicate rate calculation and CIGAR analysis
`coverageWindow`|Int|1000|Coverage window to use with mosdepth for making coverage histogram, default is 1000b


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`samstats.additionalParameters`|String?|None|Additional parameters for samtools, usually filtering flag
`samstats.jobMemory`|Int|16|Memory allocated for this job
`samstats.timeout`|Int|24|hours before task timeout
`runMosdepth.modules`|String|"mosdepth/0.2.9"|required environment modules
`runMosdepth.additionalParameters`|String?|None|Additional parameters for coverage calculation
`runMosdepth.jobMemory`|Int|16|Memory allocated for this job
`runMosdepth.timeout`|Int|4|hours before task timeout
`runBedtoolsIntersect.modules`|String|"samtools/1.16.1 bedtools/2.27"|required environment modules
`runBedtoolsIntersect.jobMemory`|Int|16|Memory allocated for this job
`runBedtoolsIntersect.timeout`|Int|12|hours before task timeout
`cumulativeDistToHistogram.modules`|String|"bam-qc-metrics/0.2.7"|required environment modules
`cumulativeDistToHistogram.coverageHistogram`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_histogram.py"|Path to script generating coverage histogram
`cumulativeDistToHistogram.outFileName`|String|"coverage_histogram.json"|Output file name, default coverage_histogram.json
`cumulativeDistToHistogram.jobMemory`|Int|8|Memory allocated for this job
`cumulativeDistToHistogram.threads`|Int|4|Requested CPU threads
`cumulativeDistToHistogram.timeout`|Int|1|hours before task timeout
`markDuplicates.additionalParameters`|String?|None|Any additional parameters for samblaster
`markDuplicates.modules`|String|"samblaster/0.1.26 samtools/1.16.1"|required environment modules
`markDuplicates.jobMemory`|Int|16|Memory allocated for this job
`markDuplicates.threads`|Int|4|Requested CPU threads
`markDuplicates.timeout`|Int|4|hours before task timeout
`getUniqueReadsCount.modules`|String|"samtools/1.16.1"|required environment modules
`getUniqueReadsCount.jobMemory`|Int|16|Memory allocated for this job
`getUniqueReadsCount.timeout`|Int|18|hours before task timeout
`downsampleBam.modules`|String|"samtools/1.16.1"|required environment modules
`downsampleBam.jobMemory`|Int|16|Memory allocated for this job
`downsampleBam.timeout`|Int|4|hours before task timeout
`runWindowedMosdepth.modules`|String|"mosdepth/0.2.9"|required environment modules
`runWindowedMosdepth.targetBed`|String?|None|Optional target bed file
`runWindowedMosdepth.jobMemory`|Int|16|Memory allocated for this job
`runWindowedMosdepth.timeout`|Int|4|hours before task timeout
`bamQCMetrics.workflowVersion`|String|"5.3.1"|Workflow version to put into report
`bamQCMetrics.modules`|String|"bam-qc-metrics/0.2.7"|required environment modules
`bamQCMetrics.bamQClite`|String|"$BAM_QC_METRICS_ROOT/bin/run_bam_qc_lite.py"|Path to bamQC lite script
`bamQCMetrics.jobMemory`|Int|8|Memory allocated for this job
`bamQCMetrics.timeout`|Int|12|hours before task timeout
`markDuplicatesMerged.additionalParameters`|String?|None|Any additional parameters for samblaster
`markDuplicatesMerged.modules`|String|"samblaster/0.1.26 samtools/1.16.1"|required environment modules
`markDuplicatesMerged.jobMemory`|Int|16|Memory allocated for this job
`markDuplicatesMerged.threads`|Int|4|Requested CPU threads
`markDuplicatesMerged.timeout`|Int|4|hours before task timeout
`mergedCoverageToHistogram.modules`|String|"bam-qc-metrics/0.2.7"|required environment modules
`mergedCoverageToHistogram.coverageMerge`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_merger.py"|Path to coverage merging script
`mergedCoverageToHistogram.jobMemory`|Int|8|Memory allocated for this job
`mergedCoverageToHistogram.threads`|Int|4|Requested CPU threads
`mergedCoverageToHistogram.timeout`|Int|1|hours before task timeout
`mergeReports.modules`|String|"bam-qc-metrics/0.2.7"|Runtime modules
`mergeReports.bamQCmerger`|String|"$BAM_QC_METRICS_ROOT/bin/bam_qc_merger.py"|Path to the merger script
`mergeReports.jobMemory`|Int|4|RAM allocated to run the merging task
`mergeReports.timeout`|Int|2|Timeout in hours for the merging task


### Outputs

Output | Type | Description | Labels
---|---|---|---
`result`|File|json file that contains metrics and meta data described in https://github.com/oicr-gsi/bam-qc-metrics/blob/master/metrics.md|vidarr_label: result


## Commands

This section lists command(s) run by bamQC lite workflow
 
* Running bamQC lite
 
bamQC lite gets most of it's metrics from external tools (except the ones it generates by running CIGAR analysis on downsampled data). 
It has fewer steps then the original bamQC but produces almost exact results.
 
### Run samstats on unaltered inputs
 
```
     samtools stats ~{bamFile} -r  ~{referenceFile} > ~{filePrefix}.stats
```
 
### Run mosdepth to get coverage metrics
 
```
     set -eo pipefail
     # ensure BAM file and index are symlinked to working directory
     ln -s ~{bamFile}
     ln -s ~{bamIndex}
     # run mosdepth
     MOSDEPTH_PRECISION=8 mosdepth -x -n ~{"--by " + targetBed} ~{additionalParameters} -t 3 ~{prefix} ~{bamFileName}
```
 
### Extract coverage histogram
 
This step creates a lane-level coverage hitogram for each of the lane-level inputs
 
```
     python3 ~{coverageHistogram} -s ~{summary} -g ~{globalDist} -o ~{outFileName}
```
 
### Merge windowed coverage, generate coverage histogram
 
Windowed coverage is used to create a call-ready coverage histogram. mosDepth produces windowed coverage in a format 
suitable for merging if we have multiple lane-level inputs
 
```
     python3 ~{coverageMerge} -f ~{sep="," coverageFiles} \
                              -o ~{outFileName} 
```
 
### Duplicate read marking with samblaster, lane-level mode
 
This task runs on a lane-level BAM and produces Duplicate Read metrics later injected into the final report (in lane-level mode)
 
```
     set -euxo pipefail
     samtools head -n ~{downsampleToReads} ~{bamFile} | \
     samtools view -h -F 2308 - | \
     samtools sort -n - | \
     samtools fixmate -m -O SAM - - | \
     samblaster --ignoreUnmated ~{additionalParameters} --output /dev/null 2> >(tee "~{filePrefix}.markDuplicates.txt")
```
 
### Duplicate read marking with samblaster, call-ready mode
 
This task accepts multiple (downsampled) BAM files and merges them on the fly, piping the results into samblaster
 
```
     set -euxo pipefail
     samtools cat ~{sep=" " bamFiles} | \
     samtools view -h -F 2308 - | \
     samtools sort -n - | \
     samtools fixmate -m -O SAM - - | \
     samblaster --ignoreUnmated ~{additionalParameters} --output /dev/null 2> >(tee "~{outputFileNamePrefix}.markDuplicates.txt")
```
 
### For targeted sequencing, count reads on target with bedtools
 
 This runs only if we have a target .bed file supplied (targeted sequencing mode)
 
```
     bedtools intersect -a ~{inputBam} -b ~{targetBed} -u | samtools view -c | perl -pe 'chomp'
```
 
### Count unique reads on downsampled data (needed for CIGAR metrics)
 
CIGAR metrics are generated using downsampled bam, this step uses the same rate of downsampling
 
```
     samtools head -n ~{downsampleToReads} ~{inputBam} | samtools view -F 256 -q 30 -c | perl -pe 'chomp'
```
 
### Run bamQC lite which aggregates metrics into lane-level json report
 
This creates a lane-level report. In lane-level mode this is going to be the final output from the workflow.
In call-ready mode these reports will be merged and the final merged metrics provisioned.
 
```
         set -euxo pipefail
         python3 ~{bamQClite} ~{"-b " + bamFile} \
         -s ~{samstatsFile} \
         -d ~{markDuplicatesStats} \
         -c ~{histogram} \
         -m ~{metadataJson} \
         -w ~{workflowVersion} ~{"-tf " + targetBed} \
         -r ~{referenceFileName} ~{"-S " + downsampleToReads} \
         -o ~{outputFileName} \
         -t ~{mosdepthSummary} ~{"-u " + uniqueReads} ~{"-ot " + readsOnTarget} 
```
 
### Run merger script which will combine reports into call-ready report if there are multiple lanes
 
This step will return lane-level report (exact copy of it's input) if there is one lane, but will
combine metrics into call-ready report if there are multiple lanes
 
```
         set -euxo pipefail
         python3 ~{bamQCmerger} -l ~{sep="," inputs} ~{"-d " + mergedDupmarkingData} ~{"-t " + mergedCoverageData} -o ~{outputFileName}
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
