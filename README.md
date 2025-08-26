# bamQC

bamQC workflow collects a number of metrics which are computed using several methods (by employing third-party software tools along with some custom code) and outputs the results in JSON format. The output also contains metadata, such as the instrument and lane names. bamQC supports downsampling for faster analysis.

## Overview

## Dependencies

* [samtools 1.16.1](https://github.com/samtools/samtools)
* [samblaster 0.1.26](https://github.com/GregoryFaust/samblaster)
* [python 3.6](https://www.python.org/downloads/)
* [bam-qc-metrics 0.2.6](https://github.com/oicr-gsi/bam-qc-metrics.git)
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
`mode`|String|running mode for the workflow, only allow value 'lane_level' and 'call_ready'
`reference`|String|Reference id, we need it to pick the right reference file


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`targetBed`|String?|None|Path to optional target bed file
`outputFileNamePrefix`|String|"bamQC"|Prefix for output files
`runBedtools`|Boolean|false|This is to run bedtools when targetBed is supplied, collect reads on target metric. Default False
`downsampleToReads`|Int|500000|Downsample to this many reads when running unique read count, duplicate rate calculation and CIGAR analysis


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`samstats.jobMemory`|Int|16|Memory allocated for this job
`samstats.timeout`|Int|24|hours before task timeout
`runMosdepth.modules`|String|"mosdepth/0.2.9"|required environment modules
`runMosdepth.jobMemory`|Int|16|Memory allocated for this job
`runMosdepth.timeout`|Int|4|hours before task timeout
`runBedtoolsIntersect.modules`|String|"samtools/1.16.1 bedtools/2.27"|required environment modules
`runBedtoolsIntersect.jobMemory`|Int|16|Memory allocated for this job
`runBedtoolsIntersect.timeout`|Int|12|hours before task timeout
`cumulativeDistToHistogram.modules`|String|"python/3.6"|required environment modules
`cumulativeDistToHistogram.jobMemory`|Int|8|Memory allocated for this job
`cumulativeDistToHistogram.threads`|Int|4|Requested CPU threads
`cumulativeDistToHistogram.timeout`|Int|1|hours before task timeout
`markDuplicates.additionalParam`|String?|None|Any additional parameters for samblaster
`markDuplicates.modules`|String|"samtools/1.16.1"|required environment modules
`markDuplicates.jobMemory`|Int|16|Memory allocated for this job
`markDuplicates.threads`|Int|4|Requested CPU threads
`markDuplicates.timeout`|Int|4|hours before task timeout
`getUniqueReadsCount.modules`|String|"samtools/1.16.1"|required environment modules
`getUniqueReadsCount.jobMemory`|Int|16|Memory allocated for this job
`getUniqueReadsCount.timeout`|Int|18|hours before task timeout
`bamQCMetrics.targetBed`|File?|None|Optional target bed
`bamQCMetrics.workflowVersion`|String|"5.3.0"|Workflow version to put into report
`bamQCMetrics.modules`|String|"bam-qc-metrics/0.2.6"|required environment modules
`bamQCMetrics.bamQClite`|String|"$BAM_QC_METRICS_ROOT/bin/run_bam_qc_lite.py"|Path to bamQC lite script
`bamQCMetrics.jobMemory`|Int|8|Memory allocated for this job
`bamQCMetrics.timeout`|Int|12|hours before task timeout
`mergeReports.modules`|String|"bam-qc-metrics/0.2.6"|Runtime modules
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
     MOSDEPTH_PRECISION=8 mosdepth -x -n -t 3 bamqc ~{bamFileName} ~{"--by " + targetBed}
```
 
### Extract coverage histogram
 
```
         python3 <<CODE
         import csv, json
         summary = open("~{summary}").readlines()
         globalDist = open("~{globalDist}").readlines()
         # read chromosome lengths from the summary
         summaryReader = csv.reader(summary, delimiter="\t")
         lengthByChr = {}
         for row in summaryReader:
           if row[0] == 'chrom' or row[0] == 'total':
             continue # skip initial header row, and final total row
           if row[0].endswith('_region'):
             continue # skip contigs from target file, if passed
           lengthByChr[row[0]] = int(row[1])
         chromosomes = sorted(lengthByChr.keys())
         # read the cumulative distribution for each chromosome
         globalReader = csv.reader(globalDist, delimiter="\t")
         cumDist = {}
         for k in chromosomes:
           cumDist[k] = {}
         for row in globalReader:
           if row[0]=="total":
             continue
           cumDist[row[0]][int(row[1])] = float(row[2])
         # convert the cumulative distributions to non-cumulative and populate histogram
         # if the input BAM is empty, chromosomes and histogram will also be empty
         histogram = {}
         for k in chromosomes:
           depths = sorted(cumDist[k].keys())
           dist = {}
           for i in range(len(depths)-1):
             depth = depths[i]
             nextDepth = depths[i+1]
             dist[depth] = cumDist[k][depth] - cumDist[k][nextDepth]
           maxDepth = max(depths)
           dist[maxDepth] = cumDist[k][maxDepth]
           # now find the number of loci at each depth of coverage to construct the histogram
           for depth in depths:
             loci = int(round(dist[depth]*lengthByChr[k], 0))
             histogram[depth] = histogram.get(depth, 0) + loci
         # if histogram is non-empty, fill in zero values for missing depths
         for i in range(max(histogram.keys(), default=0)):
           if i not in histogram:
             histogram[i] = 0
         out = open("~{outFileName}", "w")
         json.dump(histogram, out, sort_keys=True)
         out.close()
         CODE
```
 
### Extract duplicate reads metrics with samblaster
 
samblaster uses downsampled data piped into the tool. This allows having a much reduced filesystem footprint
as there are no intermediate files generated in the process
 
```
     set -euxo pipefail
     samtools head -n ~{downsampleToReads} ~{bamFile} | \
     samtools sort -n - | \
     samtools fixmate -m -O SAM - - | \
     samblaster --ignoreUnmated ~{additionalParam} --output /dev/null 2> >(tee "~{filePrefix}.markDuplicates.txt")
```
 
### For targeted sequencing, count reads on target with bedtools
 
```
     bedtools intersect -a ~{inputBam} -b ~{targetBed} -u | samtools view -c | perl -pe 'chomp'
```
 
### Count unique reads on downsampled data (needed for CIGAR metrics)
 
CIGAR metrics are generated using downsampled bam, this step uses the same rate of downsampling
 
```
     samtools head -n ~{downsampleToReads} ~{inputBam} | samtools view -F 256 -q 30 -c | perl -pe 'chomp'
```
 
### Run bamQC lite which aggregates metrics into lane-level json report
 
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
         python3 ~{bamQCmerger} -l ~{sep="," inputs} -o ~{outputFileName}
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
