# bamQC

QC metrics for BAM files

## Overview

## Dependencies

* [samtools 1.14](https://github.com/samtools/samtools)
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
`mode`|String|running mode for the workflow, only allow value 'lane_level' and 'call_ready'
`bamQCMetrics.refFasta`|String|Path to human genome FASTA reference
`bamQCMetrics.refSizesBed`|String|Path to human genome BED reference with chromosome sizes
`bamQCMetrics.workflowVersion`|String|Workflow version string


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|"bamQC"|Prefix for output files
`intervalsToParallelizeByString`|String|"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3,chr4).


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitStringToArray.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`splitStringToArray.recordSeparator`|String|"+"|Record separator - a delimiter for joining records
`splitStringToArray.jobMemory`|Int|1|Memory allocated to job (in GB).
`splitStringToArray.cores`|Int|1|The number of cores to allocate to the job.
`splitStringToArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`splitStringToArray.modules`|String|""|Environment module name and version to load (space separated) before command execution.
`getChrCoefficient.memory`|Int|2|Memory allocated for this job
`getChrCoefficient.timeout`|Int|1|Hours before task timeout
`getChrCoefficient.modules`|String|"samtools/1.14"|Names and versions of modules to load
`preFilter.filterFlags`|Int|260|Samtools filter flags to apply.
`preFilter.minMapQuality`|Int?|None|Minimum alignment quality to pass filter
`preFilter.filterAdditionalParams`|String?|None|Additional parameters to pass to samtools.
`preFilter.modules`|String|"samtools/1.14"|required environment modules
`preFilter.jobMemory`|Int|6|Memory allocated for this job
`preFilter.threads`|Int|4|Requested CPU threads
`preFilter.timeout`|Int|4|hours before task timeout
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
`filter.minQuality`|Int|30|Minimum alignment quality to pass filter
`filter.modules`|String|"samtools/1.14"|required environment modules
`filter.jobMemory`|Int|16|Memory allocated for this job
`filter.threads`|Int|4|Requested CPU threads
`filter.timeout`|Int|4|hours before task timeout
`updateMetadata.modules`|String|"python/3.6"|required environment modules
`updateMetadata.jobMemory`|Int|16|Memory allocated for this job
`updateMetadata.threads`|Int|4|Requested CPU threads
`updateMetadata.timeout`|Int|4|hours before task timeout
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
`downsample.modules`|String|"samtools/1.14"|required environment modules
`downsample.jobMemory`|Int|16|Memory allocated for this job
`downsample.threads`|Int|4|Requested CPU threads
`downsample.timeout`|Int|4|hours before task timeout
`downsampleRegion.modules`|String|"samtools/1.14"|required environment modules
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


This section lists command(s) run by WORKFLOW workflow
  
### Run getChrCoefficient
 
```
 
  CHROM_LEN=$(samtools view -H ~{bamFile} | grep ^@SQ | grep -v _ | grep -w ~{chromosome} | cut -f 3 | sed 's/LN://')
  LARGEST=$(samtools view -H ~{bamFile} | grep ^@SQ | grep -v _ | cut -f 3 | sed 's/LN://' | sort -n | tail -n 1)
  echo | awk -v chrom_len=$CHROM_LEN -v largest=$LARGEST '{print int((chrom_len/largest + 0.1) * 10)/10}'
 
```
 
### Run bamQCMetrics
   
```
  run_bam_qc.py \
   -b ~{bamFile} \
   -d ~{markDuplicates} \
   --debug \
   -i ~{normalInsertMax} \
   -o ~{resultName} \
   -r ~{refFasta} \
   -t ~{refSizesBed} \
   -T . \
   -w ~{workflowVersion} \
   ~{dsInput}
```
### Run collateResults

```
  python3 <<CODE
  import json
  data = json.loads(open("~{bamQCMetricsResult}").read())
  histogram = json.loads(open("~{histogram}").read())
  data["coverage histogram"] = histogram
  metadata = json.loads(open("~{metadata}").read())
  for key in metadata.keys():
      data[key] = metadata[key]
  out = open("~{outputFileName}", "w")
  json.dump(data, out, sort_keys=True)
  out.close()
  CODE
```

### Run cumulativeDistToHistogram

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

### Run downsample 

```
  	set -e
  	set -o pipefail
  	samtools view -b -h ~{bamFile} | \
  	~{preDownsampleCommand} ~{downsampleCommand} \
  	samtools view -b > ~{resultName}
```

### downsampleRegion

```
  	set -e
  	# ensure BAM file and index are symlinked to working directory
  	ln -s ~{bamFile}
  	ln -s ~{bamIndex}
  	samtools view -b -h ~{bamFileName} ~{region} > ~{resultName}
```

### Run filter

```
  	set -e
  	set -o pipefail
  	samtools view -h -b -F 2304 -U ~{nonPrimaryReadsFile} ~{bamFile} | \
  	samtools view -h -b -F 4 -U ~{unmappedReadsFile} | \
  	samtools view -h -b -q ~{minQuality} -U ~{lowQualityReadsFile} \
  	> ~{resultName}
  	samtools view -c ~{bamFile} > ~{totalInputReadsFile}
  	samtools view -c ~{nonPrimaryReadsFile} > ~{totalNonPrimaryReadsFile}
  	samtools view -c ~{unmappedReadsFile} > ~{totalUnmappedReadsFile}
  	samtools view -c ~{lowQualityReadsFile} > ~{totalLowQualityReadsFile}
    samtools index ~{bamFile} > ~{resultIndexName}
```

### Run findDownsampleParams

```
  python3 <<CODE
  import json, math, sys
  readsIn = ~{inputReads}
  readsTarget = ~{targetReads}
  precision = ~{precision}
  print("Input reads param =", readsIn, file=sys.stderr)
  print("Target reads param =", readsTarget, file=sys.stderr)
  minReadsAbsolute = ~{minReadsAbsolute}
  minReadsRelative = ~{minReadsRelative}
  preDownsampleMultiplier = ~{preDSMultiplier}
  if readsIn <= readsTarget:
    # absolutely no downsampling
    applyPreDownsample = False
    applyDownsample = False
    preDownsampleTarget = "no_pre_downsample"
    downSampleTarget = "no_downsample"
  elif readsIn < readsTarget * minReadsRelative or readsTarget < minReadsAbsolute:
    # no predownsampling
    applyPreDownsample = False
    applyDownsample = True
    preDownsampleTarget = "no_pre_downsample"
    downSampleTarget = str(readsTarget)
  else:
    # predownsampling and downsampling
    applyPreDownsample = True
    applyDownsample = True
    probability = (readsTarget * preDownsampleMultiplier)/readsIn
    formatString = "{:0"+str(precision)+"d}"
    preDownsampleTarget = formatString.format(int(math.floor(probability * 10**precision)))
    downSampleTarget = str(readsTarget)
  status = {
    "pre_ds": applyPreDownsample,
    "ds": applyDownsample
  }
  targets = {
    "pre_ds": preDownsampleTarget,
    "ds": downSampleTarget
  }
  statusFile = open("~{statusFile}", "w")
  json.dump(status, statusFile)
  statusFile.close()
  targetFile = open("~{targetsFile}", "w")
  json.dump(targets, targetFile)
  targetFile.close()
  CODE
```

### Run findDownsampleParamsMarkDup

```
  python3 <<CODE
  readsIn = ~{inputReads}
  threshold = ~{threshold}
  interval = ~{baseInterval}
  start = ~{intervalStart} + 1 # start of sub-chromosome window, if needed; exclude telomeres
  chromosomes = [line.strip() for line in open("~{chromosomesText}").readlines()]
  customRegions = "~{customRegions}" # overrides other chromosome/interval parameters
  ds = True # True if downsampling, false otherwise
  end = None # end of window, if needed
  if readsIn <= threshold:
      ds = False # no downsampling
  elif readsIn <= threshold*10:
      pass # default to chr12 & chr13 =~ 8% of genome
  elif readsIn <= threshold*10**2:
      end = start + interval*10**3 - 1 # default 2*15 million base window ~ 1% of genome
  elif readsIn <= threshold*10**3:
      end = start + interval*10**2 - 1
  elif readsIn <= threshold*10**4:
      end = start + interval*10 - 1
  else:
      end = start + interval - 1
  if ds:
      status = "true"
      if customRegions != "":
	  region = customRegions
      elif end == None:
	  region = " ".join(chromosomes)
      else:
	  regions = ["%s:%i-%i" % (chromosome, start, end) for chromosome in chromosomes ]
	  region = " ".join(regions)
  else:
      status = "false"
      region = ""
  outStatus = open("~{outputStatus}", "w")
  print(status, file=outStatus)
  outStatus.close()
  outRegion = open("~{outputRegion}", "w")
  print(region, file=outRegion)
  outRegion.close()
  CODE
 ```
 ### Run markDuplicates 
  ```
  	java -Xmx~{picardMaxMemMb}M \
  	-jar ${PICARD_ROOT}/picard.jar \
  	MarkDuplicates \
  	INPUT=~{bamFile} \
  	OUTPUT=~{outFileBam} \
  	VALIDATION_STRINGENCY=SILENT \
  	TMP_DIR=${PWD} \
  	METRICS_FILE=~{outFileText} \
  	OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance}
```

### Run MergeFiles
 
```
      set -euo pipefail
  
      gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
      ~{sep=" " prefix("--INPUT=", bams)} \
      --OUTPUT="~{outputFileName}~{suffix}.bam" \
      --CREATE_INDEX=true \
      --SORT_ORDER=coordinate \
      --ASSUME_SORTED=false \
      --USE_THREADING=true \
      --VALIDATION_STRINGENCY=SILENT 
  
```

### Run prefilter

```
     set -e
     set -o pipefail
     samtools view -b \
          -F ~{filterFlags} \
          ~{"-q " + minMapQuality} \
          ~{filterAdditionalParams} \
          ~{bamFile} \
          ~{sep=" " intervals} > ~{resultName}
```

### Run runMosdepth 

```
  set -eo pipefail
  # ensure BAM file and index are symlinked to working directory
  ln -s ~{bamFile}
  ln -s ~{bamIndex}
  # run mosdepth
  MOSDEPTH_PRECISION=8 mosdepth -x -n -t 3 bamqc ~{bamFileName}
```

### splitStringToArray

```
  set -euo pipefail

  echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
```

### Run updateMetadata

```
  python3 <<CODE
  import json
  metadata = json.loads(open("~{metadataJson}").read())
  metadata["total input reads meta"] = ~{totalInputReads}
  metadata["non-primary reads meta"] = ~{nonPrimaryReads}
  metadata["unmapped reads meta"] = ~{unmappedReads}
  metadata["low-quality reads meta"] = ~{lowQualityReads}
  outFile = open("~{outFileName}", "w")
  json.dump(metadata, outFile)
  outFile.close()
  CODE
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
