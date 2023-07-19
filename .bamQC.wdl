version 1.0

struct InputGroup {
  File bam
  File bamIndex
}

workflow bamQC {

  input {
    Array[InputGroup] inputGroups
    Map[String, String] metadata
    String outputFileNamePrefix = "bamQC"
    String intervalsToParallelizeByString = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
  }

  parameter_meta {
	inputGroups: "Array of objects describing sets of bams to merge together and on which to compute QC metrics"
	metadata: "JSON file containing metadata"
	outputFileNamePrefix: "Prefix for output files"
    intervalsToParallelizeByString: "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)."
    }

  call splitStringToArray {
    input:
      str = intervalsToParallelizeByString
  }
  Array[Array[String]] intervalsToParallelizeBy = splitStringToArray.out

  scatter (i in inputGroups) {
    scatter (intervals in intervalsToParallelizeBy) {
      call filter {
        input:
        bamFile = i.bam,
        bamIndex = i.bamIndex,
        intervals = intervals,
        outputFileName = outputFileNamePrefix
      }
    }

	Array[File] filteredBams = filter.filteredBam
    Array[String] totalInputReadsArray = filter.totalInputReads
    Array[String] nonPrimaryReadsArray = filter.nonPrimaryReads
    Array[String] unmappedReadsArray = filter.unmappedReads
    Array[String] lowQualityReadsArray = filter.lowQualityReads

  call mergeFiles as mergeSplitByIntervalFiles {
    input:
      bams = filteredBams,
      outputFileName = outputFileNamePrefix,
      totalInputReads = totalInputReadsArray,
      nonPrimaryReads = nonPrimaryReadsArray,
      unmappedReads = unmappedReadsArray,
      lowQualityReads = lowQualityReadsArray
    }
  }

  Array[File] processedBams = mergeSplitByIntervalFiles.mergedBam
  Array[File] sumTotalInputReads = mergeSplitByIntervalFiles.sumTotalInputReads
  Array[File] sumNonPrimaryReads = mergeSplitByIntervalFiles.sumNonPrimaryReads
  Array[File] sumUnmappedReads = mergeSplitByIntervalFiles.sumUnmappedReads
  Array[File] sumLowQualityReads = mergeSplitByIntervalFiles.sumLowQualityReads

  if (length(processedBams) > 1) {
    call mergeFiles {
			input:
			bams = processedBams,
            outputFileName = outputFileNamePrefix,
            totalInputReads = sumTotalInputReads,
            nonPrimaryReads = sumNonPrimaryReads,
            unmappedReads = sumUnmappedReads,
            lowQualityReads = sumLowQualityReads
		}
  } 
	
  call updateMetadata {
	input:
	metadata = metadata,
	outputFileNamePrefix = outputFileNamePrefix,
	totalInputReads = select_first([mergeFiles.sumTotalInputReads, mergeSplitByIntervalFiles.sumTotalInputReads[0]]),
	nonPrimaryReads = select_first([mergeFiles.sumNonPrimaryReads, mergeSplitByIntervalFiles.sumNonPrimaryReads[0]]),
	unmappedReads = select_first([mergeFiles.sumUnmappedReads, mergeSplitByIntervalFiles.sumUnmappedReads[0]]),
	lowQualityReads = select_first([mergeFiles.sumLowQualityReads, mergeSplitByIntervalFiles.sumLowQualityReads[0]])
    }

    call countInputReads {
	input:
	bamFile = select_first([mergeFiles.mergedBam, mergeSplitByIntervalFiles.mergedBam[0]])
    }

    call findDownsampleParams {
	input:
	outputFileNamePrefix = outputFileNamePrefix,
	inputReads = countInputReads.result
    }

    call findDownsampleParamsMarkDup {
	input:
	outputFileNamePrefix = outputFileNamePrefix,
	inputReads = countInputReads.result
    }

    Boolean ds = findDownsampleParams.status["ds"]
    Boolean dsMarkDup = findDownsampleParamsMarkDup.status

    if (ds) {
	call downsample {
	    input:
	    bamFile = select_first([mergeFiles.mergedBam, mergeSplitByIntervalFiles.mergedBam[0]]),
	    outputFileNamePrefix = outputFileNamePrefix,
	    downsampleStatus = findDownsampleParams.status,
	    downsampleTargets = findDownsampleParams.targets,
	}
    }

    if (dsMarkDup) {
	call downsampleRegion {
	    input:
	    bamFile = select_first([mergeFiles.mergedBam, mergeSplitByIntervalFiles.mergedBam[0]]),
	    bamIndex = select_first([mergeFiles.mergedBamIndex, mergeSplitByIntervalFiles.mergedBamIndex[0]]),
	    outputFileNamePrefix = outputFileNamePrefix,
	    region = findDownsampleParamsMarkDup.region
	}
    }

    File mergedBamFile = select_first([mergeFiles.mergedBam, mergeSplitByIntervalFiles.mergedBam[0]])
    Array[File?] markDupInputs = [downsampleRegion.result, mergedBamFile]
    call markDuplicates {
	input:
	bamFile = select_first(markDupInputs),
	outputFileNamePrefix = outputFileNamePrefix
    }

    call bamQCMetrics {
	input:
	bamFile = select_first([mergeFiles.mergedBam, mergeSplitByIntervalFiles.mergedBam[0]]),
	outputFileNamePrefix = outputFileNamePrefix,
	markDuplicates = markDuplicates.result,
	downsampled = ds,
	bamFileDownsampled = downsample.result
    }

    call runMosdepth {
	input:
	bamFile = select_first([mergeFiles.mergedBam, mergeSplitByIntervalFiles.mergedBam[0]]),
	bamIndex = select_first([mergeFiles.mergedBamIndex, mergeSplitByIntervalFiles.mergedBamIndex[0]])
    }

    call cumulativeDistToHistogram {
	input:
	globalDist = runMosdepth.globalDist,
	summary = runMosdepth.summary
    }

    call collateResults {
	input:
	bamQCMetricsResult = bamQCMetrics.result,
	metadata = updateMetadata.result,
	histogram = cumulativeDistToHistogram.histogram,
	outputFileNamePrefix = outputFileNamePrefix
    }

    output {
	File result = collateResults.result
    }

  meta {
	description: "test"
  }

}

task splitStringToArray {
  input {
    String str
    String lineSeparator = ","
    String recordSeparator = "+"

    Int jobMemory = 1
    Int cores = 1
    Int timeout = 1
    String modules = ""
  }

  command <<<
    set -euo pipefail

    echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
  >>>

  output {
    Array[Array[String]] out = read_tsv(stdout())
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    str: "Interval string to split (e.g. chr1,chr2,chr3+chr4)."
    lineSeparator: "Interval group separator - these are the intervals to split by."
    recordSeparator: "Interval interval group separator - this can be used to combine multiple intervals into one group."
    jobMemory: "Memory allocated to job (in GB)."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task filter {

    # filter out non-primary, unmapped, and low-quality aligned reads
    # count the number of reads filtered out at each step
    # return filtered read counts and the filtered BAM file

    input {
	File bamFile
    File bamIndex
    Array [String] intervals
	String outputFileName
	Int minQuality = 30
	String modules = "samtools/1.9"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
    bamIndex: "Index file for the input bam"
    intervals: "One or more genomic intervals over which to operate."
	outputFileName: "Prefix for output file"
	minQuality: "Minimum alignment quality to pass filter"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileName}.filtered.bam"
    String totalInputReadsFile = "total_input_reads.txt"
    String totalNonPrimaryReadsFile = "total_non_primary_reads.txt"
    String totalUnmappedReadsFile = "total_unmapped_reads.txt"
    String totalLowQualityReadsFile = "total_low_quality_reads.txt"
    String nonPrimaryReadsFile = "non_primary_reads.bam"
    String unmappedReadsFile = "unmapped_reads.bam"
    String lowQualityReadsFile = "low_quality_reads.bam"

    # -F 2304 excludes secondary and supplementary alignments
    # -F 4 excludes unmapped reads

    command <<<
	set -e
	set -o pipefail
	samtools view -h -b -F 2304 -U ~{nonPrimaryReadsFile} ~{bamFile} ~{sep=" " intervals} | \
	samtools view -h -b -F 4 -U ~{unmappedReadsFile} | \
	samtools view -h -b -q ~{minQuality} -U ~{lowQualityReadsFile} \
	> ~{resultName}
	samtools view -c ~{bamFile} ~{sep=" " intervals} > ~{totalInputReadsFile}
	samtools view -c ~{nonPrimaryReadsFile} > ~{totalNonPrimaryReadsFile}
	samtools view -c ~{unmappedReadsFile} > ~{totalUnmappedReadsFile}
	samtools view -c ~{lowQualityReadsFile} > ~{totalLowQualityReadsFile}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    # record read totals as String, not Int, to avoid integer overflow error
    output {
	String totalInputReads = read_string("~{totalInputReadsFile}")
	String nonPrimaryReads = read_string("~{totalNonPrimaryReadsFile}")
	String unmappedReads = read_string("~{totalUnmappedReadsFile}")
	String lowQualityReads = read_string("~{totalLowQualityReadsFile}")
	File filteredBam = "~{resultName}"
    }

    meta {
	output_meta: {
	    totalInputReads: "Total reads in original input BAM file",
	    nonPrimaryReads: "Total reads excluded as non-primary",
	    unmappedReads: "Total reads excluded as unmapped",
	    lowQualityReads: "Total reads excluded as low alignment quality",
            filteredBam: "Filtered BAM file"
	}
    }

}

task mergeFiles {
  input {
    Array[File] bams
    String outputFileName
    Array[String] totalInputReads
    Array[String] nonPrimaryReads
    Array[String] unmappedReads
    Array[String] lowQualityReads
    String suffix = ".merge"
    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 24
    String modules = "gatk/4.1.6.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
    ~{sep=" " prefix("--INPUT=", bams)} \
    --OUTPUT="~{outputFileName}~{suffix}.bam" \
    --CREATE_INDEX=true \
    --SORT_ORDER=coordinate \
    --ASSUME_SORTED=false \
    --USE_THREADING=true \
    --VALIDATION_STRINGENCY=SILENT 

    echo $(("~{sep="+" totalInputReads}")) > "sumTotalInputReads.txt"
    echo $(("~{sep="+" nonPrimaryReads}")) > "sumNonPrimaryReads.txt"
    echo $(("~{sep="+" unmappedReads}")) > "sumUnmappedReads.txt"
    echo $(("~{sep="+" lowQualityReads}")) > "sumLowQualityReads.txt"  
  >>>

  output {
    File mergedBam = "~{outputFileName}~{suffix}.bam"
    File mergedBamIndex = "~{outputFileName}~{suffix}.bai"
    String sumTotalInputReads = read_string("sumTotalInputReads.txt")
    String sumNonPrimaryReads = read_string("sumNonPrimaryReads.txt")
    String sumUnmappedReads = read_string("sumUnmappedReads.txt")
    String sumLowQualityReads = read_string("sumLowQualityReads.txt")
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }
}

task bamQCMetrics {

    input {
	File bamFile
	String outputFileNamePrefix
	File markDuplicates
	Boolean downsampled
	File? bamFileDownsampled
	String refFasta
	String refSizesBed
	String workflowVersion
	Int normalInsertMax = 1500
	String modules = "bam-qc-metrics/0.2.5"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data. Not downsampled; may be filtered."
	outputFileNamePrefix: "Prefix for output file"
	markDuplicates: "Text file output from markDuplicates task"
	downsampled: "True if downsampling has been applied"
	bamFileDownsampled: "(Optional) downsampled subset of reads from bamFile."
	refFasta: "Path to human genome FASTA reference"
	refSizesBed: "Path to human genome BED reference with chromosome sizes"
	workflowVersion: "Workflow version string"
	normalInsertMax: "Maximum of expected insert size range"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String dsInput = if downsampled then "-S ~{bamFileDownsampled}" else ""
    String resultName = "~{outputFileNamePrefix}.metrics.json"

    command <<<
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
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            output1: "JSON file with bam-qc-metrics output"
	}
  }

}

task collateResults {

    input {
	File bamQCMetricsResult
	File histogram
	File metadata
	String outputFileNamePrefix
	String modules = "python/3.6"
	Int jobMemory = 8
	Int threads = 4
	Int timeout = 1
    }

    parameter_meta {
	bamQCMetricsResult: "JSON result file from bamQCMetrics"
	histogram: "JSON file with coverage histogram"
	metadata: "JSON file with additional metadata"
	outputFileNamePrefix: "Prefix for output file"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    String outputFileName = "~{outputFileNamePrefix}.bamQC_results.json"

    command <<<
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
    >>>

    output {
	File result = "~{outputFileName}"
    }

    meta {
	output_meta: {
            output1: "JSON file of collated results"
	}
    }
}

task countInputReads {

    input {
	File bamFile
	String modules = "samtools/1.9"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned data"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    command <<<
	samtools view -c ~{bamFile}
    >>>
    
    output {
	String result = read_string(stdout())
    }

    meta {
	output_meta: {
            output1: "Number of reads in input BAM file"
	}
    }
}

task cumulativeDistToHistogram {

    input {
	File globalDist
	File summary
	String modules = "python/3.6"
	Int jobMemory = 8
	Int threads = 4
	Int timeout = 1
    }

    parameter_meta {
	globalDist: "Global coverage distribution output from mosdepth"
	summary: "Summary output from mosdepth"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String outFileName = "coverage_histogram.json"

    # mosdepth writes a global coverage distribution with 3 columns:
    # 1) Chromsome name, or "total" for overall totals
    # 2) Depth of coverage
    # 3) Probability of coverage less than or equal to (2)
    # Want to convert the above cumulative probability distribution to a histogram
    # The "total" section of the summary discards some information
    # So, we process the outputs for each chromosome to construct the histogram

    command <<<
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
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File histogram = "~{outFileName}"
    }

    meta {
	output_meta: {
	    histogram: "Coverage histogram in JSON format"
	}
    }
}

task downsample {

    # random downsampling for QC metrics (excepting MarkDuplicates)

    input {
	File bamFile
	String outputFileNamePrefix
	Map[String, Boolean] downsampleStatus
	Map[String, String] downsampleTargets
	String downsampleSuffix = "downsampled.bam"
	Int randomSeed = 42
	String modules = "samtools/1.9"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	downsampleStatus: "Map; whether to apply pre-downsampling and downsampling"
	downsampleTargets: "Map; target number of reads for pre-downsampling and downsampling"
	downsampleSuffix: "Suffix for output file"
	randomSeed: "Random seed for pre-downsampling (if any)"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileNamePrefix}.~{downsampleSuffix}"

    # unpack downsample parameters
    Boolean applyPreDownsample = downsampleStatus["pre_ds"]
    String preDownsampleTarget = downsampleTargets["pre_ds"]
    String downsampleTarget = downsampleTargets["ds"]

    # generate downsample commands
    # preDownsample = fast, random selection of approximate total with samtools view
    String preDownsample = "samtools view -h -u -s ~{randomSeed}.~{preDownsampleTarget} | "
    String preDownsampleCommand = if applyPreDownsample then "~{preDownsample}" else ""
    # downsample = slow, deterministic selection of exact total with samtools collate and sort
    # see https://github.com/samtools/samtools/issues/931
    String dsCollate = "samtools collate -O --output-fmt sam - | "
    String dsAwk = "awk '/^@/ { print; next } count < ~{downsampleTarget} || last == $1 { print; last = $1; count++ }' | "
    String dsSort = "samtools sort -T downsample_sort - | "
    String downsampleCommand = "~{dsCollate}~{dsAwk}~{dsSort}"
    
    command <<<
	set -e
	set -o pipefail
	samtools view -b -h ~{bamFile} | \
	~{preDownsampleCommand} ~{downsampleCommand} \
	samtools view -b > ~{resultName}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{resultName}"
    }
    
    meta {
	output_meta: {
            result: "BAM file downsampled to required number of reads"
	}
    }

}

task downsampleRegion {

    # downsample a specific chromosomal region for MarkDuplicates
    # this keeps a proportionate level of duplicates in the downsampled data

    input {
	File bamFile
	File bamIndex
	String outputFileNamePrefix
	String region
	String modules = "samtools/1.9"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file"
	bamIndex: "BAM index file in BAI format"
	outputFileNamePrefix: "Prefix for output file"
	region: "Region argument for samtools"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String bamFileName = basename(bamFile)
    String resultName = "~{outputFileNamePrefix}.downsampledRegion.bam"

    # need to index the (filtered) BAM file before viewing a specific chromosome

    command <<<
	set -e
	# ensure BAM file and index are symlinked to working directory
	ln -s ~{bamFile}
	ln -s ~{bamIndex}
	samtools view -b -h ~{bamFileName} ~{region} > ~{resultName}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{resultName}"
    }

    meta {
	output_meta: {
            result: "BAM file downsampled to required number of reads"
	}
    }

}

task findDownsampleParams {

    input {
	String outputFileNamePrefix
	String inputReads
	Int targetReads = 100000
	Int minReadsAbsolute = 10000
	Int minReadsRelative = 2
	Int precision = 8
	Float preDSMultiplier = 1.5
	String modules = "python/3.6"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    String statusFile = "status.json"
    String targetsFile = "targets.json"

    parameter_meta {
	outputFileNamePrefix: "Prefix for output file"
	inputReads: "Number of reads in input bamFile (represented as string to avoid integer overflow)"
	targetReads: "Desired number of reads in downsampled output"
	minReadsAbsolute: "Minimum value of targetReads to allow pre-downsampling"
	minReadsRelative: "Minimum value of (inputReads)/(targetReads) to allow pre-downsampling"
	precision: "Number of decimal places in fraction for pre-downsampling"
	preDSMultiplier: "Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative)."
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    # see comments in "task downsample" for effect of predownsampling and downsampling

    # target for predownsampling with "samtools view -s" is expressed as a probability
    # eg. to choose approximately 200 reads out of 10000, target = 0.02
    # we convert to a fixed-precision target string for easier handling in BASH
    # eg. 0.02 -> "020000"
    # subsequently, we concatenate in the form {$RANDOM_SEED}.${TARGET}, eg. "42.020000"
    # for consistency, express downsampling target (integer number of reads) as a string also
    
    command <<<
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
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	Map[String, Boolean] status = read_json("~{statusFile}")
	Map[String, String] targets = read_json("~{targetsFile}")
    }

    meta {
	output_meta: {
            status: "Boolean flags indicating whether to apply (pre)downsampling.",
            output2: "Strings representing target number of reads for (pre)downsampling."
	}
    }
}

task findDownsampleParamsMarkDup {

    # downsampling parameters for MarkDuplicates; see filter_downsample.md for details
    # choose a region of the genome instead of using random selection

    # a BAM file is *very* approximately 10M reads per GB
    # Current merged BAM files are unlikely to exceed 10**9 reads; but we scale up higher just in case

    input {
	String outputFileNamePrefix
	String inputReads
	Int threshold = 10000000
	Array[String] chromosomes = ["chr12", "chr13", "chrXII", "chrXIII"]
	Int baseInterval = 15000
	Int intervalStart = 100000
	String customRegions = ""
	String modules = "python/3.6"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	outputFileNamePrefix: "Prefix for output file"
	inputReads: "Number of reads in input bamFile"
	threshold: "Minimum number of reads to conduct downsampling"
	chromosomes: "Array of chromosome identifiers for downsampled subset"
	baseInterval: "Base width of interval in each chromosome, for very large BAMs"
	intervalStart: "Start of interval in each chromosome, for very large BAMs"
	customRegions: "Custom downsample regions; overrides chromosome and interval parameters"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String outputStatus = "~{outputFileNamePrefix}_status.txt"
    String outputRegion = "~{outputFileNamePrefix}_region.txt"
    File chromosomesText = write_lines(chromosomes)

    command <<<
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
    >>>

    output {
	Boolean status = read_boolean("~{outputStatus}")
	String region = read_string("~{outputRegion}")
    }

    meta {
	output_meta: {
	    status: "Boolean flag, indicates whether downsampling is required",
	    region: "String to specify downsampled region for samtools"
	}
    }
}

task markDuplicates {

    input {
	File bamFile
	String outputFileNamePrefix
	Int opticalDuplicatePixelDistance=100
	Int picardMaxMemMb=6000
	String modules = "picard/2.21.2"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    # See GR-899 for opticalDuplicatePixelDistance

    parameter_meta {
	bamFile: "Input BAM file, after filtering and downsampling (if any)"
	outputFileNamePrefix: "Prefix for output file"
	opticalDuplicatePixelDistance: "Maximum offset between optical duplicate clusters"
	picardMaxMemMb: "Memory requirement in MB for running Picard JAR"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String outFileBam = "~{outputFileNamePrefix}.markDuplicates.bam"
    String outFileText = "~{outputFileNamePrefix}.markDuplicates.txt"

    command <<<
	java -Xmx~{picardMaxMemMb}M \
	-jar ${PICARD_ROOT}/picard.jar \
	MarkDuplicates \
	INPUT=~{bamFile} \
	OUTPUT=~{outFileBam} \
	VALIDATION_STRINGENCY=SILENT \
	TMP_DIR=${PWD} \
	METRICS_FILE=~{outFileText} \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance}
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	File result = "~{outFileText}"
    }

    meta {
	output_meta: {
            result: "Text file with Picard markDuplicates metrics"
	}
    }

}

task runMosdepth {

    input {
	File bamFile
	File bamIndex
	String modules = "mosdepth/0.2.9"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned data"
	bamIndex: "Index file in samtools .bai format"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    String bamFileName = basename(bamFile)

    command <<<
	set -eo pipefail
	# ensure BAM file and index are symlinked to working directory
	ln -s ~{bamFile}
	ln -s ~{bamIndex}
	# run mosdepth
	MOSDEPTH_PRECISION=8 mosdepth -x -n -t 3 bamqc ~{bamFileName}
    >>>

    output {
	File globalDist = "bamqc.mosdepth.global.dist.txt"
	File summary = "bamqc.mosdepth.summary.txt"
    }

    meta {
	output_meta: {
            globalDist: "Global distribution of coverage",
	    summary: "Total bases in coverage"
	}
  }

}

task updateMetadata {

    # add extra fields to the metadata JSON file

    input {
	Map[String, String] metadata
	String outputFileNamePrefix
	String totalInputReads
	String nonPrimaryReads
	String unmappedReads
	String lowQualityReads
	String modules = "python/3.6"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	metadata: "Key/value map of input metadata"
	outputFileNamePrefix: "Prefix for output file"
	totalInputReads: "Total reads in original input BAM file"
	nonPrimaryReads: "Total reads excluded as non-primary"
	unmappedReads: "Total reads excluded as unmapped"
	lowQualityReads: "Total reads excluded as low alignment quality"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    File metadataJson = write_json(metadata)
    String outFileName = "~{outputFileNamePrefix}.updated_metadata.json"

    # Read totals are Strings in WDL to avoid integer overflow in Cromwell
    # Python3 can handle arbitrarily large integers

    command <<<
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
    >>>

    output {
	File result = "~{outFileName}"
    }

    meta {
	output_meta: {
            result: "JSON file with updated metadata"
	}
    }
}