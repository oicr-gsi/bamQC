version 1.0

struct InputGroup {
  File bam
  File bamIndex
}

workflow bamQC {

  input {
    Array[InputGroup] inputGroups
    Map[String, String] metadata
    String mode
    String outputFileNamePrefix = "bamQC"
    String intervalsToParallelizeByString = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
  }

  parameter_meta {
	inputGroups: "Array of objects describing sets of bams to merge together and on which to compute QC metrics"
	metadata: "JSON file containing metadata"
    mode: "running mode for the workflow, only allow value 'lane_level' and 'call_ready'"
	outputFileNamePrefix: "Prefix for output files"
    intervalsToParallelizeByString: "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3,chr4)."
    }

  if (( mode == "lane_level") && (length(inputGroups) ==1 )) {
    File laneLevelBam = inputGroups[0].bam
    File laneLevelBamIndex = inputGroups[0].bamIndex
  }

  if (mode == "call_ready") {
    call splitStringToArray {
        input:
            str = intervalsToParallelizeByString
        }
    
    Array[String] intervalsToParallelizeBy = flatten(splitStringToArray.out)
    scatter (i in inputGroups) {
        scatter (interval in intervalsToParallelizeBy) {
        call getChrCoefficient {
            input:
            bamFile = i.bam,
            chromosome = interval
        } 
        call preFilter {
            input:
            bamFile = i.bam,
            bamIndex = i.bamIndex,
            interval = interval,
            scaleCoefficient = getChrCoefficient.coeff,
            outputFileName = outputFileNamePrefix
        }
        }

        Array[File] filteredBams = preFilter.filteredBam

        call mergeFiles as mergeSplitByIntervalFiles {
            input:
            bams = filteredBams,
            outputFileName = outputFileNamePrefix
        }
    }
    Array[File] processedBams = mergeSplitByIntervalFiles.mergedBam

    if (length(processedBams) > 1) {
        call mergeFiles {
                input:
                bams = processedBams,
                outputFileName = outputFileNamePrefix
            }
    } 
    File mergedBam = select_first([mergeFiles.mergedBam, mergeSplitByIntervalFiles.mergedBam[0]])
    File mergedBamIndex = select_first([mergeFiles.mergedBamIndex, mergeSplitByIntervalFiles.mergedBamIndex[0]])
  }	

  File qcReadyBam = select_first([laneLevelBam, mergedBam])
  File qcReadyBamIndex = select_first([laneLevelBamIndex, mergedBamIndex])

  call filter {
    input:
    bamFile = qcReadyBam,
    outputFileNamePrefix = outputFileNamePrefix
  }

  call updateMetadata {
	input:
	metadata = metadata,
	outputFileNamePrefix = outputFileNamePrefix,
	totalInputReads = filter.totalInputReads,
	nonPrimaryReads = filter.nonPrimaryReads,
	unmappedReads = filter.unmappedReads,
	lowQualityReads = filter.lowQualityReads
    }

    call findDownsampleParams {
	input:
	outputFileNamePrefix = outputFileNamePrefix,
	inputReads = filter.totalInputReads
    }

    call findDownsampleParamsMarkDup {
	input:
	outputFileNamePrefix = outputFileNamePrefix,
	inputReads = filter.totalInputReads
    }

    Boolean ds = findDownsampleParams.status["ds"]
    Boolean dsMarkDup = findDownsampleParamsMarkDup.status

    if (ds) {
	call downsample {
	    input:
	    bamFile = filter.filteredBam,
	    outputFileNamePrefix = outputFileNamePrefix,
	    downsampleStatus = findDownsampleParams.status,
	    downsampleTargets = findDownsampleParams.targets,
	}
    }

    if (dsMarkDup) {
	call downsampleRegion {
	    input:
	    bamFile = filter.filteredBam,
        bamIndex = filter.filteredBai,
	    outputFileNamePrefix = outputFileNamePrefix,
	    region = findDownsampleParamsMarkDup.region
	}
    }

    Array[File?] markDupInputs = [downsampleRegion.result, filter.filteredBam]
    call markDuplicates {
	input:
	bamFile = select_first(markDupInputs),
	outputFileNamePrefix = outputFileNamePrefix
    }

    call bamQCMetrics {
	input:
	bamFile = filter.filteredBam,
	outputFileNamePrefix = outputFileNamePrefix,
	markDuplicates = markDuplicates.result,
	downsampled = ds,
	bamFileDownsampled = downsample.result
    }

    call runMosdepth {
	input:
	bamFile = filter.filteredBam,
	bamIndex = filter.filteredBai
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
        author: "Iain Bancarz"
        email: "ibancarz@oicr.on.ca"
        description: "bamQC workflow collects a number of metrics which are computed using several methods (by employing third-party software tools along with some custom code) and outputs the results in JSON format. The output also contains metadata, such as the instrument and lane names. Also, it contains the estimate of ribosomal rRNA contamination and other information. bamQC supports downsampling for faster analysis and splits some tasks by chromosome, which also increases speed."
        dependencies: [
            {
                name: "samtools/1.14",
                url: "https://github.com/samtools/samtools"
            },
            {
                name: "picard/2.21.2",
                url: "https://broadinstitute.github.io/picard/command-line-overview.html"
            },
            {
                name: "python/3.6",
                url: "https://www.python.org/downloads/"
            },
            {
                name: "bam-qc-metrics/0.2.5",
                url: "https://github.com/oicr-gsi/bam-qc-metrics.git"
            },
                {
                name: "mosdepth/0.2.9",
                url: "https://github.com/brentp/mosdepth"
            }
        ]
    output_meta: {
    result: {
        description: "json file that contains metrics and meta data described in https://github.com/oicr-gsi/bam-qc-metrics/blob/master/metrics.md",
        vidarr_label: "result"
    }
}
    }

}

# ================================================================
#  Scaling coefficient - use to scale RAM allocation by chromosome
# ================================================================
task getChrCoefficient {
  input {
    Int memory = 2
    Int timeout = 1
    String chromosome
    String modules = "samtools/1.14"
    File bamFile
  }

  parameter_meta {
    bamFile: ".bam file to process, we just need the header"
    timeout: "Hours before task timeout"
    chromosome: "Chromosome to check"
    memory: "Memory allocated for this job"
    modules: "Names and versions of modules to load"
  }

  command <<<
    CHROM_LEN=$(samtools view -H ~{bamFile} | grep ^@SQ | grep -v _ | grep -w ~{chromosome} | cut -f 3 | sed 's/LN://')
    LARGEST=$(samtools view -H ~{bamFile} | grep ^@SQ | grep -v _ | cut -f 3 | sed 's/LN://' | sort -n | tail -n 1)
    echo | awk -v chrom_len=$CHROM_LEN -v largest=$LARGEST '{print int((chrom_len/largest + 0.1) * 10)/10}'
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    String coeff = read_string(stdout())
  }

  meta {
    output_meta: {
      coeff: "Length ratio as relative to the largest chromosome."
    }
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
	String modules = "samtools/1.14"
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
	String modules = "samtools/1.14"
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

task filter {

    # filter out non-primary, unmapped, and low-quality aligned reads
    # count the number of reads filtered out at each step
    # return filtered read counts and the filtered BAM file

   input {
	File bamFile
	String outputFileNamePrefix
	Int minQuality = 30
	String modules = "samtools/1.14"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
	outputFileNamePrefix: "Prefix for output file"
	minQuality: "Minimum alignment quality to pass filter"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String resultName = "~{outputFileNamePrefix}.filtered.bam"
    String resultIndexName = "~{outputFileNamePrefix}.filtered.bam.bai"
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
    samtools view -h -b -F 2304 -U ~{nonPrimaryReadsFile} ~{bamFile} | \
    samtools view -h -b -F 4 -U ~{unmappedReadsFile} | \
    samtools view -h -b -q ~{minQuality} -U ~{lowQualityReadsFile} \
    > ~{resultName}
    samtools view -c ~{bamFile} > ~{totalInputReadsFile}
    samtools view -c ~{nonPrimaryReadsFile} > ~{totalNonPrimaryReadsFile}
    samtools view -c ~{unmappedReadsFile} > ~{totalUnmappedReadsFile}
    samtools view -c ~{lowQualityReadsFile} > ~{totalLowQualityReadsFile}
    samtools index ~{resultName} ~{resultIndexName}
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
    File filteredBai = "~{resultIndexName}"
    }

    meta {
	output_meta: {
    totalInputReads: "Total reads in original input BAM file",
    nonPrimaryReads: "Total reads excluded as non-primary",
    unmappedReads: "Total reads excluded as unmapped",
    lowQualityReads: "Total reads excluded as low alignment quality",
    filteredBam: "Filtered BAM file",
    filteredBai: "Filtered bam index"
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

task mergeFiles {
  input {
    Array[File] bams
    String outputFileName
    String suffix = ".merge"
    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 24
    String modules = "gatk/4.1.6.0"
  }

  parameter_meta {
    bams: "Array of bam files to merge together."
    outputFileName: "Output files will be prefixed with this."
    suffix: "suffix to use for merged bam"
    jobMemory: "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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

  >>>

  output {
    File mergedBam = "~{outputFileName}~{suffix}.bam"
    File mergedBamIndex = "~{outputFileName}~{suffix}.bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

}

task preFilter {
    input {
    File bamFile
    File bamIndex
    String interval
    String outputFileName
    Int filterFlags = 260
    Int? minMapQuality
    Float scaleCoefficient = 1.0
    String? filterAdditionalParams
    String modules = "samtools/1.14"
    Int jobMemory = 6
    Int minMemory = 2
    Int threads = 4
    Int timeout = 4
    }
    parameter_meta {
	bamFile: "Input BAM file of aligned rnaSeqQC data"
	outputFileName: "Prefix for output file"
	minMapQuality: "Minimum alignment quality to pass filter"
        filterFlags: "Samtools filter flags to apply."
        filterAdditionalParams: "Additional parameters to pass to samtools."
        interval: "Usually, a chromosome id as it would appear in input .bam header"
	modules: "required environment modules"
        scaleCoefficient: "Chromosome-dependent RAM scaling coefficient"
	jobMemory: "Memory allocated for this job"
        minMemory: "Minimum amount of RAM allocated to the task"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

   String resultName = "~{outputFileName}.filtered.bam"
   Int allocatedMemory = if minMemory > round(jobMemory * scaleCoefficient) then minMemory else round(jobMemory * scaleCoefficient)

   command <<<
   set -e
   set -o pipefail
   samtools view -b \
        -F ~{filterFlags} \
        ~{"-q " + minMapQuality} \
        ~{filterAdditionalParams} \
        ~{bamFile} \
        ~{interval} > ~{resultName}
    >>> 

    output {
    File filteredBam = "~{resultName}"
    }

    meta {
	output_meta: {
        filteredBam: "Filtered BAM file"
	}
    }

    runtime {
	modules: "~{modules}"
	memory:  "~{allocatedMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
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
    recordSeparator: "Record separator - a delimiter for joining records"
    jobMemory: "Memory allocated to job (in GB)."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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
