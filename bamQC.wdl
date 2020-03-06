version 1.0

workflow bamQC {

    input {
	File bamFile
	Map[String, String] metadata
	String outputFileNamePrefix = "bamQC"
    }

    parameter_meta {
	bamFile: "Input BAM file on which to compute QC metrics"
	metadata: "JSON file containing metadata"
	outputFileNamePrefix: "Prefix for output files"
    }

    call filter {
	input:
	bamFile = bamFile,
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

    call countInputReads {
	input:
	bamFile = filter.filteredBam,
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
	    bamFile = filter.filteredBam,
	    outputFileNamePrefix = outputFileNamePrefix,
	    downsampleStatus = findDownsampleParams.status,
	    downsampleTargets = findDownsampleParams.targets,
	}
    }

    if (dsMarkDup) {
	call downsampleChromosome {
	    input:
	    bamFile = filter.filteredBam,
	    outputFileNamePrefix = outputFileNamePrefix,
	    region = findDownsampleParamsMarkDup.region
	}
    }

    call markDuplicates {
	input:
	bamFile = filter.filteredBam,
	outputFileNamePrefix = outputFileNamePrefix,
	downsampled = dsMarkDup,
	bamFileDownsampled = downsampleChromosome.result
    }

    call bamQCMetrics {
	input:
	bamFile = filter.filteredBam,
	outputFileNamePrefix = outputFileNamePrefix,
	markDuplicates = markDuplicates.result,
	metadata = updateMetadata.result,
	downsampled = ds,
	bamFileDownsampled = downsample.result
    }

    output {
	File result = bamQCMetrics.result
    }

    meta {
	author: "Iain Bancarz"
	email: "ibancarz@oicr.on.ca"
	description: "QC metrics for BAM files"
	dependencies: [
	{
	    name: "samtools/1.9",
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
	}
	]
    }

}

task bamQCMetrics {

    input {
	File bamFile
	String outputFileNamePrefix
	File markDuplicates
	File metadata
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
	metadata: "JSON file containing metadata (including filtered read totals)"
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
	-m ~{metadata} \
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
            output1: "Placeholder text file for BAMQC output; TODO replace with JSON"
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
	bamFile: "Input BAM file of aligned rnaSeqQC data"
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
	Int result = read_int(stdout())
    }

    meta {
	output_meta: {
            output1: "Number of reads in input BAM file"
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

task downsampleChromosome {

    # downsample a specific chromosomal region for MarkDuplicates
    # this keeps a proportionate level of duplicates in the downsampled data

    input {
	File bamFile
	String outputFileNamePrefix
	String region
	String modules = "samtools/1.9"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    String resultName = "~{outputFileNamePrefix}.downsampledChromosome.bam"

    # need to index the (filtered) BAM file before viewing a specific chromosome

    command <<<
	set -e
	samtools index ~{bamFile}
	samtools view -b -h ~{bamFile} ~{region} > ~{resultName}
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
	String modules = "samtools/1.9"
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
    >>>

    runtime {
	modules: "~{modules}"
	memory:  "~{jobMemory} GB"
	cpu:     "~{threads}"
	timeout: "~{timeout}"
    }

    output {
	Int totalInputReads = read_int("~{totalInputReadsFile}")
	Int nonPrimaryReads = read_int("~{totalNonPrimaryReadsFile}")
	Int unmappedReads = read_int("~{totalUnmappedReadsFile}")
	Int lowQualityReads = read_int("~{totalLowQualityReadsFile}")
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

task findDownsampleParams {

    input {
	String outputFileNamePrefix
	Int inputReads
	Int targetReads = 100000
	Int minReadsAbsolute = 10000
	Int minReadsRelative = 2
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
	inputReads: "Number of reads in input bamFile"
	targetReads: "Desired number of reads in downsampled output"
	minReadsAbsolute: "Minimum value of targetReads to allow pre-downsampling"
	minReadsRelative: "Minimum value of (inputReads)/(targetReads) to allow pre-downsampling"
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
          precision = 6 # number of decimal places to keep
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
        print(json.dumps(status), file=statusFile)
        statusFile.close()
        targetFile = open("~{targetsFile}", "w")
        print(json.dumps(targets), file=targetFile)
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

    # downsampling parameters for MarkDuplicates
    # choose a region of the genome instead of using random selection

    # number of reads downsampled is not exact; will be:
    # - approximately 1/12 for between 10**6 and 10**8 reads
    # - approximately 1/300 for more than 10**8 reads
    # - TODO could make this smoother and lessen the "cliff edge" at 10**8

    input {
	String outputFileNamePrefix
	Int inputReads
	String chromosome = "chr1"
	Int width = 10000000
	String modules = "python/3.6"
	Int jobMemory = 16
	Int threads = 4
	Int timeout = 4
    }

    parameter_meta {
	outputFileNamePrefix: "Prefix for output file"
	inputReads: "Number of reads in input bamFile"
	chromosome: "Chromosome identifier for downsampled subset"
	width: "Width of range within chromosome (if required)"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String outputStatus = "~{outputFileNamePrefix}_status.txt"
    String outputRegion = "~{outputFileNamePrefix}_region.txt"

    command <<<
        python3 <<CODE
        readsIn = ~{inputReads}
        chromosome = "~{chromosome}"
        if readsIn <= 10**6:
            status = "false"
            region = ""
        elif readsIn <= 10**8:
            status = "true"
            region = chromosome
        else:
            status = "true"
            start = 10**6 + 1
            end = start + ~{width} - 1
            region = "%s:%i-%i" % (chromosome, start, end)
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
	Boolean downsampled
	File? bamFileDownsampled
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
	downsampled: "True if downsampling has been applied"
	opticalDuplicatePixelDistance: "Maximum offset between optical duplicate clusters"
	picardMaxMemMb: "Memory requirement in MB for running Picard JAR"
	bamFileDownsampled: "(Optional) downsampled subset of reads from bamFile"
	modules: "required environment modules"
	jobMemory: "Memory allocated for this job"
	threads: "Requested CPU threads"
	timeout: "hours before task timeout"
    }

    String outFileBam = "~{outputFileNamePrefix}.markDuplicates.bam"
    String outFileText = "~{outputFileNamePrefix}.markDuplicates.txt"

    File inputBam = if downsampled then bamFileDownsampled else bamFile

    command <<<
	java -Xmx~{picardMaxMemMb}M \
	-jar ${PICARD_ROOT}/picard.jar \
	MarkDuplicates \
	INPUT=~{inputBam} \
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

task updateMetadata {

    # add extra fields to the metadata JSON file

    input {
	Map[String, String] metadata
	String outputFileNamePrefix
	Int totalInputReads
	Int nonPrimaryReads
	Int unmappedReads
	Int lowQualityReads
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

    command <<<
        python3 <<CODE
        import json
        metadata = json.loads(open("~{metadataJson}").read())
        metadata["total input reads meta"] = ~{totalInputReads}
        metadata["non-primary reads meta"] = ~{nonPrimaryReads}
        metadata["unmapped reads meta"] = ~{unmappedReads}
        metadata["low-quality reads meta"] = ~{lowQualityReads}
        outFile = open("~{outFileName}", "w")
        print(json.dumps(metadata), file=outFile)
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
