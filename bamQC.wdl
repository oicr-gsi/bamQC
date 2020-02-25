version 1.0

workflow bamQC {

    input {
	File bamFile
	String outputFileNamePrefix = "bamQC"
    }

    parameter_meta {
	bamFile: "Input BAM file on which to compute QC metrics"
	outputFileNamePrefix: "Prefix for output files"
    }

    call countInputReads {
	input:
	bamFile = bamFile,
    }

    call findDownsampleParams {
	input:
	bamFile = bamFile,
	outputFileNamePrefix = outputFileNamePrefix,
	inputReads = countInputReads.result
    }
    
    call downsample {
	input:
	bamFile = bamFile,
	outputFileNamePrefix = outputFileNamePrefix,
	downsampleStatus = findDownsampleParams.status,
	downsampleTargets = findDownsampleParams.targets,
    }

    output {
	File result = downsample.result
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
	    name: "bam-qc-metrics/0.2.4",
	    url: "https://github.com/oicr-gsi/bam-qc-metrics.git"
	}
	]
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


    command <<<
	samtools view -c ~{bamFile}
    >>>
    
    output {
	Int result = read_int(stdout())
    }
}

task downsample {

    input {
	File bamFile
	String outputFileNamePrefix
	Map[String, Boolean] downsampleStatus
	Map[String, Int] downsampleTargets
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

    # unpack downsampleParams
    Boolean applyPreDownsample = downsampleStatus["pre_ds"]
    Boolean applyDownsample = downsampleStatus["ds"]
    Int preDownsampleTarget = downsampleTargets["pre_ds"]
    Int downsampleTarget = downsampleTargets["ds"]
    
    command <<<
	echo "Pre-downsampling status ~{applyPreDownsample}" | tee ~{resultName}
	echo "Downsampling status ~{applyDownsample}" | tee -a ~{resultName}
	echo "Pre-downsampling target ~{preDownsampleTarget}" | tee -a ~{resultName}
	echo "Downsampling target ~{downsampleTarget}" | tee -a ~{resultName}
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
	File bamFile
	String outputFileNamePrefix
	Int inputReads
	Int targetReads = 1000
    }

    String statusFile = "status.json"
    String targetsFile = "targets.json"
    
    command <<<
        python3 <<CODE
        import json, math, sys
        readsIn = ~{inputReads}
        readsTarget = ~{targetReads}
        print("Input reads param =", readsIn, file=sys.stderr)
        print("Target reads param =", readsTarget, file=sys.stderr)
        minReadsAbsolute = 10000
        minReadsRelative = 2
        preDownsampleMultiplier = 1.5
        if readsIn <= readsTarget:
          # absolutely no downsampling
          applyPreDownsample = False
          applyDownsample = False
          preDownsampleTarget = 0
          downSampleTarget = 0
        elif readsIn < readsTarget * minReadsRelative or readsTarget < minReadsAbsolute:
          # no predownsampling
          applyPreDownsample = False
          applyDownsample = True
          preDownSampleTarget = 0
          downSampleTarget = readsTarget
        else:
          # predownsampling and downsampling
          applyPreDownsample = True
          applyDownsample = True
          preDownSampleTarget = int(math.floor(readsTarget * preDownsampleMultiplier))
          downSampleTarget = readsTarget
        status = {
          "pre_ds": applyPreDownsample,
          "ds": applyDownsample
        }
        targets = {
          "pre_ds": preDownSampleTarget,
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

    output {
	Map[String, Boolean] status = read_json("~{statusFile}")
	Map[String, Int] targets = read_json("~{targetsFile}")
    }
}