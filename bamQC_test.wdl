version 1.0

struct BamAndBamIndex {
  File bam
  File bamIndex
}

struct InputGroup {
  String outputIdentifier
  Array[BamAndBamIndex]+ bamAndBamIndexInputs
}

struct CollectionGroup {
    String outputIdentifier
    String outputFileName
    Array[File] bams
    Array[File] bamIndexes
}

struct InputGroups {
  Array[InputGroup] inputGroups
}

struct CollectionGroups {
  Array[CollectionGroup] collectionGroups
}

struct OutputGroup {
  String outputIdentifier
  File bam
  File bamIndex
}

workflow bamQC {

  input {
    Array[InputGroup] inputGroups
    Boolean doFilter = true
    String outputFileNamePrefix = "bamQC"
  }

  parameter_meta {
    inputGroups: "Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name)."
    doFilter: "Enable/disable Samtools filtering."
    outputFileNamePrefix: "Prefix for output files"
  }

	scatter (i in inputGroups) {
		scatter(bamAndBamIndexInput in i.bamAndBamIndexInputs) {
			File inputGroupBam = bamAndBamIndexInput.bam
			File inputGroupBamIndex = bamAndBamIndexInput.bamIndex
		}
		Array[File] inputGroupBams = inputGroupBam
		Array[File] inputGroupBamIndexes = inputGroupBamIndex

		call preprocessBam {
			input:
			bams = inputGroupBams,
			bamIndexes = inputGroupBamIndexes,
			outputFileName = i.outputIdentifier,
			doFilter = doFilter,
		}
	}

	Array[File] preprocessedBams = preprocessBam.preprocessedBam
    Array[File] preprocessedBamIndexes = preprocessBam.preprocessedBamIndex

	call collectFilesBySample {
    input:
      inputGroups = inputGroups,
      bams = preprocessedBams,
      bamIndexes = preprocessedBamIndexes
	}

	scatter(o in collectFilesBySample.filesByOutputIdentifier.collectionGroups) {
		if(length(o.bams) > 1) {
		call mergeBams as mergeSplitByIntervalBams {
			input:
			bams = o.bams,
			outputFileName = o.outputFileName,
			suffix = "" # collectFilesBySample task generates the file name
		}
		}
		String outputprefix = o.outputIdentifier
	}

	output {
		Array[String] outputPrefix = outputprefix
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
	},
	    {
	    name: "mosdepth/0.2.9",
	    url: "https://github.com/brentp/mosdepth"
	}
	]
    }

}

task preprocessBam {
  input {
    Boolean doFilter = true
    String outputFileName

    # by default write tmp files to the current working directory (cromwell task directory)
    # $TMPDIR is set by Cromwell
    # $TMP is set by Univa
    String temporaryWorkingDir = ""

    Array[File] bams
    Array[File] bamIndexes

    # filter parameters
    String filterSuffix = ".filter"
    Int filterFlags = 260
    Int? minMapQuality
    String? filterAdditionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "samtools/1.9 gatk/4.1.6.0"
  }

  String workingDir = if temporaryWorkingDir == "" then "" else "~{temporaryWorkingDir}/"

  String baseFileName = "~{outputFileName}"

  String filteredFileName = if doFilter then
                            "~{baseFileName}.filter"
                           else
                            "~{baseFileName}"
  String filteredFilePath = "~{filteredFileName}"


  command <<<
    set -euxo pipefail
    inputBams="~{sep=" " bams}"
    inputBamIndexes="~{sep=" " bamIndexes}"

    # filter
    if [ "~{doFilter}" = true ]; then
      outputBams=()
      outputBamIndexes=()
      for inputBam in $inputBams; do
        filename="$(basename $inputBam ".bam")"
        outputBam="~{workingDir}${filename}.filtered.bam"
        outputBamIndex="~{workingDir}${filename}.filtered.bai"
        samtools view -b \
        -F ~{filterFlags} \
        ~{"-q " + minMapQuality} \
        ~{filterAdditionalParams} \
        $inputBam > $outputBam
        samtools index $outputBam $outputBamIndex
        outputBams+=("$outputBam")
        outputBamIndexes+=("$outputBamIndex")
      done
      # set inputs for next step
      inputBams=("${outputBams[@]}")
      inputBamIndexes=("${outputBamIndexes[@]}")
    else
      outputBams=()
      outputBamIndexes=()
      for inputBam in $inputBams; do
        filename="$(basename $inputBam ".bam")"
        outputBam="~{workingDir}${filename}.bam"
        outputBamIndex="~{workingDir}${filename}.bai"
        samtools view -b \
        $inputBam > $outputBam
        samtools index $outputBam $outputBamIndex
        outputBams+=("$outputBam")
        outputBamIndexes+=("$outputBamIndex")
      done
      # set inputs for next step
      inputBams=("${outputBams[@]}")
      inputBamIndexes=("${outputBamIndexes[@]}")
    fi


      gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
      ${inputBams[@]/#/--INPUT=} \
      --OUTPUT="~{filteredFileName}.bam" \
      --CREATE_INDEX=true \
      --SORT_ORDER=coordinate \
      --ASSUME_SORTED=false \
      --USE_THREADING=true \
      --VALIDATION_STRINGENCY=SILENT
  >>>

  output {
    File preprocessedBam =  if doFilter then
                            "~{filteredFilePath}.bam"
                           else "~{filteredFileName}.bam"
    File preprocessedBamIndex = if doFilter then
                                  "~{filteredFilePath}.bai"
                                else "~{filteredFileName}.bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    doFilter: "Enable/disable Samtools filtering."
    outputFileName: "Output files will be prefixed with this."
    temporaryWorkingDir: "Where to write out intermediary bam files. Only the final preprocessed bam will be written to task working directory if this is set to local tmp."
    bams: "Array of bam files to merge together."
    bamIndexes: "Array of index files for input bams."
    filterSuffix: "Suffix to use for filtered bams."
    filterFlags: "Samtools filter flags to apply."
    minMapQuality: "Samtools minimum mapping quality filter to apply."
    filterAdditionalParams: "Additional parameters to pass to samtools."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task collectFilesBySample {
  input {
    Array[InputGroup] inputGroups
    Array[File] bams
    Array[File] bamIndexes

    Int jobMemory = 1
    Int cores = 1
    Int timeout = 1
    String modules = "python/3.7"
  }

  InputGroups wrappedInputGroups = {"inputGroups": inputGroups}

  command <<<
    set -euo pipefail

    python3 <<CODE
    import json
    import os
    import re

    with open('~{write_json(wrappedInputGroups)}') as f:
        inputGroups = json.load(f)
    with open('~{write_lines(bams)}') as f:
        bamFiles = f.read().splitlines()
    with open('~{write_lines(bamIndexes)}') as f:
        bamIndexFiles = f.read().splitlines()

    filesByOutputIdentifier = []
    for outputIdentifier in [inputGroup['outputIdentifier'] for inputGroup in inputGroups['inputGroups']]:
        # select bams and bamIndexes for outputIdentifier (preprocessBam prefixes the outputIdentifier, so include that too)
        bams = [bam for bam in bamFiles if re.match("^" + outputIdentifier + "\.", os.path.basename(bam))]
        bais = [bai for bai in bamIndexFiles if re.match("^" + outputIdentifier + "\.", os.path.basename(bai))]

        fileNames = list(set([os.path.splitext(os.path.basename(f))[0] for f in bams + bais]))
        if len(fileNames) != 1:
            raise Exception("Unable to determine unique fileName from fileNames = [" + ','.join(f for f in fileNames) + "]")
        else:
            fileName = fileNames[0]

        filesByOutputIdentifier.append({
            'outputIdentifier': outputIdentifier,
            'outputFileName': fileName,
            'bams': bams,
            'bamIndexes': bais})

    # wrap the array into collectionGroups object
    wrappedFilesByOutputIdentifier = {'collectionGroups': filesByOutputIdentifier}

    with open('filesByOutputIdentifier.json', 'w') as f:
        json.dump(wrappedFilesByOutputIdentifier, f, indent=4)
    CODE
  >>>

  output {
    CollectionGroups filesByOutputIdentifier = read_json("filesByOutputIdentifier.json")
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    inputGroups: "Array of objects describing output file groups. The output file group name is used to partition input bams by name."
    bams: "Array of bams to partition by inputGroup output file name."
    bamIndexes: "Array of index files for input bams."
    jobMemory:  "Memory allocated to job (in GB)."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task mergeBams {
  input {
    Array[File] bams
    String outputFileName
    String suffix = ".merge"
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
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
    --VALIDATION_STRINGENCY=SILENT \
    ~{additionalParams}
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

  parameter_meta {
    bams: "Array of bam files to merge together."
    outputFileName: "Output files will be prefixed with this."
    additionalParams: "Additional parameters to pass to GATK MergeSamFiles."
    jobMemory: "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

