version 1.0

struct InputGroup {
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

		call preFilterBam {
			input:
			inputBam = i.bam,
			inputBamIndex = i.bamIndex,
			outputFileName = i.outputIdentifier,
			doFilter = doFilter,
		}
	}

	Array[File] preFilteredBams = preFilterBam.preFilteredBam
  Array[File] preFilteredBamIndexes = preFilterBam.preFilteredBamIndex
  String outputFileName = inputGroups[0].outputIdentifier

	call mergeBams {
			input:
			bams = preFilteredBams,
			outputFileName = outputFileName,
			suffix = "" 
		}

	output {
		String outputName = outputFileName
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

task preFilterBam {
  input {
    Boolean doFilter = true
    String outputFileName

    # by default write tmp files to the current working directory (cromwell task directory)
    # $TMPDIR is set by Cromwell
    # $TMP is set by Univa
    String temporaryWorkingDir = ""

    File inputBam
    File inputBamIndex

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
  String filteredFileName = "~{outputFileName}.filtered"


  command <<<
    set -euxo pipefail

    # filtercd
    if [ "~{doFilter}" = true ]; then
      filename="$(basename ~{inputBam} ".bam")"
      outputBam="~{workingDir}${filename}.filtered.bam"
      outputBamIndex="~{workingDir}${filename}.filtered.bai"
      samtools view -b \
      -F ~{filterFlags} \
      ~{"-q " + minMapQuality} \
      ~{filterAdditionalParams} \
      ~{inputBam} > $outputBam
      samtools index $outputBam $outputBamIndex
    else
      filename="$(basename ~{inputBam} ".bam")"
      outputBam="~{workingDir}${filename}.bam"
      outputBamIndex="~{workingDir}${filename}.bai"
      ln -s ~{inputBam} $outputBam
      ln -s ~{inputBamIndex} $outputBamIndex
    fi
  >>>

  output {
    File preFilteredBam =  "~{filteredFileName}.bam"
    File preFilteredBamIndex = "~{filteredFileName}.bai"
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
    temporaryWorkingDir: "Where to write out intermediary bam files. Only the final preFiltered bam will be written to task working directory if this is set to local tmp."
    inputBam: "bam files to filter."
    inputBamIndex: "index file for input bam."
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
