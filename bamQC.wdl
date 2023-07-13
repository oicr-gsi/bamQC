version 1.0

struct InputGroup {
  File bam
  File bamIndex
}


workflow bamQC {

  input {
    Array[InputGroup] inputGroups
    String outputFileNamePrefix = "bamQC"
    String intervalsToParallelizeByString = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
  }

  parameter_meta {
    inputGroups: "Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name)."
    doFilter: "Enable/disable Samtools filtering."
  }

  call splitStringToArray {
    input:
      str = intervalsToParallelizeByString
  }
  Array[Array[String]] intervalsToParallelizeBy = splitStringToArray.out

	scatter (intervals in intervalsToParallelizeBy) {
    scatter (i in inputGroups) {
      call preFilterBam {
        input:
        inputBam = i.bam,
        inputBamIndex = i.bamIndex,
        intervals = intervals,
        outputFileName = outputFileNamePrefix
      }
    }

	Array[File] preFilteredBams = preFilterBam.preFilteredBam
  Array[File] preFilteredBamIndexes = preFilterBam.preFilteredBamIndex

  call mergeBams as mergeSplitByIntervalBams {
    input:
      bams = preFilteredBams,
      outputFileName = outputFileNamePrefix
  }
  }

  Array[File] processedBams = mergeSplitByIntervalBams.mergedBam
  Array[File] ProcessedBamIndexes = mergeSplitByIntervalBams.mergedBamIndex

	call mergeBams {
			input:
			bams = processedBams,
      outputFileName = outputFileNamePrefix
		}

	output {
		String outputFileName = outputFileNamePrefix
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

task preFilterBam {
  input {
    File inputBam
    File inputBamIndex
    Array[String] intervals
    Int filterFlags = 260
    Int? minMapQuality
    String? filterAdditionalParams
    String outputFileName
    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 24
    String modules = "samtools/1.9 gatk/4.1.6.0"
  }
  String filteredFileName = "~{outputFileName}.filtered"
  command <<<
    set -euxo pipefail
      outputBam="~{filteredFileName}.bam"
      outputBamIndex="~{filteredFileName}.bai"
      samtools view -b \
      -F ~{filterFlags} \
      ~{"-q " + minMapQuality} \
      ~{filterAdditionalParams} \
      ~{inputBam} ~{sep=" " intervals} > $outputBam
      samtools index $outputBam $outputBamIndex
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
}


task mergeBams {
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

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
    ~{sep=" " prefix("--INPUT=", bams)} \
    --OUTPUT="~{outputFileName}~{suffix}.bam" \
    --CREATE_INDEX=true \
    --SORT_ORDER=coordinate \
    --ASSUME_SORTED=true \
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
