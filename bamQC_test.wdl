version 1.0

struct InputGroup {
  File bam
  File bamIndex
}


workflow bamQC {

  input {
    Array[InputGroup] inputGroups
    String teststring = "bamQC"
  }

  parameter_meta {
    inputGroups: "Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name)."
    doFilter: "Enable/disable Samtools filtering."
  }

	scatter (i in inputGroups) {

		call preFilterBam {
			input:
			inputBam = i.bam,
			inputBamIndex = i.bamIndex,
      outputFileName = teststring
		}
	}

	Array[File] preFilteredBams = preFilterBam.preFilteredBam
  Array[File] preFilteredBamIndexes = preFilterBam.preFilteredBamIndex

	call mergeBams {
			input:
			bams = preFilteredBams,
      outputFileName = teststring
		}

	output {
		String outputFileName = teststring
	}

  meta {
	description: "test"
  }

}

task preFilterBam {
  input {
    File inputBam
    File inputBamIndex
    Int filterFlags = 260
    Int? minMapQuality
    String? filterAdditionalParams
    String outputFileName
    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
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
      ~{inputBam} > $outputBam
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
