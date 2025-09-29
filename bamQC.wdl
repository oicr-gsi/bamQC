version 1.0

struct InputGroup {
  File bam
  File bamIndex
}

struct Resources {
  String refFasta
  String modules
}

workflow bamQC {

  input {
    Array[InputGroup] inputGroups
    Map[String, String] metadata
    String? targetBed
    String reference
    String outputFileNamePrefix = "bamQC"
    Int downsampleToReads = 500000
    Int coverageWindow = 1000
  }

  parameter_meta {
    inputGroups: "Array of objects describing sets of bams to merge together and on which to compute QC metrics"
    metadata: "JSON file containing metadata"
    targetBed: "Path to optional target bed file"
    reference: "Reference id, we need it to pick the right reference file"
    outputFileNamePrefix: "Prefix for output files"
    downsampleToReads: "Downsample to this many reads when running unique read count, duplicate rate calculation and CIGAR analysis"
    coverageWindow: "Coverage window to use with mosdepth for making coverage histogram, default is 1000b"
  }

  Map[String,Resources] resources = {
    "hg38": {
      "refFasta": "$HG38_ROOT/hg38_random.fa",
      "modules": "samtools/1.16.1 hg38/p12"
    },
    "hg19": {
      "refFasta": "$HG19_ROOT/hg19_random.fa",
      "modules": "samtools/1.16.1 hg19/p13"
    }
  }

  
  # Record mode in metadata, but process whatevever number of files we have
  scatter (i in inputGroups) {
   
        call samstats {
        input:
          bamFile = i.bam,
          modules = resources[reference].modules,
          referenceFile = resources[reference].refFasta
        }

        call runMosdepth {
        input:
          bamFile = i.bam,
          bamIndex = i.bamIndex,
          targetBed = targetBed,
          prefix = basename(i.bam, ".bam")
        }

        # Derive runs on target if there is a .bed file and metric is requested
        if (defined(targetBed)) {
          call runBedtoolsIntersect {
          input:
            inputBam = i.bam,
            targetBed = targetBed
          }
        }

        call cumulativeDistToHistogram {
        input:
          globalDist = runMosdepth.globalDist,
          summary = runMosdepth.summary
        }

        call markDuplicates {
        input:
          bamFile = i.bam,
          downsampleToReads = downsampleToReads
        }

        call getUniqueReadsCount {
         input:
           inputBam = i.bam,
           downsampleToReads = downsampleToReads
        }

        # If we have multiple lanes, downsample bam files and run windowed coverage analysis
        if (length(inputGroups) > 1) {

           call downsampleBam {
           input:
             bamFile = i.bam,
             prefix = basename(i.bam, ".bam"),
             downsampleToReads = downsampleToReads
           }

          call runMosdepth as runWindowedMosdepth {
          input:
            bamFile = i.bam,
            bamIndex = i.bamIndex,
            prefix = basename(i.bam, ".bam"),
            additionalParameters = "--by 1000"
          } 

        }
        
        # This task collates the results from prior steps but also runs CIGAR metrics collection on a downsampled file
        call bamQCMetrics {
        input:
          metadata = metadata,
          bamFile = i.bam, 
          downsampleToReads = downsampleToReads,
          readsOnTarget = runBedtoolsIntersect.readsOnTarget,
          uniqueReads = getUniqueReadsCount.uniqueReads,
          outputFileNamePrefix = outputFileNamePrefix,
          samstatsFile = samstats.statsFile,
          histogram = cumulativeDistToHistogram.histogram,
          referenceFileName = basename(resources[reference].refFasta),
          mosdepthSummary = runMosdepth.summary,
          targetBed = targetBed,
          markDuplicatesStats = markDuplicates.result
        }

  }

  # We will call markDuplicate and coverage histogram-creating functions specific for call-ready data
  # if we have multiple lanes of data
  if (length(inputGroups) > 1) {
   
    call markDuplicatesMerged {
        input:
          bamFiles = select_all(downsampleBam.downsampledBam),
          outputFileNamePrefix = outputFileNamePrefix
        }

    call mergedCoverageToHistogram {
        input:
          coverageFiles = select_all(runWindowedMosdepth.targetCoverage),
          outputFileNamePrefix = outputFileNamePrefix
    }
  }
   
  # json merger will either merge multiple inputs or return single-lane report as the final report
  call mergeReports {
  input:
    inputs = bamQCMetrics.result,
    mergedDupmarkingData = markDuplicatesMerged.result,
    mergedCoverageData = mergedCoverageToHistogram.histogram,
    prefix = outputFileNamePrefix
  }

  
  output {
    File result = mergeReports.finalReport
  }

  meta {
        author: "Iain Bancarz, Savo Lazic, Peter Ruzanov"
        email: "ibancarz@oicr.on.ca, slazic@oicr.on.ca, pruzanov@oicr.on.ca"
        description: "bamQC workflow collects a number of metrics which are computed using several methods (by employing third-party software tools along with some custom code) and outputs the results in JSON format. The output also contains metadata, such as the instrument and lane names. bamQC supports downsampling for faster analysis."
        dependencies: [
            {
                name: "samtools/1.16.1",
                url: "https://github.com/samtools/samtools"
            },
            {
                name: "samblaster/0.1.26",
                url: "https://github.com/GregoryFaust/samblaster"
            },
            {
                name: "python/3.6",
                url: "https://www.python.org/downloads/"
            },
            {
                name: "bam-qc-metrics/0.2.7",
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
#  1 of 8 : run samstats, extract selected metrics
# ================================================================
task samstats {

    # run the stats to get non-primary, unmapped, and low-quality aligned reads
    # return filtered read counts

   input {
    File bamFile
    String referenceFile
    String modules
    String? additionalParameters
    Int jobMemory = 16
    Int timeout = 24
    }

    parameter_meta {
    bamFile: "Input BAM file of aligned rnaSeqQC data"
    referenceFile: "Reference fasta file, needed for MPC stats"
    modules: "required environment modules"
    additionalParameters: "Additional parameters for samtools, usually filtering flag"
    jobMemory: "Memory allocated for this job"
    timeout: "hours before task timeout"
    }

    String filePrefix = basename(bamFile, ".bam")

    # Leaving this here in case we reconsider filtering:
    # -F 2304 excludes secondary and supplementary alignments
    # -F 4 excludes unmapped reads

    command <<<
    samtools stats ~{additionalParameters} ~{bamFile} -r  ~{referenceFile} > ~{filePrefix}.stats
    >>>

    runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    timeout: "~{timeout}"
    }

    # record read totals as String, not Int, to avoid integer overflow error
    output {
      File statsFile = "~{filePrefix}.stats"
    }

    meta {
    output_meta: {
      statsFile: "File with bam file stats"
    }
    }

}


# ================================================================
#  Optional task for downsampling
# ================================================================
task downsampleBam {
   
    input {
    File bamFile
    String modules = "samtools/1.16.1"
    String prefix
    Int downsampleToReads
    Int jobMemory = 16
    Int timeout = 4
    }

    parameter_meta {
    bamFile: "Input BAM file of aligned data"
    modules: "required environment modules"
    prefix: "prefix for output files"
    jobMemory: "Memory allocated for this job"
    timeout: "hours before task timeout"
    }

    runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    timeout: "~{timeout}"
    }

    String bamFileName = basename(bamFile)

    command <<<
    samtools head -n ~{downsampleToReads} ~{bamFile} | samtools view -hb > ~{prefix + ".downsampled.bam"}
    >>>

    output {
    File downsampledBam = "~{prefix}.downsampled.bam"
    }

    meta {
    output_meta: {
        downsampledBam: "Downsampled BAM file"
    }
  }
}

# ================================================================
#  2 of 8 : run Mosdepth, extract coverage distribution data
# ================================================================
task runMosdepth {

    input {
    File bamFile
    File bamIndex
    String modules = "mosdepth/0.2.9"
    String? additionalParameters
    String? targetBed
    String prefix
    Int jobMemory = 16
    Int timeout = 4
    }

    parameter_meta {
    bamFile: "Input BAM file of aligned data"
    bamIndex: "Index file in samtools .bai format"
    additionalParameters: "Additional parameters for coverage calculation"
    modules: "required environment modules"
    prefix: "prefix for output files"
    targetBed: "Optional target bed file"
    jobMemory: "Memory allocated for this job"
    timeout: "hours before task timeout"
    }

    runtime {
    modules: "~{modules}"
    memory:  "~{jobMemory} GB"
    timeout: "~{timeout}"
    }

    String bamFileName = basename(bamFile)

    command <<<
    set -eo pipefail
    # ensure BAM file and index are symlinked to working directory
    ln -s ~{bamFile}
    ln -s ~{bamIndex}
    # run mosdepth
    MOSDEPTH_PRECISION=8 mosdepth -x -n ~{"--by " + targetBed} ~{additionalParameters} -t 3 ~{prefix} ~{bamFileName}
    >>>

    output {
    File globalDist = "~{prefix}.mosdepth.global.dist.txt"
    File summary = "~{prefix}.mosdepth.summary.txt"
    File? targetCoverage = "~{prefix}.regions.bed.gz"
    }

    meta {
    output_meta: {
        globalDist: "Global distribution of coverage",
        summary: "Total bases in coverage",
        targetCoverage: "Optional data, targeted coverage"
    }
  }

}


# ================================================================
#  3a of 8 : format histogram data using mosDepth output
# ================================================================
task cumulativeDistToHistogram {

    input {
    File globalDist
    File summary
    String modules = "bam-qc-metrics/0.2.7"
    String coverageHistogram = "$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_histogram.py"
    String outFileName = "coverage_histogram.json"
    Int jobMemory = 8
    Int threads = 4
    Int timeout = 1
    }

    parameter_meta {
    globalDist: "Global coverage distribution output from mosdepth"
    summary: "Summary output from mosdepth"
    modules: "required environment modules"
    coverageHistogram: "Path to script generating coverage histogram"
    outFileName: "Output file name, default coverage_histogram.json"
    jobMemory: "Memory allocated for this job"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
    }

    # mosdepth writes a global coverage distribution with 3 columns:
    # 1) Chromsome name, or "total" for overall totals
    # 2) Depth of coverage
    # 3) Probability of coverage less than or equal to (2)
    # Want to convert the above cumulative probability distribution to a histogram
    # The "total" section of the summary discards some information
    # So, we process the outputs for each chromosome to construct the histogram

    command <<<
        python3 ~{coverageHistogram} -s ~{summary} -g ~{globalDist} -o ~{outFileName}
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


# ================================================================
#  3b of 8 : format histogram data using mosDepth output
# ================================================================
task mergedCoverageToHistogram {
  
    input {
    Array[File] coverageFiles
    String outputFileNamePrefix
    String modules = "bam-qc-metrics/0.2.7"
    String coverageMerge = "$BAM_QC_METRICS_ROOT/bin/bam_qc_coverage_merger.py"
    Int jobMemory = 8
    Int threads = 4
    Int timeout = 1
    }

    parameter_meta {
    coverageFiles: "Array of gzip-ed files with windowed coverage from mosDepth"
    outputFileNamePrefix: "Output file name prefix"
    modules: "required environment modules"
    coverageMerge: "Path to coverage merging script"
    jobMemory: "Memory allocated for this job"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
    }

    String outFileName = "~{outputFileNamePrefix}_coverage_histogram.json"

    command <<<
    python3 ~{coverageMerge} -f ~{sep="," coverageFiles} \
                             -o ~{outFileName} 
    >>>

    runtime {
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
      timeout: "~{timeout}"
    }

    output {
      File histogram = "~{outFileName}"
    }

    meta {
    output_meta: {
      result: "Text file with samblaster markDuplicates metrics"
    }
    }
}


# ================================================================
#  4a of 8 : run samblaster to collect read duplicates stats
# ================================================================
task markDuplicates {

    input {
    File bamFile
    String? additionalParameters
    Int downsampleToReads
    String modules = "samblaster/0.1.26 samtools/1.16.1"
    Int jobMemory = 16
    Int threads = 4
    Int timeout = 4
    }

    parameter_meta {
    bamFile: "Input BAM file, after filtering and downsampling (if any)"
    additionalParameters: "Any additional parameters for samblaster"
    modules: "required environment modules"
    jobMemory: "Memory allocated for this job"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
    }

    String filePrefix = basename(bamFile, ".bam")

    command <<<
    set -euxo pipefail
    samtools head -n ~{downsampleToReads} ~{bamFile} | \
    samtools view -h -F 2308 - | \
    samtools sort -n - | \
    samtools fixmate -m -O SAM - - | \
    samblaster --ignoreUnmated ~{additionalParameters} --output /dev/null 2> >(tee "~{filePrefix}.markDuplicates.txt")
    >>>

    runtime {
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
      timeout: "~{timeout}"
    }

    output {
      File result = "~{filePrefix}.markDuplicates.txt"
    }

    meta {
    output_meta: {
      result: "Text file with samblaster markDuplicates metrics"
    }
    }

}

# =====================================================================================
#  4b of 8 : run samblaster to collect read duplicates stats on merged lane-level data
# =====================================================================================
task markDuplicatesMerged {

    input {
    Array[File] bamFiles
    String outputFileNamePrefix
    String? additionalParameters
    String modules = "samblaster/0.1.26 samtools/1.16.1"
    Int jobMemory = 16
    Int threads = 4
    Int timeout = 4
    }

    parameter_meta {
    bamFiles: "Input BAM files, after filtering and downsampling (if any)"
    outputFileNamePrefix: "Output file name prefix"
    additionalParameters: "Any additional parameters for samblaster"
    modules: "required environment modules"
    jobMemory: "Memory allocated for this job"
    threads: "Requested CPU threads"
    timeout: "hours before task timeout"
    }

    command <<<
    set -euxo pipefail
    samtools cat ~{sep=" " bamFiles} | \
    samtools view -h -F 2308 - | \
    samtools sort -n - | \
    samtools fixmate -m -O SAM - - | \
    samblaster --ignoreUnmated ~{additionalParameters} --output /dev/null 2> >(tee "~{outputFileNamePrefix}.markDuplicates.txt")
    >>>

    runtime {
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
      timeout: "~{timeout}"
    }

    output {
      File result = "~{outputFileNamePrefix}.markDuplicates.txt"
    }

    meta {
    output_meta: {
      result: "Text file with samblaster markDuplicates metrics"
    }
    }

}

# =============================================================================
#  5 of 8: Optional bedtools intersect task to derive reads on target metric, if needed
# =============================================================================
task runBedtoolsIntersect {
    input {
    File inputBam
    File? targetBed
    String modules = "samtools/1.16.1 bedtools/2.27"
    Int jobMemory = 16
    Int timeout = 12
    }

    parameter_meta {
    inputBam: "Input BAM file, after filtering and downsampling (if any)"
    targetBed: "Target bed file"
    modules: "required environment modules"
    jobMemory: "Memory allocated for this job"
    timeout: "hours before task timeout"
    }

    command <<<
    bedtools intersect -a ~{inputBam} -b ~{targetBed} -u | samtools view -c | perl -pe 'chomp'
    >>>

    runtime {
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
      timeout: "~{timeout}"
    }
    
    output {
      Int readsOnTarget = read_int(stdout())
    }

    meta {
    output_meta: {
      readsOnTarget: "Number of reads on target"
    }
    }

}

# =============================================================================
#  6 of 8: Optional samtools countinting task to get unique reads, if needed
# =============================================================================
task getUniqueReadsCount {
    input {
    File inputBam
    Int downsampleToReads
    String modules = "samtools/1.16.1"
    Int jobMemory = 16
    Int timeout = 18
    }

    parameter_meta {
    inputBam: "Input BAM file, after filtering and downsampling (if any)"
    downsampleToReads: "Run the count on this many reads only bed file"
    modules: "required environment modules"
    jobMemory: "Memory allocated for this job"
    timeout: "hours before task timeout"
    }

    command <<<
    samtools head -n ~{downsampleToReads} ~{inputBam} | samtools view -F 256 -q 30 -c | perl -pe 'chomp'
    >>>

    runtime {
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
      timeout: "~{timeout}"
    }

    output {
      Int uniqueReads = read_int(stdout())
    }

    meta {
    output_meta: {
      uniqueReads: "Number of unique reads"
    }
    }

}

# ============================================================================
#  7 of 8 : run bamQCMetrics to assemble lane-level metrics into a json file
# ============================================================================
task bamQCMetrics {

    input {
    File? bamFile
    Map[String, String] metadata
    String outputFileNamePrefix
    File samstatsFile
    File histogram
    File markDuplicatesStats
    File mosdepthSummary
    File? targetBed
    Int? readsOnTarget
    Int? downsampleToReads
    Int? uniqueReads
    String referenceFileName
    String workflowVersion = "5.3.1"
    String modules = "bam-qc-metrics/0.2.7"
    String bamQClite = "$BAM_QC_METRICS_ROOT/bin/run_bam_qc_lite.py"
    Int jobMemory = 8
    Int timeout = 12
    }

    parameter_meta {
    bamFile: "We pass a bam file only if we want to run CIGAR analysis (long)"
    bamQClite: "Path to bamQC lite script"
    outputFileNamePrefix: "Prefix for output file"
    downsampleToReads: "Downsample to this number of reads, default is defined in bam_qc_lite (500K)"
    histogram: "JSON file(s) with coverage histogram"
    metadata: "Map object with additional metadata"
    samstatsFile: "files with metrics collected by samtools stats"
    markDuplicatesStats: "stats from markDuplicate step"
    mosdepthSummary: "Summary file from mosdepth task"
    targetBed: "Optional target bed"
    readsOnTarget: "Optional reads on target metric"
    uniqueReads: "Optional unique reads for calculating reads per start point metric"
    referenceFileName: "basename of the reference file"
    workflowVersion: "Workflow version to put into report" 
    modules: "required environment modules"
    jobMemory: "Memory allocated for this job"
    timeout: "hours before task timeout"
    }

    runtime {
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
      timeout: "~{timeout}"
    }

    String outputFileName = "~{outputFileNamePrefix}.bamQC_results.json"
    File metadataJson = write_json(metadata)
    command <<<
        python3 ~{bamQClite} ~{"-b " + bamFile} \
        -s ~{samstatsFile} \
        -d ~{markDuplicatesStats} \
        -c ~{histogram} \
        -m ~{metadataJson} \
        -w ~{workflowVersion} ~{"-tf " + targetBed} \
        -r ~{referenceFileName} ~{"-S " + downsampleToReads} \
        -o ~{outputFileName} \
        -t ~{mosdepthSummary} ~{"-u " + uniqueReads} ~{"-ot " + readsOnTarget} 
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

# ==========================================================================================
# 8 of 8: merge reports from multiple lanes or pass a single-lane report as the final report
# ==========================================================================================
task mergeReports {
    input {
    Array[File] inputs
    File? mergedDupmarkingData
    File? mergedCoverageData
    String prefix
    String modules = "bam-qc-metrics/0.2.7"
    String bamQCmerger = "$BAM_QC_METRICS_ROOT/bin/bam_qc_merger.py"
    Int jobMemory = 4
    Int timeout = 2
    }

    parameter_meta {
      inputs: "Array of .json report files"
      mergedDupmarkingData: "Optional merged (call-ready) duplicate marking report" 
      mergedCoverageData: "Optional merged (call-ready) coverage histogram"
      prefix: "Prefix for output file"
      modules: "Runtime modules"
      bamQCmerger: "Path to the merger script"
      jobMemory: "RAM allocated to run the merging task"
      timeout: "Timeout in hours for the merging task"
    }

    runtime {
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
      timeout: "~{timeout}"
    }

    String outputFileName = "~{prefix}.bamQC_results.json"
    command <<<
        # set -euxo pipefail
        # export PYTHONPATH=$PYTHONPATH:/.mounts/labs/gsi/testdata/bamqc/scripts/ 
        python3 ~{bamQCmerger} -l ~{sep="," inputs} ~{"-d " + mergedDupmarkingData} ~{"-t " + mergedCoverageData} -o ~{outputFileName}
    >>>

    output {
      File finalReport = "~{outputFileName}"
    }

    meta {
    output_meta: {
      finalReport: "JSON file with merged results"
    }
    }
}
