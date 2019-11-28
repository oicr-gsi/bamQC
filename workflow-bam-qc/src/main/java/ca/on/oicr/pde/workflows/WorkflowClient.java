package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.apache.commons.lang3.StringEscapeUtils;

public class WorkflowClient extends OicrWorkflow {

    String outdir = null;
    // private static final Logger logger = Logger.getLogger(WorkflowClient.class.getName());
    //workflow parameters
    private String queue = null;
    private String inputFile = null;
    private String markDuplicatesMetricsFile = null;
    private String sampleLevel = null;
    private String normalInsertMax = null;
    private String mapQualCut = null;
    private String targetBed = null;
    private String jsonMetadata = null;
    private String jsonMetadataFile = null;
    private String markDuplicatesTextFile = null;
    private String markDuplicatesBamFile = null;
    private String reference = null;
    //workflow directories
    private String dataDir = null;
    private String tmpDir = null;
    private Boolean manualOutput = false;
    private Boolean markDuplicates = true;
    private String workflowVersion = null;
    //module parameters
    private String bamQcMetricsModule = null;
    private String bamQcMetricsVersion = null;
    private String picardToolsModule = null;
    private String picardToolsVersion = null;
    private String picardMarkDuplicatesAdditionalParams = null;

    //Constructor - called in setupDirectory()
    private void WorkflowClient() {

        dataDir = "data/";
        tmpDir = "tmp/";
        queue = getOptionalProperty("queue", "");
        inputFile = getProperty("input_file");
        markDuplicatesMetricsFile = getOptionalProperty("input_mark_duplicates_metrics", null);
        manualOutput = Boolean.valueOf(getProperty("manual_output"));
        markDuplicates = Boolean.valueOf(getOptionalProperty("mark_duplicates", "true"));
        sampleLevel = getProperty("sample_level"); // total reads desired in sample; see bam-qc-metrics
        normalInsertMax = getProperty("normal_insert_max");
        mapQualCut = getProperty("map_qual_cut");
        reference = getProperty("reference");
        targetBed = getProperty("target_bed");
        workflowVersion = getProperty("workflow_version"); // bam-qc-metrics requires 3-part version, eg. 0.1.2

        if (hasPropertyAndNotNull("json_metadata")) {
            jsonMetadata = StringEscapeUtils.unescapeJava(getProperty("json_metadata")).replace("&#61;", "=");
            jsonMetadataFile = dataDir + "metadata.json";
        }

        bamQcMetricsModule = getProperty("bam_qc_metrics_module");
        bamQcMetricsVersion = getProperty("bam_qc_metrics_version");
        picardToolsModule = getProperty("picard_module");
        picardToolsVersion = getProperty("picard_version");
        picardMarkDuplicatesAdditionalParams = getOptionalProperty("picard_mark_duplicates_additional_params", null);
    }

    @Override
    public void setupDirectory() {
        WorkflowClient(); //Constructor call
        addDirectory(dataDir);
        addDirectory(tmpDir);
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        SqwFile file0 = this.createFile("file_in_0");
        file0.setSourcePath(inputFile);
        file0.setType("application/bam");
        file0.setIsInput(true);

        if (markDuplicatesMetricsFile != null) {
            SqwFile file1 = this.createFile("file_in_1");
            file1.setSourcePath(markDuplicatesMetricsFile);
            file1.setType("text/plain");
            file1.setIsInput(true);

            markDuplicatesTextFile = file1.getProvisionedPath();
        }

        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {
        Job job0 = null;
        if (jsonMetadata != null) {
            job0 = getWriteToFileJob(jsonMetadata, jsonMetadataFile);
            job0.setMaxMemory(getProperty("json_metadata_job_mem"));
            job0.setQueue(queue);
        }

        Job job1 = null;
        if (markDuplicates) {
            job1 = getPicardMarkDuplicatesJob();
            if (job0 != null) {
                job1.addParent(job0);
            }
        }

        Job job2 = getBamQcJob();
        job2.setMaxMemory(getProperty("bamqc_job_mem"));
        job2.setQueue(queue);
        if (job1 != null) {
            job2.addParent(job1);
        } else if (job0 != null) {
            job2.addParent(job0);
        }
    }

    private Job getWriteToFileJob(String fileContents, String outputFile) {
        // if metadata was supplied as a JSON blob in the INI, write to file
        Job job = getWorkflow().createBashJob("WriteToFile");

        List<String> c = new LinkedList<>();
        c.add("set -o pipefail;");
        c.add("set -o errexit;");

        c.add("cat << END_OF_FILE_CONTENTS");
        c.add(">");
        c.add(outputFile);
        c.add("\n");
        c.add(fileContents);
        c.add("\nEND_OF_FILE_CONTENTS\n");

        job.getCommand().setArguments(c);

        return job;
    }

    private Job getPicardMarkDuplicatesJob() {
        Integer picardMaxMemMb = Integer.parseInt(getProperty("picard_memory"));
        markDuplicatesTextFile = dataDir + "mark_duplicates.txt";
        markDuplicatesBamFile = tmpDir + "marked_duplicates.bam";
        Job job = getWorkflow().createBashJob("PicardMarkDuplicates");
        job.setMaxMemory(getProperty("picard_job_memory"));
        job.setQueue(queue);
        Command command = job.getCommand();
        command.addArgument("module load " + picardToolsModule + "/" + picardToolsVersion + ";");
        command.addArgument("java"); // java module is loaded by Picard module
        command.addArgument("-Xmx" + picardMaxMemMb + "M");
        command.addArgument("-jar " + getOptionalProperty("picard_jar", "${PICARD_ROOT}/picard.jar"));
        command.addArgument("MarkDuplicates");
        command.addArgument("INPUT=" + getFiles().get("file_in_0").getProvisionedPath());
        command.addArgument("OUTPUT=" + markDuplicatesBamFile);
        command.addArgument("VALIDATION_STRINGENCY=SILENT");
        command.addArgument("TMP_DIR=" + tmpDir);
        command.addArgument("METRICS_FILE=" + markDuplicatesTextFile);
        if (picardMarkDuplicatesAdditionalParams != null) {
            command.addArgument(picardMarkDuplicatesAdditionalParams);
        }
        return job;
    }

    private Job getBamQcJob() {
        Job job = getWorkflow().createBashJob("BamToJsonStats");
        String jsonOutputFileName = inputFile.substring(inputFile.lastIndexOf("/") + 1) + ".BamQC.json";
        String logFileName = "run_bam_qc.log";
        String inputBamFile;
        if (markDuplicatesBamFile != null) {
            inputBamFile = markDuplicatesBamFile;
        } else {
            inputBamFile = getFiles().get("file_in_0").getProvisionedPath();
        }

        Command command = job.getCommand();
        command.addArgument("module load " + bamQcMetricsModule + "/" + bamQcMetricsVersion + " && ");
        command.addArgument("run_bam_qc.py ");
        command.addArgument("-b " + inputBamFile);
        command.addArgument("--debug ");
        command.addArgument("-i " + normalInsertMax);
        command.addArgument("-l " + dataDir + logFileName);
        command.addArgument("-o " + dataDir + jsonOutputFileName);
        command.addArgument("-q " + mapQualCut);
        command.addArgument("-r " + reference);
        command.addArgument("-s " + sampleLevel);
        command.addArgument("-t " + targetBed);
        command.addArgument("-T " + tmpDir);
        command.addArgument("-w " + workflowVersion);
        if (jsonMetadataFile != null) {
            command.addArgument("-m " + jsonMetadataFile);
        }
        if (markDuplicatesTextFile != null) {
            command.addArgument("-d " + markDuplicatesTextFile);
        }

        SqwFile sqwJsonOutputFile = createOutputFile(dataDir + jsonOutputFileName, "text/json", manualOutput);

        job.addFile(sqwJsonOutputFile);
        job.setQueue(queue);
        return job;
    }

}
