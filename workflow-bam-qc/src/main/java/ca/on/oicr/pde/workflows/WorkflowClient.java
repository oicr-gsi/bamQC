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
    private String sampleRate = null;
    private String normalInsertMax = null;
    private String mapQualCut = null;
    private String targetBed = null;
    private String jsonMetadata = null;
    private String jsonMetadataFile = null;
    private String markDuplicatesTextFile = null;
    private String markDuplicatesBamFile = null;
    //workflow directories
    private String dataDir = null;
    private String tmpDir = null;
    private Boolean manualOutput = false;
    private Boolean markDuplicates = true;

    //Constructor - called in setupDirectory()
    private void WorkflowClient() {

        dataDir = "data/";
        tmpDir = "tmp/";
        queue = getOptionalProperty("queue", "");
        inputFile = getProperty("input_file");
        manualOutput = Boolean.valueOf(getProperty("manual_output"));
        markDuplicates = Boolean.valueOf(getOptionalProperty("mark_duplicates", "true"));
        sampleRate = getProperty("sample_rate");
        normalInsertMax = getProperty("normal_insert_max");
        mapQualCut = getProperty("map_qual_cut");
        targetBed = getProperty("target_bed");

        if (hasPropertyAndNotNull("json_metadata")) {
            jsonMetadata = StringEscapeUtils.unescapeJava(getProperty("json_metadata")).replace("&#61;", "=");
            jsonMetadataFile = dataDir + "metadata.json";
        }
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
        String java = getProperty("java");
        String picard = getProperty("picard");
        Integer picardMaxMemMb = Integer.parseInt(getProperty("picard_memory"));
        markDuplicatesTextFile = dataDir + "mark_duplicates.txt";
        markDuplicatesBamFile = tmpDir + "marked_duplicates.bam";
        Job job = getWorkflow().createBashJob("PicardMarkDuplicates");
        job.setMaxMemory(getProperty("picard_job_memory"));
        job.setQueue(queue);
        Command command = job.getCommand();
        command.addArgument(java);
        command.addArgument("-Xmx" + picardMaxMemMb + "M");
        command.addArgument("-jar " + picard + "/MarkDuplicates.jar");
        command.addArgument("I=" + getFiles().get("file_in_0").getProvisionedPath());
        command.addArgument("O=" + markDuplicatesBamFile);
        command.addArgument("M=" + markDuplicatesTextFile);
        return job;
    }

    private Job getBamQcJob() {
        Job job = getWorkflow().createBashJob("BamToJsonStats");
        String pythonpath = getProperty("pythonpath");
        String jsonOutputFileName = inputFile.substring(inputFile.lastIndexOf("/") + 1) + ".BamQC.json";
        String inputBamFile;
        if (markDuplicatesBamFile != null) {
            inputBamFile = markDuplicatesBamFile;
        } else {
            inputBamFile = getFiles().get("file_in_0").getProvisionedPath();
        }

        Command command = job.getCommand();
        command.addArgument("source /u/ibancarz/installed/miniconda3/bin/activate &&"); //temporary for testing
        command.addArgument("PYTHONPATH=" + pythonpath);
        command.addArgument(getProperty("metrics_script"));
        command.addArgument("-b " + inputBamFile);
        command.addArgument("-s " + sampleRate);
        command.addArgument("-i " + normalInsertMax);
        command.addArgument("-o " + dataDir + jsonOutputFileName);
        command.addArgument("-q " + mapQualCut);
        command.addArgument("-t " + targetBed);
        command.addArgument("-T " + tmpDir);
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
