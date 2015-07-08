package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

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
    private String jsonMetadataFile = null;
    //workflow directories
    private String binDir = null;
    private String dataDir = null;
    private Boolean manualOutput = false;

    //Constructor - called in setupDirectory()
    private void WorkflowClient() {

        binDir = getWorkflowBaseDir() + "/bin/";
        dataDir = "data/";
        queue = getOptionalProperty("queue", "");
        inputFile = getProperty("input_file");
        manualOutput = Boolean.valueOf(getProperty("manual_output"));
        sampleRate = getProperty("sample_rate");
        normalInsertMax = getProperty("normal_insert_max");
        mapQualCut = getProperty("map_qual_cut");
        targetBed = getProperty("target_bed");
        jsonMetadataFile = getProperty("json_metadata_file");

    }

    @Override
    public void setupDirectory() {
        WorkflowClient(); //Constructor call
        addDirectory(dataDir);
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
        Job job00 = getBamQcJob();
        job00.setMaxMemory("2000");
        job00.setQueue(queue);
    }

    private Job getBamQcJob() {
        Job job = getWorkflow().createBashJob("BamToJsonStats");

        String jsonOutputFileName = inputFile.substring(inputFile.lastIndexOf("/") + 1) + ".BamQC.json";

        Command command = job.getCommand();
        command.addArgument(binDir + "samtools-0.1.19/samtools " + "view " + getFiles().get("file_in_0").getProvisionedPath());
        command.addArgument("|"); //pipe to
        command.addArgument(binDir + "perl-5.14.1/perl");
        command.addArgument(getProperty("samstats_script"));
        command.addArgument("-s " + sampleRate);
        command.addArgument("-i " + normalInsertMax);
        command.addArgument("-q " + mapQualCut);
        command.addArgument("-r " + targetBed);
        command.addArgument("-j " + jsonMetadataFile);
        command.addArgument(">"); //redirect to
        command.addArgument(dataDir + jsonOutputFileName);

        SqwFile sqwJsonOutputFile = createOutputFile(dataDir + jsonOutputFileName, "text/json", manualOutput);

        job.addFile(sqwJsonOutputFile);
        job.setQueue(queue);
        return job;
    }

}
