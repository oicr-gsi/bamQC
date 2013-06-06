package ca.on.oicr.pde.workflows;

import java.util.Arrays;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends AbstractWorkflowDataModel {

    String outdir = null;
    private static final Logger logger = Logger.getLogger(WorkflowClient.class.getName());
    //workflow parameters
    private String queue = null;
    private String inputFile = null;
    private String outputPrefix = null;
    private String outputDir = null;
    private String outputPath = null;
    //private String jsonOutputFile = null;
    private String sampleRate = null;
    private String normalInsertMax = null;
    private String mapQualCut = null;
    private String targetBed = null;
    private String jsonMetadataFile = null;
    //workflow directories
    private String binDir = null;
    private String dataDir = null;
    private String finalOutputDir = null;

    //Constructor - called in setupDirectory()
    private void WorkflowClient() {

        binDir = getWorkflowBaseDir() + "/bin/";
        dataDir = "data/";

        try {
            
            
            queue = getProperty("queue");
            inputFile = getProperty("input_file");
            outputDir = getProperty("output_dir");
            outputPrefix = getProperty("output_prefix");
            outputPath = getProperty("output_path");
            
            if (Arrays.asList("na", "").contains(outputPath.toLowerCase().trim())) {
                finalOutputDir = outputPrefix + outputDir + "/seqware-" + getSeqware_version() + "_" + getName() + "_" + getVersion() + "/" + getRandom() + "/";
            } else {
                //make sure the path ends with a "/"
                outputPath = outputPath.lastIndexOf("/") == (outputPath.length() - 1) ? outputPath : outputPath + "/";
                finalOutputDir = outputPath;
            }
            
            sampleRate = getProperty("sample_rate");
            normalInsertMax = getProperty("normal_insert_max");
            mapQualCut = getProperty("map_qual_cut");
            targetBed = getProperty("target_bed");
            jsonMetadataFile = getProperty("json_metadata_file");
        } catch (Exception ex) {
            logger.log(Level.SEVERE, "Expected parameter missing", ex);
            System.exit(-1);
            //throw new RuntimeException(ex);
        }

    }

    @Override
    public void setupDirectory() {

        WorkflowClient(); //Constructor call
        addDirectory(dataDir);
        addDirectory(finalOutputDir);

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
        command.addArgument(getWorkflowBaseDir() + "/bin" + "/samtools-0.1.19/samtools " + "view " + inputFile);
        command.addArgument("|"); //pipe to
        command.addArgument(getWorkflowBaseDir() + "/bin" + "/perl-5.14.1/perl");
        command.addArgument(getWorkflowBaseDir() + "/bin" + "/samStats.pl");
        command.addArgument("-s " + sampleRate);
        command.addArgument("-i " + normalInsertMax);
        command.addArgument("-q " + mapQualCut);
        command.addArgument("-r " + targetBed);
        command.addArgument("-j " + jsonMetadataFile);
        command.addArgument(">"); //redirect to
        command.addArgument(dataDir + jsonOutputFileName);

        SqwFile sqwJsonOutputFile = createOutFile(dataDir + jsonOutputFileName, "text/json", finalOutputDir + jsonOutputFileName, true);
        job.addFile(sqwJsonOutputFile);

        return job;

    }

    private SqwFile createOutFile(String sourcePath, String sourceType, String outputPath, boolean forceCopy) {

        SqwFile file = new SqwFile();
        file.setSourcePath(sourcePath);
        file.setType(sourceType);
        file.setIsOutput(true);
        file.setOutputPath(outputPath);
        file.setForceCopy(forceCopy);

        return file;

    }
}