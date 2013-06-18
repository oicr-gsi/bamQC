package ca.on.oicr.pde.workflows;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
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
    private String groupId = null;
    private String groupIdDescription = null;
    private String externalName = null;
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
            groupId = getProperty("group_id");
            groupIdDescription = getProperty("group_id_description");
            externalName = getProperty("external_name");

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
        command.addArgument(getWorkflowBaseDir() + "/bin" + "/samtools-0.1.19/samtools " + "view " + getFiles().get("file_in_0").getProvisionedPath());
        command.addArgument("|"); //pipe to
        command.addArgument(getWorkflowBaseDir() + "/bin" + "/perl-5.14.1/perl");
        command.addArgument(getWorkflowBaseDir() + "/bin" + "/samStats.pl");
        command.addArgument("-s " + sampleRate);
        command.addArgument("-i " + normalInsertMax);
        command.addArgument("-q " + mapQualCut);
        command.addArgument("-r " + targetBed);
        command.addArgument("-j " + jsonMetadataFile);

        if (isPropertySet(groupId)) {
            command.addArgument("-g " + escapeStringForSeqwareShellCommand(groupId));
        }

        if (isPropertySet(groupIdDescription)) {
            command.addArgument("-d " + escapeStringForSeqwareShellCommand(groupIdDescription));
        }

        if (isPropertySet(externalName)) {
            command.addArgument("-n " + escapeStringForSeqwareShellCommand(externalName));
        }

//        command.addArgument("-w " + getName());  //workflow name of workflow that generated bam file?
//        command.addArgument("-v " + getVersion()); //workflow version of workflow that generated bam file?
//        command.addArgument("-t " + Long.toString(System.currentTimeMillis()/1000L)); // workflow run creation timestamp (number of seconds since epoch/1970)
//        command.addArgument("-b " + inputFile); //-b records bam path and last modification data
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

    private boolean isPropertySet(String property) {

        List unsetPropertyValueSet = Arrays.asList("na", "", "null");

        return !(unsetPropertyValueSet.contains(property.toLowerCase().trim()));

    }

    private String escapeStringForSeqwareShellCommand(String input) {

        //don't escape characters defined by the regex
        Pattern allowedUnescapedCharactersRegex = Pattern.compile("[0-9A-Za-z _]");

        //prefix/suffix used to create html code
        String escapedCharacterPrefix = "&#";
        String escapedCharacterSuffix = ";";

        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < input.length(); i++) {
            if (allowedUnescapedCharactersRegex.matcher(input.subSequence(i, i + 1)).find()) {
                sb.append(input.charAt(i));
            } else {
                sb.append(escapedCharacterPrefix);
                sb.append(input.codePointAt(i));
                sb.append(escapedCharacterSuffix);
            }
        }

        String result = sb.toString();
        
        //Escape all characters for bash.
        //Escape backslash for java, and then escape backslash for seqware/pegasus/condor.
        // ie, (java) \\\\X -> (seqware) \\X -> (bash) \X
        result = result.replaceAll(".", "\\\\$0");
        
        //Group together for seqware/pegasus/condor arg using single quotes
        //Note: seqware/pegasus/condor converts double quotes to single quotes
        //Note: These quotes will note be seen by bash.
        result = "'" + result + "'";
        
        return result;
    }
}