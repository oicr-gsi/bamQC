package ca.on.oicr.pde.workflows;

import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends AbstractWorkflowDataModel {

    String outdir = null;
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {
            // register an input file
            SqwFile file0 = this.createFile("file_in_0");
            file0.setSourcePath(getProperty("input_file"));
            file0.setType("text/plain");
            file0.setIsInput(true);

            return this.getFiles();
        } catch (Exception ex) {
            ex.printStackTrace();
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }

    }

    @Override
    public void setupDirectory() {

//        //String outdir = null;
//        try {
//            String outputPrefix = getProperty("output_prefix");
//            String outputDir = getProperty("output_dir");
//            if (outputPrefix != null && outputDir != null) {
//                outdir = outputPrefix + outputDir;
//                outdir = outdir.endsWith("/") ? outdir : outdir + "/";
//                this.addDirectory(outdir);
//            } else {
//                String inputPath = this.getFiles().get("file_in_0").getProvisionedPath();
//                outdir = inputPath.substring(0, outdir.lastIndexOf("/"));
//            }
//
//        } catch (Exception ex) {
//            Logger.getLogger(WorkflowClient.class.getName()).log(Level.WARNING, null, ex);
//        }
    }

    @Override
    public void buildWorkflow() {

        // Workflow parameters
        String queue = null;
        String outputPath = null;

        // load Workflow parameters
        try {
            queue = getProperty("queue");
            if (queue == null) {
                queue = "default";
            }
        } catch (Exception ex) {
            queue = "default";
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.WARNING, "queue not defined, using default queue");
        }
        
//        try {
//            outputPath = getProperty("output_path");
//
//            if (outputPath.equals("NA")) {
//                // set output path to same as input
//                this.getFiles().get("file_in_0").getProvisionedPath();
//            }
//            outputPath = outputPath.equalsIgnoreCase("NA") ? outputPath : this.getFiles().get("file_in_0").getProvisionedPath();
//        } catch (Exception ex){
//            
//        }

        // bamToJsonStats parameters
        String inputFile = null;
        String jsonOutputFile = null;
        String sampleRate = null;
        String normalInsertMax = null;
        String mapQualCut = null;
        String targetBed = null;
        String jsonMetadataFile = null;
        boolean bamToJsonStatsIsReady = false;

        // load bamToJsonStats parameters
        //TODO: Sanitize input
        try {
            inputFile = this.getFiles().get("file_in_0").getProvisionedPath();
            jsonOutputFile = inputFile + ".BamQC.json";
            sampleRate = getProperty("sample_rate");
            normalInsertMax = getProperty("normal_insert_max");
            mapQualCut = getProperty("map_qual_cut");
            targetBed = getProperty("target_bed");
            jsonMetadataFile = getProperty("json_metadata_file");
            bamToJsonStatsIsReady = true;
        } catch (Exception ex) {
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, ex);
        }

        if (bamToJsonStatsIsReady) {
            Job job00 = this.getWorkflow().createBashJob("bamToJsonStats");
            job00.setCommand(
                    getWorkflowBaseDir() + "/bin" + "/samtools-0.1.19/samtools " + "view " + inputFile
                    + " | "
                    + getWorkflowBaseDir() + "/bin" + "/perl-5.14.1/perl "
                    + getWorkflowBaseDir() + "/bin" + "/samStats.pl "
                    + "-s " + sampleRate + " "
                    + "-i " + normalInsertMax + " "
                    + "-q " + mapQualCut + " "
                    + "-r " + targetBed + " "
                    + "-j " + jsonMetadataFile
                    + " > " + jsonOutputFile);
            job00.setMaxMemory("2000");
            job00.setQueue(queue);

            SqwFile file1 = new SqwFile();
            file1.setSourcePath(jsonOutputFile);
            file1.setType("text/json");
            file1.setIsOutput(true);
            file1.setForceCopy(true);
        } else {
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, "bamToJsonStats is missing required parameters");
        }

    }
}