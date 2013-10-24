package ca.on.oicr.pde.deciders;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.common.util.maptools.MapTools;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

/**
 * @author mtaschuk@oicr.on.ca
 *
 */
public class BamQCDecider extends OicrDecider {

    private Map<String, ReturnValue> pathToAttributes = new HashMap<String, ReturnValue>();
    private String sampleRate = "1000";
    private String normalInsertMax = "1500";
    private String mapQualCut = "30";
    private String iniFile = null;
//    private String outputDir = "seqware-results";
//    private String outputPrefix = "./";
//    private String outputPath = "NA";
    private String tmp = "/tmp";
    private Random random = new Random(System.currentTimeMillis());
    private String rsconfigXmlPath = "/.mounts/labs/PDE/data/TargetedSequencingQC/rsconfig.xml";
    private boolean forceType = false;
    private String forcedResequencingType = "";
    private String forcedIntervalFile = "";

    public BamQCDecider() {
        super();
        parser.accepts("sample-rate", "Optional. Set the sampling rate in the decider "
                + "(only 1/sample_rate reads are used for memory-intensive sampling) "
                + "This needs to be set lower for mate-pair libraries. Default: 1000").withRequiredArg();
        parser.accepts("normal-insert-max", "Optional. Set the maximum insert size "
                + "to prevent skewing of insert statistics. Default:1500").withRequiredArg();
        parser.accepts("map-qual-cut", "Optional. Set the mapQ value. Default: 30").withRequiredArg();
        parser.accepts("ini-file", "Optional: an INI file with parameters to override the "
                + "installed INI file. ").withRequiredArg();
//        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to the output-path. "
//                + "Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
//        parser.accepts("output-path", "Optional: the path where the files should be copied to "
//                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("tmp", "Optional: specify the temporary directory where the JSON snippets will be stored during processing. Default: /tmp").withRequiredArg();
        parser.accepts("check-file-exists", "Optional: Flag to check whether or not the file exists before launching the workflow. WARNING! Will "
                + "not work if you are not on the same filesystem or do not have appropriate permissions!");
        parser.accepts("interval-file", "Optional: path to a file with target coordinates.").withRequiredArg();
        parser.accepts("resequencing-type", "Optional: specify resequencing type which should use the supplied interval-file. Will use the supplied interval file for this type only, "
                + "will works when interval-file is also set.").withRequiredArg();
        parser.accepts("rsconfig-file", "Optional: specify location of .xml file which should be used to configure references, "
                + "will be used if resequencing-type is different from the default."
                + "Default: " + rsconfigXmlPath).withRequiredArg();
        parser.accepts("force-type", "Optional: will process only bams without resequencing type set, "
                + " Need to have resequencing-type passed as well.").withOptionalArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setHeader(Header.FILE_SWA);
        this.setMetaType(Arrays.asList("application/bam"));


        if (options.has("sample-rate")) {
            sampleRate = options.valueOf("sample-rate").toString();
        }
        if (options.has("normal-insert-max")) {
            normalInsertMax = options.valueOf("normal-insert-max").toString();
        }
        if (options.has("map-qual-cut")) {
            mapQualCut = options.valueOf("map-qual-cut").toString();
        }
        if (options.has("ini-file")) {
            File file = new File(options.valueOf("ini-file").toString());
            if (file.exists()) {
                iniFile = file.getAbsolutePath();
                Map<String, String> iniFileMap = MapTools.iniString2Map(iniFile);
            } else {
                Log.stdout("The given INI file does not exist: " + file.getAbsolutePath());
                ret.setExitStatus(ReturnValue.INVALIDPARAMETERS);
            }

        }

        if (options.has("tmp")) {
            String temp = (String) options.valueOf("tmp");
            File tempDir = new File(temp);
            if (tempDir.exists()) {
                tmp = tempDir.getAbsolutePath();
            } else {
                Log.stdout("The temporary directory " + tempDir.getAbsolutePath() + " does not exist.");
                ret.setExitStatus(ReturnValue.INVALIDPARAMETERS);
            }
        }

        if (options.has("force-type")
                && !((options.has("resequencing-type") && options.has("interval-file"))
                || (options.has("resequencing-type")))) {
            Log.stdout("--force-type requires either --resequencing-type OR --resequencing-type and --interval-file to be set.");
            System.exit(1);
        }

        if (options.has("resequencing-type")) {
            if (!options.has("force-type")) {
                Log.stdout("Forcing --resequencing-type requires --resequencing-type to have an argument and --force-type be set.");
                System.exit(1);
            } else {
                forceType = true;
                forcedResequencingType = options.valueOf("resequencing-type").toString();
            }
        }

        if (options.has("interval-file")) {
            if (!options.has("resequencing-type") || !options.has("force-type")) {
                Log.stdout("Forcing --interval-file requires --interval-file to have an argument, --resequencing-type be set, and --force-type to be set.");
                System.exit(1);
            } else {
                forcedIntervalFile = options.valueOf("interval-file").toString();
                if (!fileExistsAndIsAccessible(forcedIntervalFile)) {
                    Log.stdout("ERROR: interval-file is not valid");
                    System.exit(1);
                }
            }
        }

        if (options.has("rsconfig-file")) {
            rsconfigXmlPath = options.valueOf("rsconfig-file").toString();
            if (!fileExistsAndIsAccessible(rsconfigXmlPath)) {
                Log.stdout("ERROR: rsconfig-file is not valid");
                System.exit(1);
            }
        }

        //allows anything defined on the command line to override the 'defaults' here.
        ret = super.init();
        return ret;

    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        Log.debug("GROUP BY ATTRIBUTE: " + getHeader().getTitle() + " " + attribute);
        return attribute;
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);

        String templateType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");

        if (forceType && !forcedResequencingType.isEmpty() && !forcedIntervalFile.isEmpty()) {
            returnValue.setAttribute("target_bed", forcedIntervalFile);
        } else if (forceType && !forcedResequencingType.isEmpty()) {

            String target_bed = null;
            try {
                target_bed = configFromParsedXML(rsconfigXmlPath, forcedResequencingType);
            } catch (Exception e) {
                Log.stdout(e.getMessage());
                System.exit(1);
            }

            if (target_bed == null) {
                //throw new Error("ERROR: FORCED TEMPLATE TYPE DOES NOT EXIST IN RSCONFIG");
                Log.stdout("ERROR: FORCED TEMPLATE TYPE DOES NOT EXIST IN RSCONFIG");
                System.exit(1);
            }

            returnValue.setAttribute("target_bed", target_bed);
        } else if ("WG".equals(templateType)) {
            returnValue.setAttribute("target_bed", "/oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed");
        } else if ("EX".equals(templateType)) {

            if ("Illumina TruSeq Exome".equals(targetResequencingType)) {
                returnValue.setAttribute("target_bed", "/oicr/data/reference/genomes/homo_sapiens_mc/TruSeq/TruSeq-Exome-Targeted-Regions-BED-file");
            } else if ("Agilent SureSelect ICGC/Sanger Exon".equals(targetResequencingType)) {
                returnValue.setAttribute("target_bed", "/oicr/data/genomes/homo_sapiens/Agilent/SureSelect_Whole_Exome_ICGC_Sanger/GRCh37hg19/sanger.exons.bed.hg19");
            } else if ("Nimblegen 2.1M Human Exome (21191)".equals(targetResequencingType)) {
                returnValue.setAttribute("target_bed", "/oicr/data/reference/genomes/homo_sapiens/Nimblegen/2.1M_Human_Exome_Annotation_21191/hg19/080904_ccds_exome_rebalfocus_hg19/processed/2.1M_Human_Exome.bed");
            } else {
                Log.stdout("ERROR: The targeted resequencing type does not have an associated BED file." + templateType + (" ") + targetResequencingType + " Modify the decider to include this type:" + targetResequencingType + " for file " + fm.getFilePath());
                return false;
            }

        } else if ("TS".equals(templateType)) {

            if (targetResequencingType == null || targetResequencingType.isEmpty()) {
                Log.stdout("ERROR: The targeted resequencing is not set.." + templateType + (" ") + "NULL" + " Set the targeted resequencing type for file " + fm.getFilePath());
                return false;
            }

            String target_bed = null;
            try {
                target_bed = configFromParsedXML(rsconfigXmlPath, targetResequencingType);
            } catch (Exception ex) {
                Log.stdout(ex.getMessage());
                System.exit(1);
            }
            if (target_bed == null) {
                Log.stdout("The targeted resequencing type does not have an associated BED file." + templateType + (" ") + targetResequencingType + " Modify the decider to include this type:" + targetResequencingType + " for file " + fm.getFilePath());
                return false;
            }
            returnValue.setAttribute("target_bed", target_bed);
        } else {
            Log.stdout("This template type is not supported for the BAM QC decider: " + templateType + " for file " + fm.getFilePath());
            return false;
        }
        if (options.has("check-file-exists")) {
            File file = new File(fm.getFilePath());
            if (!file.exists()) {
                file = new File(fm.getFilePath() + ".bak");
                if (file.exists()) {
                    Log.stdout("File does not exist! " + fm.getFilePath() + "\t .bak file exists: " + file.getAbsolutePath());
                } else {
                    Log.stdout("File does not exist! " + fm.getFilePath());
                }
                return false;
            }
        }
        pathToAttributes.put(fm.getFilePath(), returnValue);
        return super.checkFileDetails(returnValue, fm);
    }

    protected static String configFromParsedXML(String filePath, String resequencingType) throws ParserConfigurationException, SAXException, IOException, Exception {

        //TODO: call this function more efficiently

        Element eElement = null;
        File fXmlFile = new File(filePath);
        if (!fXmlFile.exists()) {
            //would be better if unchecked, but seqware catches unchecked :(
            throw new Exception("ERROR: XML FILE DOES NOT EXIST: " + fXmlFile.getCanonicalPath());
        }
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        Document doc = dBuilder.parse(fXmlFile);
        doc.getDocumentElement().normalize();

        NodeList nList = doc.getElementsByTagName("resequencing_type");

        if (nList.getLength() == 0) {
            //would be better if unchecked, but seqware catches unchecked :(
            throw new Exception("ERROR: NO RESEQUENCING TYPES LOCATED IN: " + filePath);
        }

        Map<String, String> h = new HashMap<String, String>();

        for (int i = 0; i < nList.getLength(); i++) {

            Node nNode = nList.item(i);

            if (nNode.getNodeType() == Node.ELEMENT_NODE) {

                eElement = (Element) nNode;

                String resequencingTypeId = eElement.getAttribute("id");
                String resequencingTypeBedFile = eElement.getElementsByTagName("interval_file").item(0).getTextContent();

                //validate key
                if (h.containsKey(resequencingTypeId)) {
                    //would be better if unchecked, but seqware catches unchecked :(
                    throw new Exception("ERROR: DUPLICATE RESEQUENCING_TYPE ID FOUND: " + resequencingTypeId);
                }

                //validate value
                if (resequencingTypeBedFile == null || resequencingTypeBedFile.isEmpty()) {
                    //would be better if unchecked, but seqware catches unchecked :(
                    throw new Exception("ERROR: INVALID INTERVAL FILE: " + resequencingTypeBedFile); 
                }

                h.put(resequencingTypeId, resequencingTypeBedFile);

            }

        }

        return h.get(resequencingType);

    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE: " + iniFile);
        Log.debug("Input file: " + commaSeparatedFilePaths);
        //Map<String, String> iniFileMap = new TreeMap<String, String>();
        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);

        if (iniFile != null) {
            MapTools.ini2Map(iniFile, iniFileMap);
        }

        //Remove "input_files" from ini file - BamQC workflow only accepts one BAM file at a time.
        //"input_files" is added to the iniFileMap by BasicDecider (parent of OicrDecider).
        iniFileMap.remove("input_files");

        if (commaSeparatedFilePaths.contains(",")) {
            Log.fatal("The BAM QC workflow only accepts one BAM file at a time. Try another grouping strategy, e.g. FILE_SWA. Files = " + commaSeparatedFilePaths);
            System.exit(1);
        }

        ReturnValue r = pathToAttributes.get(commaSeparatedFilePaths);

        iniFileMap.put("input_file", commaSeparatedFilePaths);
        iniFileMap.put("sample_rate", sampleRate);
        iniFileMap.put("normal_insert_max", normalInsertMax);
        iniFileMap.put("map_qual_cut", mapQualCut);
        iniFileMap.put("target_bed", r.getAttribute("target_bed"));
        iniFileMap.put("json_metadata_file", makeJsonSnippet(commaSeparatedFilePaths, r));

        return iniFileMap;
    }

    private String makeJsonSnippet(String filePath, ReturnValue r) {

        Map<String, String> atts = r.getAttributes();

        String runName = atts.get(Header.SEQUENCER_RUN_NAME.getTitle()), instrument;
        String[] tokens = runName.split("_");
        if (tokens.length >= 2) {
            instrument = tokens[1];
        } else {
            instrument = "NA";
        }

        String libraryName = atts.get(Header.SAMPLE_NAME.getTitle()), sample = "";
        tokens = libraryName.split("_");
        if (tokens.length > 2) {
            sample = libraryName.substring(0, libraryName.lastIndexOf(tokens[tokens.length - 2]) - 1);
        } else {
            sample = "NA";
        }

        String pSample = atts.get(Header.PARENT_SAMPLE_NAME.getTitle()), sampleGroup;
        tokens = pSample.split(":");
        if (tokens.length >= 1) {
            sampleGroup = tokens[tokens.length - 1];
        } else {
            sampleGroup = "NA";
        }
        long time = 1000;
        try {
            String dateString = atts.get(Header.PROCESSING_DATE.getTitle());
            java.util.Date date = (new java.text.SimpleDateFormat("yyyy-MM-dd H:mm:ss.S")).parse(dateString);
            time *= (date.getTime());
        } catch (java.text.ParseException e) {
            e.printStackTrace();
        }

        String groupId = null;
        for (String key : atts.keySet()) {
            if (key.contains("geo_group_id")) {
                groupId = escapeString(atts.get(key));
                break;
            }
        }

        String groupIdDescription = null;
        for (String key : atts.keySet()) {
            if (key.contains("geo_group_id_description")) {
                groupIdDescription = escapeString(atts.get(key));
                break;
            }
        }


        String externalName = null;
        for (String key : atts.keySet()) {
            if (key.contains("geo_tube_id")) {
                externalName = escapeString(atts.get(key));
                break;
            }
        }
        String workflowName = atts.get(Header.WORKFLOW_NAME.getTitle());
        String workflowVersion = atts.get(Header.WORKFLOW_VERSION.getTitle());

        StringBuilder sb = new StringBuilder();
        sb.append("{");

        sb.append("\"run name\":\"").append(runName).append("\",");
        sb.append("\"instrument\":\"").append(instrument).append("\",");
        sb.append("\"barcode\":\"").append(atts.get(Header.IUS_TAG.getTitle())).append("\",");
        sb.append("\"library\":\"").append(libraryName).append("\",");
        sb.append("\"sample\":\"").append(sample).append("\",");
        sb.append("\"sample group\":\"").append(sampleGroup).append("\",");
        sb.append("\"lane\":").append(atts.get(Header.LANE_NUM.getTitle())).append(",");
        sb.append("\"sequencing type\":\"").append(atts.get(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type")).append("\",");
        if (groupId != null) {
            sb.append("\"group id\":\"").append(groupId).append("\",");
        }
        if (groupIdDescription != null) {
            sb.append("\"group id description\":\"").append(groupIdDescription).append("\",");
        }
        if (externalName != null) {
            sb.append("\"external name\":\"").append(externalName).append("\",");
        }
        if (workflowName != null) {
            sb.append("\"workflow name\":\"").append(workflowName).append("\",");
        }
        if (workflowVersion != null) {
            sb.append("\"workflow version\":\"").append(workflowVersion).append("\",");
        }
        sb.append("\"last modified\":\"").append(time).append("\"");

        sb.append("}");
        Log.debug(sb.toString());
        int rand = random.nextInt();
        File file = new File(tmp + File.separator + libraryName + rand + ".json");
        writeFile(file, sb);

        Log.debug("Wrote to " + file.getAbsolutePath());

        return file.getAbsolutePath();

    }

    private void writeFile(File file, StringBuilder sb) {
        try {
            FileWriter writer = new FileWriter(file);
            writer.append(sb);
            writer.flush();
            writer.close();
            file.setReadable(true, false);

        } catch (IOException ex) {
            Log.stdout("Error writing JSON file:" + file.getAbsolutePath());
        }
    }

    public static boolean fileExistsAndIsAccessible(String filePath) {

        File file = new File(filePath);
        return (file.exists() && file.canRead() && file.isFile());

    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(BamQCDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }
}
