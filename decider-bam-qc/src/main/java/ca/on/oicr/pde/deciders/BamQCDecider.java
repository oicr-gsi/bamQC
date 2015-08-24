package ca.on.oicr.pde.deciders;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.File;
import java.util.*;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.common.util.maptools.MapTools;
import org.apache.commons.lang3.StringEscapeUtils;

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
    private String tmp = "/tmp";
    private String rsconfigXmlPath = "/.mounts/labs/PDE/data/rsconfig.xml";
    private boolean forceType = false;
    private String forcedResequencingType = "";
    private String forcedIntervalFile = "";
    private Rsconfig rs;

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
        ReturnValue rv = new ReturnValue();
        rv.setExitStatus(ReturnValue.SUCCESS);

        Log.debug("INIT");
        this.setGroupingStrategy(Header.FILE_SWA);
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
                Log.error("The given INI file does not exist: " + file.getAbsolutePath());
                rv.setExitStatus(ReturnValue.INVALIDPARAMETERS);
                return rv;
            }
        }

        if (options.has("tmp")) {
            String temp = (String) options.valueOf("tmp");
            File tempDir = new File(temp);
            if (tempDir.exists()) {
                tmp = tempDir.getAbsolutePath();
            } else {
                Log.error("The temporary directory " + tempDir.getAbsolutePath() + " does not exist.");
                rv.setExitStatus(ReturnValue.INVALIDPARAMETERS);
                return rv;
            }
        }

        if (options.has("force-type")
                && !((options.has("resequencing-type") && options.has("interval-file"))
                || (options.has("resequencing-type")))) {
            Log.error("--force-type requires either [--resequencing-type] OR [--resequencing-type and --interval-file] to be set.");
            rv.setExitStatus(ReturnValue.INVALIDPARAMETERS);
            return rv;
        }

        if (options.has("resequencing-type")) {
            if (!options.has("force-type") || !options.hasArgument("resequencing-type")) {
                Log.error("Forcing --resequencing-type requires --resequencing-type to have an argument and --force-type be set.");
                rv.setExitStatus(ReturnValue.INVALIDPARAMETERS);
                return rv;
            } else {
                forceType = true;
                forcedResequencingType = options.valueOf("resequencing-type").toString();
            }
        }

        if (options.has("interval-file")) {
            if (!options.has("resequencing-type") || !options.has("force-type")) {
                Log.error("Forcing --interval-file requires --interval-file to have an argument, --resequencing-type be set, and --force-type to be set.");
                rv.setExitStatus(ReturnValue.INVALIDPARAMETERS);
                return rv;
            }
            if (!fileExistsAndIsAccessible(options.valueOf("interval-file").toString())) {
                Log.error("The interval file is not accessible.");
                rv.setExitStatus(ReturnValue.FILENOTREADABLE);
                return rv;
            } else {
                forcedIntervalFile = options.valueOf("interval-file").toString();
            }
        }

        if (options.has("rsconfig-file")) {
            if (!options.hasArgument("rsconfig-file")) {
                Log.error("--rsconfig-file requires a file argument.");
                rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                return rv;
            }
            if (!fileExistsAndIsAccessible(options.valueOf("rsconfig-file").toString())) {
                Log.error("The rsconfig-file is not accessible.");
                rv.setExitStatus(ReturnValue.FILENOTREADABLE);
                return rv;
            } else {
                rsconfigXmlPath = options.valueOf("rsconfig-file").toString();
            }
        }

        try {
            rs = new Rsconfig(new File(rsconfigXmlPath));
        } catch (Exception e) {
            Log.error("Rsconfig file did not load properly, exeception stack trace:\n" + e.getStackTrace());
            rv.setExitStatus(ReturnValue.FAILURE);
            return rv;
        }

        if (rv.getExitStatus() != ReturnValue.SUCCESS) {
            return rv;
        }

        //allows anything defined on the command line to override the 'defaults' here.
        rv = super.init();
        return rv;

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
        String target_bed = null;

        if (forceType && !forcedResequencingType.isEmpty()) {
            Log.debug("Overriding resequencing type from [" + targetResequencingType + "] to [" + forcedResequencingType + "]");
            targetResequencingType = forcedResequencingType;
        }

        if (forceType && !forcedResequencingType.isEmpty() && !forcedIntervalFile.isEmpty()) {
            Log.debug("Overriding interval_file to [" + forcedIntervalFile + "]");
            target_bed = forcedIntervalFile;
        } else {
            target_bed = rs.get(templateType, targetResequencingType, "interval_file");
        }

        if (target_bed == null) {
            Log.error("For the file with SWID = [" + returnValue.getAttribute(Header.FILE_SWA.getTitle())
                    + "], the template type/geo_library_source_template_type = [" + templateType
                    + "] and resequencing type/geo_targeted_resequencing = [" + targetResequencingType
                    + "] could not be found in rsconfig.xml (path = [" + rsconfigXmlPath + "])");
            return false;
        }

        returnValue.setAttribute("target_bed", target_bed);

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
        try {
            iniFileMap.put("json_metadata", escapeForSeqwareIni(makeJsonSnippet(commaSeparatedFilePaths, r)));
        } catch (JsonProcessingException jpe) {
            throw new RuntimeException(jpe);
        }

        return iniFileMap;
    }

    private final ObjectMapper mapper = new ObjectMapper();

    private String makeJsonSnippet(String filePath, ReturnValue r) throws JsonProcessingException {

        FileAttributes fa = new FileAttributes(r, r.getFiles().get(0));

        Map<String, Object> j = new LinkedHashMap<>();

        String runName = fa.getSequencerRun();
        String instrument;
        String[] tokens = runName.split("_");
        if (tokens.length >= 2) {
            instrument = tokens[1];
        } else {
            instrument = "NA";
        }
        j.put("run name", runName);
        j.put("instrument", instrument);

        j.put("barcode", fa.getOtherAttribute(Header.IUS_TAG));

        String libraryName = fa.getLibrarySample();
        String sample;
        tokens = libraryName.split("_");
        if (tokens.length > 2) {
            sample = libraryName.substring(0, libraryName.lastIndexOf(tokens[tokens.length - 2]) - 1);
        } else {
            sample = "NA";
        }
        j.put("library", libraryName);
        j.put("sample", sample);

        String pSample = fa.getOtherAttribute(Header.PARENT_SAMPLE_NAME);
        String sampleGroup;
        tokens = pSample.split(":");
        if (tokens.length >= 1) {
            sampleGroup = tokens[tokens.length - 1];
        } else {
            sampleGroup = "NA";
        }
        j.put("sample group", sampleGroup);

        j.put("lane", Integer.parseInt(fa.getOtherAttribute(Header.LANE_NUM)));

        j.put("sequencing type", fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE));

        String groupId = fa.getLimsValue(Lims.GROUP_ID);
        if (groupId != null) {
            j.put("group id", escapeString(groupId));
        }

        String groupIdDescription = fa.getLimsValue(Lims.GROUP_DESC);
        if (groupIdDescription != null) {
            j.put("group id description", escapeString(groupIdDescription));
        }

        String externalName = fa.getLimsValue(Lims.TUBE_ID);
        if (externalName != null) {
            j.put("external name", escapeString(externalName));
        }

        String workflowName = fa.getOtherAttribute(Header.WORKFLOW_NAME);
        if (workflowName != null) {
            j.put("workflow name", workflowName);
        }

        String workflowVersion = fa.getOtherAttribute(Header.WORKFLOW_VERSION);
        if (workflowVersion != null) {
            j.put("workflow version", workflowVersion);
        }

        long time = 1000;
        try {
            String dateString = fa.getOtherAttribute(Header.PROCESSING_DATE);
            java.util.Date date = (new java.text.SimpleDateFormat("yyyy-MM-dd H:mm:ss.S")).parse(dateString);
            time *= (date.getTime());
        } catch (java.text.ParseException e) {
            e.printStackTrace();
        }
        j.put("last modified", Long.toString(time));

        return mapper.writeValueAsString(j);
    }

    public static String escapeForSeqwareIni(String s) {
        return StringEscapeUtils.escapeJava(StringEscapeUtils.escapeJava(s.replace("=", "&#61;")));
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
