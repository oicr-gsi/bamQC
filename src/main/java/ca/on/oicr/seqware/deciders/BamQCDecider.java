package ca.on.oicr.seqware.deciders;

import java.util.*;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.common.util.maptools.MapTools;
import net.sourceforge.seqware.pipeline.deciders.BasicDecider;

/**
 * @author mtaschuk@oicr.on.ca
 *
 */
public class BamQCDecider extends BasicDecider {

    private Map<String, ReturnValue> pathToAttributes = new HashMap<String, ReturnValue>();

    public BamQCDecider() {
        super();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setHeader(Header.IUS_SWA);
        this.setMetaType(Arrays.asList("application/bam"));


        ResourceBundle rb = PropertyResourceBundle.getBundle("decider");
        List<String> pas = Arrays.asList(rb.getString("parent-workflow-accessions").split(","));
        List<String> cwa = Arrays.asList(rb.getString("check-wf-accessions").split(","));
        this.setWorkflowAccession(rb.getString("workflow-accession"));
        this.setWorkflowAccessionsToCheck(new TreeSet(cwa));
        this.setParentWorkflowAccessions(new TreeSet(pas));

        //allows anything defined on the command line to override the 'defaults' here.
        ReturnValue val = super.init();
        return val;

    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
	Log.debug("GROUP BY ATTRIBUTE: "+getHeader().getTitle()+ " " + attribute);
        return attribute;
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
	Log.debug("CHECK FILE DETAILS:" + fm);
        pathToAttributes.put(fm.getFilePath(), returnValue);
        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
	Log.debug("INI FILE:" + commaSeparatedFilePaths);
        //Load the user-defined ini-file
        Map<String, String> iniFileMap = new TreeMap<String, String>();
        if (options.has("ini-file")) {
            MapTools.ini2Map((String) options.valueOf("ini-file"), iniFileMap, false);
        }

	StringBuilder indexFilePaths = new StringBuilder();

        for (String path : commaSeparatedFilePaths.split(",")) {
            ReturnValue rv = pathToAttributes.get(path);

            iniFileMap.put("run_name", rv.getAttribute(Header.SEQUENCER_RUN_NAME.getTitle()));
            iniFileMap.put("sample_name", rv.getAttribute(Header.SAMPLE_NAME.getTitle()));
            iniFileMap.put("sample_group", rv.getAttribute(" "));
            iniFileMap.put("library", handleGroupByAttribute(rv.getAttribute(Header.PARENT_SAMPLE_NAME.getTitle())));
            iniFileMap.put("experimental_design", getDesign(rv));
            iniFileMap.put("sequencer_run", rv.getAttribute(Header.SEQUENCER_RUN_NAME.getTitle()));
            iniFileMap.put("lane", rv.getAttribute(Header.LANE_NUM.getTitle()));
            iniFileMap.put("barcode", rv.getAttribute(Header.IUS_TAG.getTitle()));
            iniFileMap.put("last_modified", rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
	    if (indexFilePaths.length()>0) {
		indexFilePaths.append(",");
	    }
	    indexFilePaths.append(getIndexFile(rv));
        }

	iniFileMap.put("bam_indexes", indexFilePaths.toString());
	iniFileMap.put("bam_inputs", commaSeparatedFilePaths);

        return iniFileMap;
    }

    private String getIndexFile(ReturnValue rv)
    {
	List<FileMetadata> files = rv.getFiles();
	for (FileMetadata fm: files) {
	    if (fm.getMetaType().equals("application/bam-index")) {
		return fm.getFilePath();
	    }
	}	
	return "NO_INDEX_FOUND";
    }



    private String getDesign(ReturnValue ret) {
        StringBuilder design = new StringBuilder();
        for (String key : ret.getAttributes().keySet()) {
            if (key.contains("geo_library_source_template_type")) {
                if (design.length() > 0) {
                    design.append(",");
                }
                design.append(ret.getAttribute(key));
            }
        }
        return design.toString();
    }
}
