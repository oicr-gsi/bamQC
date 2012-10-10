package ca.on.oicr.seqware.deciders;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.deciders.BasicDecider;

/**
 * @author mtaschuk@oicr.on.ca
 *
 */
public class BamQCDecider extends BasicDecider {

    private Map<String, ReturnValue> pathToAttributes = new HashMap<String, ReturnValue>();
    private String sampleRate = "1000";
    private String normalInsertMax = "1500";
    private String mapQualCut = "30";

    public BamQCDecider() {
        super();
        parser.accepts("sample_rate", "Optional. Set the sampling rate in the decider "
                + "(only 1/sample_rate reads are used for memory-intensive sampling) "
                + "This needs to be set lower for mate-pair libraries. Default: 1000").withRequiredArg();
        parser.accepts("normal_insert_max", "Optional. Set the maximum insert size "
                + "to prevent skewing of insert statistics. Default:1500").withRequiredArg();
        parser.accepts("map_qual_cut", "Optional. Set the mapQ value (quality of the "
                + "alignment of the read to the reference). Default: 30").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setHeader(Header.IUS_SWA);
        this.setMetaType(Arrays.asList("application/bam"));


        if (options.has("sample_rate")) {
            sampleRate = options.valueOf("sample_rate").toString();
        }
        if (options.has("normal_insert_max")) {
            normalInsertMax = options.valueOf("normal_insert_max").toString();
        }
        if (options.has("map_qual_cut")) {
            mapQualCut = options.valueOf("map_qual_cut").toString();
        }

        //allows anything defined on the command line to override the 'defaults' here.
        ReturnValue val = super.init();
        return val;

    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        Log.debug("GROUP BY ATTRIBUTE: " + getHeader().getTitle() + " " + attribute);
        return attribute;
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);

        String templateType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX + "geo_library_source_template_type");
        if ("WG".equals(templateType)) {
            returnValue.setAttribute("target_bed", "/oicr/data/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19_random/hg19_random.genome.sizes.bed");
        } else if ("EX".equals(templateType)) {
            String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX + "sample.geo_targeted_resequencing");
            if ("Illumina TruSeq Exome".equals(targetResequencingType)) {
                returnValue.setAttribute("target_bed", "/oicr/data/reference/genomes/homo_sapiens_mc/TruSeq/TruSeq-Exome-Targeted-Regions-BED-file");
            } else if ("Agilent SureSelect ICGC/Sanger Exon".equals(targetResequencingType)) {
                returnValue.setAttribute("target_bed", "/oicr/data/genomes/homo_sapiens/Agilent/SureSelect_Whole_Exome_ICGC_Sanger/GRCh37hg19/sanger.exons.bed.hg19");
            } else if ("Nimblegen 2.1M Human Exome (21191)".equals(targetResequencingType)) {
                returnValue.setAttribute("target_bed", "/oicr/data/reference/genomes/homo_sapiens/Nimblegen/2.1M_Human_Exome_Annotation_21191/hg19/080904_ccds_exome_rebalfocus_hg19/processed/2.1M_Human_Exome.bed");
            } else {
                Log.error("The targeted resequencing type does not have an associated BED file. Modify the decider to include this type:" + targetResequencingType);
                return false;
            }
        } else {
            Log.error("This template type is not supported for the BAM QC decider: " + templateType);
            return false;
        }
        pathToAttributes.put(fm.getFilePath(), returnValue);
        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE:" + commaSeparatedFilePaths);
        Map<String, String> iniFileMap = new TreeMap<String, String>();

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


        StringBuilder sb = new StringBuilder();
        sb.append("{");

        sb.append("\n\t\"run name\":\"").append(runName).append("\"");
        sb.append("\n\t\"instrument\":\"").append(instrument).append("\"");
        sb.append("\n\t\"barcode\":\"").append(atts.get(Header.IUS_TAG.getTitle())).append("\"");
        sb.append("\n\t\"library\":\"").append(libraryName).append("\"");
        sb.append("\n\t\"sample\":\"").append(sample).append("\"");
        sb.append("\n\t\"sample group\":\"").append(sampleGroup).append("\"");
        sb.append("\n\t\"lane\":\"").append(atts.get(Header.LANE_NUM.getTitle())).append("\"");
        sb.append("\n\t\"sequencing type\":\"").append(atts.get(Header.SAMPLE_TAG_PREFIX + "geo_library_source_template_type")).append("\"");
        sb.append("\n\t\"last modified\":\"").append(atts.get(Header.PROCESSING_DATE.getTitle())).append("\"");

        sb.append("\n}");


        File file = new File(filePath.substring(0, filePath.lastIndexOf("/") + 1) + libraryName + ".json");
        try {
            FileWriter writer = new FileWriter(file);
            writer.append(sb);
            writer.flush();
            writer.close();

        } catch (IOException ex) {
            Log.error("Error writing JSON file:" + file.getAbsolutePath(), ex);
        }
        return file.getAbsolutePath();

    }
}
