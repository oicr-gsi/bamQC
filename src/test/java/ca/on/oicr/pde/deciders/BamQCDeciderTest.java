package ca.on.oicr.pde.deciders;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import javax.xml.parsers.ParserConfigurationException;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.apache.commons.io.FileUtils;
import org.xml.sax.SAXException;
import org.testng.Assert;
import org.testng.annotations.*;

/**
 *
 * @author rtahir
 */
public class BamQCDeciderTest {

    File filepath = FileUtils.toFile(this.getClass().getResource("/rsconfig.xml"));
   
    public BamQCDeciderTest() {
    }

    /**
     * Test of init method, of class BamQCDecider.
     */
    @Test(enabled = false)
    public void testInit() {
        System.out.println("init");
        BamQCDecider instance = new BamQCDecider();
        ReturnValue expResult = null;
        ReturnValue result = instance.init();
        Assert.assertEquals(result, expResult);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of handleGroupByAttribute method, of class BamQCDecider.
     */
    @Test(enabled = false)
    public void testHandleGroupByAttribute() {
        System.out.println("handleGroupByAttribute");
        String attribute = "";
        BamQCDecider instance = new BamQCDecider();
        String expResult = "";
        String result = instance.handleGroupByAttribute(attribute);
        Assert.assertEquals(result, expResult);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of checkFileDetails method, of class BamQCDecider.
     */
    @Test(enabled = false)
    public void testCheckFileDetails() {
        System.out.println("checkFileDetails");
        ReturnValue returnValue = null;
        FileMetadata fm = null;
        BamQCDecider instance = new BamQCDecider();
        boolean expResult = false;
        boolean result = instance.checkFileDetails(returnValue, fm);
        Assert.assertEquals(result, expResult);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    /**
     * Test of modifyIniFile method, of class BamQCDecider.
     */
    @Test(enabled = false)
    public void testModifyIniFile() {
        System.out.println("modifyIniFile");
        String commaSeparatedFilePaths = "";
        String commaSeparatedParentAccessions = "";
        BamQCDecider instance = new BamQCDecider();
        Map expResult = null;
        Map result = instance.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        Assert.assertEquals(result, expResult);
        // TODO review the generated test code and remove the default call to fail.
        Assert.fail("The test case is a prototype.");
    }

    @Test(expectedExceptions = Exception.class)
    public void configFromParsedXML_missingXMLFile() throws ParserConfigurationException, SAXException, IOException, Exception {

        BamQCDecider.configFromParsedXML("doesNotExist.xml", "");

    }

    @Test(expectedExceptions = Exception.class)
    public void configFromParsedXML_missingResequencingTag() throws ParserConfigurationException, SAXException, IOException, Exception {

        BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/missingTargetResequencingTag.xml")).getPath(), "");

    }

    @Test
    public void configFromParsedXML_getAllIntervalFilePaths() throws IOException, ParserConfigurationException, SAXException, Exception {

        Assert.assertEquals(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "HALT"), 
                "/.mounts/labs/PDE/data/TargetedSequencingQC/HALT/halt_coding.bed");
        Assert.assertEquals(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Illumina TruSeq Exome"), 
                "/.mounts/labs/PDE/data/TargetedSequencingQC/Illumina.TruSeq/TruSeq-Exome-Targeted-Regions-BED-file.bed");
        Assert.assertEquals(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Agilent SureSelect ICGC/Sanger Exon"), 
                "/.mounts/labs/PDE/data/TargetedSequencingQC/Agilent.SureSelect.ICGC/sanger.exons.hg19.bed");
        Assert.assertEquals(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Agilent SureSelect All Exon G3362"), 
                "/.mounts/labs/PDE/data/TargetedSequencingQC/Agilent.SureSelect.G3362/SureSelect_All_Exon_G3362_with_names.v2.bed");
        Assert.assertEquals(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Nimblegen 2.1M Human Exome (21191)"), 
                "/.mounts/labs/PDE/data/TargetedSequencingQC/Nimblegen.2.1M.Human.Exome/2.1M_Human_Exome.bed");
        Assert.assertEquals(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "TruSeq Amplicon - Cancer Panel"), 
                "/.mounts/labs/PDE/data/TargetedSequencingQC/TruSeqCancerPanel/TruSeq_Cancer_Panel_Targets_Regions_sorted.bed");

    }

    @Test(expectedExceptions = Exception.class)
    public void configFromParsedXML_duplicateTargetSequencingCheck() throws ParserConfigurationException, SAXException, IOException, Exception {

        BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig-duplicates-bottom.xml")).getPath(), "");

    }

    @Test
    public void configFromParsedXML_returnsNullWhenMissing() throws ParserConfigurationException, SAXException, IOException, Exception {

        Assert.assertNull(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "ThisDoesNotExist"));

    }
    
    @Test(expectedExceptions = Exception.class)
    public void configFromParsedXML_noneXmlFile() throws ParserConfigurationException, SAXException, IOException, Exception {

        BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/notAnXmlFile.txt")).getPath(), "");
        
    }
}