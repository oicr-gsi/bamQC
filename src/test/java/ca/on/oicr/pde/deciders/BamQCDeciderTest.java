package ca.on.oicr.pde.deciders;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import javax.xml.parsers.ParserConfigurationException;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.rules.ExpectedException;
import org.xml.sax.SAXException;

/**
 *
 * @author rtahir
 */
public class BamQCDeciderTest {

    File filepath = FileUtils.toFile(this.getClass().getResource("/rsconfig.xml"));

    public BamQCDeciderTest() {
    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of init method, of class BamQCDecider.
     */
    @Test
    @Ignore
    public void testInit() {
        System.out.println("init");
        BamQCDecider instance = new BamQCDecider();
        ReturnValue expResult = null;
        ReturnValue result = instance.init();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of handleGroupByAttribute method, of class BamQCDecider.
     */
    @Test
    @Ignore
    public void testHandleGroupByAttribute() {
        System.out.println("handleGroupByAttribute");
        String attribute = "";
        BamQCDecider instance = new BamQCDecider();
        String expResult = "";
        String result = instance.handleGroupByAttribute(attribute);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of checkFileDetails method, of class BamQCDecider.
     */
    @Test
    @Ignore
    public void testCheckFileDetails() {
        System.out.println("checkFileDetails");
        ReturnValue returnValue = null;
        FileMetadata fm = null;
        BamQCDecider instance = new BamQCDecider();
        boolean expResult = false;
        boolean result = instance.checkFileDetails(returnValue, fm);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of modifyIniFile method, of class BamQCDecider.
     */
    @Test
    @Ignore
    public void testModifyIniFile() {
        System.out.println("modifyIniFile");
        String commaSeparatedFilePaths = "";
        String commaSeparatedParentAccessions = "";
        BamQCDecider instance = new BamQCDecider();
        Map expResult = null;
        Map result = instance.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
    @Rule
    public ExpectedException exception = ExpectedException.none();

    @Test
    public void configFromParsedXML_missingXMLFile() throws ParserConfigurationException, SAXException, IOException, Exception {

        exception.expect(Exception.class);
        BamQCDecider.configFromParsedXML("doesNotExist.xml", "");

    }

    @Test
    public void configFromParsedXML_missingResequencingTag() throws ParserConfigurationException, SAXException, IOException, Exception {

        exception.expect(Exception.class);
        BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/missingTargetResequencingTag.xml")).getPath(), "");

    }

    @Test
    public void configFromParsedXML_getAllIntervalFilePaths() throws IOException, ParserConfigurationException, SAXException, Exception {

        Assert.assertEquals("/.mounts/labs/PDE/data/TargetedSequencingQC/HALT/halt_coding.bed",
                BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "HALT"));
        Assert.assertEquals("/.mounts/labs/PDE/data/TargetedSequencingQC/Illumina.TruSeq/TruSeq-Exome-Targeted-Regions-BED-file.bed",
                BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Illumina TruSeq Exome"));
        Assert.assertEquals("/.mounts/labs/PDE/data/TargetedSequencingQC/Agilent.SureSelect.ICGC/sanger.exons.hg19.bed",
                BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Agilent SureSelect ICGC/Sanger Exon"));
        Assert.assertEquals("/.mounts/labs/PDE/data/TargetedSequencingQC/Agilent.SureSelect.G3362/SureSelect_All_Exon_G3362_with_names.v2.bed",
                BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Agilent SureSelect All Exon G3362"));
        Assert.assertEquals("/.mounts/labs/PDE/data/TargetedSequencingQC/Nimblegen.2.1M.Human.Exome/2.1M_Human_Exome.bed",
                BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "Nimblegen 2.1M Human Exome (21191)"));
        Assert.assertEquals("/.mounts/labs/PDE/data/TargetedSequencingQC/TruSeqCancerPanel/TruSeq_Cancer_Panel_Targets_Regions_sorted.bed",
                BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "TruSeq Amplicon - Cancer Panel"));

    }

    @Test
    public void configFromParsedXML_duplicateTargetSequencingCheck() throws ParserConfigurationException, SAXException, IOException, Exception {

        exception.expect(Exception.class);
        BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig-duplicates-bottom.xml")).getPath(), "");

    }

    @Test
    public void configFromParsedXML_returnsNullWhenMissing() throws ParserConfigurationException, SAXException, IOException, Exception {

        Assert.assertNull(BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "ThisDoesNotExist"));

    }
    
    @Test
    public void configFromParsedXML_noneXmlFile() throws ParserConfigurationException, SAXException, IOException, Exception {
        
        exception.expect(SAXException.class);
        BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/notAnXmlFile.txt")).getPath(), "");
        
    }
}