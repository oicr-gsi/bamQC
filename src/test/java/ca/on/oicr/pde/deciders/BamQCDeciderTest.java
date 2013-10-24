/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;
import org.junit.Rule;
import org.junit.rules.ExpectedException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
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

    /**
     * Test of configFromParsedXML method, of class BamQCDecider.
     */
//    @Test
//    @Ignore
//    public void testConfigFromParsedXML() {
//        System.out.println("configFromParsedXML");
//        String fileName = "";
//        String resequencingType = "";
//        BamQCDecider instance = new BamQCDecider();
//        String expResult = "";
//        String result = instance.configFromParsedXML(fileName, resequencingType);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
    /**
     *
     */
    @Test
    @Ignore
      
    public void XMLfileExists_ConfigFromParsedXML() throws ParserConfigurationException, SAXException, IOException, DuplicatesException {
        
       exception.expect(IOException.class);
        
        String e = BamQCDecider.configFromParsedXML(
                (FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig-doesNotExist.xml")).getPath() == null) ? "" 
                : BamQCDeciderTest.class.getResource("/rsconfig-doesNotExist.xml").getPath(),
                "random");
      
    }

    @Test
    @Ignore

   
    public void targetSequencingExists_ConfigFromParsedXML() throws ParserConfigurationException, SAXException, IOException, DuplicatesException {

        exception.expect(NullPointerException.class);
         
         System.out.println((BamQCDeciderTest.class.getResource("/notargetSequencing.xml")).getPath());
         BamQCDecider.configFromParsedXML(
                 FileUtils.toFile(BamQCDeciderTest.class.getResource("/notargetSequencing.xml")).getPath(), 
                 "HALT");

        
        
//        
//        File fXmlFile = FileUtils.toFile(this.getClass().getResource("/rsconfig.xml"));
//        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
//        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
//        Document doc = dBuilder.parse(fXmlFile);
//        doc.getDocumentElement().normalize();
//        NodeList nList = doc.getElementsByTagName("resequencing_type");
//        assertTrue(nList.getLength() > 0);


    }

    @Test
      @Ignore
    public void targetSequencingMatchesExpectations_ConfigFromParsedXML() {
        
        
        
        
        try {
            File fXmlFile = FileUtils.toFile(this.getClass().getResource("/rsconfig.xml"));
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);
            doc.getDocumentElement().normalize();
            NodeList nList = doc.getElementsByTagName("resequencing_type");


            for (int temp = 0; temp < nList.getLength(); temp++) {

                Node nNode = nList.item(temp);

                if (nNode.getNodeType() == Node.ELEMENT_NODE) {

                    Element eElement = (Element) nNode;
                    assertTrue((eElement.getAttribute("id").equals("HALT"))
                            || (eElement.getAttribute("id").equals("Illumina TruSeq Exome"))
                            || (eElement.getAttribute("id").equals("Agilent SureSelect ICGC/Sanger Exon"))
                            || (eElement.getAttribute("id").equals("Agilent SureSelect All Exon G3362"))
                            || (eElement.getAttribute("id").equals("Nimblegen 2.1M Human Exome (21191)"))
                            || (eElement.getAttribute("id").equals("TruSeq Amplicon - Cancer Panel")));
                }
            }

        } catch (Exception e) {
        }

    }
    
    @Rule
	public ExpectedException exception = ExpectedException.none();

    @Test
    public void duplicateTargetSequencingCheck_configFromParsedXML() 
            throws ParserConfigurationException, SAXException, IOException, DuplicatesException {
         
        exception.expect(DuplicatesException.class);
         
         System.out.println((BamQCDeciderTest.class.getResource("/rsconfig-duplicates-bottom.xml")).getPath());
         BamQCDecider.configFromParsedXML(
                 FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig-duplicates-bottom.xml")).getPath(), 
                 "TruSeq Amplicon - Cancer Panel");
    }

    @Test
    @Ignore
    public void checksIfDesiredTargetSequencingExists_inXMLFile() throws ParserConfigurationException, SAXException, IOException, DuplicatesException {
       
        String returnVal = BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "TruSeq Amplicon - Cancer Panel");
        assertFalse(returnVal == null);
    }

    @Test
    @Ignore
      
    public void checksIftargetBEDexists_inXMLFile() throws ParserConfigurationException, SAXException, IOException, DuplicatesException {
       
        System.out.println(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath());
        String returnVal = BamQCDecider.configFromParsedXML(FileUtils.toFile(BamQCDeciderTest.class.getResource("/rsconfig.xml")).getPath(), "HALT");
        assertTrue(returnVal.equals("/.mounts/labs/PDE/data/TargetedSequencingQC/HALT/halt_coding.bed"));
    }
}