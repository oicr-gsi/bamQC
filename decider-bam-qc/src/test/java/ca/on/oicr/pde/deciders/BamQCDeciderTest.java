package ca.on.oicr.pde.deciders;

import java.io.File;
import java.util.Arrays;
import java.util.Map;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.apache.commons.io.FileUtils;
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

    private void verifyFailureOnInit(String... args) {
        BamQCDecider d = new BamQCDecider();
        d.setParams(Arrays.asList(args));

        ReturnValue parseReturnValue = d.parse_parameters();

        ReturnValue initReturnValue = new ReturnValue();
        if (parseReturnValue.getExitStatus() == ReturnValue.SUCCESS) {
            initReturnValue = d.init();
        }

        Assert.assertTrue(parseReturnValue.getExitStatus() != ReturnValue.SUCCESS || initReturnValue.getExitStatus() != ReturnValue.SUCCESS);
    }

    @Test
    public void invalidParams_missingParametersToForceType() {
        verifyFailureOnInit("--force-type");
    }

    @Test
    public void invalidParams_missingRsconfigFileArgument() {
        verifyFailureOnInit("--rsconfig-file");
    }

    @Test
    public void invalidParams_invalidRsconfigFile() {
        verifyFailureOnInit("--rsconfig-file", "/tmp/does/not/exist");
    }

    @Test
    public void invalidParams_missingResequencingType() {
        verifyFailureOnInit("--force-type", "--resequencing-type");
    }

    @Test
    public void invalidParams_missingIntervalFile() {
        verifyFailureOnInit("--force-type", "--resequencing-type", "EX", "--interval-file");
    }

    @Test
    public void invalidParams_missingForceType() {
        verifyFailureOnInit("--interval-file", "...");
        verifyFailureOnInit("--resequencing-type", "...");
        verifyFailureOnInit("--resequencing-type", "...", "--interval-file", "...");
        verifyFailureOnInit("--force-type", "--interval-file", "...");
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

}