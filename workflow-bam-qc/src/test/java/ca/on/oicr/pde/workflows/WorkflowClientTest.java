package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.testing.workflow.DryRun;
import ca.on.oicr.pde.testing.workflow.TestDefinition;
import java.io.File;
import java.io.IOException;
import org.apache.commons.io.FileUtils;
import org.testng.annotations.Test;

/**
 *
 * @author mlaszloffy
 */
public class WorkflowClientTest {

    public WorkflowClientTest() {
    }

    @Test
    public void testSomeMethod() throws IllegalAccessException, InstantiationException, IOException, Exception {
        TestDefinition td = TestDefinition.buildFromJson(FileUtils.readFileToString(new File("src/test/resources/developmentTests.json")));
        for (TestDefinition.Test t : td.getTests()) {
            DryRun d = new DryRun(System.getProperty("bundleDirectory"), t.getParameters(), WorkflowClient.class);
            d.buildWorkflowModel();
            d.validateWorkflow();
        }
    }
}
