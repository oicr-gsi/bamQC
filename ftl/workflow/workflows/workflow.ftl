<?xml version="1.0" encoding="UTF-8"?>
<adag xmlns="http://pegasus.isi.edu/schema/DAX" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://pegasus.isi.edu/schema/DAX http://pegasus.isi.edu/schema/dax-3.2.xsd" version="3.2" count="1" index="0" name="${workflow-name}">

<!-- the directory structure inside the bundle -->
<#assign workflow_name = "Workflow_Bundle_${workflow-directory-name}/${version}"/>

<!-- MACRO: to create a mkdir pre job and stage mkdir binary -->
<#macro requires_dir dir>
  <profile namespace="env" key="GRIDSTART_PREJOB">/${workflow_bundle_dir}/${workflow_name}/bin/globus/pegasus-dirmanager -c -d ${dir}</profile>
</#macro>

<!-- VARS -->
<#-- workflow and seqware versions -->
<#assign seqware_version = "${seqware-version}"/>
<#assign workflow_version = "${workflow-version}"/>
<#assign java_version = "1.6.0"/>
<#assign perl_version = "5.14.1"/>
<#assign perl_bin = "${workflow_bundle_dir}/${workflow_name}/bin/perl-5.14.1/perl" />
<#-- make sure it is a string -->
<#assign parentAccessions = "${parent_accessions}"/>
<#-- Set relative paths for files within the run-->
<#assign bin_dir = "bin"/>
<#assign data_dir = "data"/>
<#assign lib_dir = "lib"/>
<#assign dep_dir = "dependencies"/>
<#if output_path == "NA" >
 <#assign output_path="${output_dir}/seqware-${seqware_version}_${workflow_name}/${random}"/>
</#if>

<!-- BASE FILE NAME -->
<!-- Set the basename from input file name -->
<#list input_file?split("/") as tmp>
  <#assign basename = tmp/>
</#list>

<!-- EXECUTABLES INCLUDED WITH BUNDLE -->
<executable namespace="seqware" name="java" version="${java_version}" 
            arch="x86_64" os="linux" installed="true" >
  <!-- the path to the tool that actually runs a given module -->
  <pfn url="file:///${workflow_bundle_dir}/${workflow_name}/bin/jre1.6.0_29/bin/java" site="${seqware_cluster}"/>
</executable>

<executable namespace="seqware" name="perl" version="${perl_version}" 
            arch="x86_64" os="linux" installed="true" >
  <!-- the path to the tool that actually runs a given module -->
  <pfn url="file:///${workflow_bundle_dir}/${workflow_name}/bin/perl-5.14.1/perl" site="${seqware_cluster}"/>
</executable>

<executable namespace="pegasus" name="dirmanager" version="1" 
            arch="x86_64" os="linux" installed="true" >
  <!-- the path to the tool that actually runs a given module -->
  <pfn url="file:///${workflow_bundle_dir}/${workflow_name}/bin/globus/pegasus-dirmanager" site="${seqware_cluster}"/>     
</executable>


<!-- Part 1: Define all jobs -->

  <!-- Provision input file, make link if local file, copy from HTTP/S3 otherwise -->
  <job id="IDPRE1" namespace="seqware" name="java" version="${java_version}">
    <argument>
      -Xmx1000M
      -classpath ${workflow_bundle_dir}/${workflow_name}/classes:${workflow_bundle_dir}/${workflow_name}/lib/seqware-distribution-${seqware_version}-full.jar
      net.sourceforge.seqware.pipeline.runner.Runner
      --no-metadata
      --module net.sourceforge.seqware.pipeline.modules.utilities.ProvisionFiles
      --
      --input-file ${input_file}
      --output-dir ${data_dir}
    </argument>

    <profile namespace="globus" key="jobtype">condor</profile>
    <profile namespace="globus" key="count">1</profile>
     <profile namespace="globus" key="queue">${queue}</profile>
   <profile namespace="globus" key="maxmemory">2000</profile>
</job>

  <!-- The job that converts sam to bam and pipes it to samStats.pl, this uses the GenericCommandRunner to just use a tool on the command line. -->
  <#assign algo = "bamToJsonStats"/>
  <job id="ID001" namespace="seqware" name="java" version="${java_version}">
    <argument>
      -Xmx1000M
      -classpath ${workflow_bundle_dir}/${workflow_name}/bin:${workflow_bundle_dir}/${workflow_name}/lib/seqware-distribution-${seqware_version}-full.jar
      net.sourceforge.seqware.pipeline.runner.Runner
      --${metadata}
      <#list parentAccessions?split(",") as pa>
        --metadata-parent-accession ${pa}
      </#list>
      --metadata-processing-accession-file ${data_dir}/${algo}_accession
      --metadata-output-file-prefix ${output_prefix}${output_path}
      --metadata-workflow-run-accession ${workflow_run_accession}
      --module net.sourceforge.seqware.pipeline.modules.GenericCommandRunner
      --
      --gcr-algorithm ${algo}
      --gcr-command ${workflow_bundle_dir}/${workflow_name}/${bin_dir}/samtools-0.1.17/samtools view ${input_file} | 
            ${perl_bin} ${workflow_bundle_dir}/${workflow_name}/${dep_dir}/samStats.pl 
            -s ${sample_rate} 
            -i ${normal_insert_max} 
            -q ${map_qual_cut} 
            -r ${target_bed} 
            -j ${json_metadata_file} 
            > ${data_dir}/${basename}.BamQC.json
    </argument>

    <profile namespace="globus" key="jobtype">condor</profile>
    <profile namespace="globus" key="count">1</profile>
    <profile namespace="globus" key="maxmemory">2000</profile>
    <profile namespace="globus" key="queue">${queue}</profile>

  </job>


  <!-- Provision output file to either local dir or S3 -->
  <#assign parentAlgo = "bamToJsonStats" />
  <#assign algo = "provideJsonStats"/>
  <job id="IDPOST1" namespace="seqware" name="java" version="${java_version}">
    <argument>
      -Xmx1000M
      -classpath ${workflow_bundle_dir}/${workflow_name}/classes:${workflow_bundle_dir}/${workflow_name}/lib/seqware-distribution-${seqware_version}-full.jar
      net.sourceforge.seqware.pipeline.runner.Runner
      --${metadata}
      --metadata-processing-accession-file ${data_dir}/${algo}_accession
      --metadata-parent-accession-file ${data_dir}/${parentAlgo}_accession
      --metadata-output-file-prefix ${output_prefix}${output_path}
      --metadata-workflow-run-ancestor-accession ${workflow_run_accession}
      --module net.sourceforge.seqware.pipeline.modules.utilities.ProvisionFiles
      --
      --force-copy
      --input-file ${data_dir}/${basename}.BamQC.json
      --input-file-metadata ${algo}::text/json::${data_dir}/${basename}.BamQC.json
      --output-dir ${output_prefix}${output_path}
    </argument>

    <profile namespace="globus" key="jobtype">condor</profile>
    <profile namespace="globus" key="count">1</profile>
    <profile namespace="globus" key="maxmemory">2000</profile>
    <profile namespace="globus" key="queue">${queue}</profile>

  </job>

<!-- End of Job Definitions -->

<!-- Part 2: list of control-flow dependencies -->

  <!-- Define task group dependencies -->
  <child ref="ID001">
     <parent ref="IDPRE1"/>
   </child>
  <child ref="IDPOST1">
    <parent ref="ID001"/>
  </child>

<!-- End of Dependencies -->

</adag>
