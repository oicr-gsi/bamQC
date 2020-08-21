FROM modulator:latest

MAINTAINER Fenglin Chen <f73chen@uwaterloo.ca>

# packages should already be set up in modulator:latest
USER root

# move in the yaml to build modulefiles from
COPY recipes/bamQC_recipe.yaml /modulator/code/gsi/recipe.yaml

# build the modules and set folder / file permissions
RUN ./build-local-code /modulator/code/gsi/recipe.yaml --initsh /usr/share/modules/init/sh --output /modules && \
	find /modules -type d -exec chmod 777 {} \; && \
	find /modules -type f -exec chmod 777 {} \;

# add the user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 ubuntu
USER ubuntu

# copy the setup file to load the modules at startup
COPY .bashrc /home/ubuntu/.bashrc

#ENV MOSDEPTH_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/mosdepth-0.2.9"
#ENV PICARD_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/picard-2.21.2"
#ENV JAVA_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/java-8"
#ENV BAM_QC_METRICS_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/bam-qc-metrics-0.2.5"
#ENV BEDTOOLS_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/bedtools-2.27"
#ENV PYTHON_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/python-3.6"
#ENV SAMTOOLS_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/samtools-1.9"
#ENV HTSLIB_ROOT="/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9"

#ENV PATH="/modules/gsi/modulator/sw/Ubuntu18.04/mosdepth-0.2.9/bin:/modules/gsi/modulator/sw/Ubuntu18.04/java-8/bin:/modules/gsi/modulator/sw/Ubuntu18.04/bam-qc-metrics-0.2.5/bin:/modules/gsi/modulator/sw/Ubuntu18.04/bedtools-2.27/bin:/modules/gsi/modulator/sw/Ubuntu18.04/python-3.6/bin:/modules/gsi/modulator/sw/Ubuntu18.04/samtools-1.9/bin:/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
#ENV MANPATH="/modules/gsi/modulator/sw/Ubuntu18.04/java-8/man:/modules/gsi/modulator/sw/Ubuntu18.04/python-3.6/share/man:/modules/gsi/modulator/sw/Ubuntu18.04/samtools-1.9/share/man:/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/share/man"
#ENV LD_LIBRARY_PATH="/modules/gsi/modulator/sw/Ubuntu18.04/java-8/lib:/modules/gsi/modulator/sw/Ubuntu18.04/bam-qc-metrics-0.2.5/lib:/modules/gsi/modulator/sw/Ubuntu18.04/python-3.6/lib:/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/lib"
#ENV PYTHONPATH="/modules/gsi/modulator/sw/Ubuntu18.04/bam-qc-metrics-0.2.5/lib/python3.6/site-packages:/modules/gsi/modulator/sw/Ubuntu18.04/python-3.6/lib/python3.6/site-packages"
#ENV PKG_CONFIG_PATH="/modules/gsi/modulator/sw/Ubuntu18.04/python-3.6/lib/pkgconfig:/modules/gsi/modulator/sw/Ubuntu18.04/htslib-1.9/lib/pkgconfig"

CMD /bin/bash
