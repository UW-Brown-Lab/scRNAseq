#!/bin/bash
#
# brown-fastq-fetch.sh
# Retrieve fastq file via fetchngs
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****SET HOME DIR*****
#
export HOME="$_CONDOR_SCRATCH_DIR"
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****INSTALL JAVA*****
#
tar -xzf amazon-corretto-17-x64-linux-jdk.tar.gz
# add unzipped JDK folder to env
export PATH=$PWD/amazon-corretto-17.0.10.7.1-linux-x64/bin:$PATH
export JAVA_HOME=$PWD/amazon-corretto-17.0.10.7.1-linux-x64
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****INSTALL APPTAINER*****
#
curl -s https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh | \
    bash -s - apptainer_dir
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****INSTALL NEXTFLOW*****
wget -qO- https://get.nextflow.io | bash
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****RUN FETCHNGS*****
./nextflow run nf-core/fetchngs -profile apptainer --input sra_list.csv --outdir ./fastq-output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****ZIP OUTPUT*****
#
tar -czvf fastq_raw_reads.tar.gz ./fastq-output
#
#
