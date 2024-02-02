#!/bin/bash
#
# brown-nf-alignment.sh
# FASTQ processing and alignment of single-cell data using Nextflow nf-core/scrnaseq
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
# *****INSTALL NEXTFLOW*****
#
chmod +x ./nextflow.sh
./nextflow.sh
#
# Test Nextflow install:
./nextflow run hello

