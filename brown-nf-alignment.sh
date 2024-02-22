#!/bin/bash
#
# brown-nf-alignment.sh
# FASTQ processing and alignment of single-cell data using Nextflow nf-core/scrnaseq
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****UPDATE HOME DIR*****
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
# *****INSTALL NEXTFLOW*****
#
wget -qO- https://get.nextflow.io | bash
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****UNZIP FETCHNGS OUTPUT*****
#
tar -xzf fetchngs_CMs.tar.gz
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****DELETE THIS LATER*****
# *****RENAME RUN FILES*****
#
mv ./fastq-output/fastq/SRX14149591_SRR17993955_1.fastq.gz ./fastq-output/fastq/SRX14149591_SRR17993955_R1.fastq.gz
mv ./fastq-output/fastq/SRX14149591_SRR17993955_2.fastq.gz ./fastq-output/fastq/SRX14149591_SRR17993955_R2.fastq.gz
mv ./fastq-output/fastq/SRX14149589_SRR17993957_1.fastq.gz ./fastq-output/fastq/SRX14149589_SRR17993957_R1.fastq.gz
mv ./fastq-output/fastq/SRX14149589_SRR17993957_2.fastq.gz ./fastq-output/fastq/SRX14149589_SRR17993957_R2.fastq.gz
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****RUN SCRNASEQ PIPELINE*****
#
./nextflow run nf-core/scrnaseq -profile apptainer -c ./nextflow.config  -params-file params.json
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# *****ZIP OUTPUT*****
#
tar -czvf nf_output.tar.gz ./nf_output



