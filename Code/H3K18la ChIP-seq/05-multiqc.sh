#!/bin/bash
#Title: clean_multi_qc
#Author: Ao Shen
#Date: 2024-10-22
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
project=/home/yjx/projects/kla/PRJNA857271
dataset=rawdata
biosoft=multiqc
dir=${project}/4.clean_qc/2.multiqc
multiqc=/home/yjx/anaconda3/envs/chip-seq/bin/multiqc

# MultiQC
if [ ! -f ${dir}/ok.${biosoft}.status ]
then
    echo -e "MultiQC started at:" `date`
    ${multiqc} ${project}/4.clean_qc/1.fastqc/*.zip -o ${dir}/
    #echo "MultiQC has finished at:" `date`

	  if [ $? -eq 0 ]
          then
	         echo "${biosoft} succeed." `date`
                 touch ${dir}/ok.${biosoft}.status
	  else
	         echo "${biosoft} failed." `date`
	  fi
fi

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Process will sleep 10s."
sleep 10
echo "Process end at:" `date`
