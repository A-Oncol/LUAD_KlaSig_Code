#!/bin/bash
#Title: cleanfastq_qc
#Author: Ao Shen
#Date: 2024-10-22
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"


# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
dataset=3.fastp
id=`cat ${project}/sample.list | sed -n "${i}p"`
biosoft=fastqc
dir=${project}/4.clean_qc/1.fastqc
fastqc=/home/yjx/anaconda3/envs/chip-seq/bin/fastqc

# FastQC
if [ ! -f ${dir}/ok.${id}_${biosoft}.status ]
then
    echo "FastQC ${i} for ${id} started at:" `date`
    ${fastqc} --outdir ${dir} --threads 6  ${project}/${dataset}/${id}_*.fq.gz > ${dir}/${id}_${biosoft}.log 2>&1
    #echo "FastQC for ${id} has finished." `date`

    if [ $? -eq 0 ]
    then
         echo "${biosoft} succeed for ${id}." `date`
         touch ${dir}/ok.${id}_${biosoft}.status
    else
         echo "${biosoft} failed for ${id}" `date`
    fi
fi

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Process will sleep 10s"
sleep 10
echo "Process end at:" `date`
