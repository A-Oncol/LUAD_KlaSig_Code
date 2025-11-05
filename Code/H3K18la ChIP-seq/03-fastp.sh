#!/bin/bash
#Title: fastp
#Author: Ao Shen
#Date: 2024-10-22
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
id=`cat ${project}/sample.list | sed -n "${i}p"`
dataset=rawdata
biosoft=fastp
dir=${project}/3.fastp
fastp=/home/yjx/anaconda3/envs/chip-seq/bin/fastp

# fastp
if [ ! -f ${dir}/ok.${id}_${biosoft}.status ]
then
    echo "start ${i} fastp for ${id}." `date`
    ${fastp} -i ${project}/1.${dataset}/${id}_1.fastq.gz -I ${project}/1.${dataset}/${id}_2.fastq.gz \
	    -o ${dir}/${id}_1.fq.gz -O ${dir}/${id}_2.fq.gz \
	    -l 36 -q 28 -w 10 --compression=4 -R ${id} \
	    -h ${dir}/${id}_fastp.html \
	    -j ${dir}/${id}_fastp.json \
	    1>${dir}/${id}_fastp.log 2>&1

	  if [ $? -eq 0 ]
          then
	         echo "${biosoft} succeed." `date`
                 touch ${dir}/ok.${id}_${biosoft}.status
	  else
	         echo "${biosoft} failed." `date`
	  fi
fi

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Process will sleep 10s."
sleep 10
echo "Process end at:" `date`
