#!/bin/bash
#Title: bowtie2_align
#Author: Ao Shen
#Date: 2024-10-22
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
id=`cat ${project}/sample.list | sed -n "${i}p"`
biosoft=bowtie2
dir=${project}/5.align
bt2=/home/yjx/anaconda3/envs/chip-seq/bin/bowtie2
bt2_idx=/home/yjx/database/human/hg38/index/bowtie2_index/hg38

# alignment
if [ ! -f ${dir}/ok.${id}_${biosoft}.status ]
then
    echo -e "Bowtie2 alignment started at:" `date`
    ${bt2} -q -p 8 \
	    -x ${bt2_idx} \
	    -1 ${project}/3.fastp/${id}_1.fq.gz \
	    -2 ${project}/3.fastp/${id}_2.fq.gz \
	    -S ${dir}/${id}.sam \
	    > ${dir}/${id}_${biosoft}.log 2>&1
    #echo "Bowtie2 alignment has finished at:" `date`

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
