#!/bin/bash
#Title: sam2bam bamSort bamIndex
#Author: Ao Shen
#Date: 2024-10-23
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
id=`cat ${project}/sample.list | sed -n "${i}p"`
biosoft=samtools
dir=${project}/5.align
samtools=/home/yjx/anaconda3/envs/chip-seq/bin/samtools

# samtools
if [ ! -f ${dir}/ok.${id}_${biosoft}.status ]
then
    echo -e "Process started at:" `date`
    ${samtools} view -@ 8 -q 20 -S ${dir}/${id}.sam -b -o ${dir}/${id}.bam
    ${samtools} sort -@ 8 -l 5 -o ${dir}/${id}_sorted.bam ${dir}/${id}.bam
    ${samtools} index -@ 8 -b ${dir}/${id}_sorted.bam
    ${samtools} flagstat -@ 8 ${dir}/${id}_sorted.bam > ${dir}/${id}_flagstat.txt
    ${samtools} markdup -@ 8 -r ${dir}/${id}_sorted.bam ${dir}/${id}_dedup.bam
    ${samtools} flagstat -@ 8 ${dir}/${id}_dedup.bam > ${dir}/${id}_dedup_flagstat.txt
    ${samtools} index -@ 8 -b ${dir}/${id}_dedup.bam
    #echo "Process has finished at:" `date`

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
