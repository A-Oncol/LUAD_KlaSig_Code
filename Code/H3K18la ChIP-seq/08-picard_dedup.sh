#!/bin/bash
#Title: Picard MarkDuplicates
#Author: Ao Shen
#Date: 2024-10-24
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
id=`cat ${project}/sample.list | sed -n "${i}p"`
biosoft=picard
dir=${project}/5.align
picard=/home/yjx/anaconda3/envs/chip-seq/bin/picard
samtools=/home/yjx/anaconda3/envs/chip-seq/bin/samtools

# Picard MarkDuplicates
if [ ! -f ${dir}/ok.${id}_${biosoft}.status ]
then
    echo -e "Process started at:" `date`
    ${picard} AddOrReplaceReadGroups \
	      -I ${dir}/${id}_sorted.bam \
	      -O ${dir}/${id}_sorted_addRG.bam \
	      --RGID group1 \
	      --RGLB lib1 \
	      --RGPL illumina \
	      --RGPU unit1 \
	      --RGSM ${id}

    ${picard} MarkDuplicates -I ${dir}/${id}_sorted_addRG.bam \
	      -M ${dir}/${id}_picard.mat \
	      -O ${dir}/${id}_picard_dedup.bam \
	      --REMOVE_DUPLICATES true \
	      --READ_NAME_REGEX null

    ${samtools} flagstat -@ 8 ${dir}/${id}_picard_dedup.bam > ${dir}/${id}_picard_dedup_flagstat.txt
    ${samtools} index -@ 8 -b ${dir}/${id}_picard_dedup.bam
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
