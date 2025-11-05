#!/bin/bash
#Title: macs2 callpeak
#Author: Ao Shen
#Date: 2024-10-24
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
ctrl_id=`cat ${project}/data.info | sed -n "${i}p" | cut -f 2`
treat_id=`cat ${project}/data.info | sed -n "${i}p" | cut -f 3`
ctr_bam=${project}/5.align/${ctrl_id}_picard_dedup.bam
treat_bam=${project}/5.align/${treat_id}_picard_dedup.bam
outname=`cat ${project}/data.info | sed -n "${i}p" | cut -f 1`
biosoft=macs2
dir=${project}/6.macs2
macs2=/home/yjx/anaconda3/envs/chip-seq/bin/macs2

# macs2 callpeak
if [ ! -f ${dir}/ok.${outname}_${biosoft}.status ]
then
    echo -e "Process started at:" `date`
    ${macs2} callpeak -c ${ctr_bam} -t ${treat_bam} -f BAM -g hs -q 0.0001 \
	    -n ${outname} --outdir ${dir}
    #echo "Process has finished at:" `date`

	  if [ $? -eq 0 ]
          then
	         echo "${biosoft} succeed." `date`
                 touch ${dir}/ok.${outname}_${biosoft}.status
	  else
	         echo "${biosoft} failed." `date`
	  fi
fi

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Process will sleep 10s."
sleep 10
echo "Process end at:" `date`
