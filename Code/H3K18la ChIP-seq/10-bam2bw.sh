#!/bin/bash
#Title: bamCoverage bam2bw
#Author: Ao Shen
#Date: 2024-10-24
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
id=`cat ${project}/sample.list | sed -n "${i}p"`
biosoft=bamCoverage
dir=${project}/7.deepTools
bamCoverage=/home/yjx/anaconda3/envs/chip-seq/bin/bamCoverage


# bam2bw
if [ ! -f ${dir}/ok.${id}_${biosoft}.status ]
then
    echo -e "Process started at:" `date`
    ${bamCoverage} -b ${project}/5.align/${id}_picard_dedup.bam \
	    -o ${dir}/${id}.bw \
	    --outFileFormat bigwig \
	    --normalizeUsing RPKM \
	    --numberOfProcessors 8
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
