#!/bin/bash
#Title: computeMatrix
#Author: Ao Shen
#Date: 2024-10-24
echo "Process will start at:" `date`
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

# variables
i=$1
project=/home/yjx/projects/kla/PRJNA857271
ctrl_id=`cat ${project}/data.info | sed -n "${i}p" | cut -f 2`
treat_id=`cat ${project}/data.info | sed -n "${i}p" | cut -f 3`
id=${treat_id}
biosoft=computeMatrix
dir=${project}/7.deepTools
computeMatrix=/home/yjx/anaconda3/envs/chip-seq/bin/computeMatrix
input_bw=${dir}/${ctrl_id}.bw
treat_bw=${dir}/${treat_id}.bw
TSSbed=/home/yjx/database/human/hg38/hg38_gene_TSS.bed

# computeMatrix
if [ ! -f ${dir}/ok.${id}_${biosoft}.status ]
then
    echo -e "Process started at:" `date`
    ${computeMatrix} reference-point -S ${input_bw} ${treat_bw} \
	    -R ${TSSbed} \
	    --referencePoint TSS \
	    --beforeRegionStartLength 3000 \
	    --afterRegionStartLength 3000 \
	    -o ${dir}/${id}_matrix_TSS_H3K18la.gz \
	    --skipZeros -p 8
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
