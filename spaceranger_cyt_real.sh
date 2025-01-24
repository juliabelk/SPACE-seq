#!/bin/bash

ROOT= ### FILL IN
SR=${ROOT}/software/spaceranger-2.1.1 

MODE=$1 # rna or atac
FOLDER=$2 
SAMPLE=$3
GENOME=$4 
SLIDE=$5
AREA=$6
IMAGE_ALIGN=$7
IMAGE=$8

echo $IMAGE

if [ $GENOME == "mm10" ]
then
  REF=${ROOT}/software/refdata-gex-mm10-2020-A
elif [ $GENOME == "hg38" ]
then
  REF=${ROOT}/software/refdata-gex-GRCh38-2020-A
else
  echo "invalid genome: mm10 or hg38"
  exit
fi

cd $ROOT/$FOLDER/${MODE}

## put all fastqs for the sample into their own directory
if [ ! -d fastqs/${SAMPLE} ]
then
  echo "1/ organizing fastqs"
  mkdir fastqs/${SAMPLE}
  mv fastqs/${SAMPLE}* fastqs/${SAMPLE}
fi

# for atac samples, start the second mapping
if [ $MODE == "atac" ] #&& [ ! -d hisat2_${GENOME}/${SAMPLE} ]
then
  echo "submit hisat2..."
  sbatch $ROOT/scripts/atacseq/se_atacseq_R2_only.sh $FOLDER $GENOME $SAMPLE 
  #sbatch $ROOT/scripts/atacseq/se_atacseq_R2_only.sh $FOLDER GRCh38_mtMask $SAMPLE 
fi

mkdir -p spaceranger_$GENOME
cd spaceranger_$GENOME

echo $SAMPLE
if [ ! -d $SAMPLE ]
then
  echo "2/ running spaceranger"
  $SR/spaceranger count --id=$SAMPLE \
                    --description=$SAMPLE \
                    --transcriptome=$REF \
                    --fastqs=../fastqs/${SAMPLE} \
                    --image=../../pic/$IMAGE \
                    --loupe-alignment=../../pic/$IMAGE_ALIGN \
                    --unknown-slide visium-1

fi


if [ $MODE == "atac" ]
then

  source activate seurat

  mkdir -p ${SAMPLE}/outs/spaceseq
  cd ${SAMPLE}/outs/spaceseq

  TBL_O=possorted_genome_bam_tbl.tsv.gz
  if [ ! -f $TBL_O ]
  then 
    echo "3/ creating CB UMI table from bam"
    Rscript $ROOT/scripts/atacseq/bam_to_tbl.R
  fi

fi






