#!/bin/bash

ROOT= ### FILL IN

FLDR=$1
GENOME=$2
SAMPLE=$3


R1=${ROOT}/${FLDR}/atac/fastqs/${SAMPLE}/${SAMPLE}*_R1_*.gz
R2=${ROOT}/${FLDR}/atac/fastqs/${SAMPLE}/${SAMPLE}*_R2_*.gz

gref=${ROOT}/genomes
if [ $GENOME == "mm10" ]
then
  gpth=$gref/mm10/mm10
  SR_TBL=${ROOT}/${FLDR}/atac/spaceranger_${GENOME}/${SAMPLE}/outs/spaceseq/possorted_genome_bam_tbl.tsv
  ARCHR_GENOME=mm10
elif [ $GENOME == "hg38" ]
then
  gpth=$gref/hg38/hg38
  SR_TBL=${ROOT}/${FLDR}/atac/spaceranger_${GENOME}/${SAMPLE}/outs/spaceseq/possorted_genome_bam_tbl.tsv
  ARCHR_GENOME=hg38
elif [ $GENOME == "GRCh38_mtMask" ]
then
  gpth=$gref/GRCh38_mtMask_build/hisat2/GRCh38_mtMask
  SR_TBL=${ROOT}/${FLDR}/atac/spaceranger_hg38/${SAMPLE}/outs/spaceseq/possorted_genome_bam_tbl.tsv
  ARCHR_GENOME=hg38
else
  echo "invalid genome: mm10 or hg38"
  exit
fi

source activate seurat 
ml biology hisat2 samtools 

outdir=${ROOT}/${FLDR}/atac/hisat2_${GENOME}/$SAMPLE

mkdir -p ${ROOT}/${FLDR}/atac/hisat2_${GENOME}
mkdir -p $outdir
cd $outdir

r1_atrim=${SAMPLE}_R1_001.adaptertrim.fastq.gz
r2_atrim=${SAMPLE}_R2_001.adaptertrim.fastq.gz

if [ -f $R2 ] && [ ! -f $r2_atrim ]; then
    echo "1/ trimming"

    ## EXPLANATION: want to remove low quality bases at the end of R2
    ## trimming a variable number of bases messes up deduplication (see fastp github)
    ## therefore, everything is trimmed to 75bp maximum.

    fastp --in1 $R1 --in2 $R2 --out1 $r1_atrim --out2 $r2_atrim --max_len2 75 --adapter_sequence_r2 CTGTCTCTTATACACATCT -h ${SAMPLE}.adaptertrim.html -j ${SAMPLE}.adaptertrim.json

    mkdir -p fastqc
    fastqc -o fastqc *fastq.gz
fi

if [ -f $r2_atrim ] && [ ! -f $SAMPLE.bam ]; then
    echo "2/ aligning"

    hisat2 --no-spliced-alignment --very-sensitive -X 2000 -p 8 -x $gpth --summary-file ${SAMPLE}.hisat2.log --un-gz ${SAMPLE}.nalign.fastq.gz -U $r2_atrim | samtools sort - > ${SAMPLE}.bam
    samtools index ${SAMPLE}.bam
    fastqc -t 4 -o fastqc *nalign.fastq.gz
fi

if [ -f ${SAMPLE}.bam ] &&  [ ! -f ${SAMPLE}_q30.bam ]; then
    echo "3/ samtools"
    samtools view -b -F 4 -q 30 ${SAMPLE}.bam | samtools sort - -o ${SAMPLE}_q30.bam
    samtools index ${SAMPLE}_q30.bam
fi

if [ ! -f $SR_TBL ]
then
    echo "spaceranger table not found! exiting"
    exit
fi

### add CB UMI tags to the bam
if [ ! -f ${SAMPLE}_q30_augment.bam ] && [ -f $SR_TBL ] 
then
    echo "3/ adding CB to bam file"

    python $ROOT/scripts/atacseq/augment_bam.py $FLDR $GENOME $SAMPLE
    samtools view -o ${SAMPLE}_q30_augment.bam ${SAMPLE}_q30_augment.sam
    samtools index ${SAMPLE}_q30_augment.bam
    rm ${SAMPLE}_q30_augment.sam
fi

BAM_I=${SAMPLE}_q30_augment.bam 
BAM_O=${SAMPLE}_q30_augment_dedup.bam 
BAM_O2=${SAMPLE}_q30_augment_dedup2.bam 

if [ -f $BAM_I ] && [ ! -f ${SAMPLE}_q30_augment_RG.tsv ]
then
    echo "inspecting bam read groups"
    umi_tools group -I ${BAM_I} --group-out=${SAMPLE}_q30_augment_RG.tsv --log=${SAMPLE}_q30_augment_RG.log --per-cell --extract-umi-method=tag --umi-tag=UB --cell-tag=CB
fi

if [ ! -f $BAM_O ] && [ -f $BAM_I ]
then
    echo "3/ deduplicating BAM"
    umi_tools dedup --per-cell --stdin=$BAM_I --log=umi_tools_dedup.log --extract-umi-method=tag --umi-tag=UB --cell-tag=CB --output-stats=${SAMPLE}_dedup > $BAM_O
fi  

#if [ ! -f $BAM_O2 ] && [ -f $BAM_I ]
#then
#    echo "3/ deduplicating BAM option 2"
#    umi_tools dedup --per-cell --ignore-umi --stdin=$BAM_I --log=umi_tools_dedup2.log --extract-umi-method=tag --umi-tag=UB --cell-tag=CB  > $BAM_O2
#fi  


if [ ! -f ${SAMPLE}_q30_augment_dedup.tsv.gz ] && [ -f $BAM_O ]
then 
    echo "4/ creating fragments file"
    Rscript $ROOT/scripts/atacseq/bam_to_frag.R $FLDR $ARCHR_GENOME $SAMPLE ${SAMPLE}_q30_augment_dedup
fi  


cd ../ 

if [ ! -d ${SAMPLE}_ArchR ] && [ -f ${outdir}/${SAMPLE}_q30_augment_dedup_ss.tsv.gz ]
then
    echo "5/ creating Arrow file"

    mkdir -p ${SAMPLE}_ArchR
    cd ${SAMPLE}_ArchR

    Rscript $ROOT/scripts/atacseq/create_arrow.R $ARCHR_GENOME $SAMPLE ${outdir}/${SAMPLE}_q30_augment_dedup_ss.tsv.gz 
    cd ..
fi  

if [ ! -d ${SAMPLE}_ArchR_all ] && [ -f ${outdir}/${SAMPLE}_q30_augment_dedup.tsv.gz ]
then
    echo "5/ creating Arrow file of ALL fragments"

    mkdir -p ${SAMPLE}_ArchR_all
    cd ${SAMPLE}_ArchR_all

    Rscript $ROOT/scripts/atacseq/create_arrow.R $ARCHR_GENOME $SAMPLE ${outdir}/${SAMPLE}_q30_augment_dedup.tsv.gz

fi













