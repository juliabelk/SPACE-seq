# SPACE-seq

This repository includes scripts to analyze SPACE-seq data (spatial multi-omic data generated using the 10X Visium platform). Please find more details in the accompanying SPACE-seq manuscript (Huang*, Belk* et.al., _in review_).

- `spaceranger_cyt_real.sh` is the first script which runs the 10X spaceranger pipeline
- `bam_to_tbl.R` extracts the corrected cell barcodes and UMIs from the 10X spaceranger bam file for each read
- `se_atacseq_R2_only.sh` aligns the spatial ATAC-seq data using hisat2 parameters suited for ATAC-seq analysis
- `augment_bam.py` transfers the corrected CB and UB tags from the 10X spaceranger ATAC bam to the hisat2 ATAC bam
- `bam_to_frag.R` converts the hisat2 ATAC bam into a fragments file suitable for downstream analysis in ArchR
- `create_arrow.R` creates an ArchR project using the spatial ATAC fragments file
