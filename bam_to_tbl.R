suppressMessages({
    library(Rsamtools)
    library(dplyr)
})

print("entering bam to table script...")

bamFile <- BamFile("../possorted_genome_bam.bam")

sb <- scanBam(bamFile,
    param = ScanBamParam(
      tag = c("CB", "UB"),
      what = c("rname", "pos", "qname")
    ))[[1]]

df <- data.frame(
    qname = sb$qname,
    chr = sb$rname,
    start = sb$pos, 
    CB = sb$tag[["CB"]],
    UB = sb$tag[["UB"]],
    stringsAsFactors=F
)

print(paste0("total reads (M): ",nrow(df)/1000000))

write.table(df,"possorted_genome_bam_tbl.tsv",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
zipped <- bgzip("possorted_genome_bam_tbl.tsv", "possorted_genome_bam_tbl.tsv.gz")


