suppressMessages({
    library(Rsamtools)
    library(dplyr)
    library(ggplot2)
})

args = commandArgs(trailingOnly=TRUE)
fldr <- args[1]
genome <- args[2]
smpl <- args[3]
bam_name <- args[4]

root <- ### FILL IN

tissue_pos <- paste0(root,fldr,"/atac/spaceranger_",genome,"/",smpl,"/outs/spatial/tissue_positions.csv")

if (!file.exists(tissue_pos)) {
    print(tissue_pos)
    print("not found! exiting")
    q()
}

bamFile <- BamFile(paste0(bam_name,".bam"))

if (genome == "mm10") {
  valid_chr = c(1:19,"X","Y","M")
} else if (genome == "hg38") {
  valid_chr = c(1:22,"X","Y","M")
} else {
  print("genome error, exiting...")
  q()
}

sb <- scanBam(bamFile,
    param = ScanBamParam(
      tag = c("CB", "UB"),
      what = c("rname", "pos", "strand", "qwidth")
    ))[[1]]

df <- data.frame(
    chr = as.character(sb$rname),
    pos = sb$pos,
    pos2 = sb$pos + sb$qwidth,
    strand = sb$strand,
    CB = sb$tag[["CB"]],
    UB = sb$tag[["UB"]],
    stringsAsFactors=F
)
df$chr[which(df$chr == "MT")] <- "M"

# From ArchR manual: 
# @param offsetPlus The numeric offset to apply to a "+" stranded Tn5 insertion to account for the precise Tn5 binding site.
offsetPlus = 4
offsetMinus = -5

#df$start <- unlist(lapply(1:nrow(df),function(i) {
#  if (df$strand[i] == "+") return(df$pos[i] + offsetPlus)
#  if (df$strand[i] == "-") return(df$pos2[i] + offsetMinus)
#}))

df$start <- df$pos + offsetPlus
sel <- which(df$strand == "-")
df$start[sel] <- df$pos2[sel] + offsetMinus
#df$start <- as.integer(df$start)

print("number of adjusted insertions < 1:")
print(length(which(df$start == 0)))
print("setting these to 1....")
df$start[which(df$start < 1)] <- 1

print(head(df))

df <- df[which(!is.na(df$CB)),]
#df$chr <- as.character(df$chr)

if ("chr1" %in% df$chr) {
  valid_chr <- paste0("chr",valid_chr)
}

print(valid_chr)
print(nrow(df))
print('removing reads mapping to non-standard chromosomes...')
df <- df[which(df$chr %in% valid_chr),]
print(nrow(df))

if ("1" %in% df$chr) {
  print("appending chr to chromosome names...")
  df$chr <- paste0("chr",df$chr)
}
print(head(df))

#df$CB_UB <- paste0(df$CB,"_",df$UB)
#df_cu <- as.data.frame(table(df$CB_UB))
#rownames(df_cu) <- df_cu$Var1

#df <- df[which(df_cu[df$CB_UB,"Freq"] == 1),]

#print(nrow(df))

se_to_pe_frag <- function(df) {
    df <- df %>% arrange(start)
    df <- df %>% arrange(chr)
    df <- df %>% arrange(CB)

    df$odd <- (1:nrow(df)) %% 2 == 1

    df$end <- c(df$start[2:nrow(df)],NA)
    df$chr_1 <- c(df$chr[2:nrow(df)],NA)
    df$CB_1 <- c(df$CB[2:nrow(df)],NA)

    df$kp <- df$odd & df$chr == df$chr_1 & df$CB == df$CB_1
    
    pe_df <- df[which(df$kp),c("chr","start","end","CB")]
    
    pe_df <- pe_df %>% arrange(start)
    pe_df <- pe_df %>% arrange(chr)

# *make sure* the coordinates are stored as integers
# if not, if any coordinates are round numbers, R may
# convert them into scientific notation (53000000 -> 5.3e+07)
# which will not work with any of the downstream steps
    pe_df$start <- as.integer(pe_df$start)
    pe_df$end   <- as.integer(pe_df$end)
     
    return(pe_df)
}

pe_frag <- se_to_pe_frag(df)
print("number of pe fragments:")
print(nrow(pe_frag))

#spatial <- read.delim("../spatial/tissue_positions.csv",sep=",")
spatial <- read.delim(tissue_pos,sep=",") 
rownames(spatial) <- spatial$barcode 
valid_bc <- spatial$barcode[which(spatial$in_tissue == 1)]

print("number of valid barcodes:")
print(length(valid_bc))

pe_frag_ss <- pe_frag[which(pe_frag$CB %in% valid_bc),]

print(paste0("total frag (M): ",nrow(pe_frag)/1000000))
print(paste0("in tissue frag (M): ",nrow(pe_frag_ss)/1000000))

f_all <- paste0(bam_name,".tsv")
write.table(pe_frag,f_all,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
zipped <- bgzip(f_all, paste0(f_all,".gz"))
idx <- indexTabix(zipped, format="bed")

f_ss <- paste0(bam_name,"_ss.tsv")
write.table(pe_frag_ss,f_ss,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
zipped <- bgzip(f_ss, paste0(f_ss,".gz"))
idx <- indexTabix(zipped, format="bed")




