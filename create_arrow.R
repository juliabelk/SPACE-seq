suppressMessages({
    library(ArchR)
    library(dplyr)
})

args = commandArgs(trailingOnly=TRUE)

print(args)

genome <- args[1]
sample <- args[2]
frag_pth <- args[3]

if (genome == "mm10") {
  data("geneAnnoMm10")
  data("genomeAnnoMm10")
  geneAnno <- geneAnnoMm10
  genomeAnno <- genomeAnnoMm10
} else if (genome == "hg38") {
  data("geneAnnoHg38")
  data("genomeAnnoHg38")
  geneAnno <- geneAnnoHg38
  genomeAnno <- genomeAnnoHg38
} else {
  print("genome error, exiting...")
  q() 
}

if (!dir.exists("Arrow/")) dir.create("Arrow/")
setwd("Arrow")

ArrowFiles <- createSpaceArrowFiles(
    inputFiles = c(frag_pth),
    sampleNames = c(sample),
    geneAnno = geneAnno,
    genomeAnno = genomeAnno,
    TileMatParams = list("tileSize" = 5000),
    minFragSize = -1,
    maxFragSize = 10**10,
    minFrags = 100,
    minTSS = 1, #2,
    force = TRUE
)

setwd("..")


