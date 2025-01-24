
import simplesam
import sys

print(sys.argv[1:])

fldr,genome,smpl=sys.argv[1],sys.argv[2],sys.argv[3]

root = fldr + "/atac/"

bam_i = root + "hisat2_" + genome + "/" + smpl + "/" + smpl + "_q30.bam"
bam_o = root + "hisat2_" + genome + "/" + smpl + "/" + smpl + "_q30_augment.sam"

if genome == "GRCh38_mtMask":
  bam_tbl_f = root + "spaceranger_hg38/" + smpl + "/outs/spaceseq/possorted_genome_bam_tbl.tsv"
else:
  bam_tbl_f = root + "spaceranger_" + genome + "/" + smpl + "/outs/spaceseq/possorted_genome_bam_tbl.tsv"

read_info = {}

with open(bam_tbl_f) as f:
  for line in f:
    qname, chrom, pos, cb, umi = line.rstrip().split()
    read_info[qname] = (cb,umi)

with simplesam.Reader(open(bam_i)) as in_bam:
  with simplesam.Writer(open(bam_o, 'w'), in_bam.header) as out_sam:
    for read in in_bam:
      cb, umi = read_info[read.qname]
      if cb != "NA" and umi != "NA":
        read["CB"] = cb 
        read["UB"] = umi 
        out_sam.write(read)


