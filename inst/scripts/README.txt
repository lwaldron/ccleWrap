README

PART1.  PROCESS THE ACTIVITY DATA.  INSTALL THE PACKAGE WITH NO DATA.
THEN YOU HAVE nstring2vec function etc. after library(ccleWrap)


ddata = read.csv("CCLE_NP24.2009_Drug_data_2015.02.24.csv", stringsAsFactors=FALSE) # updated

# ddata = read.csv("CCLE_NP24.2009_Drug_data_2012.02.20.csv", stringsAsFactors=FALSE)


dosemat = sapply(ddata[,5], nstring2vec)
activityMeds = sapply(ddata[,"Activity.Data..median."], nstring2vec)
activitySDs = sapply(ddata[,"Activity.SD"], nstring2vec)
lines = sub("_", ",", ddata[,1])
lineName = sapply(strsplit(lines, ","), "[", 1)
lineOrg = sapply(strsplit(lines, ","), "[", 2)
#setClass("ccleExpt", representation(line="character", organ="character",
# compound="character", target="character", doses_uM="numeric",
# activityMedian="numeric", activitySD="numeric", fitType="character",
# EC50_uM="numeric", IC50_uM="numeric", Amax="numeric", ActArea="numeric"))

csvname="CCLE_NP24.2009_Drug_data_2012.02.20.csv"
csvhash.md5="b64295ef99912d1d4bead76461d0e2a1"



df = read.csv(csvname, h=TRUE, stringsAsFactors=FALSE)
nr = nrow(df)
recs = lapply(1:nr, function(x) parseCCLEline(df[x,]))
ccleRx = new("ccleSet", expts=recs, dateCreated=date(),
  csvname=csvname, csvhash.md5=csvhash.md5)
save(ccleRx, file="ccleRx.rda")


PART2.  PROCESS Expression data

exgct = read.table("CCLE_Expression_Entrez_2012-09-29.gct", skip=2, h=TRUE)

## Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : line 276 did not have 1039 elements

z = readLines("CCLE_Expression_Entrez_2012-09-29.gct")
zz = strsplit(z, "\t")
zzz = zz[-c(1:3)]
exprs = matrix(as.numeric(NA), nr=18988, nc=1037) # change to 1036
for (i in 1:nrow(exprs)) exprs[i,] = as.numeric(zzz[[i]][-c(1:2)])

esif = read.table("CCLE_Expression.Arrays.sif_2012-10-18.txt", sep="\t", h=TRUE, stringsAsFactors=FALSE)
rownames(esif) = esif$CCLE.name
exprs = exprs[, rownames(esif)]
# Error in exprs[, rownames(esif)] : no 'dimnames' attribute for array

ccleEx = ExpressionSet(exprs)
pData(ccleEx) = esif
annotation(ccleEx) = "hgu133plus2hsentrezg.db"
save(ccleEx, file="ccleEx.rda")

PART3.  PROCESS CNV Segmentation

# seg = read.delim(dir(patt="seg$"), sep="\t", h=TRUE)
# Error in file(file, "rt") : invalid 'description' argument 
seg = read.delim(file.path(datadir, "DNACopyNumber/CCLE_copynumber_2013-12-03.seg.txt"))


library(GenomicRanges)
ccleCNSeg = GRanges(paste("chr", seg$Chromosome, sep=""), IRanges(seg$Start, seg$End)) 

mcols(ccleCNSeg) = seg[,-c(2,3,4)]

ccleCNSeg$CCLE_name = as.character(ccleCNSeg$CCLE_name) # what is this?
save(ccleCNSeg, file="ccleCNSeg.rda")



byGene = read.delim(file.path(datadir, "DNACopyNumber/CCLE_copynumber_byGene_2013-12-03.txt"))

ccleCNGene = GRanges(paste("chr", byGene$CHR, sep=""), IRanges(byGene$CHRLOC, byGene$CHRLOCEND))

mcols(ccleCNGene) = byGene[,-c(3,4,5)]

## Use case to compare the seg value with the gene value

ccleCNSeg[ccleCNSeg@elementMetadata[,1] == colnames(byGene)[6],]

ccleCNGene[,colnames(byGene)[6]]

# function to check if the gene i which is fully contained in a segment has the same assigned value as that of the segment
fun <- function(i) {
print(i)
query[genes_in_segment@queryHits[i]]@elementMetadata[,1] ==
subject[genes_in_segment@subjectHits[i]]@elementMetadata[,3]
}

# function to find the values of the gene i and the values of all its overlapping segments 
find_values <- function(i) {
query_i = ccleCNGene[i,sample_name]
segments_i = findOverlaps(query_i, subject)@subjectHits
print(c("gene_id: ", i, "gene_value: ", query_i@elementMetadata[,1],  "Values of intersecting segments: ", subject[segments_i]@elementMetadata[,3]))
}

FURTHER WORK NEEDED for genotypes and targeted variants...

check_sample <- function(sample_name) {
sample_name = colnames(byGene)[6] 
#6 could be replaced by any number between 6 and ncol(byGene)#
query = ccleCNGene[,sample_name]
subject = ccleCNSeg[ccleCNSeg@elementMetadata[,1] == sample_name,] 

genes_in_segment = findOverlaps(query, subject, type = "within")


# The below code checks if the genes in segment all have the same value as segment 

lapply(1:length(genes_in_segment), fun)

# n_genes_overlap_segments = length(query) - length(genes_in_segment)

gene_id_overlap_segments = setdiff(1:length(query), genes_in_segment@queryHits) 

# the below code will report all values of genes and its overlapping segments - this is run only for genes that partially overlap - but it can be run for any gene, in general. 
lapply(gene_id_overlap_segments,find_values) 
}

lapply(6:ncol(byGene), check_sample)
## further work will be included in the biocMultiAssay package