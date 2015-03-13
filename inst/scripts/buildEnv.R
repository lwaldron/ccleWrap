#README
#
#PART1.  PROCESS THE ACTIVITY DATA.  INSTALL THE PACKAGE WITH NO DATA.
#THEN YOU HAVE nstring2vec function etc. after library(ccleWrap)
#
# for these codes to work on unix, set environment variable CCLE_FILES_PATH
# to the folder where the csv and gct and txt files are, note that 
# the segmentation file is assumed to be gzipped
#
# then start R and source this file.  it will generate 3 .rda files
#  put them in ccleWrap/data and then reinstall ccleWrap
#

prefProc = function(x) {
  fn =  paste0(Sys.getenv("CCLE_FILES_PATH"), "/", x)
  stopifnot(file.exists(fn))
  fn
  }

 nstring2vec = function(x, sep=",") {
    as.numeric(strsplit(x, sep)[[1]])
    }

library(ccleWrap)

ddata = read.csv(prefProc("CCLE_NP24.2009_Drug_data_2012.02.20.csv"), stringsAsFactors=FALSE)
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

csvname=prefProc("CCLE_NP24.2009_Drug_data_2012.02.20.csv")
csvhash.md5="b64295ef99912d1d4bead76461d0e2a1"

df = read.csv(csvname, h=TRUE, stringsAsFactors=FALSE)
nr = nrow(df)
recs = lapply(1:nr, function(x) parseCCLEline(df[x,]))
ccleRx = new("ccleSet", expts=recs, dateCreated=date(),
  csvname=csvname, csvhash.md5=csvhash.md5)
save(ccleRx, file="ccleRx.rda")


#PART2.  PROCESS Expression data

#exgct = read.table(prefProc("CCLE_Expression_Entrez_2012-09-29.gct"), skip=2, h=TRUE)
z = readLines(prefProc("CCLE_Expression_Entrez_2012-09-29.gct"))
zz = strsplit(z, "\t")
cn = zz[[3]][-(1:2)]
zzz = zz[-c(1:3)]
exprs = matrix(as.numeric(NA), nr=18988, nc=1037)
for (i in 1:nrow(exprs)) exprs[i,] = as.numeric(zzz[[i]][-c(1:2)])
colnames(exprs) = cn
fn = sapply(zz[-(1:3)], "[", 1)
rownames(exprs) = fn

esif = read.table(prefProc("CCLE_Expression.Arrays.sif_2012-10-18.txt"), sep="\t", h=TRUE, stringsAsFactors=FALSE)
rownames(esif) = esif$CCLE.name
exprs = exprs[, rownames(esif)]
ccleEx = ExpressionSet(exprs)
pData(ccleEx) = esif
annotation(ccleEx) = "hgu133plus2hsentrezg.db"
save(ccleEx, file="ccleEx.rda")

#PART3.  PROCESS CNV Segmentation

#seg = read.delim(dir(patt="seg$"), sep="\t", h=TRUE)
seg = read.delim(gzfile(prefProc("CCLE_copynumber_2012-09-29.seg.gz")), sep="\t", h=TRUE)
library(GenomicRanges)
ccleCNSeg = GRanges(paste("chr", seg$Chromosome, sep=""), IRanges(seg$Start, seg$End))
values(ccleCNSeg) = seg[,-c(2,3,4)]
ccleCNSeg$CCLE_name = as.character(ccleCNSeg$CCLE_name)
save(ccleCNSeg, file="ccleCNSeg.rda")

#FURTHER WORK NEEDED for genotypes and targeted variants...
