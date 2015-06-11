setClass("ccleExpt", representation(line="character", organ="character",
 compound="character", target="character", doses_uM="numeric",
 activityMedian="numeric", activitySD="numeric", fitType="character",
 EC50_uM="numeric", IC50_uM="numeric", Amax="numeric", ActArea="numeric"))

setMethod("show", "ccleExpt", function(object) {
 cat("Cancer Cell Line Encyclopedia experiment data for line ", object@line, 
    "\n organ ", object@organ, 
    "\n compound ",object@compound, 
    "\n target ",object@target, 
    "\n", sep="")
 })

setClass("ccleSet", representation(expts="list",
    dateCreated="character", csvname="character",
    csvhash.md5="character"))

setMethod("show", "ccleSet", function(object) {
 cat("Broad/Novartis Cancer Cell Line Encyclopedia data.\n")
 cat("There are ", length(object@expts), " lines/experiments represented.\n")
 cat("Use '[', organ(), compound(), ... to obtain more information.\n")
})

csvname="CCLE_NP24.2009_Drug_data_2012.02.20.csv"
csvhash.md5="b64295ef99912d1d4bead76461d0e2a1"

nstring2vec = function(x, sep=",") {
  as.numeric(strsplit(x, sep)[[1]])
}

parseCCLEline = function(dfline, validateNames=FALSE) {
 # parse a row of ccle CSV import
 if (nrow(dfline)!=1 || !is.data.frame(dfline)) 
    stop("input must be 1line data.frame")
 if(!identical(colnames(dfline)[1], "CCLE.Cell.Line.Name"))
   stop("First column name in dfline should be CCLE.Cell.Line.Name")
 jan2013names = c("CCLE.Cell.Line.Name", "Primary.Cell.Line.Name", "Compound", 
"Target", "Doses..uM.", "Activity.Data..median.", "Activity.SD", 
"Num.Data", "FitType", "EC50..uM.", "IC50..uM.", "Amax", "ActArea"
)
 if (validateNames && !isTRUE(all.equal(names(dfline), jan2013names)))
     warning("it appears the names of the input data.frame do not agree with expectation")
 getOrgan = function(x) {x = sub("_", "@@", x); gsub(".*@@", "", x) }

 cvec = as.character(dfline)
 new("ccleExpt",
      organ=getOrgan(cvec[1]),
      line=cvec[2], compound=cvec[3],
      target = cvec[4], doses_uM=nstring2vec(cvec[5]),
      activityMedian=nstring2vec(cvec[6]),
      activitySD=nstring2vec(cvec[7]),
      fitType=cvec[9],
      EC50_uM = as.numeric(cvec[10]),
      IC50_uM = as.numeric(cvec[11]),
      Amax = as.numeric(cvec[12]),
      ActArea = as.numeric(cvec[13]))
}

df = read.csv(csvname, h=TRUE, stringsAsFactors=FALSE)
nr = nrow(df)
recs = lapply(1:nr, function(x) parseCCLEline(df[x,]))
ccle = new("ccleSet", expts=recs, dateCreated=date(),
  csvname=csvname, csvhash.md5=csvhash.md5)
