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

ddf = function(ccleSet) {
    expts = ccleSet@expts
    dfs = lapply(expts, function(x) data.frame(logdoses=log10(x@doses_uM),
                   medact=x@activityMedian, ec50=x@EC50_uM, cell_line=x@line,
                   ic50 = x@IC50_uM,
                   organ=x@organ, compound=x@compound, stringsAsFactors=FALSE))
    df = do.call(rbind, dfs)
}


pl2 = function (x, y, ...) 
{
    dose = x@doses_uM
    medy = x@activityMedian
    df = data.frame(log10dose = log10(dose), medy = medy)
    require(ggplot2)
    ggplot(df) + geom_point(aes(x = log10dose, y = medy)) + ylab(paste("median activity in cell line ", 
        x@line)) + xlab(paste("log10 dose of ", x@compound, " in uM", 
        sep = "")) + geom_vline(aes(xintercept = log10(ic50)), 
        data = data.frame(ic50 = x@IC50_uM), linetype = 2, colour = "red") + 
        ggtitle(paste("Compound: ", x@compound, "Organ: ", tolower(x@organ), sep = "")) + 
	stat_smooth(aes(x = log10dose, y = medy)) + theme(text=element_text(size=10))
}

setMethod("plot", "ccleExpt", function(x, y, ...) {
  pl2(x, y, ...)
})
#  oldpar=par(no.readonly=TRUE)
#  on.exit(par(oldpar))
#  dose = x@doses_uM
#  medy = x@activityMedian
#  df = data.frame(dose=dose, medy=medy)
#  require(ggplot2)
#  ggplot(df) + geom_point(aes(x=dose, y=medy)) + ylab(paste(
#     "median activity in line ", x@line)) +
#     xlab( paste("dose of ", x@compound, " in uM", sep="")) + scale_x_log10() +
#     geom_vline(aes(xintercept=ec50), data=data.frame(ec50=x@EC50_uM), linetype=2, colour="red") +
#     ggtitle(paste("Organ: ", tolower(x@organ), sep=""))
#  plot(dose, medy, xlab=paste("dose of ", x@compound, " in uM", sep=""),
#    ylab=paste("median Activity in line ", x@line, "\n organ: ", x@organ, sep=""), ...)
#})

setClass("ccleSet", representation(expts="list",
    dateCreated="character", csvname="character",
    csvhash.md5="character"))

setMethod("show", "ccleSet", function(object) {
 cat("Broad/Novartis Cancer Cell Line Encyclopedia data.\n")
 cat("There are ", length(object@expts), " lines/experiments represented.\n")
 cat("Use '[', '[[', organ(), compound(), ... to obtain more information.\n")
})

setMethod("[", c("ccleSet", "ANY", "missing"), function(x, i, j, ..., drop=FALSE) {
  x@expts = x@expts[i]
  x
})
setMethod("[[", c("ccleSet", "ANY", "missing"), function(x, i, j, ..., drop=FALSE) {
  if (length(i) == 1 & is.atomic(i)) return(x@expts[[i]])
  else lapply(i, function(z) x@expts[[z]])
})

setGeneric("organ", function(x)standardGeneric("organ"))
setMethod("organ", "ccleSet", function(x) sapply(x@expts, slot, "organ"))
setGeneric("compound", function(x)standardGeneric("compound"))
setMethod("compound", "ccleSet", function(x) sapply(x@expts, slot, "compound"))
setGeneric("EC50", function(x)standardGeneric("EC50"))
setMethod("EC50", "ccleSet", function(x) {
   ans = sapply(x@expts, slot, "EC50_uM")
   names(ans) = line(x)
   ans
})
setGeneric("IC50", function(x)standardGeneric("IC50"))
setMethod("IC50", "ccleSet", function(x) {
   ans = sapply(x@expts, slot, "IC50_uM")
   names(ans) = line(x)
   ans
})
setGeneric("target", function(x)standardGeneric("target"))
setMethod("target", "ccleSet", function(x) {
   ans = sapply(x@expts, slot, "target")
   names(ans) = line(x)
   ans
})
setGeneric("line", function(x)standardGeneric("line"))
setMethod("line", "ccleSet", function(x) sapply(x@expts, slot, "line"))

plmany = function( ccleSet, ... ) {
  df = ddf(ccleSet)
  comps = unique(df$compound)
  if (length(comps)!=1) compcode = paste(comps, sep="/") 
      else compcode = comps
  p1 = ggplot(df)
  p1 + geom_point(aes(x=logdoses, y=medact, shape=cell_line)) + stat_smooth(aes(x=logdoses,
       y=medact, colour=cell_line)) + ylab("median activity") + 
          xlab(paste("log10 dose (uM) of", compcode, sep=" "))
}

setMethod("plot", "ccleSet", function(x, y, ...) plmany(x, ...))

setClass("targVar", contains="GRanges")
setMethod("show", "targVar", function(object) {
  cat("ccleWrap targVar instance with", length(object), "ranges.\n", sep=" ")
  cat("use ranges(), values(), etc. for more details.\n")
})

