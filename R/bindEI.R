
bindEI = function(eset, ccleset, fixup=TRUE) {
 snames = repairNames(ccleset)
 if (!all(snames %in% sampleNames(eset))) {
    if (!fixup) stop("can't match all sample names after repair.")
    }
 eset = eset[, intersect(sampleNames(eset), snames) ]
 # need to ensure congruence of ic50 data with expression samples
 ic50s = IC50(ccleset)
 organs = organ(ccleset)
 names(organs) = names(ic50s) = snames
 eset$IC50 = ic50s[sampleNames(eset)]
 eset$organ = organs[sampleNames(eset)]
 eset
}

