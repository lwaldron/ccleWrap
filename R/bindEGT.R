
bindEGT = function(es, chrname="chr22") {
 sm = get(load(system.file(paste("parts/", chrname, ".rda", sep=""), package="ccleWrap")))
 if (!exists("snpsif_995")) data(snpsif_995)
 rownames(sm) = snpsif_995$CCLE.name
 sml = list(sm)
 names(sml) = chrname
 make_smlSet(es, sml, harmonize=TRUE)
}

