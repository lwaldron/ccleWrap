
ieplot = function (es, pb, legx = 8, legy = .6, noleg=FALSE, ...) 
{
    sym = AnnotationDbi::get(pb, hgu133plus2hsentrezgSYMBOL)
    org = factor(es$organ)
    uorg = factor(unique(es$organ))
    plot(exprs(es)[pb, ], es$lic50, xlab = paste(sym, "expression (RMA)", 
        sep = " "), ylab = "log10 IC50", col = as.numeric(org), 
        pch = 19, ...)
    if (!noleg) legend(legx, legy, pch = 19, col = as.numeric(uorg), legend = uorg)
}

