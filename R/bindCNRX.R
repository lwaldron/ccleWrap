
bindCNRX = function( segrng, rxobj ) {
  ic50s = IC50(rxobj)
  orgs = organ(rxobj)
  comps = compound(rxobj)
  rn = repairNames(rxobj)
  names(ic50s) = names(orgs) = names(comps) = rn
  linesInSeg = segrng$CCLE_name
  names(segrng) = linesInSeg
  okn = intersect(names(ic50s), names(segrng))
  segrng = segrng[okn]
  segrng$IC50 = ic50s[okn]
  segrng$organ = orgs[okn]
  segrng$compound = comps[okn]
  segrng
}

