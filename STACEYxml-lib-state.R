

add.state <- function(A) {
  add.opennode("state", attrs=c(id="state", storeEvery="500000"))
  
  add.comment("The species or minimal clusters tree (smcTree)")
  add.opennode("tree", attrs=c(id=smcTreeID(), name="stateNode"))
  tt <- A$taxa.table$taxa.table
  mcs <- tt[,2]
  umcs <- unique(mcs)
  tx <- tt[,1] 
  add.opennode("taxonset", attrs=c(id=taxonSetOfSetsID(), spec="TaxonSet"))
  for (mc in 1:length(umcs)) {
    add.opennode("taxon", attrs=c(id=umcs[mc], spec="Taxon"))
    for (i in 1:nrow(tt)) {
      if (umcs[mc] == mcs[i]) {
        add.node("taxon", attrs=c(id=tx[i], spec="Taxon"))
      }
    }        
    add.closetag()
  }
  add.closetag() 
  add.closetag()
  
  
  add.comment("The birth-death-collapse hyper-parameters for smcTree")
  bdc <- get.bdc.model()
  value <- ifelse(bdc$growthrate$kind == "FixedValue", bdc$growthrate$par1, "100.0")
  add.node("parameter", attrs=c(id=bdcGrowthID(), name="stateNode", lower="1.0E-99", upper="1.0E99"), .children=value)
  value <- ifelse(bdc$growthrate$kind == "FixedValue", bdc$reldeath$par1, "0.5")
  add.node("parameter", attrs=c(id=bdcRelDeathID(), name="stateNode", lower="0.0", upper="1.0"), .children=value)
  value <- ifelse(bdc$w$kind == "FixedValue", bdc$w$par1, "0.5")
  add.node("parameter", attrs=c(id=bdcCollapseWtID(), name="stateNode", lower="0.0", upper="1.0"), .children=value)
  add.comment("The origin height of the smcTree. Must be estimated. Initial value ignored.")
  add.node("parameter", attrs=c(id=bdcOriginHtID(), name="stateNode", lower="1.0E-99", upper="1.0E99"), .children="1")
  
  add.comment("The hyper-parameter for the population scale for smcTree")
  psf <- get.smc.coalescent()$popSF
  value <- ifelse(psf$kind == "FixedValue", psf$par1, "0.02")
  add.node("parameter", attrs=c(id=popSFID(), name="stateNode", lower="1.0E-99", upper="1.0E99"), .children=value)
  
  add.comment("The trees for each locus")
  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    add.comment(paste0("Gene tree ", g))
    add.opennode("tree", attrs=c(id=geneTreeID.g(g), name="stateNode"))
    add.opennode("tree", attrs=c(id=geneTaxonSetID.g(g), spec="TaxonSet"))
    add.node("data", attrs=c(idref=partitiondataID.g(g), name="alignment"))
    #TODO don't understand the data refereence when trees are linked.
    add.closetag()
    add.closetag()
  }  
  add.comment("The partition clock rates, omitting first")
  # TODO this just does strict clock
  clocks <- get.clocks()
  if (length(clocks) > 1) {
    for (c in 2:length(clocks)) {
      add.node("parameter", attrs=c(id=clocks[[c]]$id, name="stateNode"), .children="1.0")
    }
  }

  siteMs <- get.siteMs() 
  add.comment("The substitution models")
  add.comment("HKY kappas and uniform frequencies")
  # TODO this just does HKY
  for (u in 1:length(siteMs)) {
    # TODO what does u index? site models? OK?
    add.node("parameter", attrs=c(id=kappaID.u(u), name="stateNode"), .children="2.0")
    attrs <- c(id=frequenciesParamID.u(u), name="stateNode", dimension="4", lower="0.0", upper="1.0")
    add.node("parameter", attrs=attrs, .children="0.25")
  }
  
  add.closetag()
}







