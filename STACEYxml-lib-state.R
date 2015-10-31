

# TODO this could be rewritten using get.all.parameters()
# in STACEYxml-libutil-find.R



add.state <- function() {
  open.xmlnode("state", attrs=c(id="state", storeEvery=get.runoption("store.every")))
  
  add.comment("The species or minimal clusters tree (smcTree)")
  open.xmlnode("tree", attrs=c(id=smcTreeID(), name="stateNode"))
  tt <- get.taxa.table()$value
  mcs <- tt[,2]
  umcs <- unique(mcs)
  tx <- tt[,1] 
  open.xmlnode("taxonset", attrs=c(id=taxonSetOfSetsID(), spec="TaxonSet"))
  for (mc in 1:length(umcs)) {
    open.xmlnode("taxon", attrs=c(id=umcs[mc], spec="TaxonSet"))
    for (i in 1:nrow(tt)) {
      if (umcs[mc] == mcs[i]) {
        add.xmlnode("taxon", attrs=c(id=tx[i], spec="Taxon"))
      }
    }        
    close.xmlnode()
  }
  close.xmlnode() 
  close.xmlnode()
  
  # TODO if fixed, not a parameter
  add.comment("The birth-death-collapse hyper-parameters for smcTree")
  bdc <- get.bdc.model()
  
  bdc.g <- get.child(bdc, "growthrate")
  if (is.estimated(bdc.g)) {
  	bdc.g.iv <- get.childvalue(bdc.g, "initval")
    attrs <- c(id=bdcGrowthID(), name="stateNode", lower="1.0E-99", upper="1.0E99")
    add.xmlnode.children("parameter", attrs=attrs, children=bdc.g.iv)
  }
  bdc.rd <- get.child(bdc, "reldeath")
  if (is.estimated(bdc.rd)) {
  	bdc.rd.iv <- get.childvalue(bdc.rd, "initval")
    attrs <- c(id=bdcRelDeathID(), name="stateNode", lower="0.0", upper="1.0")
    add.xmlnode.children("parameter", attrs=attrs, children=bdc.rd.iv)
  }
  bdc.w <- get.child(bdc, "w")
  if (is.estimated(bdc.w)) {
  	bdc.w.iv <- get.childvalue(bdc.w, "initval")
    attrs <- c(id=bdcCollapseWtID(), name="stateNode", lower="0.0", upper="1.0")
    add.xmlnode.children("parameter", attrs=attrs, children=bdc.w.iv)
  }
  add.comment("The origin height of the smcTree. Must be estimated. Initial value ignored.")
  add.xmlnode.children("parameter", attrs=c(id=bdcOriginHtID(), name="stateNode", lower="1.0E-99", upper="1.0E99"), children="1")
  
  add.comment("The hyper-parameter for the population scale for smcTree")
  psf <- get.child(get.smc.coalescent(), "popSF") 
  if (is.estimated(psf)) {
  	psf.iv <- get.childvalue(psf, "initval")
    attrs <- c(id=popSFID(), name="stateNode", lower="1.0E-99", upper="1.0E99")
    add.xmlnode.children("parameter", attrs=attrs, children=psf.iv)
  }
  
  add.comment("The trees for each locus")
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    add.comment(paste0("Gene tree ", u))
    open.xmlnode("tree", attrs=c(id=geneTreeID.u(u), name="stateNode"))
    open.xmlnode("taxonset", attrs=c(id=geneTaxonSetID.u(u), spec="TaxonSet"))
    add.xmlnode("data", attrs=c(idref=partitiondataID.u(u), name="alignment"))
    close.xmlnode()
    close.xmlnode()
  }  
  add.comment("The partition clock rates, omitting first")
  # TODO this just does strict clock
  clocks <- get.Clocks()
  if (length(clocks) > 1) {
    for (u in 2:length(clocks)) {
      prm <- get.child(clocks[[u]], "PartitionRateM")
      clkrate <- get.child(prm, "rate")
      if (is.estimated(clkrate)) {
        attrs <- c(id=clockRateID.u(u), name="stateNode")
        clkrate.iv <- get.childvalue(clkrate, "initval")
        add.xmlnode.children("parameter", attrs=attrs, children=clkrate.iv)
      } 
    }
  }
  siteMs <- get.SiteMs() 
  add.comment("The substitution models")
  for (u in 1:length(siteMs)) {
    substM <- get.child(siteMs[[u]], "SubstM")
    model <- get.child(substM, "model")
    # TODO this just does HKY
    if (model$kind == "HKY") {
      add.comment("HKY kappa and frequencies")
      kappa <- get.child(model, "kappa")
      freqs <- get.child(model, "freqs")
      if (is.numeric(kappa$value)) {        
        if (is.character(freqs$value)) {
          if (identical(freqs$value, "Empirical")) {
            # no parameters in state, nothing top do
          } else {
            # TODO frequencies "Equal"
            stop("Unknown/unimplemented option for frequencies ", freqs$value)
          }
        } else {
          stop("Unimplemented: fixed kappa, estimated freqs")
          # TODO fixed kappa, estimated freqs
        }
      } else {
        if (is.character(freqs$value)) {
          stop("Unimplemented: estimated kappa, fixed freqs")
          # TODO estimated kappa, fixed freqs
        } else {
        	# estimated kappa and freqs
        	kappa.iv <- get.child(kappa, "initval")$value
          add.xmlnode.children("parameter", attrs=c(id=kappaID.u(u), name="stateNode"), children=kappa.iv)
          kappa.iv 
          attrs <- c(id=frequenciesParamID.u(u), name="stateNode", dimension="4", lower="0.0", upper="1.0")
          add.xmlnode.children("parameter", attrs=attrs, children="0.25")
        }
      }
    }
  }
  close.xmlnode()
}













