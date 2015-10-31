


add.loggers <- function() {
  add.bigcomment("Main logger for parameters (Trace file)")
  add.main.logger(get.runoption("sampledparams.fpath"), 
  								get.runoption("params.logevery"))
  add.bigcomment("Logger for the species or minimal clusters tree")
  add.smctree.logger(get.runoption("sampledsmctrees.fpath"),
  									 get.runoption("smctree.logevery"))
  gtrees <- get.GTrees()
  add.bigcomment(paste0("Loggers for ", length(gtrees), " gene trees"))
  for (g in 1:length(gtrees)) {
    add.gtree.logger(g,
    								 get.runoption("sampledgtrees.fpathbase"),
    								 get.runoption("gtrees.logevery"))
  }
  add.bigcomment("Logger for screen")
  add.screen.logger(get.runoption("screen.logevery"))
}



add.main.logger <- function(fpath, logevery) {
  open.xmlnode("logger", attrs=c(id=mainloggerID(), fileName=fpath, 
                             logEvery=logevery, model="@posterior", sort="smart"))

  add.centredcomment("Distributions")
  add.comment("Posterior = likelihood * prior * coalescent")
  add.xmlnode("log", attrs=c(idref=posteriorID()))
  add.xmlnode("log", attrs=c(idref=likelihoodID()))
  add.xmlnode("log", attrs=c(idref=priorID()))
  add.xmlnode("log", attrs=c(idref=smcCoalescentID()))
  add.comment("birth-death-collapse model")
  add.xmlnode("log", attrs=c(idref=bdcModelID()))
  gtrees <- get.GTrees()
  add.comment("Tree likelihoods for each locus")
  for (g in 1:length(gtrees)) {
    add.xmlnode("log", attrs=c(idref=geneTreeLhoodID(g)))
  }
  
  # TODO The parameters section could be rewritten using get.all.parameters()
  # in STACEYxml-libutil-find.R
  add.centredcomment("Parameters")
  add.comment("parameters for stacey coalescent and birth-death-collapse model")
  add.xmlnode("log", attrs=c(idref=popSFID()))
  bdc <- get.bdc.model()
  bdc.g <- get.child(bdc, "growthrate")
  bdc.rd <- get.child(bdc, "reldeath")
  bdc.w <- get.child(bdc, "w")
  if (is.estimated(bdc.g)) {
    add.xmlnode("parameter", attrs=c(idref=bdcGrowthID(), name="log"))
  }
  if (is.estimated(bdc.rd)) {
    add.xmlnode("parameter", attrs=c(idref=bdcRelDeathID(), name="log")) 
  }
  if (is.estimated(bdc.w)) {
    add.xmlnode("parameter", attrs=c(idref=bdcCollapseWtID(), name="log")) 
  }
  add.xmlnode("parameter", attrs=c(idref=bdcOriginHtID(), name="log"))

  clocks <- get.Clocks()
  siteMs <- get.SiteMs()
  add.comment("clock rates (relative rates)")
  for (u in 1:length(clocks)) {
  	prm <- get.child(clocks[[u]], "PartitionRateM")
  	clkrate <- get.child(prm, "rate")
    if (is.estimated(clkrate)) {
      add.xmlnode("parameter", attrs=c(idref=clockRateID.u(u), name="log"))
    }
  }
  
  add.comment("Substitution model parameters")
  for (u in 1:length(siteMs)) {
  	substM <- get.child(siteMs[[u]], "SubstM")
  	model <- get.child(substM, "model")
    if (model$kind == "HKY") { 
    	kappa <- get.childvalue(model, "kappa")
    	freqs <- get.childvalue(model, "freqs")
      if (is.numeric(kappa)) {
        # fixed value, don't log
      } else {
        add.xmlnode("parameter", attrs=c(idref=kappaID.u(u), name="log"))
      }
      if (is.character(freqs)) {
        # fixed value, don't log
      } else {
        add.xmlnode("parameter", attrs=c(idref=frequenciesParamID.u(u), name="log"))
      }
    } else {
      # TODO
      stop()
    }
  }

  add.centredcomment("Statistics")
  add.comment("Statistic for number of clusters")
  add.xmlnode("log", attrs=c(spec="stacey.BirthDeathCollapseNClustersStatistic", 
                          bdcm=IDtoREF(bdcModelID()), smcTree=IDtoREF(smcTreeID())))
  
  add.comment("Statistic for population size")
  add.xmlnode("log", attrs=c(spec="stacey.PopSampleStatistic", 
                          popPriorScale=IDtoREF(popSFID()), piomsCoalDist=IDtoREF(smcCoalescentID())))
  
  close.xmlnode()
}



add.screen.logger <- function(logevery) {
  attrs = c(id=screenloggerID(), logEvery=logevery, model=IDtoREF(posteriorID()))
  open.xmlnode("logger", attrs=attrs)
  add.xmlnode("log", attrs=c(idref=posteriorID()))
  add.xmlnode("log", attrs=c(id="ESS.posterior", spec="util.ESS", arg=IDtoREF(posteriorID())))
  add.xmlnode("log", attrs=c(idref=likelihoodID()))
  add.xmlnode("log", attrs=c(idref=priorID()))
  add.xmlnode("log", attrs=c(idref=smcCoalescentID()))
  add.xmlnode("log", attrs=c(spec="stacey.BirthDeathCollapseNClustersStatistic",
                         bdcm=IDtoREF(bdcModelID()), smcTree=IDtoREF(smcTreeID())))
  close.xmlnode()
}




add.smctree.logger <- function(fpath, logevery) {
  attrs <- c(id=smctreeloggerID(), fileName=fpath, logEvery=logevery, mode="tree")
  open.xmlnode("logger", attrs=attrs)
  add.xmlnode("log", attrs=c(spec="beast.evolution.tree.TreeWithMetaDataLogger", tree=IDtoREF(smcTreeID())))
  close.xmlnode()
}



add.gtree.logger <- function(u, fpathbase, logevery) {
  gtrees <- get.GTrees()
  gtid <- gtrees[[u]]$id
  fpath <- paste0(fpathbase, "-", gtid, ".txt")
  attrs <- c(id=gtreeloggerID(u), fileName=fpath, logEvery=logevery, mode="tree")
  open.xmlnode("logger", attrs=attrs)
  add.xmlnode("log", attrs=c(spec="beast.evolution.tree.TreeWithMetaDataLogger", tree=IDtoREF(geneTreeID.u(u))))
  close.xmlnode()
}




