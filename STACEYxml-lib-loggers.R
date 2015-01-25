


add.loggers <- function(A) {
  add.bigcomment("Main logger for parameters (Trace file)")
  add.main.logger(A)
  add.bigcomment("Logger for the species or minimal clusters tree")
  add.smctree.logger(A)
  gtrees <- get.gtrees()
  add.bigcomment(paste0("Loggers for ", length(gtrees), " gene trees"))
  for (g in 1:length(gtrees)) {
    add.gtree.logger(A, g)
  }
  add.bigcomment("Logger for screen")
  add.screen.logger(A)
}

 


add.main.logger <- function(A) {
  add.opennode("logger", attrs=c(id=mainloggerID(), fileName=A$run.options$opts$sampledparams.fpath, 
                             logEvery=A$run.options$opts$params.logevery, model="@posterior", sort="smart"))

  add.node("log", attrs=c(idref=posteriorID()))
  add.node("log", attrs=c(idref=likelihoodID()))
  add.node("log", attrs=c(idref=priorID()))
  add.node("log", attrs=c(idref=smcCoalescentID()))
  add.node("log", attrs=c(idref=popSFID()))
  
  add.comment("birth-death-collapse model and its parameters")
  add.node("log", attrs=c(idref=bdcModelID()))
  add.node("parameter", attrs=c(idref=bdcGrowthID(), name="log"))
  add.node("parameter", attrs=c(idref=bdcRelDeathID(), name="log"))
  add.node("parameter", attrs=c(idref=bdcCollapseWtID(), name="log"))
  
  add.comment("smcTree origin height")
  add.node("parameter", attrs=c(idref=bdcOriginHtID(), name="log"))

  add.comment("Statistic for number of clusters")
  add.node("log", attrs=c(spec="stacey.BirthDeathCollapseNClustersStatistic", 
                       bdcm=IDtoREF(bdcModelID()), smcTree=IDtoREF(smcTreeID())))
  
  add.comment("Statistic for population size")
  add.node("log", attrs=c(spec="stacey.PopSampleStatistic", 
                          popPriorScale=IDtoREF(popSFID()), piomsCoalDist=IDtoREF(smcCoalescentID())))

  clocks <- get.clocks()
  gtrees <- get.gtrees()
  siteMs <- get.siteMs()
  
  add.comment("Tree likelihoods for each locus")
  for (g in 1:length(gtrees)) {
    add.node("log", attrs=c(idref=geneTreeLhoodID(g)))
  }

  add.comment("clock rates (relative rates)")
  for (c in 1:length(clocks)) {
    add.node("parameter", attrs=c(idref=clockRateID.c(c), name="log"))
  }
  
  add.comment("Substitution model parameters")
  for (u in 1:length(siteMs)) {
    add.node("parameter", attrs=c(idref=kappaID.u(u), name="log"))
    add.node("parameter", attrs=c(idref=frequenciesParamID.u(u), name="log"))
    }
  add.closetag()
}



add.screen.logger <- function(A) {
  attrs = c(id=screenloggerID(), logEvery="1000", model=IDtoREF(posteriorID()))
  add.opennode("logger", attrs=attrs)
  add.node("log", attrs=c(idref=posteriorID()))
  add.node("log", attrs=c(id="ESS.posterior", spec="util.ESS", arg=IDtoREF(posteriorID())))
  add.node("log", attrs=c(idref=likelihoodID()))
  add.node("log", attrs=c(idref=priorID()))
  add.node("log", attrs=c(idref=smcCoalescentID()))
  add.node("log", attrs=c(spec="stacey.BirthDeathCollapseNClustersStatistic",
                         bdcm=IDtoREF(bdcModelID()), smcTree=IDtoREF(smcTreeID())))
  add.closetag()
}




add.smctree.logger <- function(A) {
  fpath <- A$run.options$opts$sampledsmctrees.fpath
  logevery <- A$run.options$opts$smctree.logevery  
  attrs <- c(id=smctreeloggerID(), fileName=fpath, logEvery=logevery, mode="tree")
  add.opennode("logger", attrs=attrs)
  add.node("log", attrs=c(spec="beast.evolution.tree.TreeWithMetaDataLogger", tree=IDtoREF(smcTreeID())))
  add.closetag()
}



add.gtree.logger <- function(A, g) {
  gtrees <- get.gtrees()
  gtid <- gtrees[[g]]$id
  fpath <- paste0(A$run.options$opts$sampledgtrees.fpathbase, "-", gtid, ".txt")
  logevery <- A$run.options$opts$gtrees.logevery 
  attrs <- c(id=gtreeloggerID(g), fileName=fpath, logEvery=logevery, mode="tree")
  add.opennode("logger", attrs=attrs)
  add.node("log", attrs=c(spec="beast.evolution.tree.TreeWithMetaDataLogger", tree=IDtoREF(geneTreeID.g(g))))
  add.closetag()
}




