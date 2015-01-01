


add.loggers <- function(A) {
  add.main.logger(A)
#   add.smctree.logger(A)
#   for (g in 1:A$nloci) {
#     add.gtree.logger(A, g)
#   }
#   add.screen.logger(A)
}

 


add.main.logger <- function(A) {
  add.opennode("logger", attrs=c(id="tracelog", fileName=A$run.options$sampledparams.fpath, 
                             logEvery=A$run.options$params.logevery, model="@posterior", sort="smart"))

  add.node("log", attrs=c(idref="posterior"))
  add.node("log", attrs=c(idref="likelihood"))
  add.node("log", attrs=c(idref="prior"))
  add.node("log", attrs=c(idref=smcCoalescentID()))
  add.node("log", attrs=c(idref=popSFID()))
  
  add.comment("birth-death-collapse model and its parameters")
  add.node("log", attrs=c(idref=bdcModelID()))
  add.node("parameter", attrs=c(idref=bdcModelID(), name="log"))
  add.node("parameter", attrs=c(idref=bdcRelDeathID(), name="log"))
  add.node("parameter", attrs=c(idref=bdcCollapseWtID(), name="log"))
  
  add.comment("smcTree origin height")
  add.node("parameter", attrs=c(idref=bdcOriginHtID(), name="log"))

  add.comment("Statistic for number of clusters")
  add.node("log", atrrs=c(spec="stacey.BirthDeathCollapseNClustersStatistic", 
                       bdcm=IDtoREF(bdcModelID()), smcTree=IDtoREF(smcTreeID())))
  
  add.comment("Statistic for population size")
  add.node("log", atrrs=c(spec="stacey.PopSampleStatistic", 
                          popPriorScale=IDtoREF(popSFID()), piomsCoalDist=IDtoREF(smcCoalescentID())))

  # TODO Depends on which clocks, substs, are used
  for (g in 1:nof.alignments()) {
    add.comment(paste0("clock and substitution parameters for gene tree ", g))
    add.node("log", attrs=c(idref=geneTreeLhoodID(g)))
    add.node("parameter", attrs=c(idref=clockRateID(g), name="log"))
    add.node("parameter", attrs=c(idref=kappaID(g), name="log"))
    add.node("parameter", attrs=c(idref=frequenciesID(g), name="log"))
    }
  add.closetag()
}


 
 # smctree.logevery=1000,
#  add.smctree.logger <- function(A) {
#   #sampledspptrees.fpath
#   catln(1, "<logger",
#         id("smcTreeLogger"),
#         namevalue("fileName", sim.options$sampledspptrees.fpath),
#         namevalue("logEvery", sim.options$beastspptree.logevery),
#         namevalue("mode", "tree"), ">")
#   catln(2, "<log", 
#         id("TreeWithMetaDataLogger.smc"),
#         spec("beast.evolution.tree.TreeWithMetaDataLogger"),
#         namevalue("tree", IDtoREF(smcTreeID())), "/>")
#   catln(1, "</logger>")
#   
#   # TODO: do I need this? My smcTree is a plain tree, no demographics
#   #   <log id="SpeciesTreeLoggerX" 
#   #   popSize="@popSize" 
#   #   popSizeTop="@popSizeTop" 
#   #   spec="beast.evolution.speciation.SpeciesTreeLogger" 
#   #   speciesTreePrior="@SpeciesTreePopSize.Species"
#   #   tree="@Tree.t:Species"> 
#   #   <treetop id="treeTopFinder" spec="beast.evolution.speciation.TreeTopFinder">
#   #   <tree idref="Tree.t:29"/>
#   #   <tree idref="Tree.t:26"/>
#   #   </treetop>
#   #   </log>
#   
#   
#   
# }
# 
# 
#   gtrees.logevery=1000,
# add.gtree.logger <- function(A, g) {
#   catln(1, "<logger",
#         id(paste0("geneTreeLogger.", g)),
#         namevalue("fileName", paste0(sim.options$sampledgtrees.fpathbase, g, ".txt")),
#         namevalue("logEvery", sim.options$beastgtrees.logevery),
#         namevalue("mode", "tree"), ">")
#   catln(2, "<log", 
#         id(paste0("TreeWithMetaDataLogger.", g)),
#         spec("beast.evolution.tree.TreeWithMetaDataLogger"),
#         namevalue("tree", IDtoREF(geneTreeID(g))), "/>")
#   catln(1, "</logger>")
# }
# 
# 
#  screen.logevery=10000
# add.screen.logger <- function(A) {
#   catln(1, "<logger",
#         id("screenlog"),
#         namevalue("logEvery", sim.options$beastscreen.logevery),
#         namevalue("model", "@posterior"), ">")
#   catln(2, idref("log", "posterior", NULL))
#   catln(2, "<log",
#         id("ESS.posterior"),
#         spec("util.ESS"),
#         namevalue("arg", "@posterior"), "/>")
#   catln(2, idref("log", "likelihood", NULL))
#   catln(2, idref("log", "prior", NULL))
#   catln(2, idref("log", smcCoalescentID(), NULL))
#   
#   catln(2, "<log",
#         spec("stacey.BirthDeathCollapseNClustersStatistic"),
#         namevalue("bdcm", IDtoREF(bdcModelID())),
#         namevalue("smcTree", IDtoREF(smcTreeID())), "/>")
#   
#   catln(1, "</logger>")
#   
# }
# 








