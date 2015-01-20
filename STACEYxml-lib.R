# top level functions for user to call
#
#
#

library(XML)

source("STACEYxml-libutil-addnodes.R")
source("STACEYxml-libutil-priors.R")
source("STACEYxml-libutil-ids.R")
source("STACEYxml-libutil-find.R")

source("STACEYxml-lib-data.R")
source("STACEYxml-lib-distribution.R")
source("STACEYxml-lib-init.R")
source("STACEYxml-lib-loggers.R")
source("STACEYxml-lib-operators.R")
source("STACEYxml-lib-run.R")
source("STACEYxml-lib-state.R")



TheSTACEYxmlTree <- NULL


xmlTree.from.analysis.structure <- function(A)
{
  packages=paste(sep=":",
                 "beast.core",
                 "beast.evolution.alignment",
                 "beast.evolution.tree",
                 "beast.evolution.tree.coalescent",
                 "beast.evolution.speciation",
                 "beast.core.util",
                 "beast.evolution.nuc",
                 "beast.evolution.operators",
                 "beast.evolution.sitemodel",
                 "beast.evolution.substitutionmodel",
                 "beast.evolution.likelihood")
  attrs <- c(version="2.0", namespace=packages, beautistatus="", beautitemplate="")
  TheSTACEYxmlTree <<- suppressWarnings(xmlTree(tag="beast", attrs = attrs))
  # TODO keep getting 'empty XML document' warning
  add.data(A)
  add.node("map", attrs=c(name="Beta"),  .children=c("beast.math.distributions.Beta"))
  add.node("map", attrs=c(name="Exponential"),  .children=c("beast.math.distributions.Exponential"))
  add.node("map", attrs=c(name="InverseGamma"),  .children=c("beast.math.distributions.InverseGamma"))
  add.node("map", attrs=c(name="LogNormal"),  .children=c("beast.math.distributions.LogNormalDistributionModel"))
  add.node("map", attrs=c(name="Gamma"),  .children=c("beast.math.distributions.Gamma"))
  add.node("map", attrs=c(name="Uniform"),  .children=c("beast.math.distributions.Uniform"))
  add.node("map", attrs=c(name="LaplaceDistribution"),  .children=c("beast.math.distributions.LaplaceDistribution"))
  add.node("map", attrs=c(name="OneOnX"),  .children=c("beast.math.distributions.OneOnX"))
  add.node("map", attrs=c(name="Normal"),  .children=c("beast.math.distributions.Normal"))
  add.node("map", attrs=c(name="prior"), .children=c("beast.math.distributions.Prior"))
  add.run(A)
  add.closetag() #### TODO ????
}





