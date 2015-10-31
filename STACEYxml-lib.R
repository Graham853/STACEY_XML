# This provides one function for user to call,
# and defines two global variables.
# It is called after the analysis structure 
# has been made using functions in use-input-lib.r.


TheAlignmentData  <- NULL
TheCurrentPath <- character(0)
TheBEASTOutputConnection <- NULL


xmlFile.from.analysis.structure <- function()
{
  # Load functions that are only needed to generate XML.
  
  
  # They are done here to hide them from user code (eg user-input.r)
  
  source(paste0(RcodeSourceDirectory, "/", "analysis-structure-lib-find.r"))
  source(paste0(RcodeSourceDirectory, "/", "analysis-structure-lib-util.r"))
  
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-addnodes.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-priors.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-ids.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-loaddata.R"))
  
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-data.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-distribution.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-init.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-loggers.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-operators.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-run.R"))
  source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-state.R"))
  

  
  TheAlignmentData <<- load.data() 
  
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
  open.toplevel.xmlnode("beast", attrs)
  add.data()
  add.xmlnode.children("map", attrs=c(name="Beta"), children=c("beast.math.distributions.Beta"))
  add.xmlnode.children("map", attrs=c(name="Exponential"), children=c("beast.math.distributions.Exponential"))
  add.xmlnode.children("map", attrs=c(name="InverseGamma"), children=c("beast.math.distributions.InverseGamma"))
  add.xmlnode.children("map", attrs=c(name="LogNormal"), children=c("beast.math.distributions.LogNormalDistributionModel"))
  add.xmlnode.children("map", attrs=c(name="Gamma"), children=c("beast.math.distributions.Gamma"))
  add.xmlnode.children("map", attrs=c(name="Uniform"), children=c("beast.math.distributions.Uniform"))
  add.xmlnode.children("map", attrs=c(name="LaplaceDistribution"), children=c("beast.math.distributions.LaplaceDistribution"))
  add.xmlnode.children("map", attrs=c(name="OneOnX"), children=c("beast.math.distributions.OneOnX"))
  add.xmlnode.children("map", attrs=c(name="Normal"), children=c("beast.math.distributions.Normal"))
  add.xmlnode.children("map", attrs=c(name="prior"), children=c("beast.math.distributions.Prior"))
  add.run()
  close.xmlnode()

}



