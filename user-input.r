


RcodeSourceDirectory <- "C:/Users/Work/AAA/Programming/biology/STACEY_XML"

library(XML)
library(stringr)
library(ape)

source(paste0(RcodeSourceDirectory, "/", "user-input-lib.r"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib.r"))

source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-addnodes.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-priors.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-ids.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-libutil-find.R"))

source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-data.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-distribution.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-init.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-loggers.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-operators.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-run.R"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib-state.R"))




nloci <- 5
alignment.table <- NULL
for (i in 1:nloci) {
  file <- paste0("seqs", sprintf("%02d", i), ".nex")
  id <- paste0(i)
  partition <- list(file=file, gtree=Gtree(id), clock=Clock(id), siteM=SiteModel(id))
  alignment.table <- c(alignment.table, list(partition))
}
names(alignment.table) <- paste0("partition", 1:nloci)


taxa.table <- rbind(
  c("a01_A", "a"),
  c("a02_A", "a"),
  c("a03_A", "a"),
  c("a04_A", "a"),
  c("a05_A", "a"),
  c("b01_A", "b"),
  c("b02_A", "b"),
  c("b03_A", "b"),
  c("b04_A", "b"),
  c("b05_A", "b"),
  c("c01_A", "c"),
  c("c02_A", "c"),
  c("c03_A", "c"),
  c("c04_A", "c"),
  c("c05_A", "c")
)
colnames(taxa.table) <- c("taxon", "mincluster")


run.options <- list(
  sampledgtrees.fpathbase="gtrees",
  sampledsmctrees.fpath="smctrees.txt",
  sampledparams.fpath="params.txt",
  chainlength=10000000,
  store.every=500000,
  params.logevery=1000,
  smctree.logevery=1000,
  gtrees.logevery=1000,
  screen.logevery=10000
  )          
TheAnalysis("Test.1", data.dpath="C:/Users/Work/AAA/Programming/biology/STACEY_XML/Tests", 
            alignment.table, taxa.table, run.options)

########################################

gtrees <- get.gtrees()
for (i in 1:length(gtrees)) {
  id <- gtrees[[i]]$id
  Gtree(id, SMCtree("smctree"), BranchRM(id), 2.0)
  BranchRM(id, StrictClock(id))
  StrictClock(id, 1.0)
}


SMCtree("smctree", BDCPrior("smctree"), SMCCoalescent("staceycoal"))
BDCPrior("smctree", growthrate=LogNorm("smctree.g"), 
         reldeath=Beta("smctree.rd"), 
         w=Beta("smctree.w"), 
         ohID="origin.height",
         eps=0.00003)
LogNorm("smctree.g", 4.6, 2)
Beta("smctree.rd", 1, 1)
Beta("smctree.w", 1, 1) 
SMCCoalescent("staceycoal", invgammamix=InvGammaMix("popBV"), popSF=LogNorm("popSF"))
LogNorm("popSF", -7, 2)
InvGammaMix("popBV", weights=c(0.5,0.5), alphas=c(3,3), betas=c(1.5,2.5))


clocks <- get.clocks()
for (i in 1:length(clocks)) {
  id <- clocks[[i]]$id
  Clock(id, PartitionRate(id))
  if (i==1) {
    PartitionRate(id, 1.0)
  } else {
    clockID <- paste0("clockrate.", id)
    PartitionRate(id, LogNorm(clockID))
    LogNorm(clockID, 0.0, 1.0) 
  }
}


siteMs <- get.siteMs()
for (i in 1:length(siteMs)) {
  id <- siteMs[[i]]$id
  SiteModel(id, Subst(id), SiteHet(id))
  SiteHet(id, "None")
  Subst(id, HKY(id))
  freqsID <- paste0("HKYfreqs.", id)
  kappaID <- paste0("HKYkappa.", id)
  HKY(id, UniformUnitSimplex(freqsID), LogNorm(kappaID))
  UniformUnitSimplex(freqsID, 3)
  LogNorm(kappaID, 1.0, 1.25)
}




xmlTree.from.analysis.structure(TheAnalysisStructure)

saveXML(TheSTACEYxmlTree$value(), file=paste0(TheAnalysisStructure$data.dpath, "/test.xml"))

sink(file=paste0(TheAnalysisStructure$data.dpath, "/test.txt"))
str(TheAnalysisStructure)
sink(NULL)


