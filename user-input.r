setwd("C:/Users/Work/AAA/Programming/biology/STACEY_XML")
source("user-input-lib.r")
source("STACEYxml-lib.r")

Rprof()

alignment.table <- list(
  partition1=list(file="U1.nex", gtree=Gtree("1"), clock=Clock("1"), siteM=SiteModel("1")),
  partition2=list(file="U2.nex", gtree=Gtree("2"), clock=Clock("2"), siteM=SiteModel("2")),
  partition3=list(file="U3.nex", gtree=Gtree("3"), clock=Clock("3"), siteM=SiteModel("3")),
  partition4=list(file="U4.nex", gtree=Gtree("4"), clock=Clock("4"), siteM=SiteModel("4")),
  partition5=list(file="U5.nex", gtree=Gtree("5"), clock=Clock("5"), siteM=SiteModel("5"))
)

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
TheAnalysis("Test.1", data.dpath="C:/...", alignment.table, taxa.table, run.options)

########################################

gtrees <- get.gtrees()
for (i in 1:length(gtrees)) {
  id <- gtrees[[i]]$id
  Gtree(id, SMCtree("smctree"), BranchRM(id))
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
saveXML(TheSTACEYxmlTree$value(), file="test.xml")

Rprof(NULL)

sink(file="test.txt")
str(TheAnalysisStructure)
sink(NULL)

#str(TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$prior)


