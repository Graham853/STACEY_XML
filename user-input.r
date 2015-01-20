setwd("C:/Users/Work/AAA/Programming/biology/STACEY_XML")
source("user-input-lib.r")
source("STACEYxml-lib.r")

alignment.table <- list(
  partition1=list(file="xxx.nex", gtree=Gtree("1"), clock=Clock("1"),     siteM=SiteModel("1")),
  partition2=list(file="yyy.nex", gtree=Gtree("2"), clock=Clock("2and3"), siteM=SiteModel("2and3")),
  partition3=list(file="zzz.nex", gtree=Gtree("3"), clock=Clock("2and3"), siteM=SiteModel("2and3"))
)

taxa.table <- rbind(
c("a1", "a"),
c("a2", "a"),
c("a3", "a"),
c("b1", "b"),
c("b2", "b")
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
         reldeath=Uniform("smctree.rd"), 
         w=Uniform("smctree.w"), 
         ohID="origin.height",
         eps=0.00003)
LogNorm("smctree.g", 4.6, 2)
Uniform("smctree.rd", 0, 1)
Uniform("smctree.w", 0, 1) 
SMCCoalescent("staceycoal", invgammamix=InvGammaMix("popBV"), popSF=LogNorm("popSF"))
LogNorm("popSF", -7, 2)
InvGammaMix("popBV", weights=c(0.5,0.5), alphas=c(3,3), betas=c(1.5,2.5))


Clock("1", PartitionRate("1"))
PartitionRate("1", 1.0)
Clock("2and3", PartitionRate("2and3"))
PartitionRate("2and3", LogNorm("clock.2and3"))
LogNorm("clock.2and3",0,1)

SiteModel("1", Subst("1"), SiteHet("1"))
SiteModel("2and3", Subst("2and3"), SiteHet("2and3"))
SiteHet("1", "None")
SiteHet("2and3")
SiteHet("2and3", "None")
Subst("1", HKY("1"))
HKY("1", UniformUnitSimplex("hky.1.freqs"), LogNorm("hky.kappa"))
UniformUnitSimplex("hky.1.freqs", 3)
Subst("2and3", HKY("2and3"))
HKY("2and3", UniformUnitSimplex("hky.2and3.freqs"), LogNorm("hky.kappa"))
UniformUnitSimplex("hky.2and3.freqs", 3)
LogNorm("hky.kappa",1,1.25)


xmlTree.from.analysis.structure(TheAnalysisStructure)
saveXML(TheSTACEYxmlTree$value(), file="test.xml")


sink(file="test.txt")
str(TheAnalysisStructure)
sink(NULL)

str(TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$prior)


