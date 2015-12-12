


library(ape)
library(hash)

RcodeSourceDirectory <- "C:/Users/Work/AAA/Programming/biology/STACEY_XML"

source(paste0(RcodeSourceDirectory, "/", "analysis-structure-lib.r"))
source(paste0(RcodeSourceDirectory, "/", "STACEYxml-lib.r"))

cat(date(),"\n")

###########################################################################

nloci <- 3
data.dpath <- "C:/Users/Work/AAA/Programming/biology/STACEY_XML/Tests"
beastxml.fpath  <- paste0(data.dpath, "/beast.xml")
data.fnames <- paste0("seqs", sprintf("%02d", 1:nloci), ".nex")
names(data.fnames) <- paste0("seq", sprintf("%02d", 1:nloci))

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


run.options <- c(
  sampledgtrees.fpathbase="gtrees",
  sampledsmctrees.fpath="smctrees.txt",
  sampledparams.fpath="params.txt",
  chainlength="10000000",
  store.every="500000",
  params.logevery="1000",
  smctree.logevery="1000",
  gtrees.logevery="1000",
  screen.logevery="10000"
  )          
# note in current implementation, always one alignment per data file
TheAnalysis("Test.1", data.dpath=data.dpath, alignment.table=data.fnames, taxa.table=taxa.table, run.options=run.options)

########################################################################



gtree.priors <- get.GTreePriors()
for (i in 1:length(gtree.priors)) {
  GTreePrior(gtree.priors[[i]]$id, SMCtree("smctree"))
}

SMCtree("smctree", BDCPrior("smctree"), SMCCoalescent("staceycoal"))
BDCPrior("smctree", growthrate=LogNorm("smctree.g"), 
         reldeath=Beta("smctree.rd"), 
         w=Beta("smctree.w"), 
         oh="origin.height",
         eps=0.00003)
LogNorm("smctree.g", 4.6, 2)
Beta("smctree.rd", 1, 1)
Beta("smctree.w", 1, 1) 
SMCCoalescent("staceycoal", invgammamix=InvGammaMix("popBV"), popSF=LogNorm("popSF"))
LogNorm("popSF", -7, 2)
InvGammaMix("popBV", weights=c(0.5,0.5), alphas=c(3,3), betas=c(1.5,2.5))


gtree.brms <- get.gtree.BranchRMs()
for (i in 1:length(gtree.brms)) {
  id <- gtree.brms[[i]]$id
  BranchRM(id, StrictClock(id))
  StrictClock(id, 1.0)
}

gtree.ploidys <- get.gtree.Ploidys()
for (i in 1:length(gtree.ploidys)) {
  Ploidy(gtree.ploidys[[i]]$id, 2)
}



gtree.prms <- get.gtree.PartitionRateMs()
for (i in 1:length(gtree.prms)) {
  id <- gtree.prms[[i]]$id
  if (i==1) {
    PartitionRateM(id, 1.0)
  } else {
    PartitionRateM(id, LogNorm(id))
    LogNorm(id, 0.0, 1.0) 
  }
}


gtree.substs <- get.gtree.SubstMs()
for (i in 1:length(gtree.substs)) {
  id <- gtree.substs[[i]]$id
  SubstM(id, HKY(id))
  freqsID <- paste0("HKYfreqs.", id)
  kappaID <- paste0("HKYkappa.", id)
  HKY(id, UniformUnitSimplex(freqsID), LogNorm(kappaID))
  UniformUnitSimplex(freqsID, 4)
  LogNorm(kappaID, 1.0, 1.25)
}

gtree.sitehets <- get.gtree.SiteHets()
for (i in 1:length(gtree.sitehets)) {
  SiteHet(gtree.sitehets[[i]]$id, "None")
}


###################################################################

cat(date(),"\n")

cat.structure("C:/Users/Work/Desktop/TheAnalysisStructure.txt")


xmlFile.from.analysis.structure(beastxml.fpath)




