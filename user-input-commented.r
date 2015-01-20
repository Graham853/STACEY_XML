setwd("C:/Users/Work/AAA/Programming/biology/STACEY_XML")
source("user-input-lib.r")
source("STACEYxml-lib.r")



##### DATA 

# This is like the first Beauti tab:
alignment.table <- list(
  partition1=list(file="locus1.nex", gtree=Gtree("1"), clock=Clock("1"),     siteM=SiteModel("1")),
  partition2=list(file="locus2.nex", gtree=Gtree("2"), clock=Clock("2and3"), siteM=SiteModel("2and3")),
  partition3=list(file="locus3.nex", gtree=Gtree("3"), clock=Clock("2and3"), siteM=SiteModel("2and3"))
)

# This is like the second Beauti tab, for a *BEAST/DISSECT/STACEY analysis:
taxa.table <- rbind(
c("a1", "a"),
c("a2", "a"),
c("a3", "a"),
c("b1", "b"),
c("b2", "b")
)
colnames(taxa.table) <- c("taxon", "mincluster")


##### RUN OPTIONS

# This is like the last Beauti tab:
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


##### INTIALISE THE ANALYSIS

TheAnalysis("Test.1", data.dpath="C:/...", alignment.table, taxa.table, run.options)
# This makes empty structures for the gene trees, clocks and site models.
# These now need to be filled in recursively.

# The general mechanism is: 
# 1. make empty sub-structures by calling functions with 
# a single argument (which is an ID):
#      Object(ID, SubObject1(ID1), SubObject2(ID2))
# 2. fill in the sub-structures in simailar fashion:
#      SubObject1(ID1, ...)
#      SubObject2(ID2, ...)
# 3. The recursion bottoms out with actual values, usually real 
# numbers or vectors of real numbers, sometimes character strings 
# or TRUE/FALSE.

# Note that every creation function is called twice, once with just
# an ID, later with the same ID and one or more subobjects or values.
# The first calls to Gtree(), Clock(), SiteModel() have been done.

# The structure is stored in the global variable TheAnalysisStructure.
# It is a list of list of lists... In other words it is a tree, with values
# at the tips. More precisely, it is a multi-labeled tree. The same values
# can appear at more than one tip, and quite big subtrees can appear more than once.
# for example, in the multispecies coalescent model, all gene trees share the same prior,
# so the structure representing the multispecies coalescent appears for each gene tree.

# An equivalent way of representing the structure mathematically is as a directed acyclic
# graph, or DAG. A DAG would be more compact but trickier to manipulate. 

######  MODEL

gtrees <- get.gtrees()
for (i in 1:length(gtrees)) {
  id <- gtrees[[i]]$id
  Gtree(id, SMCtree("smctree"), BranchRM(id))
  BranchRM(id, StrictClock(id))
  StrictClock(id, 1.0)
}
# get.gtrees() is a convenience function. Here it is used to get the 
# unique IDs of the gene trees in alignment.table. ("1", "2", "3" here.)
# Empty gene trees have been made by the call to TheAnalysis()
# and now their details need to be filled in. 

# A gene tree needs a prior and a branch rate model. 
# The call Gtree(id, SMCtree("smctree"), BranchRM(id))
# makes a place for these: it creates an empty STACEY 
# multispecies coalescent model (same ID for all genes trees)
# and an empty branch rate model (unique ID for each tree). Now we
# have to fill in the branch rate models and the prior.

# The branch rate models are done in the same loop, first creating 
# empty strict clock model for each, then filling it in with a 
# rate (1.0 here). The SMC-tree is filled in next.

SMCtree("smctree", BDCPrior("smctree"), SMCCoalescent("staceycoal"))
# The SMC-tree needs a prior, a birth-death-collapse model, and a 
# coalescent model. The above makes empty versions.

BDCPrior("smctree", growthrate=LogNorm("smctree.g"), 
         reldeath=Uniform("smctree.rd"), 
         w=Uniform("smctree.w"), 
         eps=0.00003)
LogNorm("smctree.g", 4.6, 2)
Uniform("smctree.rd", 0, 1)
Uniform("smctree.w", 0, 1) 
# The above fills in the birth-death-collapse model.

SMCCoalescent("staceycoal", invgammamix=InvGammaMix("popBV"), popSF=LogNorm("popSF"))
LogNorm("popSF", -7, 2)
InvGammaMix("popBV", weights=c(0.5,0.5), alphas=c(3,3), betas=c(1.5,2.5))
# The above fills in the coalescent model. "popBV" is the ID for the 
# distribution representing branch-to-branch varibility of the population
# size. "popSF" the ID for the overall population scaling factor.
# The gene trees are now done. Clocks and site models remain.

Clock("1", PartitionRate("1"))
PartitionRate("1", 1.0)
Clock("2and3", PartitionRate("2and3"))
PartitionRate("2and3", LogNorm("clock.2and3"))
LogNorm("clock.2and3",0,1)
# Clocks deal with the variability of rates between partitions.
# The above sets the first one to 1.0, and the other is estimated
# and has a log-normal prior.

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
# each site model needs a substitution models and a site rate heterogeneity model.
# The above sets no site rate heterogeneity and HKY substitution models.
# The UniformUnitSimplex provides the standard prior for the base frequencies.

#####  MAKE AND SAVE THE XML 

xmlTree.from.analysis.structure(TheAnalysisStructure)
saveXML(TheSTACEYxmlTree$value(), file="test.xml")


##### TESTING, VIEWING

sink(file="test.txt")
str(TheAnalysisStructure)
sink(NULL)

ki <- make.ki(TheAnalysisStructure,  character(0),  character(0))
cbind(ki$kinds, ki$ids)

TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$id

