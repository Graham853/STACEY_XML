
setwd("C:/Users/Work/AAA/Programming/biology/STACEY_XML")
source("user-input-lib.r")
source("STACEYxml-lib.r")




#################################### DATA  #####################################



alignment.table <- list(
locus1=list(file="xxx.nex", gtree=Gtree("1"), clock=Clock("1"),     subst=Subst("1"),     site=SiteHet("all")),
locus2=list(file="yyy.nex", gtree=Gtree("2"), clock=Clock("2and3"), subst=Subst("2and3"), site=SiteHet("all")),
locus3=list(file="zzz.nex", gtree=Gtree("3"), clock=Clock("2and3"), subst=Subst("2and3"), site=SiteHet("all"))
)


taxa.table <- rbind(
c("a1", "a"),
c("a2", "a"),
c("a3", "a"),
c("b1", "b"),
c("b2", "b")
)
colnames(taxa.table) <- c("taxon", "mincluster")

###################################### RUN #####################################

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
          

################################ ANALYSIS STRUCTURE ############################


TheAnalysis("Test.1", data.dpath="C:/...", alignment.table, taxa.table, run.options)


####################################  MODEL  ###################################

SMCtree("smctree", BDCPrior("smctree"), SMCCoalescent("staceycoal"))
BDCPrior("smctree", growthrate=LogNorm("smctree.g"), 
         reldeath=Uniform("smctree.rd"), 
         w=Uniform("smctree.w"), 
         eps=0.00003)
LogNorm("smctree.g", 4.6, 2)
Uniform("smctree.rd", 0, 1)
Uniform("smctree.w", 0, 1) 
SMCCoalescent("staceycoal", invgammamix=InvGammaMix("popBV"), popSF=LogNorm("popSF"))
LogNorm("popSF", -7, 2)
InvGammaMix("popBV", weights=c(0.5,0.5), alphas=c(3,3), betas=c(1.5,2.5))


Clock("1", StrictClock("1"))
StrictClock("1", 1.0)
Clock("2and3", StrictClock("2and3"))
StrictClock("2and3", LogNorm("clock.2and3"))
LogNorm("clock.2and3",0,1)

Subst("1", HKY("1"))
HKY("1", Dirichlet("hky.1.freqs"), LogNorm("hky.kappa"))
Dirichlet("hky.1.freqs", rep(0.25,4), 1.0)
Subst("2and3", HKY("2and3"))
HKY("2and3", Dirichlet("hky.2and3.freqs"), LogNorm("hky.kappa"))
Dirichlet("hky.2and3.freqs", rep(0.25,4), 1.0)
LogNorm("hky.kappa",1,1.25)


sink(file="test.txt")
str(TheAnalysisStructure)
sink(NULL)


xmlTree.from.analysis.structure(TheAnalysisStructure)

saveXML(TheSTACEYxmlTree$value(), file="test.xml")

ki <- make.ki(TheAnalysisStructure,  character(0),  character(0))
cbind(ki$kinds, ki$ids)


