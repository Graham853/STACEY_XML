

add.init <- function(A) {
  # TODO ids as functions
  
  add.comment("Initialisation. First makes a random smcTree of fixed height. It  must be fixed ")
  add.comment("so that BirthDeathCollapseModel can set a bigger initial originHeight.")
  st.attrs <- c(id="randomSpeciesTree", spec="RandomTree", estimate="false", 
             taxonset="@taxonSetOfSets", initial="@smcTree", rootHeight="0.05")
  add.opennode("init", attrs=st.attrs)
  sm.attrs <- c(id="ConstantPopulation.smc", spec="ConstantPopulation")
  add.opennode("populationModel", attrs=sm.attrs)
  sp.attrs <- c(id="randomPopSize.smc", name="popSize")
  add.node("parameter", attrs=sp.attrs, .children="0.05")
  add.closetag()
  add.comment("Next make initial gene trees with random topologies and nodes deeper than smcTree root.")

  for (g in 1:nof.alignments()) {
    gt.attrs <- c(id=paste0("randomGeneTree.", g), spec="RandomGeneTree", estimate="false", 
               taxonset=IDtoREF(geneTaxonSetID(g)), initial=IDtoREF(geneTreeID(g)),
               speciesTree=IDtoREF(smcTreeID()))
    add.opennode("init", attrs=gt.attrs)
    gm.attrs <- c(id=paste0("RGTPopulationModel.", g), spec="ConstantPopulation")
    add.opennode("populationModel", attrs=gm.attrs)
    gp.attrs <- c(id=paste0("RGTPopSize.", g), name="popSize")
    add.node("parameter", attrs=gp.attrs, .children="0.05")
    add.closetag()
    add.closetag() 
  }
  add.closetag()
}

