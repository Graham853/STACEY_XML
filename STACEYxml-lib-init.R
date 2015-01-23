

add.init <- function(A) {
  add.comment("Initialisation. First makes a random smcTree of fixed height. It  must be fixed ")
  add.comment("so that BirthDeathCollapseModel can set a bigger initial originHeight.")
  st.attrs <- c(id=initialsmctreeID(), spec="RandomTree", estimate="false", 
             taxonset="@taxonSetOfSets", initial="@smcTree", rootHeight="0.05")
  add.opennode("init", attrs=st.attrs)
  sm.attrs <- c(id=initialsmctreepopmodelID(), spec="ConstantPopulation")
  add.opennode("populationModel", attrs=sm.attrs)
  sp.attrs <- c(id=initialsmctreepopsizeID(), name="popSize")
  add.node("parameter", attrs=sp.attrs, .children="0.05")
  add.closetag()
  add.comment("Next make initial locus trees with random topologies and nodes deeper than smcTree root.")

  gtrees <- get.gtrees()
  for (g in 1:length(gtrees)) {
    gt.attrs <- c(id=initialgtreeID(g), spec="RandomGeneTree", estimate="false", 
               taxonset=IDtoREF(geneTaxonSetID.g(g)), initial=IDtoREF(geneTreeID.g(g)),
               speciesTree=IDtoREF(smcTreeID()))
    add.opennode("init", attrs=gt.attrs)
    gm.attrs <- c(id=initialgtreepopmodelID(g), spec="ConstantPopulation")
    add.opennode("populationModel", attrs=gm.attrs)
    gp.attrs <- c(id=initialgtreepopsizeID(g), name="popSize")
    add.node("parameter", attrs=gp.attrs, .children="0.05")
    add.closetag()
    add.closetag() 
  }
  add.closetag()
}



