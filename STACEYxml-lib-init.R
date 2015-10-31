

add.init <- function(A) {
  add.comment("Initialisation. First makes a random smcTree of fixed height. It  must be fixed ")
  add.comment("so that BirthDeathCollapseModel can set a bigger initial originHeight.")
  st.attrs <- c(id=initialsmctreeID(), spec="RandomTree", estimate="false", 
             taxonset="@taxonSetOfSets", initial=IDtoREF(smcTreeID()), rootHeight="0.05")
  open.xmlnode("init", attrs=st.attrs)
  sm.attrs <- c(id=initialsmctreepopmodelID(), spec="ConstantPopulation")
  open.xmlnode("populationModel", attrs=sm.attrs)
  sp.attrs <- c(id=initialsmctreepopsizeID(), name="popSize")
  add.xmlnode.children("parameter", attrs=sp.attrs, children="0.05")
  close.xmlnode()
  close.xmlnode()
  
  add.comment("Next make initial locus trees with random topologies and nodes deeper than smcTree root.")
  gtrees <- get.GTrees()
  for (u in 1:length(gtrees)) {
    gt.attrs <- c(id=initialgtreeID(u), spec="RandomGeneTree", estimate="false", 
               taxonset=IDtoREF(geneTaxonSetID.u(u)), initial=IDtoREF(geneTreeID.u(u)),
               speciesTree=IDtoREF(smcTreeID()))
    open.xmlnode("init", attrs=gt.attrs)
    gm.attrs <- c(id=initialgtreepopmodelID(u), spec="ConstantPopulation")
    open.xmlnode("populationModel", attrs=gm.attrs)
    gp.attrs <- c(id=initialgtreepopsizeID(u), name="popSize")
    add.xmlnode.children("parameter", attrs=gp.attrs, children="0.05")
    close.xmlnode()
    close.xmlnode() 
  }
}



