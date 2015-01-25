# top level functions for user to call
#
#
#




TheAlignmentData  <- NULL
TheSTACEYxmlTree <- NULL


xmlTree.from.analysis.structure <- function(A)
{
  TheAlignmentData <<- load.data(A) 
  
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
  TheSTACEYxmlTree <<- suppressWarnings(xmlTree(tag="beast", attrs = attrs))
  # I suppress warnings because I keep getting 'empty XML document' warning.
  # Online search suggests this is a bug
  add.data(A)
  add.node.children("map", attrs=c(name="Beta"), children=c("beast.math.distributions.Beta"))
  add.node.children("map", attrs=c(name="Exponential"), children=c("beast.math.distributions.Exponential"))
  add.node.children("map", attrs=c(name="InverseGamma"), children=c("beast.math.distributions.InverseGamma"))
  add.node.children("map", attrs=c(name="LogNormal"), children=c("beast.math.distributions.LogNormalDistributionModel"))
  add.node.children("map", attrs=c(name="Gamma"), children=c("beast.math.distributions.Gamma"))
  add.node.children("map", attrs=c(name="Uniform"), children=c("beast.math.distributions.Uniform"))
  add.node.children("map", attrs=c(name="LaplaceDistribution"), children=c("beast.math.distributions.LaplaceDistribution"))
  add.node.children("map", attrs=c(name="OneOnX"), children=c("beast.math.distributions.OneOnX"))
  add.node.children("map", attrs=c(name="Normal"), children=c("beast.math.distributions.Normal"))
  add.node.children("map", attrs=c(name="prior"), children=c("beast.math.distributions.Prior"))
  add.run(A)
  add.closetag() #### TODO ????
}




load.data <- function(A) {
  alignments <- A$alignment.table$alignments
  all.sequences <- NULL
  all.taxonnames <- NULL
  for (a in 1:length(alignments)) {
    filetype <- guess.filetype(a)
    fpath <- paste0(A$data.dpath, "/", alignments[[a]]$file)
    if (filetype == "nex") {
      nexdata <- read.nexus.data(fpath)
      sequences <- rep("", length(nexdata))
      taxonnames <- rep("", length(nexdata)) 
      lengthfirstseq <- nchar(paste(nexdata[[1]], collapse=""))
      for (q in 1:length(nexdata)) {
        sequences[q] <- paste(nexdata[[q]], collapse="")
        stopifnot(lengthfirstseq==nchar(sequences[q]))
        taxonnames[q] <- names(nexdata)[q]
      }
    } else if (filetype == "phy") {
      phydata <- as.character(read.dna(fpath, format="sequential"))
      sequences <- rep("", nrow(phydata))
      taxonnames <- rep("", nrow(phydata))
      for (q in 1:nrow(phydata)) {
        sequences[q] <- paste(phydata[q,], collapse="")
        taxonnames[q] <- rownames(phydata)[q]
      }
    } else {
      stop()
      #TODO fasta...
    }
    all.sequences <- c(all.sequences, list(sequences))
    all.taxonnames <- c(all.taxonnames, list(taxonnames))
  }
  list(all.taxonnames=all.taxonnames, all.sequences=all.sequences)
}



guess.filetype <- function(a) {
  fname <- TheAnalysisStructure$alignment.table$alignments[[a]]$file
  parts <- str_split(fname, fixed("."))[[1]]
  nparts <- length(parts)
  if (nparts == 1) {
    type <- ""
  } else {
    type <- parts[nparts]
  }
  if (str_detect(type, ignore.case("nex"))) {
    return ("nex")
  } else if (str_detect(type, ignore.case("fas")) || 
               str_detect(type, ignore.case("fasta"))) {
    return ("fas")
  } else {
    return ("phy")
  }
}
