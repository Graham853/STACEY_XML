
nof.alignments <- function() {
  length(TheAnalysisStructure$alignment.table$alignments)
}



get.taxa.table <- function() {
  TheAnalysisStructure$taxa.table$taxa.table
}


get.minclusters <- function() {
  tt <- get.taxa.table()
  mcs <- tt[,2]
  unique(mcs)
}


get.bdc.model <- function() {
  # TODO not good that access first gtree for all
  TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$prior$branching
}


get.bdcm.growthrate.fixed <- function() {
  (get.bdc.model()$growthrate$kind == "FixedValue")
}

get.bdcm.reldeath.fixed <- function() {
  (get.bdc.model()$reldeath$kind == "FixedValue")
}


get.bdcm.w.fixed <- function() {
  (get.bdc.model()$w$kind == "FixedValue")
}



get.smc.coalescent <- function() {
  # TODO not good that access first gtree for all
  TheAnalysisStructure$alignment.table$alignments[[1]]$gtree$prior$smccoal
}

get.smct.popsf.fixed <- function() {
  (get.smc.coalescent()$popSF$kind == "FixedValue")
}





get.gtrees <- function() {
  gtrees <- NULL
  gtreeids <- character(0)
  for (a in 1:nof.alignments()) {
    gtree <- TheAnalysisStructure$alignment.table$alignments[[a]]$gtree
    if (length(which(gtreeids==gtree$id)) == 0) {
      gtreeids <- c(gtreeids, gtree$id)
      gtree.aug <- c(gtree, list(first.partition=a))
      gtrees <- c(gtrees, list(gtree.aug))
    }
  }
  gtrees
}




get.clocks <- function() {
  clocks <- NULL
  clockids <- character(0)
  for (a in 1:nof.alignments()) {
    clock <- TheAnalysisStructure$alignment.table$alignments[[a]]$clock
    if (length(which(clockids==clock$id)) == 0) {
      clockids <- c(clockids, clock$id)
      clocks <- c(clocks, list(clock))
    }
  }
  clocks
}



get.siteMs <- function() {
  siteMs <- NULL
  siteMids <- character(0)
  for (a in 1:nof.alignments()) {
    siteM <- TheAnalysisStructure$alignment.table$alignments[[a]]$siteM
    if (length(which(siteMids==siteM$id)) == 0) {
      siteMids <- c(siteMids, siteM$id)
      siteM.aug <- c(siteM, list(first.partition=a))
      siteMs <- c(siteMs, list(siteM.aug))        
    }
  }
  siteMs
}



