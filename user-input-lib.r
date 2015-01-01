

TheAnalysisStructure <- NULL


TheAnalysis <- function(id, data.dpath, alignment.table, taxa.table, run.options) {
  stopifnot(is.character(data.dpath))
  stopifnot(length(data.dpath)==1)
  stopifnot(is.list(alignment.table))
  stopifnot(is.list(run.options))
  taxa.table <- as.matrix(taxa.table)
  stopifnot(ncol(taxa.table)==2)
  
  aug.alignments <- list()
  for (a in 1:length(alignment.table)) {
    row <- alignment.table[[a]]
    stopifnot(is.list(row))
    stopifnot(length(row)==5)
    stopifnot(names(row)[[1]]=="file")
    stopifnot(names(row)[[2]]=="gtree")
    stopifnot(names(row)[[3]]=="clock")
    stopifnot(names(row)[[4]]=="subst")
    stopifnot(names(row)[[5]]=="site")

    locus <- c(kind="Locus", id=names(alignment.table)[a], row)
    aug.alignments <- c(aug.alignments, list(locus))
  }
  TheAnalysisStructure <<- list(kind="TheAnalysisStructure",
                                id=id, 
                                data.dpath=data.dpath,
                                alignment.table=list(kind="AlignmentTable", id="alignment.table", alignments=aug.alignments),
                                taxa.table=list(kind="TaxaTable", id="taxa.table", taxa.table=taxa.table),
                                run.options=list(kind="RunOptions", id="run.options", run.options=run.options)
                                )
  
  for (i in 1:nof.alignments()) {
    Gtree(paste(i), SMCtree("smctree"))
  }
}




sub.locate <- function(kind, id, want.unfilled, lst, indices) {
  if (!is.list(lst)) {
    #cat("1.", indices, "\n")
    return (c(indices,-1))
    }
  n <- length(lst)
  if (n == 0) {
    #cat("2.", indices, "\n")
    return (c(indices,-1))
    }
  got <- identical(lst$kind, kind)  &&  identical(lst$id, id)
  if (want.unfilled) {
    filled <- TRUE
    for (i in length(lst)) {
      if (is.null(lst[[i]])) { filled <- FALSE } 
      }
    if (filled) {
      got <- FALSE
      }
    } 
  if (got) {
    #cat("3.", indices, "\n")
    return (indices)
    }
  for (i in 1:n) {
    indices <- c(indices,i)
    indices <- sub.locate(kind, id, want.unfilled, lst[[i]], indices)
    #cat("4.", indices, "\n")
    m <- length(indices)
    if (m > 0  &&  indices[m] > 0) {
      #cat("5.", indices, "\n")
      return (indices)
      }
    indices <- indices[c(-(m-1),-m)]   
    }
  #cat("6.", indices, "\n")
  return (c(indices,-1))
}


locate.first.unfilled <- function(kind, id) {
  sub.locate(kind, id, TRUE, TheAnalysisStructure$alignment.table, integer(0))
}






kind.id <- function(lst) {
  paste0(lst$kind,"(",lst$id,")")
  }



check.unused <- function(lst) {
  stopifnot(is.list(lst))
  n <- length(lst)
  lst.names <- names(lst)
  for (i in 1:n) {
    if (!identical(lst.names[i],"kind")  && !identical(lst.names[i],"id")) {
      if (!is.null(lst[[i]])) {
        stop(paste0(kind.id(lst), " already filled."))
        }
      }
    }
  }
 

    
insert.first.unfilled  <- function(x) {
  i <- locate.first.unfilled(x$kind, x$id)
  m <- length(i)
  if (m <= 0  ||  i[m] < 0) {
    return (FALSE)
    }
   
  # I don't know how to do this properly. It shows.
  stopifnot(m <= 10)
  if (m == 1) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]])
    TheAnalysisStructure$alignment.table[[i[1]]] <<- x
    }
  if (m == 2) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]] <<- x
    }  
  if (m == 3) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]] <<- x
    }  
  if (m == 4) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]] <<- x
    }  
  if (m == 5) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]] <<- x
    }  
  if (m == 6) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]] <<- x
    }  
  if (m == 7) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]] <<- x
    }  
  if (m == 8) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]][[i[8]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]][[i[8]]] <<- x
    }  
  if (m == 9) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]][[i[8]]][[i[9]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]][[i[8]]][[i[9]]] <<- x
    }  
  if (m == 10) {
    check.unused(TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]][[i[8]]][[i[9]]][[i[10]]])
    TheAnalysisStructure$alignment.table[[i[1]]][[i[2]]][[i[3]]][[i[4]]][[i[5]]][[i[6]]][[i[7]]][[i[8]]][[i[9]]][[i[10]]] <<- x
    } 
  TRUE
}


insert  <- function(x) {
  ok <- insert.first.unfilled(x)
  if (!ok) {
    stop(paste0(kind.id(x), " has nowhere to go"), call.=FALSE)
    }  
  while (ok) {
    ok <- insert.first.unfilled(x)
    }
  }
      

real.kids <- function(x) {
  stopifnot(is.list(x))
  stopifnot(length(x) > 2)
  rk <- FALSE
  rk2 <- TRUE
  for (i in 3:length(x)) {
    if (is.null(x[[i]])) {
      rk2 <- FALSE
    } else {
      rk <- TRUE
    }
  }
  if (rk != rk2) {
    stop(paste0(kind.id(x), " is incomplete."))
  }
  rk
}


######################################################################


make.ki <- function(A, kinds, ids) {
  if (is.list(A)) {
    if (is.character(A$kind)  &&  is.character(A$id)) {
      wk <- which(kinds == A$kind)
      wi <- integer(0)
      if (length(wk) > 0) {
        wi <- which(ids[wk] == A$id)
      }
      if (length(wi) == 0) {
        kinds <- c(kinds, A$kind)
        ids <- c(ids, A$id)
      }
    }
    for (i in 1:length(A)) {   
      ki <- make.ki(A[[i]], kinds, ids)
      kinds <- ki$kinds
      ids <- ki$ids      
    }
  } 
  list(kinds=kinds, ids=ids)
}




######################################################################



Gtree <- function(id, prior=NULL) {
  x <- list(kind="Gtree", id=id, prior=prior)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}


SMCtree <- function(id, branching=NULL, smccoal=NULL) {
  x <- list(kind="SMCtree", id=id, branching=branching, smccoal=smccoal)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}



BDCPrior <- function(id, growthrate=NULL, reldeath=NULL, w=NULL, eps=NULL) {
  x <- list(kind="BDCPrior", id=id, growthrate=growthrate, reldeath=reldeath, w=w, eps=eps)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}  



SMCCoalescent <- function(id, invgammamix=NULL, popSF=NULL) {
  x <- list(kind="SMCCoalescent", id=id, invgammamix=invgammamix, popSF=popSF)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}    



Clock <- function(id, clock=NULL) {
  x <- list(kind="Clock", id=id, clock=clock)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
  }



StrictClock <- function(id, rate=NULL) {
  x <- list(kind="StrictClock", id=id, rate=rate)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
  }



Subst <- function(id, model=NULL) {
  x <- list(kind="Subst", id=id, model=model)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}
  
  
HKY <- function(id, freqs=NULL, kappa=NULL) {
  x <- list(kind="HKY", id=id, freqs=freqs, kappa=kappa)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}
  
  
SiteHet <- function(id, model=NULL) {
  x <- list(kind="SiteHet", id=id, model=model)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}




InvGammaMix <- function(id, weights=NULL, alphas=NULL, betas=NULL) {
  x <- list(kind="InvGammaMix", id=id, weights=weights, alphas=alphas, betas=betas)
  if (real.kids(x)) {
    n <- length(weights)
    if (n != length(alphas)  ||  n != length(betas)) {
      stop("Vectors weights, alphas, and betas must have same dimension in InvGammaMix.")
    }
    insert(x)  
  }
  invisible(x)
}  




MeanOneGamma <- function(id, alpha=NULL) {
  x <- list(kind="MeanOneGamma", id=id, alpha=alpha)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
} 



LogNorm <- function(id, meanlog=NULL, sdlog=NULL) {
  x <- list(kind="LogNorm", id=id, meanlog=meanlog, sdlog=sdlog)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}

Uniform <- function(id, lower=NULL, upper=NULL) {
  x <- list(kind="Uniform", id=id, lower=lower, upper=upper)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}


Dirichlet <- function(id, mean=NULL, alpha=NULL) {
  x <- list(kind="Dirichlet", id=id, mean=mean, alpha=alpha)
  if (real.kids(x)) {
    insert(x)  
  }
  invisible(x)
}


