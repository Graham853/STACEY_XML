



load.data <- function() {
  all.sequences <- NULL
  all.taxonnames <- NULL
  fpaths <- get.data.fpaths()
  for (a in 1:length(fpaths)) {
    filetype <- guess.filetype(fpaths[a])
    if (filetype == "nex") {
      nexdata <- read.nexus.data(fpaths[a])
      sequences <- rep("", length(nexdata))
      taxonnames <- rep("", length(nexdata)) 
      lengthfirstseq <- nchar(paste(nexdata[[1]], collapse=""))
      for (q in 1:length(nexdata)) {
        sequences[q] <- paste(nexdata[[q]], collapse="")
        stopifnot(lengthfirstseq==nchar(sequences[q]))
        taxonnames[q] <- names(nexdata)[q]
      }
    } else if (filetype == "phy") {
      phydata <- as.character(read.dna(fpaths[a], format="sequential"))
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



guess.filetype <- function(fpath) {
  parts <- str_split(fpath, fixed("."))[[1]]
  nparts <- length(parts)
  if (nparts == 1) {
    type <- ""
  } else {
    type <- parts[nparts]
  }
  if (str_detect(type, fixed("nex", ignore_case = TRUE))) {
    return ("nex")
  } else if (str_detect(type, fixed("nex", ignore_case = TRUE)) || 
               str_detect(type, fixed("nex", ignore_case = TRUE))) {
    return ("fas")
  } else {
    return ("phy")
  }
}

