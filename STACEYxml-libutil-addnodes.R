
add.node <- function(...) {
  TheSTACEYxmlTree$addNode(...)
}


add.opennode <- function(...) {
  TheSTACEYxmlTree$addNode(..., close=FALSE)
}

add.closetag <- function() {
  TheSTACEYxmlTree$closeTag()
}




add.comment <- function(txt) {
  TheSTACEYxmlTree$addComment(txt)
}



add.bigcomment <- function(txt) {
  fullline <- paste0(rep("~", 45), collapse="") 
  halftextlen <- as.integer(round(nchar(txt)/2))
  if (halftextlen > 22) { halftextlen <- 22}
  partline <- paste0(rep("~", 23 - halftextlen), collapse="")
  TheSTACEYxmlTree$addComment(fullline)
  TheSTACEYxmlTree$addComment(paste0(partline, " ", txt, " ", partline))
  TheSTACEYxmlTree$addComment(fullline)
}


add.hugecomment <- function(txt) {
  xline <- paste0(rep("~", 45), collapse="")
  emptyline <- paste0(rep("~", 45), collapse="")
  partline <- paste0(rep("~", 23 - as.integer(round(nchar(txt)/2))), collapse="")
  TheSTACEYxmlTree$addComment(xline)
  TheSTACEYxmlTree$addComment(emptyline)
  TheSTACEYxmlTree$addComment(emptyline)
  TheSTACEYxmlTree$addComment(paste0(partline, " ", txt, " ", partline))
  TheSTACEYxmlTree$addComment(emptyline)
  TheSTACEYxmlTree$addComment(emptyline)
  TheSTACEYxmlTree$addComment(xline)
}




