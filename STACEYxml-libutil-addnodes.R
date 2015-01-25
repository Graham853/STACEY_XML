
add.node <- function(name, attrs) {
  TheSTACEYxmlTree$addNode(name=name, attrs=attrs)
}

add.node.children <- function(name, attrs, children) {
  TheSTACEYxmlTree$addNode(name=name, attrs=attrs, .children=children)
}


add.opennode <- function(name, attrs) {
  TheSTACEYxmlTree$addNode(name=name, attrs=attrs, close=FALSE)
}

add.closetag <- function() {
  TheSTACEYxmlTree$closeTag()
}




add.comment <- function(txt) {
  TheSTACEYxmlTree$addComment(txt)
}



add.bigcomment <- function(txt) {
  add.surrounded.comment(txt, 1)
}


add.hugecomment <- function(txt) {
  add.surrounded.comment(txt, 3)
}


add.surrounded.comment <- function(txt, n) {
  fullline <- "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"
  halftextlen <- as.integer(floor(nchar(txt)/2))
  if (halftextlen > 32) { halftextlen <- 32}
  partline <- paste0(rep("~ ", (36 - halftextlen)/2), collapse="")
  for (i in 1:n) {
    TheSTACEYxmlTree$addComment(fullline)
  }
  TheSTACEYxmlTree$addComment(paste0(partline, " ", txt, " ", partline, "~"))
  for (i in 1:n) {
    TheSTACEYxmlTree$addComment(fullline)
  } 
}


