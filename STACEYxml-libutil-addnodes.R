
# commenting out the line in cat0() means 30min not 40min on 500 loci

cat0 <- function(...) {
 cat(..., file=TheBEASTOutputConnection, sep="")
}


cat.indent <- function(n) {
  cat0(paste0(rep("    ", n), collapse=""))
}


add.xmlnode <- function(name, attrs) {
  cat.indent(length(TheCurrentPath))
  cat0("<", name, " ")
  attrnames <- names(attrs)
  for (a in 1:length(attrs)) {
    cat0(attrnames[a], "=\"", attrs[a], "\" ")
  }
  cat0("/>\n")
}

add.xmlnode.children <- function(name, attrs, children) {
  cat.indent(length(TheCurrentPath))
  cat0("<", name, " ")
  attrnames <- names(attrs)
  for (a in 1:length(attrs)) {
    cat0(attrnames[a], "=\"", attrs[a], "\"")
    if (a < length(attrs)) { cat0(" ") }
  }
  cat0(">", children, "</", name, ">\n")
}


open.xmlnode <- function(name, attrs) {
  n <- length(TheCurrentPath)
  cat.indent(n)
  cat0("<", name, " ")
  attrnames <- names(attrs)
  for (a in 1:length(attrs)) {
    cat0(attrnames[a], "=\"", attrs[a], "\"")
    if (a < length(attrs)) { cat0(" ") }
  }
  cat0(">\n")
  
  TheCurrentPath <<- c(TheCurrentPath, name)
}


open.toplevel.xmlnode <- function(name, attrs) {
  n <- length(TheCurrentPath)
  stopifnot(n == 0)
  cat0("<?xml version=\"1.0\"?>\n")
  cat0("<", name, " ")
  attrnames <- names(attrs)
  for (a in 1:length(attrs)) {
    cat0(attrnames[a], "=\"", attrs[a], "\" ")
  }
  cat0(">\n")
  
  TheCurrentPath <<- c(TheCurrentPath, name)
}





close.xmlnode <- function() {
  n <- length(TheCurrentPath)
  cat.indent(n-1)
  cat0("</", TheCurrentPath[n], ">\n")
  TheCurrentPath <<- TheCurrentPath[1:(n-1)]
}




add.comment <- function(txt) {
  n <- length(TheCurrentPath)
  cat.indent(n)
  cat0("<!--", txt, "-->\n")
}


add.centredcomment <- function(txt) {
  add.surrounded.comment(txt, 0)
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
  
  add.n.surrounds <- function(n) {
    if (n >= 1) {
      for (i in 1:n) {
        add.comment(fullline)
      }
    }
  }
  
  if (halftextlen > 32) { halftextlen <- 32}
  partline <- paste0(rep("~ ", (36 - halftextlen)/2), collapse="")
  add.n.surrounds(n)
  add.comment(paste0(partline, " ", txt, " ", partline, "~"))
  add.n.surrounds(n)
  
}


