




add.data <- function() {
  for (a in 1:get.num.alignments()) {
    open.xmlnode("data", attrs=c(name="alignment", id=partitiondataID.a(a)))
    for (q in 1:length(TheAlignmentData$all.sequences[[a]])) {
      txnam <- TheAlignmentData$all.taxonnames[[a]][q]
      sq <- TheAlignmentData$all.sequences[[a]][q]
      attrs <- c(id=sequenceID.a(a, txnam), taxon=txnam, totalcount="4", value=sq)
      add.xmlnode("sequence", attrs=attrs)
    }
    close.xmlnode()
  }
}


