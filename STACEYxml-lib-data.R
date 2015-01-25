




add.data <- function(A) {
  alignments <- A$alignment.table$alignments
  for (a in 1:length(alignments)) {
    add.opennode("data", attrs=c(name="alignment", id=partitiondataID.a(a)))
    for (q in 1:length(TheAlignmentData$all.sequences[[a]])) {
      txnam <- TheAlignmentData$all.taxonnames[[a]][q]
      sq <- TheAlignmentData$all.sequences[[a]][q]
      attrs <- c(id=sequenceID.a(a, txnam), taxon=txnam, totalcount="4", value=sq)
      add.node("sequence", attrs=attrs)
    }
    add.closetag()
  }
}


