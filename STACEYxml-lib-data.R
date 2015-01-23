
#<sequence id="seq_b05_A1" taxon="b05_A" totalcount="4" value="aaa

add.data <- function(A) {
  # TODO read the actual nex files. taxons from files
  alignments <- A$alignment.table$alignments
  for (a in 1:length(alignments)) {
    xxx <- read.nexus.data(alignments[[a]]$file)
    taxonnames <- names(xxx)
    add.opennode("data", attrs=c(name="alignment", id=partitiondataID.a(a)))
    for (q in 1:length(xxx)) {
      seq <- paste(xxx[[q]], collapse="")
      taxonname <- taxonnames[q]
      add.node("sequence", attrs=c(id=sequenceID.a(a, taxonname), taxon=taxonname, totalcount="4", value=seq))
    }
    add.closetag()
  }

}


