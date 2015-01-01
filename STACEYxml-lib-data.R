

add.data <- function(A) {
  # TODO read from nex and A
  add.opennode("data", attrs=c(name="alignment", id = "gene.1"))
  add.node("sequence", attrs=c(id="b01.g1", value="?", totalcount="4", taxon="b01"))
  add.node("sequence", attrs=c(id="b02.g1", value="?", totalcount="4", taxon="b02"))
  add.node("sequence", attrs=c(id="b03.g1", value="?", totalcount="4", taxon="b03"))
  add.closetag()
}





