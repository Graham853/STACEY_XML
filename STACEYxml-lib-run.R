

add.run <- function(A) {
  add.opennode("run", attrs=c(id="mcmc", spec="MCMC", chainLength="1000000", storeEvery="500000"))
  add.hugecomment("state")
  add.state(A)
  add.hugecomment("init")
  add.init(A)
  add.hugecomment("distribution")
  add.distribution(A)
  add.hugecomment("operators")
  add.operators(A)
  add.hugecomment("loggers")
  add.loggers(A)
  add.closetag()
}

