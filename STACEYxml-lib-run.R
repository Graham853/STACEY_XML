

add.run <- function() {
  open.xmlnode("run", attrs=c(id="mcmc", spec="MCMC", 
  														chainLength=get.runoption("chainlength"), 
  														storeEvery=get.runoption("store.every")))
  cat("done run\n")
  add.hugecomment("state")
  add.state()
  cat("done state\n")
  add.hugecomment("init")
  add.init()
  cat("done init\n")
  add.hugecomment("distribution")
  add.distribution()
  cat("done distribution\n")
  add.hugecomment("operators")
  add.operators()
  cat("done operators\n")
  add.hugecomment("loggers")
  add.loggers()
  cat("done loggers\n")
  close.xmlnode()
}

