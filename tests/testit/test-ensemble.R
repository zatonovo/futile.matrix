# :vim set filetype=R
assert('test.max_eigen', {
  en <- Ensemble(50, WishartModel(50, 200))
  mx <- max_eigen(en)
  #theo <- qmp(1, svr=4, var=1)
  #cat("\ntheoretical cutoff:",theo,"\n")
  #print(mx)
  #checkEquals(rep(theo, length(mx)), mx, tolerance=0.1)
})
