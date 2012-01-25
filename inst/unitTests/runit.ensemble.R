test.max_eigen <- function()
{
  en <- create(Ensemble, 50, create(WishartModel, 50, 200))
  mx <- max_eigen(en)
  theo <- qmp(1, svr=4, var=1)
  cat("\ntheoretical cutoff:",theo,"\n")
  print(mx)
  checkEquals(rep(theo, length(mx)), mx, tolerance=0.1)
}
