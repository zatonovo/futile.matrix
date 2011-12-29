density_lim %when% (model %isa% WignerModel)
density_lim %as% function(model) density_lim(-2,2,100, model)

density_lim %when% (model %isa% WignerModel)
density_lim %as% function(min, max, steps, model)
{
  x <- seq(min, max, length.out=steps)
  sqrt(4 - x^2) / (2 * pi)
}

density_lim %when% (model %isa% WishartModel)
density_lim %as% function(model) density_lim(-3,3,100, model)

density_lim %when% (isa(WishartModel,model))
density_lim %as% function(min, max, steps, model)
{
  x <- seq(min, max, length.out=steps)
  c <- model$n / model$m
  b.neg <- (1 - sqrt(c))^2
  b.pos <- (1 + sqrt(c))^2
  ind <- ifelse(b.neg < x & x < b.pos, 1, 0)
  sqrt(ind * (x - b.neg) * (b.pos - x)) / (2*pi*x*c)
}

density_lim %when% (model %isa% JacobiModel)
density_lim %as% function(model) density_lim(-1,1,100, model)

density_lim %when% (model %isa% JacobiModel)
density_lim %as% function(min, max, steps, model)
{
  lg <- getLogger("futile.matrix")
  lg(WARN, "This function is incomplete")
  x <- seq(min, max, length.out=steps)
  c1 <- model$n / model$m1
  c2 <- model$n / model$m2
  c0 <- c1*x + x^3*c1 - 2*c1*x^2 - c2*x^3 + c2*x^2

  b0 <- c1*x - c2*x - c1 + 2
  b1 <- -2*c2*x^2 + 2*x - 3*c1*x + c1 + c2*x - 1 + 2*c1*x^2
  b2 <- c0

  #num <- (1 - c1)^2
  #b.neg <- num / (c1^2 - c1 + 2 + c2 - c1 * c2 + 2 * sqrt(c1 + c2 - c1 * c2))
  #b.pos <- num / (c1^2 - c1 + 2 + c2 - c1 * c2 - 2 * sqrt(c1 + c2 - c1 * c2))

  #ind <- ifelse(b.neg < x & x < b.pos, 1, 0)
  #sqrt(ind * (x - b.neg) * (x + b.neg)) / (2 * pi * c0)
  sqrt(4*b2*b0 - b1^2) / (2*pi*b2)
}

