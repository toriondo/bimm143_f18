#new function
add <- function(x, y=1) {
  x+y
}

rescale <- function(x) {
  rng <-range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}
