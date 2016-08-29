library(microbenchmark)

f = function(n) {
  H = matrix(1,n,n)
  H[upper.tri(H)] = 0
  Hinv = solve(H)
  Hinv
}

g = function(n) {
  Ginv = diag(n)
  idx = row(Ginv)==(col(Ginv)+1)
  Ginv[idx] = -1
  Ginv
}

microbenchmark(f(50), g(50))
