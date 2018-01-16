# for diagonal matrix normal elementwise exp seems to work?

Alpha <- matrix(rnorm(100^2), 100, 100)

res <- list(exp(Alpha),
            expm::expm(Alpha),
            expm::expm(Alpha, method="Higham08"),
            expm::expm(Alpha, method="AlMohy-Hi09"),
            expm::expm(Alpha, method="Ward77"),
            expm::expm(Alpha, method="PadeRBS"),
            expm::expm(Alpha, method="Pade"),
            expm::expm(Alpha, method="PadeO"),
            expm::expm(Alpha, method="TaylorO"),
            expm::expm(Alpha, method="R_Eigen"),
            expm::expm(Alpha, method="R_Pade"),
            expm::expm(Alpha, method="R_Ward77"),
            expm::expm(Alpha, method="hybrid_Eigen_Ward"))

checks <- sapply(1:length(res), function(i) {
  max(abs(res[[i]]-res[[2]]))
})

names(checks) <- c("exp", "Higham08.b", "Higham08", "AlMohy-Hi09", "Ward77", "PadeRBS", "Pade", "PadeO", "TaylorO", "R_Eigen", "R_Pade", "R_Ward77", "hybrid_Eigen_Ward")
microbenchmark::microbenchmark(
  exp(Alpha),
  expm::expm(Alpha),
  expm::expm(Alpha, method="Higham08"),
  expm::expm(Alpha, method="AlMohy-Hi09"),
  expm::expm(Alpha, method="Ward77"),
  expm::expm(Alpha, method="PadeRBS"),
  expm::expm(Alpha, method="Pade"),
  expm::expm(Alpha, method="PadeO"),
  expm::expm(Alpha, method="TaylorO"),
  expm::expm(Alpha, method="R_Eigen"),
  expm::expm(Alpha, method="R_Pade"),
  expm::expm(Alpha, method="R_Ward77"),
  expm::expm(Alpha, method="hybrid_Eigen_Ward"),
  expm(Alpha))

