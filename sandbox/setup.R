for (p in list.files("src", pattern = "[.]cpp$", full.names = TRUE)) {
  Rcpp::sourceCpp(p)
}
for (p in setdiff(list.files("R", pattern = "[.]R$", full.names = TRUE), "R/em_t.R")) {
  source(p, keep.source = TRUE)
}
source("R/em_t.R", keep.source = TRUE)
