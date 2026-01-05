#!/usr/bin/env Rscript
# Spatiotemporal Poisson model with INLA:
#   BYM (spatial) + RW1 (temporal) + IID(time) + IID(space-time interaction)
#
# NOTE: This script is "model-only" and intentionally does NOT read restricted outcome data.
# Provide your own prepared objects locally (not committed to GitHub):
#   - callyears: data.frame with columns: anotias, expected, year, ID (or similar)
#   - matt: adjacency matrix (NxN) matching the spatial units used in callyears

suppressPackageStartupMessages({
  library(INLA)
  library(dplyr)
})

options(mc.cores = max(1, parallel::detectCores() - 1))

# -----------------------------
# User-provided inputs (LOCAL)
# -----------------------------
# You should create these objects upstream in a private script:
#   callyears <- ...
#   matt <- ...
stopifnot(exists("callyears"), exists("matt"))

required_cols <- c("anotias", "expected", "year", "ID")
missing_cols <- setdiff(required_cols, names(callyears))
if (length(missing_cols) > 0) {
  stop("callyears is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# Prepare indexing for INLA
# -----------------------------
callyears <- callyears %>%
  mutate(
    year = as.integer(year),
    ID   = as.integer(ID)
  ) %>%
  arrange(year, ID)

years <- sort(unique(callyears$year))
T <- length(years)

callyears <- callyears %>%
  mutate(
    year3     = year - min(years) + 1L,     # 1..T
    area.year = seq_len(n()),
    re_u      = as.integer(factor(ID)),     # 1..N (spatial index)
    re_v      = re_u,
    re_y      = year3,
    re_y2     = year3
  )

# Basic checks: adjacency dimensions must match spatial units
N <- length(unique(callyears$re_u))
if (!is.matrix(matt) || nrow(matt) != ncol(matt)) stop("matt must be a square matrix.")
if (nrow(matt) != N) {
  stop("Adjacency matrix dimension (", nrow(matt), ") does not match number of spatial units in callyears (", N, ").")
}

# -----------------------------
# Model: BYM + RW1 + IID(time) + IID(space-time)
# -----------------------------
pc_prec <- list(theta = list(
  prec = list(prior = "pc.prec", param = c(0.2 / 0.31, 0.01), initial = 5)
))

formula <- anotias ~ 1 +
  f(re_u, model = "bym", graph = matt, hyper = pc_prec) +
  f(re_y, model = "rw1", hyper = pc_prec) +
  f(re_y2, model = "iid", hyper = pc_prec) +
  f(area.year, model = "iid", hyper = pc_prec)

# Optional: linear combinations for year-specific temporal effect summaries
lcs <- inla.make.lincombs(
  re_y  = diag(T),
  re_y2 = diag(T)
)

fit <- inla(
  formula,
  family = "poisson",
  data = callyears,
  E = expected,
  control.predictor = list(compute = TRUE, link = 1),
  control.compute   = list(dic = TRUE, waic = TRUE),
  lincomb = lcs
)

# -----------------------------
# Save + minimal outputs
# -----------------------------
saveRDS(fit, file = "inla_spacetime_fit.rds")

cat("\nModel fit saved to: inla_spacetime_fit.rds\n")
cat("DIC:  ", fit$dic$dic, "\n")
cat("WAIC: ", fit$waic$waic, "\n")

# Posterior mean of fitted relative risk (E=expected => fitted is RR)
callyears$RR_mean <- fit$summary.fitted.values[, "mean"]
saveRDS(callyears[, c("ID", "year", "RR_mean")], file = "rr_fitted_means.rds")
cat("Fitted RR means saved to: rr_fitted_means.rds\n")
