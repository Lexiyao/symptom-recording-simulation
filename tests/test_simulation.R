# Invariants that make the simulation interpretable.
# Run: Rscript tests/test_simulation.R

source("R/simulate_recording_bias.R")

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg, call. = FALSE)
  cat("pass:", msg, "\n")
}

set.seed(1)

# --- structural: the observed code can never exceed the true symptom ---------
d <- simulate_cohort(sim_params(n = 50000, mnar_strength = 1.5))
ok(all(d$S_obs <= d$S), "no symptom is coded unless it is truly present (specificity = 1)")
ok(mean(d$S_obs[d$S == 1]) < 1, "some truly present symptoms go uncoded (sensitivity < 1)")

# --- recording is informative about Y, and only U explains why ---------------
# Y is generated from S and U only, with no arrow from R. But U drives both R
# and Y, so conditional on S alone, R IS predictive of Y -- that is precisely
# what makes an absent code informative. Once U is held fixed the association
# must vanish; if it does not, R has leaked into the outcome and the DAG in the
# README is wrong.
m_noU <- glm(Y ~ S + R + age_z,     family = binomial, data = d)
m_U   <- glm(Y ~ S + R + age_z + U, family = binomial, data = d)
ok(summary(m_noU)$coefficients["R", "Pr(>|z|)"] < 0.01,
   "without U, recording predicts Y -- an absent code carries information")
ok(summary(m_U)$coefficients["R", "Pr(>|z|)"] > 0.01,
   "conditioning on U removes it -- U is the only route from R to Y")

# under MAR the association must be absent even without conditioning on U
m_mar <- glm(Y ~ S + R + age_z, family = binomial,
             data = simulate_cohort(sim_params(n = 50000, mnar_strength = 0)))
ok(summary(m_mar)$coefficients["R", "Pr(>|z|)"] > 0.01,
   "at mnar_strength = 0 recording carries no information about Y")

# --- MNAR strength must actually make coding depend on severity -------------
d0 <- simulate_cohort(sim_params(n = 50000, mnar_strength = 0))
d2 <- simulate_cohort(sim_params(n = 50000, mnar_strength = 2))
gap <- function(x) mean(x$U[x$S == 1 & x$S_obs == 1]) - mean(x$U[x$S == 1 & x$S_obs == 0])
ok(abs(gap(d0)) < 0.05, "at mnar_strength = 0, coded and uncoded cases have equal severity (MAR)")
ok(gap(d2) > 0.5, "at mnar_strength = 2, coded cases are substantially more severe (MNAR)")

# --- the headline claim: lower coding probability attenuates more ------------
set.seed(20260726)
a <- run_sweep(mnar_grid = 0, base_coding_prob = 0.85, nsim = 40)
v <- run_sweep(mnar_grid = 0, base_coding_prob = 0.45, nsim = 40)
ok(v$or_standard < a$or_standard,
   "under-coding attenuates the odds ratio even with coding at random")
ok(a$or_standard < exp(2.34),
   "even the well-coded regime is attenuated relative to the declared true OR")

# --- reproducibility --------------------------------------------------------
r1 <- run_sweep(mnar_grid = c(0, 1), nsim = 20, seed = 99)
r2 <- run_sweep(mnar_grid = c(0, 1), nsim = 20, seed = 99)
ok(isTRUE(all.equal(r1, r2)), "same seed gives identical results")

cat("\nAll checks passed.\n")
