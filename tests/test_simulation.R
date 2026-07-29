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
# Cancer is rare here (~0.8%), so this needs a large cohort to have power:
# n = 50k yields only ~350 cases and the test is underpowered at p < 0.01.
d_big <- simulate_cohort(sim_params(n = 400000, mnar_strength = 1.5))
m_noU <- glm(Y ~ S + R + age_z,     family = binomial, data = d_big)
m_U   <- glm(Y ~ S + R + age_z + U, family = binomial, data = d_big)
ok(coef(m_noU)["R"] > 0 && summary(m_noU)$coefficients["R", "Pr(>|z|)"] < 0.001,
   "without U, recording predicts Y -- an absent code carries information")
ok(summary(m_U)$coefficients["R", "Pr(>|z|)"] > 0.01,
   "conditioning on U removes it -- U is the only route from R to Y")

# under MAR the association must be absent even without conditioning on U
m_mar <- glm(Y ~ S + R + age_z, family = binomial,
             data = simulate_cohort(sim_params(n = 400000, mnar_strength = 0)))
ok(summary(m_mar)$coefficients["R", "Pr(>|z|)"] > 0.01,
   "at mnar_strength = 0 recording carries no information about Y")

# --- MNAR strength must actually make coding depend on severity -------------
# Compare the severity gap against its own standard error rather than a fixed
# tolerance: with n = 400k only ~32k patients have the symptom, and of those the
# uncoded group is small, so the SE is not negligible.
gap_t <- function(x) {
  a <- x$U[x$S == 1 & x$S_obs == 1]; b <- x$U[x$S == 1 & x$S_obs == 0]
  (mean(a) - mean(b)) / sqrt(var(a)/length(a) + var(b)/length(b))
}
d0 <- simulate_cohort(sim_params(n = 400000, mnar_strength = 0))
d2 <- simulate_cohort(sim_params(n = 400000, mnar_strength = 2))
ok(abs(gap_t(d0)) < 3,
   sprintf("at mnar_strength = 0, coded and uncoded cases have equal severity (t = %.2f)", gap_t(d0)))
ok(gap_t(d2) > 10,
   sprintf("at mnar_strength = 2, coded cases are far more severe (t = %.1f)", gap_t(d2)))

# --- the sweep must vary informativeness ONLY, not coding frequency ----------
# plogis(qlogis(p) + m*U) does not average to p once m > 0, so without solving
# for the intercept the marginal coding rate drifts along the x-axis and the two
# mechanisms become confounded. These checks pin that down.
ok(isTRUE(all.equal(solve_coding_intercept(0.45, 0), qlogis(0.45))),
   "intercept solve reduces to qlogis(target) when mnar_strength = 0")
for (tgt in c(0.45, 0.85)) {
  rates <- vapply(c(0, 1, 2), function(m) {
    dd <- simulate_cohort(sim_params(n = 200000, mnar_strength = m, base_coding_prob = tgt))
    mean(dd$S_obs[dd$S == 1])
  }, numeric(1))
  ok(max(abs(rates - tgt)) < 0.01,
     sprintf("coding rate stays at %.2f across the MNAR sweep (observed %s)",
             tgt, paste(sprintf("%.3f", rates), collapse = ", ")))
}

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
