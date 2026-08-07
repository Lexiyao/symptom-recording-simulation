# ---------------------------------------------------------------------------
# Two-group branch: does a group difference in symptom CODING produce a group
# difference in the CALIBRATION of one model applied to everybody?
#
# Extends R/simulate_recording_bias.R. The generative model is unchanged except
# for one thing: base_coding_prob differs between two groups. S, U, Y and the
# true symptom-cancer effect are IDENTICAL in both groups by construction, so
# any calibration difference that appears is caused by recording alone.
#
# Everything is declared. Nothing here is estimated from patient data.
# ---------------------------------------------------------------------------

source("simulate_recording_bias.R")   # sim_params, solve_coding_intercept, simulate_cohort

#' Parameters for the two-group branch
#'
#' @param coding_or_group odds ratio for being coded, group B vs group A.
#'   DECLARED, not estimated. The default is the scale d'Elia et al. (2025)
#'   report for "any symptom coded" in CPRD Aurum (OR 1.17 Black, 1.16 South
#'   Asian vs White) rounded up to 1.3 so the mechanism is visible; individual
#'   symptoms in that paper run to 2.02. Sweep it with sensitivity_coding_or().
#' @param prop_b share of the cohort in the better-coded group. 0.5 maximises
#'   precision; it is not a claim about any real population.
two_group_params <- function(n = 100000,
                             base_coding_prob = 0.45,   # vague-symptom regime
                             coding_or_group  = 1.30,
                             prop_b           = 0.5,
                             mnar_strength    = 0,
                             ...) {
  p <- sim_params(n = n, base_coding_prob = base_coding_prob,
                  mnar_strength = mnar_strength, ...)
  p$coding_or_group <- coding_or_group
  p$prop_b <- prop_b
  p$coding_prob_b <- plogis(qlogis(base_coding_prob) + log(coding_or_group))
  p
}

#' One cohort with two groups that differ ONLY in coding probability
simulate_two_group_cohort <- function(p) {
  G     <- rbinom(p$n, 1, p$prop_b)          # 0 = group A, 1 = group B
  age_z <- rnorm(p$n)
  S     <- rbinom(p$n, 1, p$prev_symptom)    # same true prevalence in both
  U     <- rnorm(p$n)                        # same severity distribution

  # the ONLY thing that differs by group
  a_A <- solve_coding_intercept(p$base_coding_prob, p$mnar_strength)
  a_B <- solve_coding_intercept(p$coding_prob_b,   p$mnar_strength)
  coding_logit <- ifelse(G == 1, a_B, a_A) + p$mnar_strength * U
  R     <- rbinom(p$n, 1, plogis(coding_logit))
  S_obs <- as.integer(S & R)

  # outcome: identical structural equation in both groups. G does not enter.
  eta <- qlogis(p$baseline_risk) + p$log_or_symptom * S + 0.35 * U + 0.25 * age_z
  Y   <- rbinom(p$n, 1, plogis(eta))

  data.frame(G, age_z, S, U, R, S_obs, Y)
}

#' Calibration slope of ONE pooled model, evaluated separately in each group
#'
#' A deployed risk score is a single model applied to everybody. So the model is
#' fitted once on the pooled coded data, and its linear predictor is then taken
#' to each group and the recalibration slope estimated there:
#'     logit P(Y) = a_g + b_g * LP
#' b_g = 1 is perfect calibration in group g. The reported quantity is b_A - b_B.
summarise_two_group <- function(d, p) {
  fit <- glm(Y ~ S_obs + age_z, family = binomial, data = d)
  lp  <- predict(fit, type = "link")

  slope <- function(idx) {
    if (sum(d$Y[idx]) < 5) return(NA_real_)
    unname(coef(glm(d$Y[idx] ~ lp[idx], family = binomial))[2])
  }
  bA <- slope(d$G == 0)
  bB <- slope(d$G == 1)

  # calibration intercept: logit P(Y) = a_g + LP, slope fixed at 1.
  # a_g < 0 means the model over-predicts in group g; a_g > 0 under-predicts.
  cal_int <- function(idx) {
    if (sum(d$Y[idx]) < 5) return(NA_real_)
    unname(coef(glm(d$Y[idx] ~ offset(lp[idx]), family = binomial))[1])
  }
  aA <- cal_int(d$G == 0); aB <- cal_int(d$G == 1)

  pr <- predict(fit, type = "response")
  A <- d$G == 0; B <- d$G == 1
  # calibration-in-the-large, computed the way it can be computed in CPRD:
  # observed events over model-predicted events, whole group, no need to know S
  oe <- function(idx) mean(d$Y[idx]) / mean(pr[idx])
  # and the same restricted to patients the model puts near the referral line
  near <- pr > 0.02 & pr < 0.05
  oe_near <- function(idx) if (sum(idx & near) < 50) NA_real_ else
                            mean(d$Y[idx & near]) / mean(pr[idx & near])

  data.frame(
    mnar_strength = p$mnar_strength,
    or_pooled  = unname(exp(coef(fit)["S_obs"])),
    or_oracle  = unname(exp(coef(glm(Y ~ S + age_z, family = binomial, data = d))["S"])),
    slope_A    = bA,
    slope_B    = bB,
    slope_gap  = bA - bB,
    oe_A       = oe(A),
    oe_B       = oe(B),
    oe_gap     = oe(A) - oe(B),
    cal_int_A  = aA,
    cal_int_B  = aB,
    cal_int_gap = aA - aB,
    oe_near_A  = oe_near(A),
    oe_near_B  = oe_near(B),
    coded_rate_A = mean(d$R[A]),
    coded_rate_B = mean(d$R[B]),
    events       = sum(d$Y)
  )
}

run_two_group_sweep <- function(mnar_grid = seq(0, 2, by = 0.25),
                                nsim = 40, seed = 20260726, ...) {
  set.seed(seed)
  do.call(rbind, lapply(mnar_grid, function(m) {
    p <- two_group_params(mnar_strength = m, ...)
    reps <- do.call(rbind, lapply(seq_len(nsim),
             function(i) summarise_two_group(simulate_two_group_cohort(p), p)))
    num  <- setdiff(names(reps), "mnar_strength")
    est  <- vapply(reps[num], function(x) mean(x, na.rm = TRUE), numeric(1))
    mcse <- vapply(reps[num], function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))), numeric(1))
    data.frame(mnar_strength = m, nsim = nsim,
               t(est), t(setNames(mcse, paste0(num, "_mcse"))))
  }))
}

# --- invariants: if any of these fail the result is not interpretable --------
check_invariants <- function(res, p) {
  ok <- c(
    "coding rate in B exceeds A at every grid point" =
      all(res$coded_rate_B > res$coded_rate_A),
    "coding rates stay flat along the sweep (marginal rate held fixed)" =
      diff(range(res$coded_rate_A)) < 0.02 && diff(range(res$coded_rate_B)) < 0.02,
    "oracle OR recovers the declared truth" =
      all(abs(res$or_oracle - exp(p$log_or_symptom)) < 1.5),
    "pooled OR is attenuated against the oracle" =
      all(res$or_pooled < res$or_oracle),
    "enough events for a slope" = all(res$events > 300)
  )
  ok
}
