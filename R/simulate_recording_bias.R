# ---------------------------------------------------------------------------
# Informative (MNAR) symptom under-recording and cancer risk prediction
#
# Question: if a symptom is present but never coded, and whether it gets coded
# depends on severity the analyst never observes, what happens to a QCancer-style
# risk model fitted on the coded data?
#
# The generative model is deliberately small so every assumption is visible.
# All parameters are declared in sim_params() -- none are hidden in the body.
# ---------------------------------------------------------------------------

#' Default simulation parameters
#'
#' @param n patients per replicate
#' @param prev_symptom P(symptom truly present)
#' @param baseline_risk 2-year cancer risk with no symptom and mean covariates
#' @param log_or_symptom log odds ratio for a truly-present symptom.
#'   NOTE: this is a DECLARED parameter, not an estimate from any published
#'   study. It is set to the order of magnitude reported for strong alarm
#'   symptoms (e.g. rectal bleeding in colorectal cancer) so that attenuation
#'   is visible on the plot. exp(2.34) ~ 10.4. Vary it via this argument;
#'   nothing downstream depends on the specific value.
#' @param base_coding_prob target P(coded | symptom present), held FIXED across
#'   the whole mnar_strength sweep by solving for the intercept (see
#'   solve_coding_intercept). Without that solve, raising mnar_strength would
#'   also shift the marginal coding rate, and the two mechanisms this simulation
#'   is meant to separate -- how OFTEN symptoms are coded, and how INFORMATIVE
#'   the coding is -- would be confounded along the x-axis.
#'   Two regimes: an alarm symptom that is reliably coded, and a vague symptom
#'   (fatigue, abdominal discomfort) that often is not.
#' @param mnar_strength dependence of coding on unobserved severity, on the
#'   log-odds scale. 0 = coding at random given covariates (MAR).
#' @param threshold NICE NG12 referral threshold used as the reporting anchor
sim_params <- function(n = 20000,
                       prev_symptom = 0.08,
                       baseline_risk = 0.004,
                       log_or_symptom = 2.34,
                       base_coding_prob = 0.85,
                       mnar_strength = 0,
                       threshold = 0.03) {
  stopifnot(n > 0, prev_symptom > 0, prev_symptom < 1,
            base_coding_prob > 0, base_coding_prob <= 1,
            mnar_strength >= 0, threshold > 0, threshold < 1)
  as.list(environment())
}

#' Generate one cohort
#'
#' Structure (see DAG in README):
#'   U  unobserved severity / clinical suspicion
#'   S  symptom truly present
#'   R  symptom recorded, depends on U when mnar_strength > 0
#'   S* observed code = S & R
#'   Y  cancer within 2 years, depends on S and U -- never on R
#'
#' Y does not depend on R. Any association between the observed code and Y
#' beyond the causal path therefore comes from the recording process, which is
#' the whole point of the exercise.
#' Intercept that holds the marginal coding rate at `target` for a given
#' mnar_strength, assuming U ~ N(0, 1). E[plogis(a + m*U)] is not plogis(a)
#' once m > 0, so the intercept has to be solved for rather than set to
#' qlogis(target). Returns qlogis(target) exactly when m = 0.
solve_coding_intercept <- function(target, mnar_strength, n_quad = 4001) {
  if (mnar_strength == 0) return(qlogis(target))
  u <- seq(-6, 6, length.out = n_quad)
  w <- dnorm(u); w <- w / sum(w)
  f <- function(a) sum(w * plogis(a + mnar_strength * u)) - target
  stats::uniroot(f, interval = c(-30, 30), tol = 1e-10)$root
}

simulate_cohort <- function(p) {
  age_z <- rnorm(p$n)
  S     <- rbinom(p$n, 1, p$prev_symptom)
  U     <- rnorm(p$n)                        # severity, observed by nobody

  # coding depends on U among those who have the symptom, with the intercept
  # chosen so that the marginal coding rate stays at base_coding_prob
  coding_logit <- solve_coding_intercept(p$base_coding_prob, p$mnar_strength) +
                  p$mnar_strength * U
  R     <- rbinom(p$n, 1, plogis(coding_logit))
  S_obs <- as.integer(S & R)

  # outcome: symptom and severity raise risk; the code itself does not
  eta <- qlogis(p$baseline_risk) + p$log_or_symptom * S + 0.35 * U + 0.25 * age_z
  Y   <- rbinom(p$n, 1, plogis(eta))

  data.frame(age_z, S, U, R, S_obs, Y)
}

#' Fit the standard analysis and the unattainable reference
#'
#' standard: uncoded treated as absent -- what a real model can be fitted on
#' oracle:   true symptom status -- unavailable in EHR data, shown for contrast
fit_models <- function(d) {
  list(
    standard = glm(Y ~ S_obs + age_z, family = binomial, data = d),
    oracle   = glm(Y ~ S      + age_z, family = binomial, data = d)
  )
}

#' Summarise one replicate
#'
#' Reported quantities:
#'   or_standard  estimated symptom-cancer OR under 'uncoded = absent'
#'   pct_below    % of TRULY symptomatic patients whose predicted risk from the
#'                standard model falls below the referral threshold
#'   cal_oe       observed / expected risk among truly symptomatic patients
#'                (< 1 means predicted risk exceeds observed; > 1 means the
#'                model under-predicts for this group)
summarise_replicate <- function(d, p) {
  m  <- fit_models(d)
  pr <- predict(m$standard, type = "response")
  sym <- d$S == 1
  data.frame(
    mnar_strength = p$mnar_strength,
    base_coding_prob = p$base_coding_prob,
    or_standard = unname(exp(coef(m$standard)["S_obs"])),
    or_oracle   = unname(exp(coef(m$oracle)["S"])),
    pct_below   = 100 * mean(pr[sym] < p$threshold),
    cal_oe      = mean(d$Y[sym]) / mean(pr[sym])
  )
}

#' Sweep MNAR strength for one coding regime
#'
#' @param nsim replicates per grid point. Monte Carlo SE of every reported
#'   quantity is returned so readers can judge whether nsim is adequate.
run_sweep <- function(mnar_grid = seq(0, 2, by = 0.25),
                      base_coding_prob = 0.85,
                      nsim = 200,
                      seed = 20260726,
                      ...) {
  set.seed(seed)
  out <- do.call(rbind, lapply(mnar_grid, function(m) {
    p <- sim_params(mnar_strength = m, base_coding_prob = base_coding_prob, ...)
    reps <- do.call(rbind, lapply(seq_len(nsim),
                                 function(i) summarise_replicate(simulate_cohort(p), p)))
    num <- c("or_standard", "or_oracle", "pct_below", "cal_oe")
    est <- vapply(reps[num], mean, numeric(1))
    mcse <- vapply(reps[num], function(x) sd(x) / sqrt(length(x)), numeric(1))
    data.frame(mnar_strength = m, base_coding_prob = base_coding_prob,
               nsim = nsim, t(est), t(setNames(mcse, paste0(num, "_mcse"))))
  }))
  rownames(out) <- NULL
  out
}
