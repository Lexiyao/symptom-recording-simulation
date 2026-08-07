source("simulate_two_group.R")
# How many patients per 1,000 does the ONE model place on the wrong side of the
# 3% referral line, and does that differ by group? This is the Van Calster
# translation: a calibration difference expressed as patients selected.
set.seed(20260726)
p <- two_group_params(n = 400000, mnar_strength = 1)
agg <- NULL
for (i in 1:25) {
  d   <- simulate_two_group_cohort(p)
  fit <- glm(Y ~ S_obs + age_z, family = binomial, data = d)
  pr  <- predict(fit, type = "response")
  # true risk for each patient under the generative model (known only in simulation)
  true <- plogis(qlogis(p$baseline_risk) + p$log_or_symptom*d$S + 0.35*d$U + 0.25*d$age_z)
  sel  <- pr   >= p$threshold          # what the model does
  need <- true >= p$threshold          # what the truth says
  f <- function(idx) c(
    selected   = 1000*mean(sel[idx]),
    should     = 1000*mean(need[idx]),
    missed     = 1000*mean(need[idx] & !sel[idx]),   # truly over the line, scored under
    overtreat  = 1000*mean(sel[idx] & !need[idx]))
  agg <- rbind(agg, c(A = f(d$G==0), B = f(d$G==1)))
}
m <- colMeans(agg); s <- apply(agg,2,sd)/sqrt(nrow(agg))
cat("per 1,000 patients, at the 3% referral line (25 reps of 400,000):\n\n")
cat(sprintf("%-34s %-18s %-18s %s\n","","worse-coded (A)","better-coded (B)","difference"))
for (k in c("selected","should","missed","overtreat")) {
  a <- m[paste0("A.",k)]; b <- m[paste0("B.",k)]
  cat(sprintf("%-34s %8.1f          %8.1f          %+7.1f\n", k, a, b, a-b))
}
cat(sprintf("\nMonte Carlo SE on the 'missed' difference: %.2f per 1,000\n",
    sqrt(s["A.missed"]^2 + s["B.missed"]^2)))
