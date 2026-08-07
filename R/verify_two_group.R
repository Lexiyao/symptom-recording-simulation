source("simulate_two_group.R")
pass <- function(name, ok, detail) cat(sprintf("  [%s] %s\n        %s\n", if (ok) "PASS" else "FAIL", name, detail))

cat("TEST 1  Null: identical coding in both groups must give a zero gap\n")
n0 <- run_two_group_sweep(mnar_grid = 1, nsim = 40, n = 400000, coding_or_group = 1.0)
z  <- n0$cal_int_gap / n0$cal_int_gap_mcse
pass("no group difference in coding -> no calibration gap",
     abs(z) < 2,
     sprintf("gap %+.5f, MCSE %.5f, z = %+.2f", n0$cal_int_gap, n0$cal_int_gap_mcse, z))

cat("\nTEST 2  Label swap: reversing which group codes better must flip the sign\n")
fwd <- run_two_group_sweep(mnar_grid = 1, nsim = 30, n = 400000, coding_or_group = 3.0)
rev <- run_two_group_sweep(mnar_grid = 1, nsim = 30, n = 400000, coding_or_group = 1/3)
pass("sign flips, magnitude preserved",
     sign(fwd$cal_int_gap) == -sign(rev$cal_int_gap) &&
       abs(abs(fwd$cal_int_gap) - abs(rev$cal_int_gap)) < 3*fwd$cal_int_gap_mcse,
     sprintf("OR 3.0 -> %+.4f ; OR 1/3 -> %+.4f", fwd$cal_int_gap, rev$cal_int_gap))

cat("\nTEST 3  Perfect coding: no misclassification -> no attenuation, no gap\n")
pc <- run_two_group_sweep(mnar_grid = 0, nsim = 30, n = 400000,
                          base_coding_prob = 0.999, coding_or_group = 1.0)
pass("oracle and pooled OR agree when nothing is missed",
     abs(pc$or_pooled - pc$or_oracle) < 0.5 && abs(pc$cal_int_gap) < 3*pc$cal_int_gap_mcse,
     sprintf("pooled %.2f vs oracle %.2f ; gap %+.4f", pc$or_pooled, pc$or_oracle, pc$cal_int_gap))

cat("\nTEST 4  Closed form: attenuation at random coding, imperfect sensitivity\n")
cf <- run_two_group_sweep(mnar_grid = 0, nsim = 30, n = 400000,
                          base_coding_prob = 0.45, coding_or_group = 1.0)
# with perfect specificity and sensitivity q, the observed OR relates to the true OR
# through the mixing of true positives into the unexposed group
q <- 0.45; prev <- 0.08; OR <- exp(2.34); base <- 0.004
p1 <- plogis(qlogis(base) + 2.34)              # risk | symptom
p0 <- base                                      # risk | no symptom
# observed 'exposed' = truly symptomatic AND coded ; observed 'unexposed' = the rest
pe <- prev*q
pu <- 1 - pe
re <- p1
ru <- (prev*(1-q)*p1 + (1-prev)*p0) / pu
or_expected <- (re/(1-re)) / (ru/(1-ru))
pass("simulated OR matches the closed form for sensitivity 0.45",
     abs(cf$or_pooled - or_expected) < 0.6,
     sprintf("simulated %.2f vs closed form %.2f (declared true %.2f)", cf$or_pooled, or_expected, OR))
