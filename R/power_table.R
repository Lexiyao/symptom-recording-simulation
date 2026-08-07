source("simulate_two_group.R")
K <- 1.339   # SD(cal_int_gap) * sqrt(events per group), from the n=400,000 run
cat("coding OR | coded A / B | cal-int gap (log-odds) | events/group for 80% power | patients/group @1.5%\n")
for (cor in c(1.3, 2.0, 3.0, 5.0)) {
  s <- run_two_group_sweep(mnar_grid = 1, nsim = 30, n = 400000, coding_or_group = cor)
  g <- s$cal_int_gap
  need <- (2.802 * K / g)^2
  cat(sprintf("%9.1f | %.2f / %.2f | %+21.4f | %25.0f | %19.0f\n",
      cor, s$coded_rate_A, s$coded_rate_B, g, need, need/0.015))
}
