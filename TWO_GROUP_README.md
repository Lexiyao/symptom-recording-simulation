# Two-group branch — what it does, what it found, and what you must do before the interview

## Before anything else

The deck now says, on slide 4 and in Appendix W, that this branch has been run and that the code
is public. **Neither is true until you push it.** Right now `R/simulate_two_group.R` does not exist
in `github.com/Lexiyao/symptom-recording-simulation`, and slide 4 still carries
`Own work, code public`. Push the three files below, or change those lines, before the deck leaves
your machine.

You should also run it yourself once, end to end, before you present it. If a panel member asks
"what does `solve_coding_intercept` do at `mnar_strength = 0`?" the answer has to be yours.

## Files

| file | what it is |
|---|---|
| `simulate_two_group.R` | the branch. Sources your existing `simulate_recording_bias.R` and changes exactly one thing |
| `two_group_sweep.csv` | the sweep output behind the slide-4 figure: 9 grid points, 40 replicates of 200,000 |
| `two_group_fig.png` | the slide-4 figure |
| `power_table.R` | the events-per-group calculation behind the numbers on slide 5 |

## What the branch changes

One line. `base_coding_prob` differs between two groups:

```
group A: coded 45% of the time when the symptom is present
group B: odds(0.45) x 1.30  ->  coded 51.5% of the time
```

`S`, `U`, `Y`, the true prevalence, the true log odds ratio and the outcome equation are **identical**
in both groups. `G` does not enter the outcome model. So any calibration difference that appears is
produced by recording and nothing else. That is the whole point, and it is the first thing to say
if anyone asks whether the result is baked in.

The 1.30 is declared, not estimated. d'Elia reports 1.16–1.17 for "any symptom coded" and up to 2.02
for individual symptoms, so 1.30 sits inside the range they measured. `power_table.R` sweeps it.

## What it found

**1. The association recovers; the calibration does not, and the groups never converge.**
Pooled odds ratio rises 7.45 → 9.76 against a declared truth of 10.4. Observed/predicted sits at
1.015 in the worse-coded group and 0.986 in the better-coded, at every one of the nine grid points.

**2. It falsified the estimand.** You expected the gap in the calibration *slope*. It is not there:
across nine grid points the slope gap ran −0.018 to +0.003 against a Monte Carlo SE of 0.008. Even
at a coding difference of 0.45 vs 0.80 the slope gap is 0.03. The gap lands on the calibration
*level*, not the gradient.

**3. It changed which contrast is feasible.** Events per group for 80% power:

| coding OR between groups | coded A / B | gap (log-odds) | events/group | patients/group @1.5% |
|---|---|---|---|---|
| 1.30 (d'Elia between-ethnicity) | 0.45 / 0.52 | 0.030 | 15,960 | 1,064,000 |
| 2.00 | 0.45 / 0.62 | 0.078 | 2,292 | 152,795 |
| 3.00 | 0.45 / 0.71 | 0.119 | 996 | 66,429 |
| 5.00 (practice extremes) | 0.45 / 0.80 | 0.160 | 551 | 36,706 |

The between-practice contrast is powered. The between-ethnic-group contrast, at the effect size
d'Elia actually measured, is not — it would need about a million patients per group. Aim 3's primary
contrast should be practice coding propensity, with the ethnic-group contrast prespecified and
reported as underpowered. Your Appendix R, T and U already make practice coding propensity the
exposure, so the simulation supports the design you had and kills the version you were about to put
on the slide.

## Limitations you must state, not wait to be asked

**Two predictors, one of them age.** `age_z` carries most of the variance of the linear predictor,
so the calibration slope is dominated by an age gradient that coding cannot touch. In a model where
symptoms carry more of the linear predictor the slope may well move. The claim that this generalises
to QCancer is an argument, not something this simulation demonstrated. Appendix W says so.

**The SD of the gap does not fall as 1/sqrt(n).** `SD x sqrt(events)` came out 1.777 at 752
events/group and 1.339 at 3,001. The power table uses 1.339, the value from the larger run, which is
the less conservative of the two. If you want the conservative figure, multiply the events/group
column by (1.777/1.339)^2 = 1.76.

**`prop_b = 0.5`.** Equal group sizes maximise precision. It is not a claim about any real
population, and a realistic minority share makes the ethnic-group contrast worse, not better.

## Reproducing

```r
source("simulate_two_group.R")
res <- run_two_group_sweep(mnar_grid = seq(0, 2, by = 0.25), nsim = 40, n = 200000)
check_invariants(res, two_group_params())   # all five must be TRUE
```

Runtime about 8 minutes. The five invariants check that the coding rates differ in the intended
direction, that the marginal coding rate stays fixed along the sweep, that the oracle recovers the
declared truth, that the pooled estimate is attenuated against it, and that there are enough events
for a slope. If any is FALSE the figure means nothing.
