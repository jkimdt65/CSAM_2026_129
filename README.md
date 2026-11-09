# Code for "Bayesian shared frailty models for interval-censored survival data in Stan"

Supplementary code. All results in the paper and the Supplementary Material
are produced by the scripts below. Each R script states at its top what it
does, how it is run, and what it writes.

## Software

| Software | Version |
|---|---|
| R | 4.3.3 |
| CmdStan | 2.37.0 |
| cmdstanr | 0.9.0 |
| loo, posterior | 2.10.1, 1.7.0 |
| splines2 | 0.5.4 |
| survivalMPL (breast cosmesis data, `bcos2`) | 0.2-4 |
| UdderQuarterInfectionData (udder quarter data) | 1.0.0 |

The scripts call `set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")`; change this
line if CmdStan is installed elsewhere. Base R packages `splines` and
`parallel` are also used. Figures are drawn with base graphics.

Seeds are fixed in every script (2026 for the data analyses; `seed0 + rep`
for each simulation replication), so a rerun with the same software versions
and number of cores reproduces the reported numbers.

## Layout

```
stan/      Stan programs (12), shared by all analyses
bcos/      breast cosmesis data (Sections 1, 2.2 and 5)
udder/     udder quarter infection data (Sections 3 and 5, Tables 2-4, Figure 2,
           Supplementary Tables S1 and S2)
sim/       simulation study (Section 4, Table 1, Supplementary Table S4)
figures/   Figure 1
```

Every script is run from its own directory and reads Stan programs from
`../stan/`.

## Run order

`run_all.sh` runs every data analysis below in order and records `sessionInfo()`:

```
setsid nohup bash run_all.sh > run_all_log.txt 2>&1 < /dev/null &
```

The steps can also be run one by one, as follows.

### Breast cosmesis data (`bcos/`)

| Step | Script | Output | Used in |
|---|---|---|---|
| 1 | `01_data.R` (sourced by 2 and 3) | Stan data | |
| 2 | `02_fit_baselines.R` | `fit_bcos_{pwe,isp,rp}.rds`, `bcos_baselines.md` | Section 5, Supplementary Table S3 |
| 3 | `03_sens_prior.R` | `bcos_sens_prior.md` | Section 5 |

### Udder quarter infection data (`udder/`)

| Step | Script | Output | Used in |
|---|---|---|---|
| 1 | `01_data_weibull.R` (sourced) | Stan data, Weibull baseline | |
| 2 | `02_fit_frailty.R` | `fit_gamma.rds`, `fit_lognormal.rds`, `loo_frailty.rds`, `frailty_comparison.md` | Tables 2 and 3, Section 5 |
| 3 | `fig2_frailty.R` | `fig2_frailty.pdf` | Figure 2 |
| 4 | `sens_prior_theta.R` | `sens_prior_theta.md` | Supplementary Table S2 |
| 5 | `sens_prior_sigma.R` | `sens_prior_sigma.md` | Section 5 |
| 6 | `03_data_baselines.R` (sourced) | Stan data, flexible baselines | |
| 7 | `04_fit_baselines.R` | `fit_{wb_A,pwe_A,isp_A,wb_B}.rds`, `baselines.md` | Table 4 |
| 8 | `05_tables_95.R` | `tables_95.md` | Tables 2, 3 and 4 (95% intervals) |
| 9 | `sens_resolution.R` | `sens_resolution.md` | Supplementary Table S1 |

Steps 3-5 need step 2; steps 8 and 9 need steps 2 and 7.

### Simulation (`sim/`)

`sim_functions.R` holds the data generator and the construction of the Stan
data; it is sourced by the two drivers.

| Script | Output | Used in |
|---|---|---|
| `sim_table1.R` | `sim_table1_*.rds`, `sim_table1_results.md` | Table 1 |
| `sim_shapes.R` | `sim_shapes_*.rds`, `sim_shapes_results.md` | Supplementary Table S4 |

Both drivers save each configuration separately and skip configurations
already saved, so an interrupted run can be resumed. `SIM_TEST=1` runs three
replications only. Long runs are best started detached, e.g.

```
setsid nohup Rscript sim_table1.R > sim_table1_log.txt 2>&1 < /dev/null &
```

### Figure 1 (`figures/`)

`fig1_baselines.R` writes `fig1_baselines.pdf`; it needs no data.

## Running times

On a 12-core Linux server (Ubuntu 24.04), shared with other users:

| Step | Time |
|---|---|
| `udder/02_fit_frailty.R` | 21 min |
| `udder/fig2_frailty.R` | 14 s |
| `udder/sens_prior_theta.R` | 22 min |
| `udder/sens_prior_sigma.R` | 10 min |
| `udder/04_fit_baselines.R` | 56 min |
| `udder/05_tables_95.R` | 18 min |
| `udder/sens_resolution.R` | 42 min |
| `bcos/02_fit_baselines.R`, `bcos/03_sens_prior.R` | under 1 min |
| `figures/fig1_baselines.R` | under 1 s |
| `run_all.sh` in all | about 2 h 50 min |
| `sim/sim_table1.R` | about 18 h |
| `sim/sim_shapes.R` | about 13 h |

## Stan programs

| File | Model |
|---|---|
| `ph_ic_pwe_az.stan` | PH, piecewise exponential baseline, no frailty |
| `ph_ic_ispline_az.stan` | PH, I-spline baseline, no frailty |
| `ph_ic_ncs_az.stan` | PH, Royston-Parmar baseline, no frailty |
| `weibull_gamma_frailty_ic.stan` | Weibull baseline, shared gamma frailty |
| `weibull_lognormal_frailty_ic.stan` | Weibull baseline, shared log-normal frailty |
| `weibull_gamma_frailty_X.stan` | Weibull baseline, gamma frailty, general covariate matrix |
| `ph_ic_pwe_gamma_frailty.stan` | piecewise exponential baseline, gamma frailty |
| `ph_ic_ispline_gamma_frailty.stan` | I-spline baseline, gamma frailty |
| `weibull_gamma_priorvar.stan` | as `weibull_gamma_frailty_ic.stan`, prior on theta selected by `prior_id` |
| `weibull_lognormal_priorvar.stan` | as `weibull_lognormal_frailty_ic.stan`, prior on sigma_b selected by `prior_id` |
| `sim_pwe_tdc.stan` | simulation: piecewise exponential, time-varying covariate, delayed entry |
| `sim_isp_tdc.stan` | simulation: I-spline, time-varying covariate, delayed entry |
