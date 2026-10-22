# IceCubeLIV

M1 Physics internship project by Debasish Chabanya at IAP, with LPNHE.

This repository contains the GRB–neutrino matching code, a seven-event baseline calculation, and a thirteen-event extension. The updated scripts sample energy over an assumed ±30% range and sample missing redshifts. This file explains how to reproduce the saved results and what we can do for verification.

## 1. What is in each folder

| Folder or file | What it contains |
|---|---|
| `data/` | The input catalogs: `Neutrino.npy` and `GRB.txt`. |
| `src/` | The catalog readers and redshift sampler. |
| `analysis/` | The four scripts described below. |
| `results/` | Seven uploaded plots, with `_new` in their names. |
| `CSV (for values)/` | Saved numerical results from the updated runs. These can be opened in a spreadsheet or read with code. |
| `old report (full latex zip)/` | The original report and LaTeX source in a ZIP. Some methods and results are now out of date. |
| `Project.toml` | The direct Julia package dependencies: Plots and QuadGK. |
| `Manifest.toml` | The package versions from the tested Julia 1.10.12 environment. |
| `Plots_overview.png` | The four main plots together, for a quick look. |
| `Validation.json` | Notes from the earlier test run. Its script hashes refer to the delivered files before later comment and console-text edits. It is not a scientific validation. |
| `README.md` | updated instructions. |

## 2. What each code file does

| File | Purpose |
|---|---|
| `src/Neutrino.jl` | Reads the supplied IceCube `.npy` file. The reader is written for that file's binary layout; not to assume it can read a different catalog unchanged. |
| `src/GRB.jl` | Reads the GRB names, times, positions, position errors, durations and redshifts. |
| `src/GRB_RedshiftSampler.jl` | Keeps measured redshifts and draws missing ones from separate long- and short-GRB curves. |
| `analysis/Correlation.jl` | Searches the two catalogs for candidate pairs using energy, time and angular cuts. |
| `analysis/Check_RedshiftSampler.jl` | Compares sampled redshifts with the curves implemented in the sampler. This checks sampling against those curves, not whether the curves are the right physical model. |
| `analysis/Replication.jl` | Runs the seven-event fixed-proxy baseline, energy sampling with fixed proxies, and joint energy/redshift sampling. |
| `analysis/Extension.jl` | Runs joint energy/redshift sampling for the thirteen selected events. |

The two fitting scripts have their event values written directly into arrays. They do not automatically read the output of `Correlation.jl`. The extension uses energies and delays rounded to two decimal places in TeV and days.

## 3. Run the code

Use Julia 1.10.12 to match the tested environment. Download or clone the whole repository, then open a terminal in its main folder.

Install the packages once:

```sh
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Run these as separate commands, in this order:

```sh
julia --project=. analysis/Correlation.jl
julia --project=. analysis/Check_RedshiftSampler.jl
julia --project=. analysis/Replication.jl
julia --project=. analysis/Extension.jl
```

The first script should report 164 neutrinos, 9,150 GRBs and 15 candidate pairs. The earlier matching check found 13 distinct neutrinos. Neutrino 107 has three candidate GRBs; the extension retains GRB170805C, which had the smallest angular separation. This is a selection rule we took.

To inspect the candidates, start `julia --project=.` from the main folder and enter:

```julia
include("analysis/Correlation.jl")
for m in matches
    println((neutrino=m.nu.id, grb=m.grb.id,
             energy_TeV=m.E_TeV, delay_days=m.delta_t_days,
             z=m.grb.redshift, t90=m.grb.t90,
             separation_deg=m.d_psi))
end
```

Compare these rows with the arrays in `Extension.jl`. The script does not save this candidate table itself.

### Where the new outputs appear

- `Replication.jl` and `Extension.jl` write their PNGs and CSVs into `results/`.
- Their generated PNG names do **not** include `_new`. For example, `1_Replication.png` corresponds to the uploaded `1_Replication_new.png`.
- Compare the new CSVs in `results/` with the saved copies in `CSV (for values)/`. The scripts do not update that second folder.
- `Check_RedshiftSampler.jl` writes `InverseCDF_Visual.png` in the folder from which Julia was started. This corresponds to uploaded plot 7. It has no fixed random seed, so its histogram will change slightly between runs.
- `Plots_overview.png` is a prepared overview; the four scripts do not recreate it.

## 4. Plots and CSV files

| Uploaded plot in `results/` | What it shows |
|---|---|
| `1_Replication_new.png` | Seven-event fit through the origin using fixed redshift proxies and ordinary least squares. |
| `2_Energy_uncertainty_new.png` | Fixed redshift proxies, with energy sampling and the mean-ratio estimate. |
| `3_Redshift_and_energy_new.png` | Seven events with energy and missing-redshift sampling. |
| `4_Extended_13_Events_new.png` | The thirteen-event extension. |
| `5_Eta_distribution_7_events_new.png` | The distribution of eta across the seven-event joint simulations. |
| `6_Eta_distribution_13_events_new.png` | The corresponding thirteen-event distribution. |
| `7_Redshift_sampler_check_new.png` | Long- and short-GRB samples compared with the implemented curves. |

Plots 3 and 4 show one example redshift draw at the recorded energies. The horizontal bars show the assumed energy range. Their mean and median lines summarize all simulations; they are not fits to the particular points displayed. The bars do not show the full redshift uncertainty.

| CSV in `CSV (for values)/` | Contents |
|---|---|
| `replication_baseline.csv` | Baseline eta and its ordinary-least-squares standard error. |
| `replication_summary.csv` | Summaries for energy-only sampling and joint energy/redshift sampling. |
| `extension_summary.csv` | Summary for the thirteen-event joint calculation. |
| `replication_mc_draws.csv` | All 10,000 eta values from the seven-event joint calculation. |
| `extension_mc_draws.csv` | All 10,000 eta values from the thirteen-event joint calculation. |
| `replication_cutoff_sensitivity.csv` | Seven-event results for four minimum redshifts. |
| `extension_cutoff_sensitivity.csv` | Thirteen-event results for the same four minimum redshifts. |

`mean` is the average, `sd` is the sample standard deviation, and `median` is the middle value. `p025` and `p975` mark the central 95% of simulated values. They are not automatically a confidence interval for the physical LIV parameter. `draw` is the simulation number, not an event ID. The draw files contain eta values, not the individual sampled energies and redshifts.

At the default minimum redshift of 0.0001, the saved results are:

| Calculation | Eta mean | Spread | Median | Central 95% of draws |
|---|---:|---:|---:|---|
| Fixed-proxy baseline | 21.7144 | 4.5090, OLS standard error | — | — |
| Energy sampling, fixed proxies | 23.9848 | 1.9071, MC standard deviation | 23.9113 | 20.5376–27.8700 |
| Seven events, energy and redshift | 29.2371 | 445.3353, MC standard deviation | 15.3022 | 10.9358–58.0199 |
| Thirteen events, energy and redshift | 51.2536 | 966.3341, MC standard deviation | 21.0031 | 15.0678–87.7561 |

## 5. What changed in the update

The old code drew horizontal ±30% bars but kept energies fixed in the calculations. The update samples each energy independently and uniformly between `0.7E` and `1.3E` on every simulation. Energy draws are also independent of the redshift draws.

**The uniform distribution was an added choice.** Saying “±30%” alone does not specify a distribution. Here it is a bounded assumed range, not a measured event-by-event detector response and not a Gaussian one-sigma error. The 30% value was chosen for this study after discussion with other researchers associated with IceCube, not estimated from these fits.

Other additions were fixed random seeds, CSV outputs, median/percentile summaries, distribution plots, and minimum-redshift comparisons. The seeds are 2025 for energy-only sampling, 2026 for joint sampling, and 2027 for example scatter plots. The main calculations use 10,000 draws. Fixed seeds make a run repeatable; they do not show that its tail has converged.

The energy-only error shown now is the spread across simulations. The earlier code used the standard error of the seven individual ratios. Those are different quantities and should not be described as the same error calculation.

The event arrays, cosmology, catalog readers, catalog matching, redshift curves and mean-ratio estimator were retained. The default clamp of 0.0001 was already in the old code. The extra checks use 0.000001, 0.001 and 0.01, reusing the same energy and redshift draws for comparison.

## 6. Assumptions and checks

### Fitting method

The baseline uses an ordinary least-squares line through the origin. Its redshift array contains six values of 1.38 and one of 0.60; it does not give all seven events the same redshift.

The later calculations use:

```text
K = E × D(z)/D(1)
eta = mean(delta_t / K) / (D(1)/M_P)
```

This equals a fit through the origin with relative weights proportional to `1/K²`. It gives each event's ratio the same weight. It is not a full fit using a detector likelihood for uncertain energy.

Moving from plot 1 to plot 2 changes both the estimator and the energy treatment. With fixed energies, the mean-ratio estimator already gives eta about 23.2437, compared with the OLS value 21.7144.

### Redshift model and the clamp

The sampler normalizes its long- and short-GRB curves over `0 ≤ z ≤ 6` and uses them directly as probability distributions. These are working assumptions. The long-GRB rate shape and the short-GRB shape is from Wanderman & Piran (2010) (2015). This is a distribution we chose, there are others as well.

There is a numerical observation: the cumulative sum starts with a positive value at `z = 0`. The inverse sampler therefore returns exactly zero for a small, finite fraction of draws. The analysis replaces all draws below `Z_MIN` with `Z_MIN`; it does not reject and redraw them. This creates a pile-up at the cutoff. Since eta involves division by `D(z)`, very small redshifts produce very large values.

For the thirteen-event run, changing the cutoff from 0.0001 to 0.000001 changes the mean/SD from about 51.3/966.3 to 1510.5/90723.7. The median stays about 21.0 in these saved draws. This is strong sensitivity to the lower boundary, including the sampler artifact. We can try a zero-starting cumulative integral and a grid/convergence check before these extreme tails are used as a result.

### Event selection and interpretation

- Matching uses recorded energies from 60 to 500 TeV, neutrinos arriving 0–3 days after a GRB, and separation within three combined position-error widths.
- The matching code assumes 15 degrees for showers/cascades and 1 degree for other reconstruction types. Missing GRB position errors are replaced with the catalog median.
- Energy sampling happens after selection. The code does not repeat the matching cuts when a sampled energy crosses 60 or 500 TeV.
- Measured redshifts, event times and cosmological parameters are held fixed. The model does not fit an intrinsic source emission delay or an intercept.
- The seven-event duration values mainly encode long/short classifications; they are not all measured durations.
- The thirteen-event arrays should be tied to a saved table of neutrino IDs and GRB names, with the duplicate-match rule stated clearly.
- A large Monte Carlo spread shows instability under the chosen assumptions. Our result does not by itself prove that LIV exists, that it does not exist.

We can check the redshift model and sampler, justify the energy model and estimator, verify the candidate table, and repeat the calculations with more draws.

References for checking the rate curves:

- Wanderman & Piran (2010), *The luminosity function and the rate of Swift's gamma-ray bursts*: https://arxiv.org/abs/0912.0709 (especially equations 2 and 18).
- Wanderman & Piran (2015), *The rate, luminosity function and time delay of non-Collapsar short GRBs*: https://arxiv.org/abs/1405.5878 (especially equation 9).

## 7. What was checked for this README

The GitHub files were checked on 22 October 2026. All seven uploaded PNGs and all seven CSVs match the delivered files. The two updated scripts have the same calculations as the delivered versions.
