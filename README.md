# IceCubeLIV

GRB–neutrino Lorentz invariance violation (LIV) study at IAP, in collaboration with LPNHE.

The repository contains the GRB–neutrino matching code, the seven-event AC2023 replication and the thirteen-event extension. It also includes the redshift sampling code, Monte Carlo calculations, plots and CSV files with the numerical results.

## 1. Repository structure

| File / Folder | Contents |
|---|---|
| `data/` | IceCube neutrino and GRB catalogues |
| `src/` |  Code reading the data and redshift sampling code |
| `analysis/` | Julia scripts for matching, sampling and fitting |
| `results/` | Plots from the analysis |
| `CSV (for values)/` | Saved numerical results and Monte Carlo draws |
| `Project.toml` | Julia package dependencies |
| `Manifest.toml` | Package versions used in the Julia environment |
| `Validation.json` | Supplementary run-check information |
| `old report (full latex zip)/` | Internship report and LaTeX files |

### Julia files

| File | What it does |
|---|---|
| `src/Neutrino.jl` | Reads `Neutrino.npy` |
| `src/GRB.jl` | Reads `GRB.txt`, including positions, durations and redshifts |
| `src/GRB_RedshiftSampler.jl` | Samples missing redshifts using the long- and short-GRB distributions |
| `analysis/Correlation.jl` | Finds candidate GRB–neutrino associations |
| `analysis/Check_RedshiftSampler.jl` | Compares sampled redshifts against the distributions used |
| `analysis/Replication.jl` | Seven-event calculations, including the AC2023 baseline and Monte Carlo sampling |
| `analysis/Extension.jl` | Thirteen-event calculation with energy and redshift sampling |

The event arrays used for fit are written in `Replication.jl` and `Extension.jl`. The fit scripts do not directly take their inputs from `Correlation.jl`.

## 2. Running the code

The project uses Julia. The package versions are recorded in `Manifest.toml`, with Julia 1.10.12 as the reference environment.

Clone or download the repository, then open a terminal in the project directory.

Install the required packages:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Run the analysis scripts:

```bash
julia --project=. analysis/Correlation.jl
julia --project=. analysis/Check_RedshiftSampler.jl
julia --project=. analysis/Replication.jl
julia --project=. analysis/Extension.jl
```

Each script can also be run separately.

`Replication.jl` and `Extension.jl` save their generated plots and CSVs in `results/`. The numerical CSVs included in `CSV (for values)/` are saved copies and are not automatically overwritten by these scripts.

`Check_RedshiftSampler.jl` saves `InverseCDF_Visual.png` in the working directory. The corresponding plot in the repository is named `7_Redshift_sampler_check.png`.

## 3. GRB–neutrino matching

The matching is done in `Correlation.jl` using the IceCube neutrino catalogue and the GRB catalogue.

The selection criteria are:

- Neutrino energy between 60 and 500 TeV.
- Neutrino arrival between 0 and 3 days after the GRB.
- Angular separation within three times the combined positional uncertainty.

The angular separation is calculated using spherical geometry.

For neutrinos, the code uses an angular uncertainty of 15° for shower/cascade events and 1° for other reconstruction types. Missing GRB positional uncertainties are replaced using the median of the valid catalogue values.

The catalogue contains 164 neutrinos and 9,150 GRBs.

The matching gives 15 candidate pairs involving 13 distinct neutrinos. Neutrino 107 has three possible GRB associations. For the thirteen-event analysis, GRB170805C is retained because it has the smallest angular separation among those candidates.

These are candidate associations selected using the stated criteria, not confirmed physical GRB–neutrino associations.

To print the candidate information, start Julia from the main repository folder:

```julia
include("analysis/Correlation.jl")

for m in matches
    println((
        neutrino = m.nu.id,
        grb = m.grb.id,
        energy_TeV = m.E_TeV,
        delay_days = m.delta_t_days,
        z = m.grb.redshift,
        t90 = m.grb.t90,
        separation_deg = m.d_psi
    ))
end
```

The event values in the fitting scripts can be compared with this output.

## 4. Redshift sampling

The redshift sampling is handled by `src/GRB_RedshiftSampler.jl`.

For GRBs with known redshifts, the measured values are retained. Missing redshifts are sampled using separate distributions for long and short GRBs.

The classification uses the GRB duration, T90:

- Long GRBs: T90 > 2 seconds.
- Short GRBs: T90 ≤ 2 seconds.
- Missing or invalid T90 values are treated as long GRBs.

The distributions are based on Wanderman & Piran (2010, 2015).

### Long GRBs

The function used in the code is proportional to:

\[
p(z) \propto
\begin{cases}
(1+z)^{2.1}, & 0 \leq z < 3 \\
128(1+z)^{-1.4}, & 3 \leq z \leq 6
\end{cases}
\]

### Short GRBs

\[
p(z) \propto
\begin{cases}
\exp((z-0.9)/0.39), & 0 \leq z \leq 0.9 \\
\exp(-(z-0.9)/0.26), & 0.9 < z \leq 6
\end{cases}
\]

Both distributions are normalized over the range \(0 \leq z \leq 6\).

The sampler builds the cumulative distribution function using trapezoidal integration on a grid of 10,000 redshift values. It then draws a uniform random number and uses inverse-CDF sampling to obtain a redshift.

The module provides two functions:

- `get_final_physical_redshift()` — retains a measured redshift or randomly samples a missing value.
- `get_final_ac2023_redshift()` — retains a measured redshift or uses the fixed proxies of 1.38 for long GRBs and 0.60 for short GRBs. (As the AC2023 did)

The Monte Carlo calculations in `Replication.jl` and `Extension.jl` use `get_final_physical_redshift()`.

For the thirteen-event analysis, 2 redshifts are measured (1.38 and 0.645), while the remaining 11 are sampled.

A minimum redshift of \(10^{-4}\) is used in the main analysis calculation.

## 5. Energy sampling and fitting

The energy sampling assumes an uncertainty range of ±30%. (As per said by Dr. Pablo) (Can be taken as 30 and 11 depending on cascade or track, but other explanation needed for that)

For each simulation, the energy of each event is sampled independently and uniformly between:

\[
0.7E \leq E_{\mathrm{sampled}} \leq 1.3E
\]

This is an assumed bounded uncertainty range, not a Gaussian measurement error.

The sampled energies and redshifts are used to calculate the kinematic factor:

\[
K = E \frac{D(z)}{D(1)}
\]

where

\[
D(z)=\frac{1}{H_0}\int_0^z
\frac{1+z'}{\sqrt{\Omega_m(1+z')^3+\Omega_\Lambda}}
\,dz'
\]

The cosmological parameters used are:

- \(H_0 = 67.4\) km/s/Mpc
- \(\Omega_m = 0.315\)
- \(\Omega_\Lambda = 0.685\)
- \(M_P = 1.2209\times10^{16}\) TeV

The AC2023 baseline uses an ordinary least-squares fit through the origin.

For the Monte Carlo calculations, the estimator used is:

\[
\eta =
\frac{\operatorname{mean}(\Delta t_i/K_i)}
{D(1)/M_P}
\]

This is equivalent to a through-origin weighted fit with relative weights proportional to \(1/K_i^2\).

Each event contributes equally to the average of \(\Delta t_i/K_i\).

The analysis performs 10,000 Monte Carlo draws. In each draw, the energies are sampled and the missing GRB redshifts are sampled independently.

The random seeds used are:

| Calculation | Seed |
|---|---|
| Energy-only sampling | 2025 |
| Energy and redshift sampling | 2026 |
| Example redshift scatter plots | 2027 |

The seeds make the calculations reproducible with the same code and environment.

## 6. Results

The saved values from the current CSV files are listed below.

| Calculation | Mean η | Standard deviation | Median η | Central 95% of draws |
|---|---:|---:|---:|---|
| AC2023 baseline | 21.7144 | 4.5090* | — | — |
| Energy sampling, fixed redshifts | 23.9848 | 1.9071 | 23.9113 | 20.5376–27.8700 |
| Seven events, energy and redshift sampling | 23.7583 | 128.8837 | 15.2995 | 10.9351–57.7762 |
| Thirteen events, energy and redshift sampling | 36.3571 | 317.0267 | 20.9983 | 15.0655–86.6631 |

*The 4.5090 value for the AC2023 baseline is the OLS standard error, not a Monte Carlo standard deviation.*

For the Monte Carlo results, the standard deviation describes the spread of the simulated eta values.

The 2.5th and 97.5th percentiles give the central 95% of the Monte Carlo draws. Not to be interpreted as confidence intervals for the physical LIV parameter.

The thirteen-event analysis uses 13 selected candidates, with their energies and time delays entered directly in `Extension.jl`. The time delays are recorded in days and converted to seconds for the calculation.

## 7. Plots

The `results/` folder contains seven plots.

| Plot | Contents |
|---|---|
| `1_Replication.png` | AC2023 seven-event baseline with an OLS fit |
| `2_Energy_uncertainty.png` | Energy sampling with fixed redshift proxies |
| `3_Redshift_and_energy.png` | Seven-event calculation with sampled energies and redshifts |
| `4_Extended_13_Events.png` | Thirteen-event analysis |
| `5_Eta_distribution_7_events.png` | Distribution of eta for the seven-event Monte Carlo |
| `6_Eta_distribution_13_events.png` | Distribution of eta for the thirteen-event Monte Carlo |
| `7_Redshift_sampler_check.png` | Sampled long- and short-GRB distributions compared with the implemented curves |

Plots 3 and 4 show one example draw of the redshifts, with the measured energies and assumed ±30% energy ranges.

The mean and median lines in those plots come from the full Monte Carlo calculations, not from fitting only the points shown in the particular example draw.

Plots 5 and 6 show the distribution of eta values across the 10,000 simulations, using logarithmic axes.

Plot 7 uses 100,000 samples for each GRB class to compare the generated redshift distributions with the distributions implemented in the sampler.

`Plots_overview.png` gives a combined view of the four main analysis plots.

## 8. CSV files

The `CSV (for values)/` folder contains the numerical results.

| File | Contents |
|---|---|
| `replication_baseline.csv` | AC2023 baseline eta and OLS standard error |
| `replication_summary.csv` | Energy-only and joint-sampling summaries for seven events |
| `extension_summary.csv` | Thirteen-event Monte Carlo summary |
| `replication_mc_draws.csv` | 10,000 eta values from the seven-event analysis |
| `extension_mc_draws.csv` | 10,000 eta values from the thirteen-event analysis |
| `replication_cutoff_sensitivity.csv` | Seven-event results for different minimum redshifts |
| `extension_cutoff_sensitivity.csv` | Thirteen-event results for different minimum redshifts |

In the summary files:

- `mean` is the average eta.
- `sd` is the standard deviation.
- `median` is the median eta.
- `p025` and `p975` are the 2.5th and 97.5th percentiles.

The individual draw files contain the simulation number and the calculated eta. They do not store the individual sampled redshifts or energies.

The cutoff files contain the results for minimum redshifts of \(10^{-6}\), \(10^{-4}\), \(10^{-3}\) and \(10^{-2}\). The same sampled values are used when comparing these cutoffs.

## 9. Notes on the analysis

A few assumptions are useful to keep in mind when working with the code:

- The redshift distributions are models for the GRB populations. They are not measurements of the redshifts of individual GRBs.
- The energy uncertainty is assumed to be uniformly distributed within ±30%.
- The matching criteria are applied using the catalogue values. Events are not matched again after energy sampling.
- The fitting scripts use selected event arrays rather than automatically reading the full matching output.
- The calculations keep the measured redshifts, arrival times and cosmological parameters fixed.
- The estimator assumes the fitted relation passes through the origin and does not separately model intrinsic source emission delays.
- The results depend on the selected associations, redshift models, energy assumptions and fitting method.

The purpose of the analysis is to examine the values of eta obtained under these assumptions and how they behave when energy and unknown redshift values are included in the sampling.

## 10. References

**Wanderman & Piran (2010)**  
*The luminosity function and the rate of Swift's gamma-ray bursts*  
https://arxiv.org/abs/0912.0709

**Wanderman & Piran (2015)**  
*The rate, luminosity function and time delay of non-Collapsar short GRBs*  
https://arxiv.org/abs/1405.5878