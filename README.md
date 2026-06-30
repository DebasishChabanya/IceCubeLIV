\#IceCube LIV Correlation

This repository contains the computational pipeline for testing Lorentz Invariance Violation (LIV) using the 12-year IceCube High-Energy Starting Event (HESE) catalog and Gamma-Ray Burst (GRB) triggers.



This work was conducted as part of an M1 Physics internship project at the Institut d'Astrophysique de Paris (IAP) (Sorbonne Université / Université Paris Cité) in collaboration with Laboratoire de physique nucléaire et des hautes énergies (LPNHE).





\# Project Overview

The analysis replicates and re-evaluates the baseline time-delay correlations claimed by Amelino-Camelia et al. (2023). It extends the methodology by introducing:

\* Strict Haversine spherical geometry for spatiotemporal cross-matching.

\* Realistic instrumental energy uncertainties ($\\pm 30\\%$).

\* Dynamic Monte Carlo redshift sampling using Wanderman \& Piran (2015) cosmological density models.





\# Repository Structure

\* `data/`: Contains the IceCube neutrino catalog (`.npy`) and GRB trigger dataset (`.txt`).

\* `src/`: Core Julia modules (`Neutrino.jl`, `GRB.jl`, `GRB\\\\\\\_RedshiftSampler.jl`) for data ingestion and cosmological kinematic constraints.

\* `analysis/`: Julia scripts (`Correlation.jl`, `Replication.jl`, `Extension.jl`, `Check\\\\\\\_RedshiftSampler.jl`) to execute the cross-matching pipeline and perform regressions for both AC2023 data and new dataset, along with verifying the redshift sampling code matching it with theory.

\* `results/`: Generated Monte Carlo plots, Inverse Transform Sampling validations, and correlation outputs.





\# Toolstack

\* Language: Julia

\* Methodology: Spatiotemporal engine, Weighted Least Squares (WLS), Inverse Transform Sampling.

