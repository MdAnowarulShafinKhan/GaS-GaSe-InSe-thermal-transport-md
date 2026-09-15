# Thermal-conductivity validation benchmarks

## Current finite-size extrapolation

The curated shorter-length dataset gives approximately:

| Material | kappa_inf (W m^-1 K^-1) | effective characteristic length (nm) | reciprocal-fit R^2 |
|---|---:|---:|---:|
| GaS | 14.05 | 34.19 | 0.991 |
| GaSe | 10.19 | 35.55 | 0.988 |
| InSe | 3.30 | 23.57 | 0.970 |

The detailed numerical points are in `results/validation/finite_size_summary.csv`.

## Independent first-principles benchmark

Wang et al., *Journal of Applied Physics* 125, 245104 (2019), DOI 10.1063/1.5094663, report 300 K DFT-BTE lattice thermal conductivities:

| Material | Wang et al. 2019 (W m^-1 K^-1) | This repository kappa_inf (W m^-1 K^-1) |
|---|---:|---:|
| GaS | 14.26 | 14.05 |
| GaSe | 13.35 | 10.19 |
| InSe | 2.69 | 3.30 |

The comparison is normalization-compatible because both use 8.712 A, 8.875 A, and 8.775 A for GaS, GaSe, and InSe, respectively.

This is an **independent physical benchmark**, not a claim that classical NEMD should numerically reproduce DFT-BTE exactly.

## EIAP reference

Karaaslan, *J. Phys. D: Appl. Phys.* 58, 015306 (2025), DOI 10.1088/1361-6463/ad7c58, is the key source for the GaS/GaSe/InSe Tersoff-type EIAP workflow and provides the closer same-potential/same-method context. Differences from that paper should be discussed separately from agreement with the DFT-BTE benchmark.
