# Source provenance of the curated package

The repository was curated from the project's working Google Drive folder containing GaS/GaSe/InSe strain, vacancy, VDOS, harmonic-phonon, and finite-size-validation branches.

## Preserved directly from the working archive

- original Excel result/validation workbooks;
- Python analysis/generation scripts;
- representative LAMMPS inputs;
- Tersoff potential files used by the MD/VDOS and phonon branches;
- pristine and representative defect/strain structures;
- raw/smoothed VDOS tables and VDOS metadata;
- harmonic phonon dispersion/PDOS tables and figures;
- representative raw cumulative thermostat-energy files and diagnostic PNGs.

## Added during GitHub curation

The following are **derived convenience artifacts**, not additional simulations:

- `results/results_summary.csv`;
- `results/strain/strain_summary.csv`;
- `results/vacancy/vacancy_summary.csv`;
- `results/validation/finite_size_summary.csv`;
- `results/validation/literature_benchmark.csv`;
- standardized strain/vacancy/finite-size PNG summary plots;
- repository documentation (`README.md`, `docs/*`, citation/licensing notices, manifest/checksums).

The CSV summaries were derived from the original Excel workbooks. Finite-size fit values in the derived CSV are recomputed with an unweighted linear regression of the workbook's displayed `1/kappa` versus `1/L` points; therefore very small rounding differences from text values embedded in a workbook are possible. The original workbooks remain the provenance record.

## Literature

The Wang et al. (2019) DFT-BTE values included in `literature_benchmark.csv` are literature values, not project-generated data. See `references/REFERENCES.md`.
