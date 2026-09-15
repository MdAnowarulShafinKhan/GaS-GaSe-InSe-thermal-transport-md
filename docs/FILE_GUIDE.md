# File guide

## scripts

- `scripts/nemd/thermal_conductivity_postprocess.py` - representative Fourier-law post-processing from a LAMMPS temperature profile plus cumulative thermostat energy. Geometry values are hard-coded for the source case; edit them before reuse on another geometry.
- `scripts/nemd/average_profile.py` - averages `profile.langevin` bin temperatures.
- `scripts/vacancy/make_vacancies.py` - batch generator for metal vacancies, chalcogen vacancies, and bonded metal-chalcogen divacancies.
- `scripts/vdos/make_vdos_structures.py` - generates fully periodic VDOS structures, including symmetry-balanced defect structures.
- `scripts/vdos/pdos_publication.py` - MD-derived, mass-weighted VDOS/PDOS with COM removal, Welch/Hann processing, species projections, and common normalization.
- `scripts/phonon/*_phonon_dispersion_pdos.py` - harmonic phonon dispersion/projected DOS workflows using LAMMPS + phonoLAMMPS + phonopy.

## results

- `results/results_summary.csv` - consolidated machine-readable strain, vacancy, and finite-size values.
- `results/strain/` - original Excel summary workbooks, normalized CSV summary, and GitHub-ready figures.
- `results/vacancy/` - original Excel summary workbooks, normalized CSV summary, and figures.
- `results/validation/` - original validation workbooks, re-derived finite-size CSV, fit figures, and Wang-2019 comparison table.
- `results/vdos/pristine/` - pristine MD VDOS for all three materials: total/components figures, raw and smoothed spectra, and metadata.
- `results/vdos/GaS_defects/` - detailed 5% GaS S-vacancy, Ga-vacancy, and bonded-divacancy examples.
- `results/vdos/overview/` - source overview figures retained from the Drive archive.
- `results/phonon/` - harmonic dispersion and projected DOS tables/figures for all three materials.

## structures

- `structures/pristine/` - compact pristine EIAP z-fixed cells used by the VDOS branch.
- `structures/strain_representative/` - representative equilibrated strained transport structure.
- `structures/vacancy_representative/` - representative defect structures from both transport and VDOS branches.
- `structures/phonon_cells/` - initial/relaxed four-atom cells used in the harmonic phonon branch.

## validation

`validation/representative_nemd/` contains compact raw examples. Large `profile.langevin` trajectories are not included; the raw cumulative thermostat-energy file, LAMMPS input, and final diagnostic figures are retained.
