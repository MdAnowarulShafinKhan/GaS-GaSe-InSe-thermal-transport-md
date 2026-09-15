# Thermal transport engineering in monolayer GaS, GaSe, and InSe

This repository is a curated research package for molecular-dynamics and lattice-vibration calculations on monolayer **GaS, GaSe, and InSe**. It collects the project files needed to inspect the main thermal-transport results without reproducing the much larger working Google Drive archive.

The study combines:

- non-equilibrium molecular dynamics (NEMD) thermal transport;
- finite-length thermal-conductivity analysis;
- armchair/zigzag tensile-strain studies;
- metal, chalcogen, and bonded metal-chalcogen vacancy studies;
- MD-derived mass-weighted VDOS/PDOS;
- harmonic phonon dispersion and projected DOS calculations using LAMMPS + phonoLAMMPS/phonopy;
- validation against published GaS/GaSe/InSe literature.

## Main numerical results contained here

The summary tables are in [`results/`](results/). The present finite-size analysis gives approximately:

| Material | Extrapolated kappa_inf (W m^-1 K^-1) | Effective characteristic length (nm) |
|---|---:|---:|
| GaS | 14.05 | 34.19 |
| GaSe | 10.19 | 35.55 |
| InSe | 3.30 | 23.57 |

These values come from the shorter-length finite-size dataset and should be read together with [`docs/REPRODUCIBILITY_NOTES.md`](docs/REPRODUCIBILITY_NOTES.md).

At 6% tensile strain, the fixed-length armchair conductivities increase to about 12.95, 8.67, and 3.18 W m^-1 K^-1 for GaS, GaSe, and InSe, while the zigzag response is substantially weaker. Vacancy creation strongly suppresses thermal conductivity for all three materials.

## Repository map

- [`scripts/`](scripts/) - vacancy generation, VDOS/PDOS processing, NEMD post-processing, and harmonic phonon scripts.
- [`inputs/`](inputs/) - representative LAMMPS NEMD, VDOS, and phonon inputs.
- [`potentials/`](potentials/) - Tersoff EIAP files used in the archived workflows. Two variants are retained because the working project used differently formatted files in MD/VDOS and harmonic-phonon branches.
- [`structures/`](structures/) - pristine cells and representative strain/vacancy/phonon structures.
- [`results/strain/`](results/strain/) - Excel summaries, CSV table, and strain figures.
- [`results/vacancy/`](results/vacancy/) - Excel summaries, CSV table, and vacancy figures.
- [`results/vdos/`](results/vdos/) - pristine GaS/GaSe/InSe VDOS and representative GaS defect spectra, including raw/smoothed spectral tables and metadata.
- [`results/phonon/`](results/phonon/) - harmonic phonon dispersion/PDOS data and figures.
- [`results/validation/`](results/validation/) - finite-size workbooks, derived CSV, figures, and a literature benchmark table.
- [`validation/representative_nemd/`](validation/representative_nemd/) - compact raw NEMD examples containing the LAMMPS input, cumulative thermostat energy, and processed figures.
- [`docs/`](docs/) - file guide, exclusions, and important reproducibility notes.
- [`references/REFERENCES.md`](references/REFERENCES.md) - bibliography/DOIs instead of redistributed publisher PDFs.

## Effective thickness convention

The conductivity post-processing follows effective thicknesses of **8.712 A (GaS), 8.875 A (GaSe), and 8.775 A (InSe)**. These are the same thicknesses used by Wang et al., *J. Appl. Phys.* 125, 245104 (2019), DOI: 10.1063/1.5094663. For 2D materials, reported W m^-1 K^-1 values depend directly on the adopted thickness convention.

## Software

The scripts require Python plus common scientific packages. See [`requirements.txt`](requirements.txt). Harmonic phonon scripts additionally require a working LAMMPS installation with its Python module, plus phonoLAMMPS and phonopy.

A typical Python environment can be prepared with:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

LAMMPS itself should be installed separately with the pair styles/packages required by the provided input files.

## Reproducibility scope

This is a **curated GitHub package**, not a byte-for-byte mirror of the full simulation workspace. Very large raw trajectories (`profile.langevin`, `pdos.lammpstrj`), repeated LAMMPS logs, duplicate potential/script copies, and publisher PDFs are intentionally excluded. Representative NEMD raw energy files and processed VDOS/phonon tables are retained.

Known methodological/statistical limitations are documented openly in [`docs/REPRODUCIBILITY_NOTES.md`](docs/REPRODUCIBILITY_NOTES.md). In particular, the archived NEMD equilibration scripts include a 3D `fix npt ... iso` stage that should be re-evaluated for a vacuum-separated monolayer before treating the current package as the final publication archive.

## Citation

GitHub can read [`CITATION.cff`](CITATION.cff). Update the citation metadata with the associated manuscript DOI and final author/contributor list when available.

## License

No open-source or open-data license has been assigned automatically. See [`LICENSE_NOTICE.md`](LICENSE_NOTICE.md) before granting reuse rights, especially for interatomic-potential files or any third-party-derived assets.
