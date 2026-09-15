# Scripts

The scripts were curated from the working research archive. They are preserved close to their source form; case-specific constants and legacy comments have not been silently rewritten.

- `nemd/`: profile averaging and representative thermal-conductivity post-processing.
- `vacancy/make_vacancies.py`: reproducible vacancy generator with multiple defect classes and RNG seeds.
- `vdos/make_vdos_structures.py`: VDOS-cell/defect generator.
- `vdos/pdos_publication.py`: publication-oriented MD VDOS/PDOS processing.
- `phonon/`: material-specific harmonic phonon/PDOS calculations.

Read `../docs/REPRODUCIBILITY_NOTES.md` before reusing the transport input or interpreting case-specific constants.

### Vacancy generator configuration

`make_vacancies.py` is a batch script rather than an argparse CLI. Its source snapshot is configured for `INPUT_DATA_FILE = "GaSe_114p4nm.lmp"`, `OUTPUT_DIR = "GaSe_structures"`, and seeds 101/202/303. Edit the configuration block at the top before running it for a different source structure/material. The exact long-strip source filename referenced by that snapshot is not duplicated in this curated GitHub package; representative generated defect structures are included instead, while the full generated structure ensemble remains in the working archive.
