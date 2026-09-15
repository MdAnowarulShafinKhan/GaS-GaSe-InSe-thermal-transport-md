# Representative NEMD inputs

These are representative production inputs copied from the working archive:

- `GaS_strain_6pct_armchair.lmp`
- `GaS_D1_5pct_vacancy.lmp`

The source inputs expect the relevant structure and Tersoff file to be present under the filenames referenced inside the script. Paths have **not** been rewritten automatically; adjust them for your local directory layout.

Important: the archived equilibration protocol contains the 3D `fix npt ... iso` issue described in `../../docs/REPRODUCIBILITY_NOTES.md`. Treat these as provenance/reproducibility records, not as an endorsement of that barostat choice for future 2D-vacuum simulations.
