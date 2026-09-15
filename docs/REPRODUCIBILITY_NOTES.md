# Reproducibility notes and known limitations

This file records issues identified during the curation/audit of the working project. They are documented rather than silently corrected so that the public repository accurately represents the simulation history.

## 1. 2D NEMD barostatting

Representative transport scripts use a 3D simulation with `fix npt ... iso` during equilibration. In LAMMPS, `iso` couples pressure control in x, y, and z. For a vacuum-separated monolayer, pressure-controlling the artificial vacuum direction is not the preferred 2D equilibration protocol. The strain workflow also contains a stage in which z is pressure controlled while x is deformed.

Before treating the current transport results as a final publication archive, validate a corrected 2D protocol on representative cases (in-plane relaxation only; fixed vacuum z) and quantify whether the conductivity changes materially.

## 2. Statistical replication

`make_vacancies.py` supports seeds 101, 202, and 303. The production vacancy results inspected during curation are primarily single-realization (`seed101`) values. Accordingly, the vacancy curves in the supplied Excel workbooks should not be presented as multi-realization mean +/- SD unless the missing production simulations are completed.

The strain summary likewise represents single production values per strain/direction in the inspected archive. Small zigzag changes should therefore be interpreted cautiously until seed-to-seed uncertainty is quantified.

## 3. Vacancy concentration convention

The vacancy generator defines concentration as the percentage of **all atoms removed**, not the percentage of a single sublattice. For a stoichiometric MX sheet, a 1% total-atom metal-vacancy concentration corresponds to roughly 2% of metal sites.

Protected NEMD end/bath zones are excluded from vacancy placement while the percentage is referenced to the full atom count, so the local defect density in the eligible interior region is slightly higher than the global quoted percentage.

For D3, a bonded metal-chalcogen divacancy removes two atoms per defect center. Comparing D3 with monovacancies at equal removed-atom percentage therefore also changes the number density of scattering centers.

## 4. Transport vs VDOS defect ensembles

The VDOS generator enforces approximate 50/50 vacancy balance between equivalent upper/lower sublayers to avoid artificial flexural asymmetry. The transport defect generator does not impose exactly the same constraint. VDOS defect configurations should therefore be described as symmetry-balanced representative configurations, not exact replicas of every transport realization.

## 5. VDOS processing

The archive contains more than one Welch processing recipe. Some output filenames contain `4p096ps` even when the metadata shows a 16384-frame block at an effective 1 fs sampling interval (16.384 ps block duration). Use the metadata files, not the legacy filename, as the authoritative processing record.

For direct cross-condition comparison, use a common effective time step, block length, smoothing width, normalization, and frequency range. The raw `pdos.lammpstrj` trajectories are excluded from this GitHub package because of size.

## 6. VDOS duration comments

Some working LAMMPS comments refer to approximately 500 ps NVE while inspected metadata corresponds to about 100 ps of stored trajectory. Treat the metadata/actual step counts as authoritative.

## 7. NEMD duration/comments

The inspected strain production run uses 12,000,000 steps at 0.5 fs = 6 ns, despite stale comments elsewhere referring to 8 ns. Some comments also refer to 315/285 K reservoirs while the executable variables use 330/270 K. The code values are authoritative.

## 8. Temperature-gradient fit comment

The representative post-processing script describes a "middle 70%" fit, but its numerical limits are close to 0.10L to 0.90L (about 80% of the cell). Check the intended fit window before reusing the script.

## 9. Finite-size characteristic length

The reciprocal finite-size fit provides an **effective characteristic length**. It should not automatically be interpreted as a uniquely measured phonon mean free path, especially under defects or strain where heat capacity, velocities, dispersion, and scattering can all change.

## 10. Effective thickness

Conductivity is normalized with effective thicknesses 8.712 A (GaS), 8.875 A (GaSe), and 8.775 A (InSe). These match Wang et al. (2019). Because 2D W m^-1 K^-1 values scale inversely with the chosen thickness, state the convention explicitly in any manuscript or comparison.
