# Files intentionally excluded from the GitHub package

The working Google Drive contains substantially more data than is appropriate for a normal Git repository. The following categories are intentionally not copied into this curated package:

- publisher/literature PDFs - replaced by DOI references in `references/REFERENCES.md`;
- the research-planning/blueprint PDF - useful internally, but not a generated research result;
- large raw NEMD `profile.langevin` files (often tens of MB per case);
- raw `pdos.lammpstrj` velocity trajectories;
- repeated `log.lammps` files;
- repeated identical/near-identical copies of scripts and Tersoff files in each production directory;
- duplicate PDF exports of figures when the PNG and underlying data table are already included;
- unrelated files outside the GaS/GaSe/InSe project.

The full working archive should be preserved separately. If a journal or reviewer requires raw trajectories, deposit them in an appropriate research-data repository rather than inflating the Git history.
