#!/usr/bin/env python3
# =============================================================================
#  make_vacancies.py
#
#  Creates realistic vacancy defects in a LAMMPS data file (built with ATOMSK)
#  for monolayer GaS / GaSe / InSe NEMD thermal-conductivity studies.
#
#  Defect types (matching the paper blueprint):
#    D1 = chalcogen monovacancy  (remove single S / Se atoms)
#    D2 = metal    monovacancy   (remove single Ga / In atoms)
#    D3 = divacancy              (remove a bonded metal + chalcogen PAIR)
#
#  "Industry-realistic" features built in:
#    1. MINIMUM DISTANCE between vacancies  -> no artificial clustering;
#       mimics the dilute, randomly distributed vacancies seen in real
#       as-grown / irradiated samples.
#    2. PROTECTED ZONES along the heat-flow (x) direction -> never deletes
#       atoms inside the fixed walls or Langevin bath regions, so the
#       thermostats act on perfect lattice (required for clean NEMD).
#    3. EXACT atom-count bookkeeping -> concentration is defined as
#       (% of TOTAL atoms in the file), identical for D1/D2/D3, so the three
#       defect types are compared at equal missing-atom count.
#    4. FULL REPRODUCIBILITY -> every file is generated from a named RNG seed,
#       and a CSV log records every deleted atom (id, type, x, y, z).
#    5. BATCH MODE -> generates all concentrations x realizations x defect
#       types for one material in a single run.
#
#  Usage:
#    1. Edit the CONFIG block below.
#    2. python3 make_vacancies.py
#
#  Requires: numpy  (pip install numpy)
# =============================================================================

import os
import sys
import csv
import numpy as np

# =============================================================================
# CONFIG  -- edit this block only
# =============================================================================

INPUT_DATA_FILE = "GaSe_114p4nm.lmp"   # pristine ATOMSK-built LAMMPS data file
OUTPUT_DIR      = "GaSe_structures"         # created if it does not exist
MATERIAL_TAG    = "GaSe"                       # used in output file names

# --- Atom-type IDs in YOUR data file (check the Masses section!) -------------
METAL_TYPES     = [1]     # e.g. Ga (or In for InSe)
CHALCOGEN_TYPES = [2]     # e.g. S  (or Se for GaSe/InSe)

# --- What to generate --------------------------------------------------------
DEFECT_TYPES    = ["D1", "D2", "D3"]                       # any subset
CONCENTRATIONS  = [0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0]     # % of TOTAL atoms
REALIZATION_SEEDS = [101, 202, 303]                        # 3 independent structures

# --- Realism constraints -----------------------------------------------------
MIN_VACANCY_DIST = 8.0   # Angstrom. Minimum distance between vacancy sites
                         # (for D3, measured between pair midpoints).
                         # ~2x the potential cutoff avoids overlapping strain
                         # fields / artificial clustering at low concentration.
AUTO_RELAX_DIST  = True  # If the target count cannot be placed (high conc.),
                         # shrink MIN_VACANCY_DIST in 0.5 A steps down to
                         # MIN_DIST_FLOOR and warn. If False -> hard error.
MIN_DIST_FLOOR   = 3.5   # Never go below roughly one bond length.

BOND_CUTOFF      = 3.2   # Angstrom. Max metal-chalcogen distance that counts
                         # as a "bond" when building D3 pairs (Ga-S ~2.3 A,
                         # Ga-Se ~2.4 A, In-Se ~2.6 A, so 3.2 is safe).

# --- Protected zones along the heat-flow direction (x) -----------------------
# Atoms with  x < xlo + X_EXCLUDE_LO   or   x > xhi - X_EXCLUDE_HI
# are NEVER deleted. Set these to cover: fixed wall + Langevin bath + a small
# buffer, on each end. Example: 4 A wall + 20 A bath + 6 A buffer = 30 A.
X_EXCLUDE_LO = 30.0      # Angstrom, measured from the low-x box edge
X_EXCLUDE_HI = 30.0      # Angstrom, measured from the high-x box edge

# --- Periodicity for the min-distance check ----------------------------------
PERIODIC_X = False       # NEMD strip: x runs wall-to-wall  -> not periodic
PERIODIC_Y = True        # width direction is periodic      -> True

# =============================================================================
# END CONFIG
# =============================================================================


# ---------------------------------------------------------------- file parsing
def read_lammps_data(path):
    """Parse a LAMMPS data file (atom_style atomic or charge)."""
    with open(path, "r") as f:
        lines = f.readlines()

    header = {"natoms": None, "ntypes": None,
              "xlo": None, "xhi": None, "ylo": None, "yhi": None,
              "zlo": None, "zhi": None, "tilt": None}
    masses_lines, atoms_raw = [], []
    atom_style_hint = None

    section = "header"
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.split("#")[0].strip()

        # Section switches
        first_word = line.strip().split("#")[0].strip()
        if first_word in ("Masses", "Atoms", "Velocities",
                          "Bonds", "Angles", "Pair Coeffs", "PairIJ Coeffs"):
            section = first_word
            if first_word == "Atoms" and "#" in line:
                atom_style_hint = line.split("#")[1].strip()
            i += 2  # skip the blank line after the section keyword
            continue

        if section == "header" and stripped:
            toks = stripped.split()
            if "atoms" in toks and len(toks) == 2:
                header["natoms"] = int(toks[0])
            elif toks[-2:] == ["atom", "types"]:
                header["ntypes"] = int(toks[0])
            elif toks[-2:] == ["xlo", "xhi"]:
                header["xlo"], header["xhi"] = float(toks[0]), float(toks[1])
            elif toks[-2:] == ["ylo", "yhi"]:
                header["ylo"], header["yhi"] = float(toks[0]), float(toks[1])
            elif toks[-2:] == ["zlo", "zhi"]:
                header["zlo"], header["zhi"] = float(toks[0]), float(toks[1])
            elif toks[-3:] == ["xy", "xz", "yz"]:
                header["tilt"] = stripped
        elif section == "Masses" and stripped:
            masses_lines.append(line.rstrip("\n"))
        elif section == "Atoms" and stripped:
            atoms_raw.append(stripped.split())
        # Velocities of a pristine 0 K ATOMSK file are not needed; LAMMPS will
        # re-create them with `velocity all create`. We deliberately drop them.
        i += 1

    if not atoms_raw:
        sys.exit("ERROR: no Atoms section found in " + path)

    # Detect column layout: atomic = id type x y z [img], charge = id type q x y z [img]
    ncol = len(atoms_raw[0])
    if atom_style_hint == "charge" or (atom_style_hint is None and ncol in (6, 9)
                                       and _looks_like_charge(atoms_raw)):
        style = "charge"
        xcol = 3
    else:
        style = "atomic"
        xcol = 2

    ids   = np.array([int(a[0]) for a in atoms_raw])
    types = np.array([int(a[1]) for a in atoms_raw])
    xyz   = np.array([[float(a[xcol]), float(a[xcol+1]), float(a[xcol+2])]
                      for a in atoms_raw])
    extra = [a[2] if style == "charge" else None for a in atoms_raw]  # charges

    if header["natoms"] != len(ids):
        print("WARNING: header atom count != parsed atoms; using parsed count.")
        header["natoms"] = len(ids)

    return header, masses_lines, ids, types, xyz, extra, style


def _looks_like_charge(atoms_raw):
    """Heuristic: 3rd column non-integer-like -> charge style."""
    try:
        v = float(atoms_raw[0][2])
        return abs(v - round(v)) > 1e-9 or abs(v) < 1e-9 and len(atoms_raw[0]) == 6
    except ValueError:
        return False


def write_lammps_data(path, header, masses_lines, ids, types, xyz, extra, style,
                      comment):
    with open(path, "w") as f:
        f.write("# " + comment + "\n\n")
        f.write("%d atoms\n" % len(ids))
        f.write("%d atom types\n\n" % header["ntypes"])
        f.write("%.8f %.8f xlo xhi\n" % (header["xlo"], header["xhi"]))
        f.write("%.8f %.8f ylo yhi\n" % (header["ylo"], header["yhi"]))
        f.write("%.8f %.8f zlo zhi\n" % (header["zlo"], header["zhi"]))
        if header["tilt"]:
            f.write(header["tilt"] + "\n")
        f.write("\nMasses\n\n")
        for m in masses_lines:
            f.write(m + "\n")
        f.write("\nAtoms # %s\n\n" % style)
        for k in range(len(ids)):
            if style == "charge":
                f.write("%d %d %s %.8f %.8f %.8f\n" %
                        (k + 1, types[k], extra[k], xyz[k, 0], xyz[k, 1], xyz[k, 2]))
            else:
                f.write("%d %d %.8f %.8f %.8f\n" %
                        (k + 1, types[k], xyz[k, 0], xyz[k, 1], xyz[k, 2]))


# ---------------------------------------------------------- geometry utilities
def min_image_dist(p, q, Lx, Ly):
    """In-plane distance with minimum-image convention where periodic."""
    d = p - q
    if PERIODIC_X:
        d[0] -= Lx * round(d[0] / Lx)
    if PERIODIC_Y:
        d[1] -= Ly * round(d[1] / Ly)
    return np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)


def far_from_all(site, accepted, dmin, Lx, Ly):
    for s in accepted:
        if min_image_dist(site, s, Lx, Ly) < dmin:
            return False
    return True


# ------------------------------------------------------------- site selection
def select_monovacancies(rng, cand_idx, xyz, n_target, dmin, Lx, Ly):
    """Greedy random selection of single sites with a min-distance constraint."""
    order = rng.permutation(cand_idx)
    accepted_sites, accepted_idx = [], []
    for i in order:
        if len(accepted_idx) == n_target:
            break
        if far_from_all(xyz[i], accepted_sites, dmin, Lx, Ly):
            accepted_idx.append(i)
            accepted_sites.append(xyz[i])
    return accepted_idx


def select_divacancies(rng, metal_idx, chalc_idx, xyz, n_pairs, dmin,
                       Lx, Ly):
    """Pick bonded metal+chalcogen pairs; min distance between pair midpoints."""
    chalc_pos = xyz[chalc_idx]
    order = rng.permutation(metal_idx)
    accepted_mid, pairs = [], []
    used_chalc = set()
    for im in order:
        if len(pairs) == n_pairs:
            break
        # nearest chalcogen neighbour of this metal atom (min-image in y)
        d = chalc_pos - xyz[im]
        if PERIODIC_Y:
            d[:, 1] -= Ly * np.round(d[:, 1] / Ly)
        if PERIODIC_X:
            d[:, 0] -= Lx * np.round(d[:, 0] / Lx)
        dist = np.sqrt((d ** 2).sum(axis=1))
        j = int(np.argmin(dist))
        ic = chalc_idx[j]
        if dist[j] > BOND_CUTOFF or ic in used_chalc:
            continue
        mid = 0.5 * (xyz[im] + xyz[ic])
        if far_from_all(mid, accepted_mid, dmin, Lx, Ly):
            pairs.append((im, ic))
            accepted_mid.append(mid)
            used_chalc.add(ic)
    return pairs


# --------------------------------------------------------------- one structure
def build_one(defect, conc, seed, header, masses_lines, ids, types, xyz, extra,
              style, logwriter):
    rng = np.random.RandomState(seed)
    N   = len(ids)
    Lx  = header["xhi"] - header["xlo"]
    Ly  = header["yhi"] - header["ylo"]

    n_remove = int(round(conc / 100.0 * N))
    if defect == "D3":
        n_pairs  = max(1, n_remove // 2)
        n_remove = 2 * n_pairs   # keep it even; report the true value

    # candidate atoms: right sublattice AND inside the allowed x-window
    xlo_ok = header["xlo"] + X_EXCLUDE_LO
    xhi_ok = header["xhi"] - X_EXCLUDE_HI
    in_window = (xyz[:, 0] >= xlo_ok) & (xyz[:, 0] <= xhi_ok)

    metal_mask = np.isin(types, METAL_TYPES) & in_window
    chalc_mask = np.isin(types, CHALCOGEN_TYPES) & in_window
    metal_idx  = np.where(metal_mask)[0]
    chalc_idx  = np.where(chalc_mask)[0]

    dmin = MIN_VACANCY_DIST
    while True:
        if defect == "D1":
            removed = select_monovacancies(rng, chalc_idx, xyz, n_remove,
                                           dmin, Lx, Ly)
        elif defect == "D2":
            removed = select_monovacancies(rng, metal_idx, xyz, n_remove,
                                           dmin, Lx, Ly)
        elif defect == "D3":
            pairs = select_divacancies(rng, metal_idx, chalc_idx, xyz,
                                       n_remove // 2, dmin, Lx, Ly)
            removed = [i for p in pairs for i in p]
        else:
            sys.exit("Unknown defect type: " + defect)

        if len(removed) == n_remove:
            break
        if AUTO_RELAX_DIST and dmin - 0.5 >= MIN_DIST_FLOOR:
            dmin -= 0.5
            print("    [info] could only place %d/%d sites at d_min=%.1f A; "
                  "retrying with d_min=%.1f A"
                  % (len(removed), n_remove, dmin + 0.5, dmin))
            rng = np.random.RandomState(seed)   # restart deterministically
        else:
            sys.exit("ERROR: cannot place %d vacancies with d_min >= %.1f A "
                     "for %s at %.2f%%. Lower MIN_VACANCY_DIST or enable "
                     "AUTO_RELAX_DIST." % (n_remove, dmin, defect, conc))

    keep = np.ones(N, dtype=bool)
    keep[removed] = False

    # log every deleted atom for the Supporting Information
    for i in removed:
        logwriter.writerow([defect, conc, seed, int(ids[i]), int(types[i]),
                            "%.5f" % xyz[i, 0], "%.5f" % xyz[i, 1],
                            "%.5f" % xyz[i, 2], "%.2f" % dmin])

    actual_pct = 100.0 * len(removed) / N
    fname = "%s_%s_%.2fpct_seed%d.data" % (MATERIAL_TAG, defect, conc, seed)
    fpath = os.path.join(OUTPUT_DIR, fname)
    comment = ("%s | %s | target %.2f%% | actual %.4f%% (%d/%d atoms) | "
               "seed %d | d_min %.1f A | x-window [%.1f, %.1f] A"
               % (MATERIAL_TAG, defect, conc, actual_pct, len(removed), N,
                  seed, dmin, xlo_ok, xhi_ok))
    write_lammps_data(fpath, header, masses_lines,
                      ids[keep], types[keep], xyz[keep],
                      [e for e, k in zip(extra, keep) if k],
                      style, comment)
    print("  wrote %-45s  (%d atoms removed, actual %.4f%%, d_min %.1f A)"
          % (fname, len(removed), actual_pct, dmin))


# ----------------------------------------------------------------------- main
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    header, masses_lines, ids, types, xyz, extra, style = \
        read_lammps_data(INPUT_DATA_FILE)

    print("Read %s: %d atoms, style '%s', box %.1f x %.1f x %.1f A"
          % (INPUT_DATA_FILE, len(ids), style,
             header["xhi"] - header["xlo"],
             header["yhi"] - header["ylo"],
             header["zhi"] - header["zlo"]))
    for t in METAL_TYPES + CHALCOGEN_TYPES:
        print("  type %d: %d atoms" % (t, int((types == t).sum())))

    logpath = os.path.join(OUTPUT_DIR, MATERIAL_TAG + "_deleted_atoms_log.csv")
    with open(logpath, "w", newline="") as lf:
        logwriter = csv.writer(lf)
        logwriter.writerow(["defect", "target_pct", "seed", "atom_id",
                            "atom_type", "x", "y", "z", "d_min_used_A"])
        for defect in DEFECT_TYPES:
            for conc in CONCENTRATIONS:
                for seed in REALIZATION_SEEDS:
                    build_one(defect, conc, seed, header, masses_lines,
                              ids, types, xyz, extra, style, logwriter)
    print("\nDeletion log written to " + logpath)
    print("Done: %d structures."
          % (len(DEFECT_TYPES) * len(CONCENTRATIONS) * len(REALIZATION_SEEDS)))


if __name__ == "__main__":
    main()
