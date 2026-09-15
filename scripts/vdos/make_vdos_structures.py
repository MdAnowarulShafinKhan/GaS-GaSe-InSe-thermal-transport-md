#!/usr/bin/env python3
"""
make_vdos_structures.py   (v2 - standards-aligned with make_vacancies.py)
=========================================================================
Generate the VDOS "story" structures for the GaS / GaSe / InSe monolayer
paper from ONE pristine ~10 nm x 10 nm fully periodic sheet (LAMMPS data
file, atom_style atomic OR charge).

This version is formally aligned with the NEMD defect generator
(make_vacancies.py) so that the VDOS defect ensembles are statistically
the same as the kappa defect ensembles they are meant to explain:

  ALIGNED WITH make_vacancies.py
  * MIN_VACANCY_DIST default 8.0 A, auto-relaxed in -0.5 A steps down to
    a floor of 3.5 A (same schedule, same floor). The RNG is restarted
    deterministically at each relaxation step, so the final structure is
    a pure function of (seed, final d_min) - exactly like the old script.
  * Multi-seed batch mode: defect cases are generated for EVERY seed in
    --seeds (default 101 202 303, the NEMD realization seeds).
  * Per-atom CSV deletion log with the same columns:
      defect, target_pct, seed, atom_id, atom_type, x, y, z, d_min_used_A
  * atom_style 'charge' is parsed and preserved (charges carried through).
  * Concentration = % of TOTAL atoms; D3 midpoint spacing; bond cutoff
    3.2 A; velocities deliberately dropped.

  DIFFERENT ON PURPOSE (correct for the VDOS geometry)
  * Fully periodic in x AND y; NO protected zones (no walls/baths here).
  * Vacancies are ENFORCED 50/50 between the two equivalent sublayers
    (top/bottom chalcogen decks for D1, upper/lower metal decks for D2;
    D3 seed-metal decks alternate). This is stricter than the old
    script's on-average balance: an asymmetric sheet acquires a bending
    bias that contaminates the low-frequency flexural (ZA) region of the
    spectrum, which is exactly where the D1-vs-D2 story lives.
  * Removal count is forced EVEN for all defect types, so D1, D2 and D3
    remove *exactly* the same number of atoms (fixes the old script's
    odd-count mismatch where D3 could remove one atom fewer).
  * Pristine and strain cases are seed-independent -> generated once.

Output per material (with default seeds 101 202 303): 12 files
    <MAT>_pristine.data                          (VDOS case 1)
    <MAT>_D1_V<X>_5pct_seed<S>.data   x3 seeds   (VDOS case 2)
    <MAT>_D2_V<M>_5pct_seed<S>.data   x3 seeds   (VDOS case 3)
    <MAT>_D3_V<MX>_5pct_seed<S>.data  x3 seeds   (VDOS case 4)
    <MAT>_Spar_6pct.data                         (VDOS case 5, main text)
    <MAT>_Sperp_6pct.data                        (VDOS case 6, SI)
plus <MAT>_deleted_atoms_log.csv and <MAT>_generation_log.txt

Usage:
    python make_vdos_structures.py GaS_10x10.data  --material GaS
    python make_vdos_structures.py GaSe_10x10.data --material GaSe
    python make_vdos_structures.py InSe_10x10.data --material InSe \
           --conc 5.0 --strain 6.0 --min-sep 8.0 --seeds 101 202 303

IMPORTANT for statistical parity with your NEMD structures: check the
'd_min_used_A' column of the OLD script's CSV for the 5 % cases. If it
relaxed to some value (e.g. 6.5 A), pass that value here as --min-sep so
both generators sample the same defect-spacing distribution.

VDOS run reminders (also written into every file header / log):
  * defected structures: energy-minimize (or short damped dynamics) first,
    then NVT-equilibrate at 300 K (~100 ps), then pure NVE for the
    velocity dump (every 5 fs, ~200 ps). No thermostat during the dump.
  * strained structures: NVT/NVE with FIXED box only - NPT would relax
    the strain away and you would measure a pristine spectrum.
"""

import argparse
import csv
import os
import sys
from datetime import datetime

import numpy as np

# ----------------------------------------------------------------------------
# Element identification by mass (amu)
# ----------------------------------------------------------------------------
KNOWN = {"Ga": 69.723, "In": 114.818, "S": 32.06, "Se": 78.971}
METALS = {"Ga", "In"}
CHALCOGENS = {"S", "Se"}
MASS_TOL = 2.0  # amu

MIN_DIST_FLOOR = 3.5      # A  (same floor as make_vacancies.py)
RELAX_STEP = 0.5          # A  (same -0.5 A schedule)


def identify_element(mass):
    best, best_d = None, 1e9
    for el, m in KNOWN.items():
        d = abs(mass - m)
        if d < best_d:
            best, best_d = el, d
    if best_d > MASS_TOL:
        raise ValueError(
            f"Mass {mass:.3f} amu does not match Ga/In/S/Se within "
            f"{MASS_TOL} amu. Is this really a GaS/GaSe/InSe data file?")
    return best


# ----------------------------------------------------------------------------
# LAMMPS data I/O  (atomic or charge, orthogonal box)
# ----------------------------------------------------------------------------
class Structure:
    def __init__(self):
        self.box = {}          # {'x': (lo, hi), 'y': ..., 'z': ...}
        self.masses = {}       # {type: mass}
        self.ids = None        # (N,) int   original atom IDs
        self.types = None      # (N,) int
        self.xyz = None        # (N, 3) float
        self.charges = None    # (N,) float or None
        self.style = "atomic"

    @property
    def natoms(self):
        return len(self.ids)

    def L(self, axis):
        lo, hi = self.box[axis]
        return hi - lo

    def clone_frame(self):
        """Copy box/masses/style only (for strain builder)."""
        s2 = Structure()
        s2.box = dict(self.box)
        s2.masses = dict(self.masses)
        s2.style = self.style
        s2.ids = self.ids.copy()
        s2.types = self.types.copy()
        s2.xyz = self.xyz.copy()
        s2.charges = None if self.charges is None else self.charges.copy()
        return s2


def _looks_like_charge(row):
    """Heuristic when no '# charge' hint: 3rd column non-integer-like."""
    try:
        v = float(row[2])
        return abs(v - round(v)) > 1e-9
    except (ValueError, IndexError):
        return False


def read_data(path):
    s = Structure()
    with open(path) as f:
        lines = f.readlines()

    n_declared = None
    style_hint = None
    section = None
    atoms_rows = []

    i = 0
    while i < len(lines):
        raw = lines[i]
        ln = raw.split("#")[0].strip()
        first = ln.split()[0] if ln else ""
        if first in ("Masses", "Atoms", "Velocities", "Bonds", "Angles",
                     "Pair", "PairIJ"):
            section = first
            if first == "Atoms" and "#" in raw:
                style_hint = raw.split("#")[1].strip()
            i += 1
            continue
        if not ln:
            i += 1
            continue
        low = ln.lower()
        if section is None:
            toks = ln.split()
            if low.endswith(" atoms") or low == "atoms" or \
               (len(toks) == 2 and toks[1] == "atoms"):
                n_declared = int(toks[0])
            elif low.endswith("xlo xhi"):
                s.box["x"] = (float(toks[0]), float(toks[1]))
            elif low.endswith("ylo yhi"):
                s.box["y"] = (float(toks[0]), float(toks[1]))
            elif low.endswith("zlo zhi"):
                s.box["z"] = (float(toks[0]), float(toks[1]))
            elif low.endswith("xy xz yz"):
                if any(abs(float(v)) > 1e-8 for v in toks[:3]):
                    raise ValueError(
                        "Triclinic (tilted) box detected. Rebuild the "
                        "pristine sheet with an orthogonal cell.")
        elif section == "Masses":
            p = ln.split()
            s.masses[int(p[0])] = float(p[1])
        elif section == "Atoms":
            atoms_rows.append(ln.split())
        # Velocities etc.: deliberately skipped (re-created by LAMMPS)
        i += 1

    if not atoms_rows:
        raise ValueError("No Atoms section found - is this a LAMMPS data file?")

    ncol = len(atoms_rows[0])
    if style_hint == "charge" or (style_hint is None and ncol in (6, 9)
                                  and _looks_like_charge(atoms_rows[0])):
        s.style = "charge"
        xcol = 3
    else:
        s.style = "atomic"
        xcol = 2

    atoms_rows.sort(key=lambda r: int(r[0]))
    s.ids = np.array([int(r[0]) for r in atoms_rows], dtype=int)
    s.types = np.array([int(r[1]) for r in atoms_rows], dtype=int)
    s.xyz = np.array([[float(r[xcol]), float(r[xcol + 1]),
                       float(r[xcol + 2])] for r in atoms_rows], dtype=float)
    if s.style == "charge":
        s.charges = np.array([float(r[2]) for r in atoms_rows], dtype=float)

    if n_declared is not None and n_declared != len(atoms_rows):
        print(f"  [warn] header says {n_declared} atoms, Atoms section has "
              f"{len(atoms_rows)}; using the section.")
    for axis in ("x", "y", "z"):
        if axis not in s.box:
            raise ValueError(f"Missing {axis} box bounds in header.")
    return s


def write_data(struct, keep_mask, path, comment):
    """Write kept atoms, renumbering IDs 1..N; preserves charge style."""
    types = struct.types[keep_mask]
    xyz = struct.xyz[keep_mask]
    q = None if struct.charges is None else struct.charges[keep_mask]
    ntypes = max(struct.masses)
    with open(path, "w") as f:
        f.write(f"# {comment}\n")
        f.write(f"# generated by make_vdos_structures.py v2 on "
                f"{datetime.now():%Y-%m-%d %H:%M}\n\n")
        f.write(f"{len(types)} atoms\n{ntypes} atom types\n\n")
        for axis in ("x", "y", "z"):
            lo, hi = struct.box[axis]
            f.write(f"{lo:.10f} {hi:.10f} {axis}lo {axis}hi\n")
        f.write("\nMasses\n\n")
        for t in sorted(struct.masses):
            f.write(f"{t} {struct.masses[t]:.4f}\n")
        f.write(f"\nAtoms # {struct.style}\n\n")
        for k in range(len(types)):
            if struct.style == "charge":
                f.write(f"{k + 1} {types[k]} {q[k]:.6f} "
                        f"{xyz[k, 0]:.10f} {xyz[k, 1]:.10f} {xyz[k, 2]:.10f}\n")
            else:
                f.write(f"{k + 1} {types[k]} "
                        f"{xyz[k, 0]:.10f} {xyz[k, 1]:.10f} {xyz[k, 2]:.10f}\n")


# ----------------------------------------------------------------------------
# Geometry helpers (fully periodic in x and y; z is the vacuum direction)
# ----------------------------------------------------------------------------
def min_image_dxy(p, q, Lx, Ly):
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dx -= Lx * round(dx / Lx)
    dy -= Ly * round(dy / Ly)
    return (dx * dx + dy * dy) ** 0.5


def min_image_d3(p, q, Lx, Ly):
    d = p - q
    d[0] -= Lx * round(d[0] / Lx)
    d[1] -= Ly * round(d[1] / Ly)
    return float(np.linalg.norm(d))


def far_from_all(site, accepted, dmin, Lx, Ly):
    for s in accepted:
        if min_image_dxy(site, s, Lx, Ly) < dmin:
            return False
    return True


# ----------------------------------------------------------------------------
# Site selection with the make_vacancies.py relaxation contract:
# try at d_min; if the full count cannot be placed, restart the RNG from the
# seed and retry at d_min - 0.5, down to the 3.5 A floor; else hard error.
# ----------------------------------------------------------------------------
def _mono_attempt(rng, struct, top, bot, n_top, n_bot, dmin, Lx, Ly):
    """One deterministic attempt at fixed d_min. Returns list or None."""
    accepted_idx, accepted_xy = [], []
    for pool, n_want in ((rng.permutation(top), n_top),
                         (rng.permutation(bot), n_bot)):
        placed = 0
        for idx in pool:
            if placed == n_want:
                break
            p = struct.xyz[idx]
            if far_from_all(p, accepted_xy, dmin, Lx, Ly):
                accepted_idx.append(int(idx))
                accepted_xy.append(p)
                placed += 1
        if placed < n_want:
            return None
    return accepted_idx


def select_monovacancies(struct, sublattice_index, n_remove, dmin0, seed,
                         label):
    """Enforced 50/50 top/bottom split + old-script relaxation schedule."""
    Lx, Ly = struct.L("x"), struct.L("y")
    z_mid = np.median(struct.xyz[:, 2])
    z = struct.xyz[sublattice_index, 2]
    top = sublattice_index[z >= z_mid]
    bot = sublattice_index[z < z_mid]
    n_top = n_remove // 2 + (n_remove % 2)
    n_bot = n_remove // 2
    if len(top) < n_top or len(bot) < n_bot:
        sys.exit(f"ERROR: {label}: sublayer too small for requested count.")

    dmin = dmin0
    while True:
        rng = np.random.default_rng(seed)      # deterministic restart
        removed = _mono_attempt(rng, struct, top, bot, n_top, n_bot,
                                dmin, Lx, Ly)
        if removed is not None:
            return removed, dmin
        if dmin - RELAX_STEP >= MIN_DIST_FLOOR:
            print(f"    [info] {label}: could not place {n_remove} sites at "
                  f"d_min={dmin:.1f} A; retrying with d_min="
                  f"{dmin - RELAX_STEP:.1f} A")
            dmin -= RELAX_STEP
        else:
            sys.exit(f"ERROR: cannot place {n_remove} vacancies with "
                     f"d_min >= {MIN_DIST_FLOOR:.1f} A for {label}. "
                     f"Lower the concentration or the floor.")


def _div_attempt(rng, struct, m_top, m_bot, chalc_index, chalc_xyz,
                 n_pairs, dmin, bond_cut, Lx, Ly):
    removed = set()
    centres = []
    pools = [rng.permutation(m_top), rng.permutation(m_bot)]
    ptrs = [0, 0]
    side = 0
    while len(centres) < n_pairs:
        if ptrs[0] >= len(pools[0]) and ptrs[1] >= len(pools[1]):
            return None
        if ptrs[side] >= len(pools[side]):
            side = 1 - side
            continue
        m_idx = int(pools[side][ptrs[side]])
        ptrs[side] += 1
        if m_idx in removed:
            continue
        # nearest not-yet-removed chalcogen (3D minimum image)
        p = struct.xyz[m_idx]
        d = chalc_xyz - p
        d[:, 0] -= Lx * np.round(d[:, 0] / Lx)
        d[:, 1] -= Ly * np.round(d[:, 1] / Ly)
        dist = np.sqrt((d ** 2).sum(axis=1))
        order = np.argsort(dist)
        c_idx, c_d = -1, 1e9
        for j in order[:8]:
            cand = int(chalc_index[j])
            if cand not in removed:
                c_idx, c_d = cand, float(dist[j])
                break
        if c_idx < 0 or c_d > bond_cut:
            continue
        centre = 0.5 * (struct.xyz[m_idx] + struct.xyz[c_idx])
        if far_from_all(centre, centres, dmin, Lx, Ly):
            removed.update((m_idx, c_idx))
            centres.append(centre)
            side = 1 - side                    # alternate decks
    return sorted(removed)


def select_divacancies(struct, metal_index, chalc_index, n_pairs, dmin0,
                       seed, bond_cut, label):
    Lx, Ly = struct.L("x"), struct.L("y")
    z_mid = np.median(struct.xyz[:, 2])
    mz = struct.xyz[metal_index, 2]
    m_top = metal_index[mz >= z_mid]
    m_bot = metal_index[mz < z_mid]
    chalc_xyz = struct.xyz[chalc_index].copy()

    dmin = dmin0
    while True:
        rng = np.random.default_rng(seed)
        removed = _div_attempt(rng, struct, m_top, m_bot, chalc_index,
                               chalc_xyz, n_pairs, dmin, bond_cut, Lx, Ly)
        if removed is not None:
            return removed, dmin
        if dmin - RELAX_STEP >= MIN_DIST_FLOOR:
            print(f"    [info] {label}: could not place {n_pairs} pairs at "
                  f"d_min={dmin:.1f} A; retrying with d_min="
                  f"{dmin - RELAX_STEP:.1f} A")
            dmin -= RELAX_STEP
        else:
            sys.exit(f"ERROR: cannot place {n_pairs} divacancy pairs with "
                     f"d_min >= {MIN_DIST_FLOOR:.1f} A for {label}.")


def build_strained(struct, axis, strain):
    s2 = struct.clone_frame()
    ax = {"x": 0, "y": 1}[axis]
    lo, hi = struct.box[axis]
    s2.box[axis] = (lo, lo + (hi - lo) * (1.0 + strain))
    s2.xyz[:, ax] = lo + (s2.xyz[:, ax] - lo) * (1.0 + strain)
    return s2


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate the VDOS-case structures (pristine, D1/D2/D3 "
                    "at one concentration x several seeds, S_par, S_perp) "
                    "from one pristine GaS/GaSe/InSe sheet.")
    ap.add_argument("input", help="pristine LAMMPS data file "
                                  "(atom_style atomic or charge)")
    ap.add_argument("--material", required=True,
                    choices=["GaS", "GaSe", "InSe"])
    ap.add_argument("--conc", type=float, default=5.0,
                    help="defect concentration in %% of TOTAL atoms "
                         "(default 5.0; same convention as make_vacancies.py)")
    ap.add_argument("--strain", type=float, default=6.0,
                    help="uniaxial tensile strain in %% (default 6.0)")
    ap.add_argument("--min-sep", type=float, default=8.0,
                    help="minimum in-plane separation between vacancy sites "
                         "(D3: pair midpoints), Angstrom. Default 8.0 = the "
                         "make_vacancies.py default. SET THIS TO THE "
                         "d_min_used_A VALUE FROM YOUR NEMD 5%% CSV LOG for "
                         "statistical parity. Auto-relaxes in 0.5 A steps to "
                         "a 3.5 A floor, restarting the RNG each step.")
    ap.add_argument("--bond-cut", type=float, default=3.2,
                    help="max M-X distance counted as a bond for D3 "
                         "(default 3.2 A, same as make_vacancies.py)")
    ap.add_argument("--seeds", type=int, nargs="+", default=[101, 202, 303],
                    help="realization seeds for the defect cases "
                         "(default: 101 202 303, the NEMD seeds)")
    ap.add_argument("--outdir", default="vdos_structures")
    args = ap.parse_args()

    conc_frac = args.conc / 100.0
    strain_frac = args.strain / 100.0
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Reading {args.input} ...")
    s = read_data(args.input)
    Lx, Ly = s.L("x"), s.L("y")
    print(f"  {s.natoms} atoms, style '{s.style}', "
          f"box {Lx:.2f} x {Ly:.2f} A in-plane")

    # ---- classify species ----
    type_elem = {t: identify_element(m) for t, m in s.masses.items()}
    metal_types = [t for t, e in type_elem.items() if e in METALS]
    chalc_types = [t for t, e in type_elem.items() if e in CHALCOGENS]
    if not metal_types or not chalc_types:
        sys.exit("ERROR: could not find both a metal and a chalcogen type.")
    m_el = type_elem[metal_types[0]]
    x_el = type_elem[chalc_types[0]]
    expected = {"GaS": ("Ga", "S"), "GaSe": ("Ga", "Se"), "InSe": ("In", "Se")}
    if (m_el, x_el) != expected[args.material]:
        print(f"  [warn] --material {args.material} but file contains "
              f"{m_el}/{x_el}; check you loaded the right file!")

    metal_index = np.where(np.isin(s.types, metal_types))[0]
    chalc_index = np.where(np.isin(s.types, chalc_types))[0]
    print(f"  {len(metal_index)} {m_el} (metal) + "
          f"{len(chalc_index)} {x_el} (chalcogen)")

    n_remove = int(round(conc_frac * s.natoms))
    if n_remove % 2:
        n_remove += 1        # identical atom count for D1/D2/D3
    n_pairs = n_remove // 2
    print(f"  {args.conc:.2f}% of total atoms -> {n_remove} atoms removed "
          f"per defect case ({n_remove / len(chalc_index):.1%} of the "
          f"{x_el} sublattice for D1; {n_pairs} bonded pairs for D3)")

    mat = args.material
    summary = [
        f"make_vdos_structures.py v2 summary ({datetime.now():%Y-%m-%d %H:%M})",
        f"input file        : {args.input}   (style {s.style})",
        f"material          : {mat} ({m_el}/{x_el})",
        f"total atoms       : {s.natoms}",
        f"box (A)           : {Lx:.3f} x {Ly:.3f}",
        f"concentration     : {args.conc:.2f} % of total atoms "
        f"({n_remove} atoms / defect case)",
        f"strain            : {args.strain:.2f} %",
        f"min separation    : {args.min_sep:.2f} A requested "
        f"(relax -0.5 A steps, floor {MIN_DIST_FLOOR} A)",
        f"seeds             : {args.seeds}",
        "",
    ]

    def out(name):
        return os.path.join(args.outdir, name)

    keep_all = np.ones(s.natoms, dtype=bool)

    # ---- per-atom deletion log (same columns as make_vacancies.py) ----
    csv_path = out(f"{mat}_deleted_atoms_log.csv")
    csv_f = open(csv_path, "w", newline="")
    logw = csv.writer(csv_f)
    logw.writerow(["defect", "target_pct", "seed", "atom_id", "atom_type",
                   "x", "y", "z", "d_min_used_A"])

    # ---- 1. pristine ----
    f1 = out(f"{mat}_pristine.data")
    write_data(s, keep_all, f1, f"{mat} pristine reference, VDOS case 1")
    summary.append(f"[1] pristine          -> {os.path.basename(f1)}  "
                   f"({s.natoms} atoms)")

    # ---- 2-4. defect cases, one file per seed ----
    defect_specs = [
        ("D1", f"V{x_el}", "chalcogen monovacancy"),
        ("D2", f"V{m_el}", "metal monovacancy"),
        ("D3", f"V{m_el}{x_el}", "adjacent M-X divacancy"),
    ]
    for dcode, dtag, dname in defect_specs:
        for seed in args.seeds:
            label = f"{dcode} {dtag} seed {seed}"
            if dcode in ("D1", "D2"):
                sub = chalc_index if dcode == "D1" else metal_index
                removed, dused = select_monovacancies(
                    s, sub, n_remove, args.min_sep, seed, label)
            else:
                removed, dused = select_divacancies(
                    s, metal_index, chalc_index, n_pairs, args.min_sep,
                    seed, args.bond_cut, label)

            for i in removed:
                logw.writerow([dcode, f"{args.conc:.2f}", seed,
                               int(s.ids[i]), int(s.types[i]),
                               f"{s.xyz[i, 0]:.5f}", f"{s.xyz[i, 1]:.5f}",
                               f"{s.xyz[i, 2]:.5f}", f"{dused:.2f}"])

            keep = np.ones(s.natoms, dtype=bool)
            keep[np.array(removed, dtype=int)] = False
            actual_pct = 100.0 * len(removed) / s.natoms
            fname = f"{mat}_{dcode}_{dtag}_{args.conc:.0f}pct_seed{seed}.data"
            write_data(s, keep, out(fname),
                       f"{mat} {dcode} ({dname}): {len(removed)} atoms "
                       f"removed = {actual_pct:.4f}% of total | seed {seed} "
                       f"| d_min {dused:.2f} A | 50/50 deck split enforced "
                       f"| minimize + NVT-equilibrate before the NVE "
                       f"velocity dump")
            summary.append(f"[{dcode}] seed {seed:<4}      -> {fname}  "
                           f"({keep.sum()} atoms, d_min used {dused:.2f} A)")

    csv_f.close()

    # ---- 5. strain parallel (x) ----
    s_par = build_strained(s, "x", strain_frac)
    f5 = out(f"{mat}_Spar_{args.strain:.0f}pct.data")
    write_data(s_par, keep_all, f5,
               f"{mat} S-parallel: +{args.strain:.1f}% affine strain along "
               f"x. RUN NVT/NVE WITH FIXED BOX - NPT would relax the strain.")
    summary.append(f"[5] S_par  +{args.strain:.0f}%       -> "
                   f"{os.path.basename(f5)}  "
                   f"(box x: {Lx:.2f} -> {s_par.L('x'):.2f} A)")

    # ---- 6. strain perpendicular (y), SI ----
    s_perp = build_strained(s, "y", strain_frac)
    f6 = out(f"{mat}_Sperp_{args.strain:.0f}pct.data")
    write_data(s_perp, keep_all, f6,
               f"{mat} S-perpendicular (SI): +{args.strain:.1f}% affine "
               f"strain along y. RUN NVT/NVE WITH FIXED BOX.")
    summary.append(f"[6] S_perp +{args.strain:.0f}%       -> "
                   f"{os.path.basename(f6)}  "
                   f"(box y: {Ly:.2f} -> {s_perp.L('y'):.2f} A)")

    summary += [
        "",
        f"per-atom deletion log : {os.path.basename(csv_path)}",
        "",
        "VDOS run reminders:",
        "  * defected structures: minimize (or short damped run), then",
        "    NVT-equilibrate at 300 K ~100 ps, then pure NVE for the dump",
        "  * dump all velocities every 5 fs for ~200 ps, no thermostat",
        "  * strained structures: never NPT / box relax",
        "  * for parity with NEMD: --min-sep should equal the d_min_used_A",
        "    of the old make_vacancies.py CSV at the same concentration",
    ]
    log = out(f"{mat}_generation_log.txt")
    with open(log, "w") as f:
        f.write("\n".join(summary) + "\n")
    print("\n".join(summary))
    print(f"\nAll files in ./{args.outdir}/")


if __name__ == "__main__":
    main()
