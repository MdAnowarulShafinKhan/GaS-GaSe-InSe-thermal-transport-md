#!/usr/bin/env python3
# =============================================================================
#  Monolayer GaS - harmonic phonon dispersion + total/projected phonon DOS
# =============================================================================
#
#  Method:
#    1. Build the four-atom X-M-M-X primitive cell.
#    2. Relax the in-plane cell and atomic coordinates with the GaS Tersoff EIAP.
#    3. Obtain second-order force constants by finite displacement through phonoLAMMPS.
#    4. Calculate Gamma-M-K-Gamma phonon dispersion and total/element-projected DOS.
#
#  Structural source:
#    Y. Karaaslan, J. Phys. D: Appl. Phys. 58 (2025) 015306,
#    Supporting Information, Table S6.
#
#  Paper values (Angstrom):
#                    a0       c        h       d_M-M    d_M-X
#    DFT            3.640   25.000    4.657    2.475    2.368
#    Tersoff EIAP   3.835   23.926    4.136    2.237    2.409
#
#  The starting bonded geometry below uses the EIAP values because the calculation
#  is performed with that EIAP. The z cell length is kept at 25 A because c is only
#  the vacuum-containing simulation-box height, not a bonded monolayer dimension.
#
#  Required input file next to this script:
#      GaS.tersoff
#
#  Run in WSL after activating the environment containing your existing LAMMPS
#  Python interface, phonopy and phonoLAMMPS:
#      python gas_phonon_dispersion_pdos.py
# =============================================================================

from pathlib import Path
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lammps import lammps
from phonolammps import Phonolammps
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections


# =============================================================================
# 1. USER SETTINGS
# =============================================================================
SYSTEM = "GaS"
M_ELEMENT = "Ga"
X_ELEMENT = "S"
POTENTIAL = "GaS.tersoff"
ELEMENTS = [M_ELEMENT, X_ELEMENT]  # LAMMPS type 1 = M, type 2 = X
MASSES = {1: 69.723000, 2: 32.060000}  # atomic mass units

# EIAP equilibrium geometry from Supporting Information Table S6.
A_GUESS = 3.835000       # in-plane lattice constant a0 (A)
M_M_GUESS = 2.237000    # vertical metal-metal distance d_M-M (A)
THICK_GUESS = 4.136000   # outer X-X monolayer thickness h (A)
D_MX_REFERENCE = 2.409000  # nearest metal-chalcogen bond d_M-X (A)

# The paper reports c = 23.926 A for its relaxed EIAP cell. Here 25 A is used
# as a standardized vacuum-containing box height; it does not change isolated-layer
# forces because the Tersoff cutoff is much shorter than the vacuum gap.
VACUUM_C = 25.0

# Finite-displacement supercell. The shortest in-plane repeat is > 20 A, much
# larger than twice the largest Tersoff cutoff used by this potential.
SUPERCELL = [[5, 0, 0],
             [0, 5, 0],
             [0, 0, 1]]

# Hexagonal reciprocal-space path for direct vectors
# a1=(a,0,0), a2=(a/2,sqrt(3)a/2,0).
# For this 60-degree convention, K=(2/3,1/3,0), not (1/3,1/3,0).
BAND_PATH = [[[0.0, 0.0, 0.0],
              [0.5, 0.0, 0.0],
              [2.0/3.0, 1.0/3.0, 0.0],
              [0.0, 0.0, 0.0]]]
BAND_LABELS = [r"$\Gamma$", "M", "K", r"$\Gamma$"]
BAND_NPTS = 101

DOS_MESH = [40, 40, 1]
DOS_SIGMA = 0.10       # THz
DOS_FREQ_PITCH = 0.02  # THz
DOS_MARGIN = 0.50      # THz beyond the mesh extrema
MASS_MATCH_TOL = 0.20  # amu
IMAGINARY_TOL = -0.10  # THz; more negative values trigger a warning

BASE_DIR = Path(__file__).resolve().parent


# =============================================================================
# 2. GEOMETRY AND FILE HELPERS
# =============================================================================
def expected_mx_bond(a, d_mm, thickness):
    """Nearest M-X bond implied by the X-M-M-X geometry."""
    in_plane = a / np.sqrt(3.0)
    vertical = 0.5 * (thickness - d_mm)
    return float(np.sqrt(in_plane**2 + vertical**2))


def build_primitive():
    """Build the four-atom S-Ga-Ga-S primitive cell."""
    a = A_GUESS
    z_m = M_M_GUESS / 2.0
    z_x = THICK_GUESS / 2.0
    z_center = VACUUM_C / 2.0

    lattice = np.array([
        [a,       0.0,                  0.0],
        [a / 2.0, a * np.sqrt(3.0)/2.0, 0.0],
        [0.0,     0.0,                  VACUUM_C],
    ], dtype=float)

    # X atoms occupy fractional (0,0); M atoms occupy fractional (1/3,1/3).
    cart = np.array([
        [0.0,     0.0,                  z_center + z_x],
        [0.0,     0.0,                  z_center - z_x],
        [a / 2.0, a * np.sqrt(3.0)/6.0, z_center + z_m],
        [a / 2.0, a * np.sqrt(3.0)/6.0, z_center - z_m],
    ], dtype=float)

    types = np.array([2, 2, 1, 1], dtype=int)  # X, X, M, M
    return lattice, cart, types


def write_lammps_data(path, lattice, cart, types):
    """Write a triclinic LAMMPS data file with atom_style atomic."""
    path = Path(path)
    lx = float(lattice[0, 0])
    xy = float(lattice[1, 0])
    ly = float(lattice[1, 1])
    lz = float(lattice[2, 2])

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{SYSTEM} primitive cell\n\n")
        handle.write(f"{len(cart)} atoms\n2 atom types\n\n")
        handle.write(f"0 {lx:.10f} xlo xhi\n")
        handle.write(f"0 {ly:.10f} ylo yhi\n")
        handle.write(f"0 {lz:.10f} zlo zhi\n")
        handle.write(f"{xy:.10f} 0 0 xy xz yz\n\n")
        handle.write("Masses\n\n")
        for atom_type, mass in MASSES.items():
            handle.write(f"{atom_type} {mass:.8f}\n")
        handle.write("\nAtoms # atomic\n\n")
        for atom_id, (atom_type, pos) in enumerate(zip(types, cart), start=1):
            handle.write(
                f"{atom_id} {int(atom_type)} "
                f"{pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n"
            )


def write_lammps_input(path, data_file):
    """Write the force-engine input consumed by phonoLAMMPS."""
    path = Path(path)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("units metal\n")
        handle.write("boundary p p p\n")
        handle.write("box tilt large\n")
        handle.write("atom_style atomic\n")
        handle.write(f"read_data {data_file}\n")
        handle.write("pair_style tersoff\n")
        handle.write(f"pair_coeff * * {POTENTIAL} {' '.join(ELEMENTS)}\n")
        handle.write("neighbor 2.0 bin\n")
        handle.write("neigh_modify delay 0 every 1 check yes\n")


def minimum_pair_distance(lattice, cart, types, type_a, type_b):
    """Minimum distance between two atom types, including nearby x-y images."""
    indices_a = np.flatnonzero(types == type_a)
    indices_b = np.flatnonzero(types == type_b)
    best = np.inf

    shifts = [
        i * lattice[0] + j * lattice[1]
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
    ]

    for i in indices_a:
        for j in indices_b:
            for shift in shifts:
                if type_a == type_b and i == j and np.allclose(shift, 0.0):
                    continue
                distance = np.linalg.norm((cart[j] + shift) - cart[i])
                best = min(best, float(distance))
    return best


def structural_metrics(lattice, cart, types):
    """Return a0, cell angle, thickness, d_M-M and nearest d_M-X."""
    a1 = lattice[0]
    a2 = lattice[1]
    a1_len = float(np.linalg.norm(a1))
    a2_len = float(np.linalg.norm(a2))
    cos_angle = np.dot(a1, a2) / (a1_len * a2_len)
    angle = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))

    x_indices = np.flatnonzero(types == 2)
    thickness = float(np.ptp(cart[x_indices, 2]))
    d_mm = minimum_pair_distance(lattice, cart, types, 1, 1)
    d_mx = minimum_pair_distance(lattice, cart, types, 1, 2)

    return {
        "a1": a1_len,
        "a2": a2_len,
        "angle": angle,
        "h": thickness,
        "d_mm": d_mm,
        "d_mx": d_mx,
    }


def print_metrics(label, metrics):
    print(
        f"  {label}: a1={metrics['a1']:.4f} A, a2={metrics['a2']:.4f} A, "
        f"angle={metrics['angle']:.3f} deg, h={metrics['h']:.4f} A, "
        f"d_M-M={metrics['d_mm']:.4f} A, d_M-X={metrics['d_mx']:.4f} A"
    )


# =============================================================================
# 3. RELAXATION
# =============================================================================
def relax_cell(lattice, cart, types):
    """Relax x, y, xy and all atomic coordinates while holding vacuum height fixed."""
    data_file = f"{SYSTEM}_unit_initial.data"
    write_lammps_data(data_file, lattice, cart, types)

    lmp = lammps(cmdargs=["-log", "none", "-screen", "none"])
    try:
        lmp.commands_string(f"""
            units metal
            boundary p p p
            atom_style atomic
            box tilt large
            read_data {data_file}
            pair_style tersoff
            pair_coeff * * {POTENTIAL} {' '.join(ELEMENTS)}
            neighbor 2.0 bin
            neigh_modify delay 0 every 1 check yes
            min_style cg
            fix br all box/relax x 0.0 y 0.0 xy 0.0 vmax 0.001
            minimize 1.0e-12 1.0e-12 20000 200000
        """)

        xlo = float(lmp.extract_global("boxxlo"))
        xhi = float(lmp.extract_global("boxxhi"))
        ylo = float(lmp.extract_global("boxylo"))
        yhi = float(lmp.extract_global("boxyhi"))
        zlo = float(lmp.extract_global("boxzlo"))
        xy = float(lmp.extract_global("xy"))

        rel_cart = np.asarray(lmp.gather_atoms("x", 1, 3), dtype=float).reshape(-1, 3)
        # Rebase coordinates because the regenerated data file uses lower bounds of zero.
        rel_cart -= np.array([xlo, ylo, zlo], dtype=float)
        rel_types = np.asarray(lmp.gather_atoms("type", 0, 1), dtype=int)
        forces = np.asarray(lmp.gather_atoms("f", 1, 3), dtype=float).reshape(-1, 3)

        e_atom = float(lmp.get_thermo("pe")) / int(lmp.get_natoms())
        max_force = float(np.linalg.norm(forces, axis=1).max())
    finally:
        lmp.close()

    rel_lattice = np.array([
        [xhi - xlo, 0.0,       0.0],
        [xy,        yhi - ylo, 0.0],
        [0.0,       0.0,       VACUUM_C],
    ], dtype=float)

    print(f"  relaxed energy/atom = {e_atom:.8f} eV")
    print(f"  maximum residual force = {max_force:.3e} eV/A")
    return rel_lattice, rel_cart, rel_types


# =============================================================================
# 4. FORCE CONSTANTS
# =============================================================================
def compute_phonon(rel_lattice, rel_cart, rel_types):
    data_file = f"{SYSTEM}_unit_relaxed.data"
    input_file = f"in.{SYSTEM}_phonon.lammps"
    write_lammps_data(data_file, rel_lattice, rel_cart, rel_types)
    write_lammps_input(input_file, data_file)

    bridge = Phonolammps(
        input_file,
        supercell_matrix=SUPERCELL,
        primitive_matrix=np.eye(3),
        show_log=False,
    )
    return bridge.get_phonopy_phonon()


# =============================================================================
# 5. DISPERSION
# =============================================================================
def compute_dispersion(phonon):
    qpoints, connections = get_band_qpoints_and_path_connections(
        BAND_PATH, npoints=BAND_NPTS
    )
    phonon.run_band_structure(
        qpoints,
        path_connections=connections,
        labels=BAND_LABELS,
    )
    band = phonon.band_structure
    frequencies = np.vstack(band.frequencies)
    distances = np.concatenate(band.distances)
    ticks = [segment[0] for segment in band.distances] + [band.distances[-1][-1]]
    return frequencies, distances, ticks


def save_dispersion_data(frequencies, distances):
    columns = [distances] + [frequencies[:, branch] for branch in range(frequencies.shape[1])]
    array = np.column_stack(columns)
    header = "q_distance " + " ".join(
        f"branch_{branch + 1}_THz" for branch in range(frequencies.shape[1])
    )
    np.savetxt(f"{SYSTEM}_dispersion.dat", array, header=header)


def plot_dispersion(frequencies, distances, ticks):
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for branch in range(frequencies.shape[1]):
        ax.plot(distances, frequencies[:, branch], color="#1a1a1a", lw=1.25)
    for tick in ticks[1:-1]:
        ax.axvline(tick, color="0.72", lw=0.8)
    ax.axhline(0.0, color="0.75", lw=0.8)
    ax.set_xticks(ticks)
    ax.set_xticklabels(BAND_LABELS, fontsize=13)
    ax.set_xlim(distances[0], distances[-1])
    ax.set_ylim(bottom=min(-0.2, float(frequencies.min()) - 0.2))
    ax.set_ylabel("Frequency (THz)", fontsize=13, fontweight="bold")
    ax.set_title(f"Monolayer {SYSTEM} - phonon dispersion", fontsize=13)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    fig.savefig(f"{SYSTEM}_dispersion.png", dpi=1200, bbox_inches="tight")
    fig.savefig(f"{SYSTEM}_dispersion.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {SYSTEM}_dispersion.png/.pdf and {SYSTEM}_dispersion.dat")


# =============================================================================
# 6. TOTAL DOS AND ELEMENT-PROJECTED DOS
# =============================================================================
def _auto_frequency_limits(mesh_frequencies):
    mesh_min = float(np.nanmin(mesh_frequencies))
    mesh_max = float(np.nanmax(mesh_frequencies))
    freq_min = min(-0.5, np.floor((mesh_min - DOS_MARGIN) * 2.0) / 2.0)
    freq_max = max(1.0, np.ceil((mesh_max + DOS_MARGIN) * 2.0) / 2.0)
    return float(freq_min), float(freq_max)


def compute_dos(phonon):
    phonon.run_mesh(DOS_MESH, with_eigenvectors=True, is_mesh_symmetry=False)
    mesh_dict = phonon.get_mesh_dict()
    freq_min, freq_max = _auto_frequency_limits(np.asarray(mesh_dict["frequencies"]))

    phonon.run_total_dos(
        freq_min=freq_min,
        freq_max=freq_max,
        freq_pitch=DOS_FREQ_PITCH,
        sigma=DOS_SIGMA,
    )
    total_dict = phonon.get_total_dos_dict()
    frequency = np.asarray(total_dict["frequency_points"])
    total = np.asarray(total_dict["total_dos"])

    phonon.run_projected_dos(
        freq_min=freq_min,
        freq_max=freq_max,
        freq_pitch=DOS_FREQ_PITCH,
        sigma=DOS_SIGMA,
    )
    projected_dict = phonon.get_projected_dos_dict()
    projected = np.asarray(projected_dict["projected_dos"])

    # Match atoms to the two user-specified masses. This matches each species directly by its specified atomic mass and therefore
    # avoids fragile mass-threshold classification.
    primitive_masses = np.asarray(phonon.primitive.masses, dtype=float)
    m_indices = np.flatnonzero(
        np.isclose(primitive_masses, MASSES[1], rtol=0.0, atol=MASS_MATCH_TOL)
    )
    x_indices = np.flatnonzero(
        np.isclose(primitive_masses, MASSES[2], rtol=0.0, atol=MASS_MATCH_TOL)
    )

    if len(m_indices) == 0 or len(x_indices) == 0:
        raise RuntimeError(
            "Could not assign primitive atoms to element masses. "
            f"Primitive masses: {primitive_masses}; expected {MASSES}"
        )
    if len(m_indices) + len(x_indices) != len(primitive_masses):
        raise RuntimeError(
            "At least one primitive atom was assigned ambiguously or not assigned. "
            f"Primitive masses: {primitive_masses}"
        )

    m_pdos = projected[m_indices].sum(axis=0)
    x_pdos = projected[x_indices].sum(axis=0)

    # Normalize by the total number of atoms in the primitive cell. Therefore,
    # the two element curves add to the total curve and the total integral is 3.
    atom_count = len(primitive_masses)
    return (
        frequency,
        total / atom_count,
        m_pdos / atom_count,
        x_pdos / atom_count,
        freq_min,
        freq_max,
    )


def save_dos_data(frequency, total, m_pdos, x_pdos):
    array = np.column_stack([frequency, total, m_pdos, x_pdos])
    header = (
        f"frequency_THz total_states_per_THz_per_atom "
        f"{M_ELEMENT}_contribution {X_ELEMENT}_contribution"
    )
    np.savetxt(f"{SYSTEM}_PDOS.dat", array, header=header)


def plot_dos(frequency, total, m_pdos, x_pdos, freq_min, freq_max):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(frequency, total, color="#1a1a1a", lw=1.9, label="Total")
    ax.plot(frequency, m_pdos, color="#c0392b", lw=1.7, label=M_ELEMENT)
    ax.plot(frequency, x_pdos, color="#2471a3", lw=1.7, label=X_ELEMENT)
    ax.axvline(0.0, color="0.75", lw=0.8)
    ax.set_xlim(freq_min, freq_max)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel("Frequency (THz)", fontsize=13, fontweight="bold")
    ax.set_ylabel("PDOS (states / THz / atom)", fontsize=13, fontweight="bold")
    ax.set_title(f"Monolayer {SYSTEM} - projected phonon DOS", fontsize=13)
    ax.legend(fontsize=11, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    fig.savefig(f"{SYSTEM}_PDOS.png", dpi=1200, bbox_inches="tight")
    fig.savefig(f"{SYSTEM}_PDOS.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {SYSTEM}_PDOS.png/.pdf and {SYSTEM}_PDOS.dat")


# =============================================================================
# 7. MAIN
# =============================================================================
def main():
    os.chdir(BASE_DIR)

    potential_path = BASE_DIR / POTENTIAL
    if not potential_path.is_file():
        raise FileNotFoundError(
            f"Required potential file not found: {potential_path}\n"
            f"Put {POTENTIAL} in the same directory as this script."
        )

    implied_mx = expected_mx_bond(A_GUESS, M_M_GUESS, THICK_GUESS)
    if not np.isclose(implied_mx, D_MX_REFERENCE, atol=0.003, rtol=0.0):
        raise ValueError(
            f"Geometry inconsistency: a0, h and d_M-M imply d_M-X={implied_mx:.6f} A, "
            f"but the paper value is {D_MX_REFERENCE:.6f} A."
        )

    print(f"{SYSTEM} harmonic phonon calculation")
    print(f"  potential: {potential_path}")
    print(
        f"  EIAP reference: a0={A_GUESS:.3f} A, h={THICK_GUESS:.3f} A, "
        f"d_M-M={M_M_GUESS:.3f} A, d_M-X={D_MX_REFERENCE:.3f} A"
    )

    print("[1/4] building and relaxing the primitive cell ...")
    lattice, cart, types = build_primitive()
    print_metrics("initial", structural_metrics(lattice, cart, types))
    rel_lattice, rel_cart, rel_types = relax_cell(lattice, cart, types)
    relaxed_metrics = structural_metrics(rel_lattice, rel_cart, rel_types)
    print_metrics("relaxed", relaxed_metrics)
    write_lammps_data(f"{SYSTEM}_unit_relaxed.data", rel_lattice, rel_cart, rel_types)

    print("[2/4] computing second-order force constants through phonoLAMMPS ...")
    phonon = compute_phonon(rel_lattice, rel_cart, rel_types)

    print("[3/4] computing Gamma-M-K-Gamma dispersion ...")
    frequencies, distances, ticks = compute_dispersion(phonon)
    save_dispersion_data(frequencies, distances)
    plot_dispersion(frequencies, distances, ticks)
    minimum_frequency = float(frequencies.min())
    maximum_frequency = float(frequencies.max())
    print(
        f"  branches={frequencies.shape[1]}, min={minimum_frequency:.4f} THz, "
        f"max={maximum_frequency:.4f} THz"
    )
    if minimum_frequency < IMAGINARY_TOL:
        print(
            "  WARNING: a significant imaginary branch is present. Do not treat the "
            "structure as dynamically stable until relaxation, geometry, supercell "
            "and potential mapping have been checked."
        )

    print(f"[4/4] computing total DOS and {M_ELEMENT}/{X_ELEMENT}-projected DOS ...")
    frequency, total, m_pdos, x_pdos, freq_min, freq_max = compute_dos(phonon)
    save_dos_data(frequency, total, m_pdos, x_pdos)
    plot_dos(frequency, total, m_pdos, x_pdos, freq_min, freq_max)

    print("done.")


if __name__ == "__main__":
    main()
