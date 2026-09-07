#!/usr/bin/env python3
"""
build_topology.py
=================
Reads the Packmol output PDB (stage1_packed.pdb) and writes:
  - stage1.data     : LAMMPS full-style data file with complete topology
  - mg_ow_table.dat : Tabulated 12-6-4 potential for Mg2+-OW (and Mg-O_anion)

Force fields:
  - Water     (SOL): SPC/E
  - Mg2+      (MG):  Li-Merz 12-6-4   Li & Merz, JCTC 2014, doi:10.1021/ct400526u
  - SO4 2-    (SO4): Cannon et al.     J. Phys. Chem. B 118, 2014
  - CO3 2-    (CO3): Raiteri et al.    J. Phys. Chem. C 114, 2010
  - HCO3-     (HCO): Raiteri et al.
  - NH4+      (NH4): OPLS-AA
  - NH3       (NH3): OPLS-AA

WARNING: Parameters marked [APPROX] are literature-sourced estimates.
Verify all values against the original papers before production runs.
"""

from pathlib import Path

import math
import sys
from collections import defaultdict
# =============================================================================
# FORCE FIELD PARAMETERS
# atom_types[name] = (mass, charge, epsilon_kcal_per_mol, sigma_angstrom)
# Convention: LAMMPS lj/cut/coul/long  (real units)
# NOTE: Mg LJ here covers Mg-Mg and Mg-non-O interactions only.
#       Mg-O pairs (water, carbonate, sulfate) are overridden by the
#       tabulated 12-6-4 potential in LAMMPS via pair_coeff table.
# =============================================================================

from forcefield import ATOM_TYPES, TYPE_NAMES, TYPE_INDEX, BOND_TYPES, BOND_INDEX, ANGLE_TYPES, ANGLE_INDEX

# =============================================================================
# RESIDUE -> ATOM TYPE MAP
# (residue_name, atom_name_in_PDB) -> atom_type_string
# =============================================================================

RES_ATOM_TYPE = {
    ("SOL", "OW"):  "OW",
    ("SOL", "HW1"): "HW",
    ("SOL", "HW2"): "HW",
    ("MG",  "MG"):  "MG",
    ("SO4", "S"):   "S_so4",
    ("SO4", "O1"):  "O_so4",
    ("SO4", "O2"):  "O_so4",
    ("SO4", "O3"):  "O_so4",
    ("SO4", "O4"):  "O_so4",
    ("NH4", "N"):   "N_nh4",
    ("NH4", "H1"):  "H_nh4",
    ("NH4", "H2"):  "H_nh4",
    ("NH4", "H3"):  "H_nh4",
    ("NH4", "H4"):  "H_nh4",
    ("CO3", "C"):   "C_co3",
    ("CO3", "O1"):  "O_co3",
    ("CO3", "O2"):  "O_co3",
    ("CO3", "O3"):  "O_co3",
    ("HCO", "C"):   "C_hco3",
    ("HCO", "O1"):  "O_hco3_do",
    ("HCO", "O2"):  "O_hco3_so",
    ("HCO", "O3"):  "O_hco3_oh",
    ("HCO", "H1"):  "H_hco3",
    ("NH3", "N"):   "N_nh3",
    ("NH3", "H1"):  "H_nh3",
    ("NH3", "H2"):  "H_nh3",
    ("NH3", "H3"):  "H_nh3",
}

# =============================================================================
# INTRAMOLECULAR CONNECTIVITY
# =============================================================================

RES_BONDS = {
    "SOL": [("OW","HW1"), ("OW","HW2")],
    "MG":  [],
    "SO4": [("S","O1"), ("S","O2"), ("S","O3"), ("S","O4")],
    "NH4": [("N","H1"), ("N","H2"), ("N","H3"), ("N","H4")],
    "CO3": [("C","O1"), ("C","O2"), ("C","O3")],
    "HCO": [("C","O1"), ("C","O2"), ("C","O3"), ("O3","H1")],
    "NH3": [("N","H1"), ("N","H2"), ("N","H3")],
}

RES_ANGLES = {
    "SOL": [("HW1","OW","HW2")],
    "MG":  [],
    "SO4": [("O1","S","O2"), ("O1","S","O3"), ("O1","S","O4"),
            ("O2","S","O3"), ("O2","S","O4"), ("O3","S","O4")],
    "NH4": [("H1","N","H2"), ("H1","N","H3"), ("H1","N","H4"),
            ("H2","N","H3"), ("H2","N","H4"), ("H3","N","H4")],
    "CO3": [("O1","C","O2"), ("O1","C","O3"), ("O2","C","O3")],
    "HCO": [("O1","C","O2"), ("O1","C","O3"), ("O2","C","O3"),
            ("C","O3","H1")],
    "NH3": [("H1","N","H2"), ("H1","N","H3"), ("H2","N","H3")],
}

# =============================================================================
# HELPERS
# =============================================================================

def _bond_key(at1, at2):
    """Canonical bond type string (try both orderings)."""
    for k in (f"{at1}-{at2}", f"{at2}-{at1}"):
        if k in BOND_INDEX:
            return k
    raise KeyError(f"Unknown bond type: {at1}-{at2}")


def _angle_key(at1, at2, at3):
    """Canonical angle type string (central atom = at2).
    Collapses O_hco3_* subtypes to O_hco3 for the shared angle lookup."""
    def _gen(s):
        return (s.replace("O_hco3_do", "O_hco3")
                 .replace("O_hco3_so", "O_hco3")
                 .replace("O_hco3_oh", "O_hco3"))
    for k in (f"{at1}-{at2}-{at3}", f"{at3}-{at2}-{at1}",
              _gen(f"{at1}-{at2}-{at3}"), _gen(f"{at3}-{at2}-{at1}")):
        if k in ANGLE_INDEX:
            return k
    raise KeyError(f"Unknown angle type: {at1}-{at2}-{at3}")


# =============================================================================
# PDB PARSER
# =============================================================================

def parse_pdb(filename):
    atoms = []
    with open(filename) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            atoms.append({
                "serial":  int(line[6:11]),
                "name":    line[12:16].strip(),
                "resname": line[17:20].strip(),
                "resseq":  int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
    return atoms


# =============================================================================
# BUILD LAMMPS DATA FILE
# =============================================================================

def build_lammps_data(pdb_file, out_file, box_lo=0.0, box_hi=50.0):
    atoms_raw = parse_pdb(pdb_file)

    # Group atoms by molecule (resname + resseq)
    molecules = defaultdict(list)
    mol_order = []
    for a in atoms_raw:
        key = (a["resname"], a["resseq"])
        if key not in molecules:
            mol_order.append(key)
        molecules[key].append(a)

    atoms_out  = []
    bonds_out  = []
    angles_out = []

    global_id = 1

    for mol_id, key in enumerate(mol_order, start=1):
        resname, _ = key
        res_atoms = molecules[key]
        local_map = {}  # atom_name -> global_id

        for a in res_atoms:
            atype_str = RES_ATOM_TYPE.get((resname, a["name"]))
            if atype_str is None:
                raise ValueError(f"Unrecognized atom: resname={resname!r} name={a['name']!r}")
            charge = ATOM_TYPES[atype_str][1]
            atoms_out.append((global_id, mol_id, TYPE_INDEX[atype_str],
                               charge, a["x"], a["y"], a["z"]))
            local_map[a["name"]] = global_id
            global_id += 1

        for n1, n2 in RES_BONDS.get(resname, []):
            at1 = RES_ATOM_TYPE[(resname, n1)]
            at2 = RES_ATOM_TYPE[(resname, n2)]
            try:
                bkey = _bond_key(at1, at2)
            except KeyError as e:
                print(f"  WARNING: {e} — skipping")
                continue
            bonds_out.append((BOND_INDEX[bkey], local_map[n1], local_map[n2]))

        for n1, n2, n3 in RES_ANGLES.get(resname, []):
            at1 = RES_ATOM_TYPE[(resname, n1)]
            at2 = RES_ATOM_TYPE[(resname, n2)]
            at3 = RES_ATOM_TYPE[(resname, n3)]
            try:
                akey = _angle_key(at1, at2, at3)
            except KeyError as e:
                print(f"  WARNING: {e} — skipping")
                continue
            angles_out.append((ANGLE_INDEX[akey],
                                local_map[n1], local_map[n2], local_map[n3]))

    # ---- verify total charge ----
    total_q = sum(a[3] for a in atoms_out)
    print(f"  Total system charge: {total_q:+.4f} e  (should be ~0.0)")
    if abs(total_q) > 0.1:
        print("  WARNING: non-zero total charge — check ion counts or FF charges")

    # ---- write ----
    with open(out_file, "w") as fh:
        fh.write("LAMMPS data file — Stage 1 (70 C, pH 9.0) Mg-carbonate-sulfate\n\n")
        fh.write(f"{len(atoms_out)} atoms\n")
        fh.write(f"{len(bonds_out)} bonds\n")
        fh.write(f"{len(angles_out)} angles\n")
        fh.write("0 dihedrals\n0 impropers\n\n")
        fh.write(f"{len(TYPE_NAMES)} atom types\n")
        fh.write(f"{len(BOND_TYPES)} bond types\n")
        fh.write(f"{len(ANGLE_TYPES)} angle types\n\n")
        fh.write(f"{box_lo:.4f} {box_hi:.4f} xlo xhi\n")
        fh.write(f"{box_lo:.4f} {box_hi:.4f} ylo yhi\n")
        fh.write(f"{box_lo:.4f} {box_hi:.4f} zlo zhi\n\n")

        fh.write("Masses\n\n")
        for i, name in enumerate(TYPE_NAMES, 1):
            fh.write(f"  {i:2d}  {ATOM_TYPES[name][0]:8.3f}  # {name}\n")
        fh.write("\n")

        # NOTE: deliberately NO "Pair Coeffs" section. Pairs are defined in stage1.in
        # via pair_style hybrid/overlay, which needs sub-style names the single-style
        # data-file format can't express. Emitting it here caused read_data to fail with
        # "Must define pair_style before Pair Coeffs".

        fh.write("Bond Coeffs  # harmonic\n\n")
        for i, (name, k, r0) in enumerate(BOND_TYPES, 1):
            fh.write(f"  {i:2d}  {k:.2f}  {r0:.4f}  # {name}\n")
        fh.write("\n")

        fh.write("Angle Coeffs  # harmonic\n\n")
        for i, (name, k, th) in enumerate(ANGLE_TYPES, 1):
            fh.write(f"  {i:2d}  {k:.2f}  {th:.2f}  # {name}\n")
        fh.write("\n")

        fh.write("Atoms  # full\n\n")
        for (aid, mid, atype, q, x, y, z) in atoms_out:
            fh.write(f"  {aid:7d}  {mid:6d}  {atype:2d}  {q:8.5f}"
                     f"  {x:11.4f}  {y:11.4f}  {z:11.4f}\n")
        fh.write("\n")

        if bonds_out:
            fh.write("Bonds\n\n")
            for i, (bt, a1, a2) in enumerate(bonds_out, 1):
                fh.write(f"  {i:8d}  {bt:2d}  {a1:7d}  {a2:7d}\n")
            fh.write("\n")

        if angles_out:
            fh.write("Angles\n\n")
            for i, (at, a1, a2, a3) in enumerate(angles_out, 1):
                fh.write(f"  {i:8d}  {at:2d}  {a1:7d}  {a2:7d}  {a3:7d}\n")
            fh.write("\n")

    print(f"  Written: {out_file}")
    print(f"  {len(atoms_out)} atoms | {len(bonds_out)} bonds | {len(angles_out)} angles")


# =============================================================================
# GENERATE 12-6-4 TABULATED POTENTIAL  (Mg2+ — oxygen interactions)
# =============================================================================

def write_mg_o_table(out_file, r_min=1.5, r_max=9.0, n_points=3000):
    """
    Tabulated 12-6-4 potential for Mg2+—O pairs.
    Written for three O types: OW, O_so4, O_co3/O_hco3_*

    U(r) = eps_ij * [(Rmin_ij/r)^12 - 2*(Rmin_ij/r)^6] - C4/r^4
    F(r) = -dU/dr

    Li-Merz 12-6-4 parameters for Mg2+ (SPC/E, [APPROX]):
      Rmin/2_Mg = 1.4120 A,  eps_Mg = 0.00935 kcal/mol,  C4 = 22.43 kcal.A4/mol
    Verify against: Li & Merz, JCTC 2014, Table S6 (doi:10.1021/ct400526u)

    Mixing rules (Lorentz-Berthelot for Rmin, geometric for eps):
      Rmin_ij = Rmin/2_Mg + Rmin/2_O
      eps_ij  = sqrt(eps_Mg * eps_O)
      C4 is atom-specific (on Mg only, not mixed)
    """

    # Li-Merz Mg parameters [APPROX — verify from paper]
    Rmin2_Mg = 1.4120   # A
    eps_Mg   = 0.00935  # kcal/mol
    C4       = 22.43    # kcal.A4/mol

    # O type definitions: (label, Rmin/2_O, eps_O)
    o_types = [
        ("MG_OW",    1.76830, 0.15535),  # SPC/E   OW  (Rmin/2 = sigma*2^(1/6)/2)
        ("MG_OSO4",  1.57500, 0.12800),  # Cannon  O_so4  [APPROX]
        ("MG_OCO3",  1.58000, 0.03700),  # Raiteri O_co3  [APPROX]
    ]

    dr = (r_max - r_min) / (n_points - 1)

    with open(out_file, "w") as fh:
        fh.write("# Tabulated 12-6-4 potentials: Mg2+ -- O pairs\n")
        fh.write("# Li-Merz JCTC 2014 [APPROX] — verify before production\n")
        fh.write("# U(r) = eps_ij[(Rmin_ij/r)^12 - 2(Rmin_ij/r)^6] - C4/r^4\n")
        fh.write("# Units: real (kcal/mol, Angstrom)\n\n")

        for label, Rmin2_O, eps_O in o_types:
            Rmin_ij = Rmin2_Mg + Rmin2_O
            eps_ij  = math.sqrt(eps_Mg * eps_O)

            fh.write(f"{label}\n")
            fh.write(f"N {n_points} R {r_min:.4f} {r_max:.4f}\n\n")

            for idx in range(n_points):
                r   = r_min + idx * dr
                x   = Rmin_ij / r
                U   = eps_ij * (x**12 - 2.0*x**6) - C4 / r**4
                dU  = eps_ij * (-12.0*x**12/r + 12.0*x**6/r) + 4.0*C4/r**5
                F   = -dU
                fh.write(f"{idx+1:6d}  {r:.6f}  {U:.8f}  {F:.8f}\n")
            fh.write("\n")

    print(f"  Written: {out_file}  ({n_points} pts × {len(o_types)} pair types)")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Stage 1 topology builder")
    print("=" * 60)

    here = Path(__file__).parent
    pdb_in    = here / "stage1_packed.pdb"
    data_out  = here / "stage1.data"
    table_out = here / "mg_o_table.dat"

    print(f"\n[1/2] Building LAMMPS data file from {pdb_in}...")
    build_lammps_data(pdb_in, data_out)

    print(f"\n[2/2] Generating Mg-O tabulated potentials...")
    write_mg_o_table(table_out)

    print("\nDone. Wrote stage1.data (no Pair Coeffs section) + mg_o_table.dat.")
    print("Workflow: packmol < packmol.inp  ->  python build_topology.py  ->  sbatch submit_stage1.sh")
