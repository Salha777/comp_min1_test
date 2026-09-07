ATOM_TYPES = {
    # --- SPC/E water ---
    "OW":        (16.000,  -0.84760,  0.15535,  3.16560),  # SPC/E
    "HW":        ( 1.008,  +0.42380,  0.00000,  0.00001),  # SPC/E (dummy sigma>0)

    # --- Mg2+  Li-Merz 12-6-4 (SPC/E parametrization) ---
    # [APPROX]  Rmin/2=1.4120 A, eps=0.0093 kcal/mol; C4=22.43 kcal.A4/mol
    # Background 12-6 params below; Mg-O pairs are replaced by table.
    "MG":        (24.305,  +2.00000,  0.00935,  2.51440),  # [APPROX] Li-Merz 2014

    # --- SO4 2-  Cannon et al. 2014 ---
    # [APPROX]  q(S)=+1.284, q(O)=-0.821 (net=-2.0 ✓)
    "S_so4":     (32.060,  +1.28400,  0.09480,  3.55000),  # [APPROX]
    "O_so4":     (16.000,  -0.82100,  0.12800,  3.15000),  # [APPROX]

    # --- CO3 2-  Raiteri et al. 2010 ---
    # [APPROX]  q(C)=+1.123, q(O)=-1.041 (net=-2.0 ✓)
    "C_co3":     (12.011,  +1.12300,  0.00000,  3.50000),  # [APPROX]
    "O_co3":     (16.000,  -1.04100,  0.03700,  3.16000),  # [APPROX]

    # --- HCO3-  Raiteri et al. (adapted) ---
    # [APPROX]  net = +0.908 -0.932 -0.932 -0.610 +0.400 = -1.166... adjust below
    # q(C)=+0.964  q(O_do)=-0.879  q(O_so)=-0.879  q(O_oh)=-0.606  q(H)=+0.400
    # net = 0.964 - 0.879 - 0.879 - 0.606 + 0.400 = -1.000 ✓
    "C_hco3":    (12.011,  +0.96400,  0.00000,  3.50000),  # [APPROX]
    "O_hco3_do": (16.000,  -0.87900,  0.03700,  3.16000),  # C=O   [APPROX]
    "O_hco3_so": (16.000,  -0.87900,  0.03700,  3.16000),  # C-O-  [APPROX]
    "O_hco3_oh": (16.000,  -0.60600,  0.03700,  3.16000),  # C-OH  [APPROX]
    "H_hco3":    ( 1.008,  +0.40000,  0.00000,  0.00001),  # O-H   [APPROX]

    # --- NH4+  OPLS-AA ---
    # q(N)=-0.30, q(H)=+0.325  ->  net = -0.30 + 4*0.325 = +1.0 ✓
    "N_nh4":     (14.007,  -0.30000,  0.17000,  3.25000),  # OPLS-AA [APPROX]
    "H_nh4":     ( 1.008,  +0.32500,  0.00000,  0.00001),  # OPLS-AA [APPROX]

    # --- NH3  OPLS-AA ---
    # q(N)=-1.02, q(H)=+0.34   ->  net = -1.02 + 3*0.34 = 0.0 ✓
    "N_nh3":     (14.007,  -1.02000,  0.17000,  3.25000),  # OPLS-AA [APPROX]
    "H_nh3":     ( 1.008,  +0.34000,  0.00000,  0.00001),  # OPLS-AA [APPROX]
}

TYPE_NAMES = list(ATOM_TYPES.keys())
TYPE_INDEX = {name: i+1 for i, name in enumerate(TYPE_NAMES)}

# =============================================================================
# BOND TYPES  (harmonic: k [kcal/mol/A^2], r0 [A])
# Water bonds are constrained by SHAKE; large k is nominal only.
# =============================================================================

BOND_TYPES = [
    ("OW-HW",              1000.00,  1.0000),  # SPC/E  (SHAKE)
    ("S_so4-O_so4",         400.00,  1.4900),  # [APPROX] Cannon
    ("N_nh4-H_nh4",         400.00,  1.0120),  # OPLS-AA [APPROX]
    ("C_co3-O_co3",         900.00,  1.2900),  # Raiteri [APPROX]
    ("C_hco3-O_hco3_do",    900.00,  1.2400),  # HCO3 C=O  [APPROX]
    ("C_hco3-O_hco3_so",    700.00,  1.3400),  # HCO3 C-O- [APPROX]
    ("C_hco3-O_hco3_oh",    700.00,  1.3400),  # HCO3 C-OH [APPROX]
    ("O_hco3_oh-H_hco3",    553.00,  0.9700),  # HCO3 O-H  [APPROX]
    ("N_nh3-H_nh3",         400.00,  1.0120),  # OPLS-AA [APPROX]
]

BOND_INDEX = {b[0]: i+1 for i, b in enumerate(BOND_TYPES)}

# =============================================================================
# ANGLE TYPES  (harmonic: k [kcal/mol/rad^2], theta0 [degrees])
# =============================================================================

ANGLE_TYPES = [
    ("HW-OW-HW",                    55.00,  109.47),  # SPC/E  (SHAKE)
    ("O_so4-S_so4-O_so4",          100.00,  109.47),  # [APPROX] Cannon
    ("H_nh4-N_nh4-H_nh4",           40.00,  109.47),  # OPLS-AA [APPROX]
    ("O_co3-C_co3-O_co3",          100.00,  120.00),  # Raiteri [APPROX]
    ("O_hco3-C_hco3-O_hco3",       100.00,  120.00),  # HCO3  [APPROX]
    ("C_hco3-O_hco3_oh-H_hco3",     35.00,  110.00),  # HCO3 C-O-H [APPROX]
    ("H_nh3-N_nh3-H_nh3",           40.00,  106.70),  # OPLS-AA [APPROX]
]

ANGLE_INDEX = {a[0]: i+1 for i, a in enumerate(ANGLE_TYPES)}
