import femm
# Initialize
femm.openfemm()
femm.newdocument(0)
femm.mi_probdef(0, 'millimeters', 'planar', 1.e-8)

# custom BH curve points (H in A/m, B in T)
bh_data = [
    [0, 0.00],
    [107, 0.03],
    [194, 0.07],
    [344, 0.13],
    [501, 0.21],
    [944, 0.41],
    [1425, 0.60],
    [2031, 0.77],
    [2941, 0.94],
    [4563, 1.12],
    [7828, 1.29]
]

# Add material with non-linear properties
# Parameters: Name, mux muy Hc HcAng Nlam LamFill LamType CrossSection
femm.mi_addmaterial('CustomIron', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)