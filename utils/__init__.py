AXICEL_SIGS = {
    "AXICEL_good": ["KLRK1", "CCL22", "CD45RA", "SOX11", "CD19", "SIGLEC5"],
    "AXICEL_bad": [
        "IL18R1",
        "GPC4",
        "KIR3DL2",
        "ITGB8",
        "PSMB5",
        "RPS6KB1",
        "BCL2",
        "TNFSF4",
        "SERPINA9",
        "DUSP5",
        "NBN",
        "GLUD1",
        "ESR1",
        "CD45RO",
        "ARID1A",
        "KLRB1",
        "SLC16A1",
    ],
}

b_cell_score_genes = [
    "BLK",
    "CD19",
    "MS4A1",
    "TNFRSF17",
    "FCRL2",
    "FAM30A",
    "PNOC",
    "SPIB",
    "TCL1A",
]

axicel_signatures = {
    "Linear-Good": AXICEL_SIGS["AXICEL_good"],
    "Linear-Bad": AXICEL_SIGS["AXICEL_bad"],
    "Log2-Good": ['CD19', 'CYBB', 'SIGLEC5'],
    "Log2-Bad": ['AKT1', 'ARID1A', 'C7', 'CD45RB',
                 'DUSP5', 'MET', 'NCAM1', 'PDGFA', 'RPS6KB1']
}