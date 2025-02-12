# Colors for heatmap
palette_length <- 100
my_color <- colorRampPalette(c("Darkblue", "white", "brown"))(palette_length)

# Disease colors
tfs_colors <- c(
    "ESR1" = "#1f77b4",
    "FOS" = "#ff7f0e",
    "JUN" = "#2ca02c",
    "MYC" = "#e7ba52",
    "SP1" = "#843c39",
    "STAT1" = "#9467bd",
    "AP1" = "#a55194",
    "AR" = "#8c6d31",
    "HIF1A" = "#5254a3",
    "NFKB" = "#ad494a",
    "NFKB1" = "#9ecae1",
    "RELA" = "#973aa8",
    "SP3" = "#fff75e",
    "SPI1" = "#264653",
    "TP53" = "#941c2f"
    )

disease_colors <- c(
    "healthy" = "#808080",

    "RA" = "#264653",
    "PS" = "#287271",
    "PSA" = "#2a9d8f",
    "CD" = "#e76f51",
    "UC" = "#e9c46a",
    "SLE" = "#941c2f",
    "MS" = "#8ab17d",

    "asthma" = "#ea698b",
    "COPD" = "#c05299",
    "cirrhosis" = "#973aa8",

    "sepsis" = "#ef233c",

    "HIV" = "#e7ecef",
    "HBV" = "#a3cef1",
    "COVID" = "#6096ba",
    "flu" = "#274c77",

    "BRCA" = "#fff75e",
    "NPC" = "#fdb833",
    "HNSCC" = "#d9981a",
    "CRC" = "#9e7524"
)

cell_colors <- c(
    "B" = "#7bc6d6",
    "Plasma" = "#025566",

    "pDC" = "#a7c957",
    "DC" = "#6a994e",
    "Mono" = "#386641",

    "T_CD4_Naive" = "#fff3b0",
    "T_CD4_NonNaive" = "#e09f3e",
    "T_CD8_Naive" = "#9e2a2b",
    "T_CD8_NonNaive" = "#540b0e",

    "UTC" = "#88657f",
    "ILC" = "#67253a",

    "Cycling_cells" = "#d4a373",
    "Progenitors" = "#ccd5ae",

    "Platelets" = "#808080",  # To remove
    "RBC" = "#000000"         # To remove
)
