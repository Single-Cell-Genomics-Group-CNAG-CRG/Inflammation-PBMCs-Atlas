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


# Create the named vector
cell_level2_colors <- c(
  # B cells -- x7 (violet)
  'B_Transitional' = "#563C5C",
  'B_Naive' = "#6C3082",
  'B_Naive_activated' = "#4B0082",
  'B_Memory_ITGAX' = "#880085",
  'B_Memory_switched' = "#8B004B",
  'B_Memory_unswitched' = "#B284BE",
  'B_IFNresponder' = "#856088",

  # Plasma cells -- x3+1 (pink)
  'Plasma_IGHA' = "#FF91A4",
  'Plasma_IGHG' = "#FC0FC0",
  'Plasma_XBP1' = "#F4C2C2",
  'Plasma_Proliferative' = "#FF004F",
  
  # pDC -- x1 (yellow)
  'pDC' = "#FFFF00",

  # DC -- x6+1 (orange)
  'cDC1' = "#FF9966",
  'cDC2' = "#ED9121",
  'cDC3' = "#FF4F00",
  'DC4' = "#FF8C00",
  'DC5' = "#C46210",
  'DC_CCR7' = "#EDC9AF",
  'DC_Proliferative' = "#C04000",
  
  # Monocytes -- x5 (red)
  'Mono_classical' = "#960018",
  'Mono_nonClassical' = "#E23D28",
  'Mono_inflammatory' = "#FF0800",
  'Mono_IFNresponse' = "#C40234",
  'Mono_regulatory' = "#AB4E52",

  # T cells -- x11 CD4, x7 CD8, x1 (blue + turquoise)
  'T_CD4_Naive' = "#246BCE",
  'T_CD4_CM' = "#6CB4EE",
  'T_CD4_CM_ribo' = "#89CFF0",
  'T_CD4_EM' = "#1F75FE",
  'T_CD4_EMRA' = "#16166B",
  'T_CD4_eff' = "#00416A",
  'Th0' = "#0087BD",
  'Th1' = "#7C9ED9",
  'Th2' = "#00BFFF",
  'Tregs' = "#F0F8FF",
  'Tregs_activated' = "#a4d4ff",
  
  'T_CD8_Naive' = "#0FFFFF",
  'T_CD8_CM' = "#29AB87",
  'T_CD8_CM_stem' = "#20B2AA",
  'T_CD8_EM_CX3CR1high' = "#99FFFF",
  'T_CD8_EM_CX3CR1int' = "#81D8D0",
  'T_CD8_eff_HOBIT' = "#007A74",
  'T_CD8_IFNresponse' = "#3EB489", 
  'T_CD8_Mem_cytotoxic' = "#37a17b",
  'T_CD8_activated' = "#5ec7a1",
  'T_CD8_arrested' = "#5ec76d",

  'T_Proliferative' = "#008B8B",
  
  # UTC --- x4 ()
  'gdT_V1' = "#C0FF00",
  'gdT_V2_Vγ9' = "#E3F988",
  'MAIT' = "#D0F0C0",
  'MAIT_17' = "#effae9",
  
  # ILC cells -- x7 (green)
  'NK_CD16high' = "#7BA05B",
  'NK_CD56dimCD16' = "#00563B",
  'NK_CD56high' = "#ACE1AF",
  'NK_IFN1response' = "#03C03C",
  'NK_adaptive' = "#50C878",
  'NK_lowRibocontent' = "#80FF00",
  'NK_Proliferative' = "#00A86B",
  
  # Progenitors -- x5 (brown)
  'HSC_LMP' = "#81613E",
  'HSC_MEMP' = "#4A412A",
  'HSC_MMP' = "#592720",
  'B_Progenitors' = "#914034",
  'T_Progenitors' = "#b0885b",
  
  # To remove
  'Platelet' = "#808080",
  'HBcell' = "#000000"
)

