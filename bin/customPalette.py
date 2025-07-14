############################## TECHNICAL COUNFOUNDERS ##############################

chemistry_palette =  {
    '3_GEX_V3' : "#9d0208", 
    '3_GEX_V2' : "#dc2f02", 
    '5_GEX_V2' : "#f48c06", 
    '5_GEX_V1' : "#ffba08"
}


studyID_palette = {
    'SCGT03': "#8dd3c7",
    'SCGT00': "#b3de69",
    'SCGT02': "#bebada",
    'COMBAT2022': "#fb8072",
    'Ren2021': "#80b1d3",
    'SCGT05': "#fdb462",
    'Zhang2022': "#ffffb3",
    'Wang2020': "#fccde5",
    'Cillo2020': "#d9d9d9",
    'Schafflick2020': "#bc80bd",
    'Liu2021': "#ccebc5",
    'Perez2022': "#ffed6f",
    'SCGT04': "#1f78b4",
    'Reyes2020': "#33a02c",
    'SCGT01': "#e31a1c",
    'Terekhova2023': "#ff7f00",
    '10XGenomics':'#a0522d',
    'Jiang2020': '#5bd6ef',
    'Martin2019': '#4c001f',
    'Mistry2019': '#044761',
    'Palshikar2022': '#936d35',
    'Ramachandran2019': '#057636',
    'SCGT00val': '#85f80a',
    'SCGT06': '#62cbab',
    'Savage2021': '#d4a9b3'  
}

 

############################## CLINICAL COUNFOUNDERS ##############################

sex_palette = {
    'male': "#2a9d8f", 
    'female': "#f4a261",
    'na': '#ffffff'
}

binned_age_palette = {
    '<18': '#ade8f4',
    '18-30': '#90e0ef',
    '31-40': '#48cae4',
    '41-50': '#00b4d8',
    '51-60': '#0096c7',
    '61-70': '#0077b6',
    '71-80': '#023e8a',
    '>80': '#03045e',
    'unknown': '#ffffff'
}

############################## DISEASES ##############################

disease_sortednames = ['RA', "PS", 'PSA', 'CD', 'UC', "SLE", 'MS','sepsis', 'asthma', 'COPD', 'cirrhosis','HIV', 'HBV',"COVID", "flu",'BRCA', 'NPC', 'HNSCC', 'CRC', 'healthy']
disease_sortednames_filtered = ['RA', "PS", 'PSA', 'CD', 'UC', "SLE", 'MS', 'sepsis', 'asthma', 'COPD', 'cirrhosis','HIV', 'HBV',"COVID", "flu", 'BRCA', 'NPC', 'HNSCC', 'CRC']

diseaseGroup_sortednames = ['IMIDs', 'acute_inflammation', 'chronic_inflammation', 'infection', 'solid_tumor', 'healthy']

############################## DISEASES CATEGORIES ##############################

diseaseCategories = {
    'healthy':'healthy',
    'RA':'IMIDs',
    'CD':'IMIDs',
    'UC':'IMIDs',
    'PS':'IMIDs',
    'SLE':'IMIDs',
    'PSA':'IMIDs',
    'MS':'IMIDs',
    'asthma':'chronic_inflammation',
    'COPD':'chronic_inflammation',
    'cirrhosis':'chronic_inflammation',
    'sepsis':'acute_inflammation',
    'HIV':'infection',
    'HBV':'infection',
    'COVID':'infection',
    'flu':'infection',
    'BRCA':'solid_tumor',
    'NPC':'solid_tumor',
    'HNSCC':'solid_tumor',
    'CRC':'solid_tumor'
}


############################## PALETTE ##############################

###### DISEASES

diseases_palette = {
  'healthy': "#808080",
    
  'RA': '#264653',
  'PS': '#287271',
  'PSA': '#2a9d8f',
  'CD': '#e76f51',
  'UC': '#e9c46a',
  'SLE': '#941c2f', 
  'MS': '#8ab17d',
    
  'asthma': '#ea698b',
  'COPD': '#c05299',
  'cirrhosis': '#973aa8',
    
  'sepsis': '#ef233c',
    
  'HIV': '#e7ecef',
  'HBV': '#a3cef1',
  'COVID': '#6096ba', 
  'flu': '#274c77', 
    
  'BRCA': '#fff75e',
  'NPC': '#fdb833',
  'HNSCC': '#d9981a',
  'CRC': '#9e7524'
}

###### DISEASE GROUPS
diseasesGroup_palette = {
  'IMIDs' : '#2a9d8f',
  'solid_tumor' : '#e3a52d',
  'chronic_inflammation' : '#ffafcc',
  'acute_inflammation' : '#ef233c',
  'infection' : '#abc4ff',
  'healthy' : '#808080'
}


COVID_severity_palette = {
    'healthy': "#808080",
    'COVID_MILD' : "#FFD6FF",
    'COVID_SEV' : "#C8B6FF",
    'COVID_CRIT' : "#E7C6FF",
    'sepsis': '#ef233c', 
    'flu': '#274c77'
}

###### Level 1
annotation_Level1_palette = {
      'B': '#7bc6d6',
      'Plasma': '#025566',
        
      'pDC': '#a7c957',
      'DC': '#6a994e',
      'Mono': '#386641',
        
      'T_CD4_Naive': '#fff3b0',
      'T_CD4_NonNaive': '#e09f3e',
      'T_CD8_Naive': '#9e2a2b',
      'T_CD8_NonNaive': '#540b0e',
        
      'UTC': '#88657f',
      'ILC': '#67253a',
        
      'Cycling_cells': '#d4a373',
      'Progenitors': '#ccd5ae',
    
    'Platelets': '#808080',  # To remove
    'RBC': '#000000'         # To remove
}


############################## OTHER ##############################

disease_marker = {
    'RA': '^',
 'CD': '^',
 'UC': '^',
 'healthy': 'o',
 'PS': '^',
 'SLE': '^',
 'PSA': '^',
 'asthma': 's',
 'COPD': 's',
 'BRCA': '*',
 'sepsis': 'D',
 'HNSCC': '*',
 'cirrhosis': 's',
 'HBV': 'h',
 'flu': 'h',
 'COVID': 'h',
 'MS': '^',
 'NPC': '*',
 'HIV': 'h',
 'CRC': '*'
}

diseaseGroup_marker = dict({
    'healthy': "o",
      'infection': "X",
      'chronic_inflammation': "^",
      'acute_inflammation': "s",
      'solid_tumor': "P",
      'IMIDs': "D",
})




############################## ANNOTATION LEVEL 1 ##############################

annotation_Level1_LowQFilt_sortednames = ["B", "Plasma", "pDC", "DC", "Mono", "T_CD4_Naive", "T_CD4_NonNaive", "T_CD8_Naive", "T_CD8_NonNaive", "UTC", "ILC",
                                          "Cycling_cells", "Progenitors", "Platelets", "RBC"]
annotation_Level1_NonUsedFilt_sortednames = ["B", "Plasma", "pDC", "DC", "Mono", "T_CD4_Naive", "T_CD4_NonNaive", "T_CD8_Naive", "T_CD8_NonNaive", "UTC", "ILC",
                                             "Cycling_cells", "Progenitors"]

annotation_Level1_NonUsedFilt_palette = {
    'B': '#6491c9',    
    'Plasma': '#fb7afa', 
    'pDC': '#844533',  
    'DC': '#7acbad',   
    'Mono': '#59c120',
    'T_CD4_Naive': '#ce1c5e',   
    'T_CD4_NonNaive': '#3225c0',  
    'T_CD8_Naive': '#652d88',  
    'T_CD8_NonNaive': '#326800', 
    'UTC': '#d48080',  
    'ILC': '#d65733',
    'Cycling_cells':'#c9aa54',   
    'Progenitors': '#b171b3'
}


############################## ANNOTATION LEVEL 2 ##############################

annotation_Level2_LowQFilt_sortednames = [
    'B_Transitional', 'B_Naive', 'B_Naive_activated', 'B_Memory_ITGAX', 'B_Memory_unswitched', 'B_Memory_switched','B_IFNresponder',
    'Plasma_XBP1', 'Plasma_IGHA', 'Plasma_IGHG', 
    'pDC',
    'cDC1', 'cDC2', 'cDC3', 'DC4', 'DC5', 'DC_CCR7',
    'Mono_classical', 'Mono_nonClassical', 'Mono_inflammatory', 'Mono_IFNresponse', 'Mono_regulatory',
    'T_CD4_Naive',
    'T_CD4_CM', 'T_CD4_CM_ribo', 'T_CD4_EM', 'T_CD4_EMRA', 'T_CD4_eff', 'Th0','Th1', 'Th2', 'Tregs', 'Tregs_activated',
    'T_CD8_Naive',
    'T_CD8_CM', 'T_CD8_CM_stem', 'T_CD8_EM_CX3CR1high', 'T_CD8_EM_CX3CR1int',  'T_CD8_eff_HOBIT', 'T_CD8_IFNresponse', 'T_CD8_Mem_cytotoxic', 'T_CD8_activated', 'T_CD8_arrested',
    'gdT_V1', 'gdT_V2_Vγ9',  'MAIT', 'MAIT_17',
    'NK_CD16high', 'NK_CD56dimCD16', 'NK_CD56high', 'NK_IFN1response', 'NK_adaptive', 'NK_lowRibocontent',
    'HSC_LMP', 'HSC_MEMP', 'HSC_MMP','B_Progenitors',  'T_Progenitors',
    'Plasma_Proliferative', 'DC_Proliferative', 'T_Proliferative', 'NK_Proliferative', 
    'RBC',
    'Platelets'
]

annotation_Level2_NonUsedFilt_sortednames = [
    'B_Transitional', 'B_Naive', 'B_Naive_activated', 'B_Memory_ITGAX', 'B_Memory_unswitched', 'B_Memory_switched','B_IFNresponder',
    'Plasma_XBP1', 'Plasma_IGHA', 'Plasma_IGHG', 
    'pDC',
    'cDC1', 'cDC2', 'cDC3', 'DC4', 'DC5', 'DC_CCR7',
    'Mono_classical', 'Mono_nonClassical', 'Mono_inflammatory', 'Mono_IFNresponse', 'Mono_regulatory',
    'T_CD4_Naive',
    'T_CD4_CM', 'T_CD4_CM_ribo', 'T_CD4_EM', 'T_CD4_EMRA', 'T_CD4_eff', 'Th0','Th1', 'Th2', 'Tregs', 'Tregs_activated',
    'T_CD8_Naive',
    'T_CD8_CM', 'T_CD8_CM_stem', 'T_CD8_EM_CX3CR1high', 'T_CD8_EM_CX3CR1int',  'T_CD8_eff_HOBIT', 'T_CD8_IFNresponse', 'T_CD8_Mem_cytotoxic', 'T_CD8_activated', 'T_CD8_arrested',
    'gdT_V1', 'gdT_V2_Vγ9',  'MAIT', 'MAIT_17',
    'NK_CD16high', 'NK_CD56dimCD16', 'NK_CD56high', 'NK_IFN1response', 'NK_adaptive', 'NK_lowRibocontent',
    'HSC_LMP', 'HSC_MEMP', 'HSC_MMP','B_Progenitors',  'T_Progenitors',
    'Plasma_Proliferative', 'DC_Proliferative', 'T_Proliferative', 'NK_Proliferative'
]


annotation_Level2_palette = {
    
    # B cells -- x7 (violet)
    'B_Transitional' : "#563C5C",
    'B_Naive' : "#6C3082",
    'B_Naive_activated' : "#4B0082",
    'B_Memory_ITGAX' : "#880085",
    'B_Memory_switched' : "#8B004B",
    'B_Memory_unswitched' : "#B284BE",
    'B_IFNresponder' : "#856088",

    # Plasma cells -- x3+1 (pink)
    'Plasma_IGHA' : "#FF91A4",
    'Plasma_IGHG' : "#FC0FC0",
    'Plasma_XBP1' : "#F4C2C2",
    'Plasma_Proliferative' : "#FF004F",
    
    # pDC -- x1 (yellow)
    'pDC' : "#FFFF00",

    # DC -- x6+1 (orange)
    'cDC1' : "#FF9966",
    'cDC2' : "#ED9121",
    'cDC3' : "#FF4F00",
    'DC4' : "#FF8C00",
    'DC5' : "#C46210",
    'DC_CCR7' : "#EDC9AF",
    'DC_Proliferative' : "#C04000",
    
    # Monocytes -- x5 (red)
    'Mono_classical' : "#960018",
    'Mono_nonClassical' : "#E23D28",
    'Mono_inflammatory' : "#FF0800",
    'Mono_IFNresponse' : "#C40234",
    'Mono_regulatory' : "#AB4E52",

    # T cells -- x11 CD4, x7 CD8, x1 (blue + turquoise)
    'T_CD4_Naive' : "#246BCE",
    'T_CD4_CM' : "#6CB4EE",
    'T_CD4_CM_ribo' : "#89CFF0",
    'T_CD4_EM' : "#1F75FE",
    'T_CD4_EMRA' : "#16166B",
    'T_CD4_eff' : "#00416A",
    'Th0' : "#0087BD",
    'Th1' : "#7C9ED9",
    'Th2' : "#00BFFF",
    'Tregs' : "#F0F8FF",
    'Tregs_activated' : "#a4d4ff",
    
    'T_CD8_Naive' : "#0FFFFF",
    'T_CD8_CM' : "#29AB87",
    'T_CD8_CM_stem' : "#20B2AA",
    'T_CD8_EM_CX3CR1high' : "#99FFFF",
    'T_CD8_EM_CX3CR1int' : "#81D8D0",
    'T_CD8_eff_HOBIT' : "#007A74",
    'T_CD8_IFNresponse' : "#3EB489", 
    'T_CD8_Mem_cytotoxic' : "#37a17b",
    'T_CD8_activated': "#5ec7a1",
    'T_CD8_arrested': "#5ec76d",

    'T_Proliferative' : "#008B8B",
    
    # UTC --- x4 ()
    'gdT_V1' : "#C0FF00",
    'gdT_V2_Vγ9' : "#E3F988",
    'MAIT' : "#D0F0C0",
    'MAIT_17' : "#effae9",
    
    # ILC cells -- x7 (green)
    'NK_CD16high' : "#7BA05B",
    'NK_CD56dimCD16' : "#00563B",
    'NK_CD56high' : "#ACE1AF",
    'NK_IFN1response' : "#03C03C",
    'NK_adaptive' : "#50C878",
    'NK_lowRibocontent' : "#80FF00",
    'NK_Proliferative' : "#00A86B",
    
    # Progenitors -- x5 (brown)
    'HSC_LMP' : "#81613E",
    'HSC_MEMP' : "#4A412A",
    'HSC_MMP' : "#592720",
    'B_Progenitors' : "#914034",
    'T_Progenitors': "#b0885b",
    
    # To remove
    'Platelet' : "#808080",
    'HBcell' : "#000000"
}



annotation_Level1Level2_NonUsedFilt_sortednames_dict = {
    "B" : ['B_Transitional', 'B_Naive', 'B_Naive_activated', 'B_Memory_ITGAX', 'B_Memory_unswitched', 'B_Memory_switched','B_IFNresponder'],
    "Plasma": ['Plasma_XBP1', 'Plasma_IGHA', 'Plasma_IGHG'], 
    "pDC" : ["pDC"],
     "DC" : ['cDC1', 'cDC2', 'cDC3', 'DC4', 'DC5', 'DC_CCR7'],
    "Mono" : ['Mono_classical', 'Mono_nonClassical', 'Mono_inflammatory', 'Mono_IFNresponse', 'Mono_regulatory'],
    'T_CD4_Naive' : ['T_CD4_Naive'],
    "T_CD4_NonNaive": ['T_CD4_CM', 'T_CD4_CM_ribo', 'T_CD4_EM', 'T_CD4_EMRA', 'T_CD4_eff', 'Th0','Th1', 'Th2', 'Tregs', 'Tregs_activated'],
    'T_CD8_Naive': ['T_CD8_Naive'],
    "T_CD8_NonNaive" : ['T_CD8_CM', 'T_CD8_CM_stem', 'T_CD8_EM_CX3CR1high', 'T_CD8_EM_CX3CR1int',  'T_CD8_eff_HOBIT', 'T_CD8_IFNresponse', 'T_CD8_Mem_cytotoxic', 'T_CD8_activated', 'T_CD8_arrested'],
    "UTC": ['gdT_V1', 'gdT_V2_Vγ9',  'MAIT', 'MAIT_17'],
    "ILC" : ['NK_CD16high', 'NK_CD56dimCD16', 'NK_CD56high', 'NK_IFN1response', 'NK_adaptive', 'NK_lowRibocontent'],
    "Progenitors": ['HSC_LMP', 'HSC_MEMP', 'HSC_MMP','B_Progenitors',  'T_Progenitors'],
    "Cycling_cells": ['Plasma_Proliferative', 'DC_Proliferative', 'T_Proliferative', 'NK_Proliferative']
}


############################## GENE MARKERS at LEVEL 1 & 2 ##############################

geneMarkers_level1 = {
    "B":  ["MS4A1", "CD19", "CD79A", "CD79B"],
    "Plasma" :  ["MZB1", "JCHAIN", "XBP1", "IRF4", "PRDM1", "FKBP11"],
    'pDC': ["IL3RA", "IRF7", "LILRA4", "IRF8", "GZMB"],
    'DC': ["CD1C", "FCER1A", "CLEC9A", "SERPINA1"],
    'Mono': ["CD14", "FCGR3A", "LYZ", "VCAN"],
    'T' : ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"],
    'T Naive': ["SELL", "LEF1", "CCR7", "NOSIP", "CD27"],
    'T NonNaive': ["IL7R", "PASK", "CCL5"],
    'T Unconventional': ["KLRB1"],
    'ILC': ["NCAM1", "FCGR3A"],
    'Proliferation': ["MKI67", "TOP2A"],
    'Progenitors': ["CD34", "SOX4"],
    'Platelets': ["PPBP"],  # To remove
    'RBC': ["HBA1", "HBB"]         # To remove
}

geneMarkers_level2 = {

    "B" : {
        'B' : ["MS4A1", "CD19", "CD79A", "CD79B"],
       'B_Transitional' : ["VPREB3"],
       'B_Naive' : ["IGHD", "IGHM", "SELL"],
       'B_Naive_activated' : ["CD69", "CD86"],
       'B_Memory_ITGAX' : ["ITGAX", "ITGB2", "NEAT1", "FCRL5"],
       'B_Memory_unswitched' : ["CD37", "CD40", "MIF"],
       'B_Memory_switched' : ["CD27", "ITGB1", "S100A10", "GPR183"],
       'B_IFNresponder' : ["IFITM1", "IFITM2", "IFNGR2", "ISG15"],
    },

    "Plasma" : {
        "Plasma" :  ["MZB1", "JCHAIN", "IRF4", "PRDM1", "FKBP11"],
       'Plasma_XBP1' : ["CD38", "CXCR3", "XBP1"],
       'Plasma_IGHA' : ["IGHA1", "IGHA2", "IGLC1"],
       'Plasma_IGHG' : ["IGHG1", "IGHG2", "IGHG3", "IGHG4"],
       },
        
        'pDC' : ["IL3RA", "IRF7", "LILRA4", "IRF8", "GZMB"],

    "DC" : {      
       'cDC1' : ["CLEC9A", "C1orf54", "IDO1", "CLNK"],
       'cDC2' : ["HLA-DPB1", "CD1C", "CLEC10A", "FCER1A", "FCGR2B"],
       'cDC3' : ["VCAN", "S100A9", "S100A8", "LYZ"],
       'DC4' : ["SERPINA1", "LILRA2", "LILRB2", "CSF1R"],
       'DC5' : ["AXL", "SIGLEC6", "TMSB4X", "LST1"],
       'DC_CCR7' : ["CCR7", "ITGAX", "CXCL16"]
    },

    "Mono" : {        
       'Mono_classical' : ["CD14", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1"],
       'Mono_nonClassical' : ["FCGR3A", "LST1", "IFITM2", "NAP1L1"],
       'Mono_inflammatory' : ["IL1B", "CXCL8", "G0S2"], #"DSUSP6", 
       'Mono_IFNresponse' : ["ISG15", "IFI6", "LY6E", "IFI44L"],
       'Mono_regulatory' : ["PHLDA1", "THBD", "THBS1", "HAVCR2"]
    },

   
    'T_CD4_NonNaive' : {
        'T' : ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"],
       'T_CD4_CM' : ["PASK", "IL7R", "CCR7", "CD28"],
       'T_CD4_CM_ribo' : ["RPL32", "RPL30", "PASK"],
       'T_CD4_EM' : ["GZMK", "KLRB1", "CCL5"],
       'T_CD4_EMRA' : ["GZMB", "CCL4", "GZMA", "PRF1", "CST7"],
       'T_CD4_eff' : ["ITGB1", "PFN1", "ITM2A", "TCF7"],
       'Th0' : ["FOS", "JUNB", "DUSP1", "ANXA1"],
       'Th1' : ["ISG15", "STAT1", "IFI6"], #, "ANXA3"
       'Th2' : ["GATA3", "IL4R", "TNFSF10"], #"TNFRSF42",
       'Tregs' : ["FOXP3", "CTLA4", "TIGIT", "IL2RA"],
       'Tregs_activated' : ["HLA-DPA1", "HLA-DRA", "FOXP3"]
    },

    'T_CD8_NonNaive' : {
        'T' : ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"],
       'T_CD8_CM' : ["LTB", "TCF7", "CD27"], # "IL7R2", 
       'T_CD8_CM_stem' : ["IL6ST",  "SELL", "IL7R"],
       'T_CD8_EM_CX3CR1high' : ["CX3CR1","GNLY", "GZMH", "FGFBP2", "NKG7"],
       'T_CD8_EM_CX3CR1int' : ["KLRF1", "FCGR3A", "KLRC2"],
       'T_CD8_eff_HOBIT' : ["ZNF683", "TYROBP", "KLRC3"],
       'T_CD8_IFNresponse' : ["IFIT2", "IFIT3", "OASL", "ISG15"],
       'T_CD8_Mem_cytotoxic' : ["KLRB1", "GZMK", "GZMA", "KLRG1"],
       'T_CD8_activated': ["HLA-DRA", "HLA-DRB1", "HLA-DMA"],
       'T_CD8_arrested': ["CXCR6", "DUSP5", "NR4A3", "DUSP4"]
    },

    'UTC' : {
        'T' : ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"],
       'gdT_V1' : ["TCF7", "TRDC", "TRGC1", "TRGC2"],
       'gdT_V2_Vγ9' : ["TRDV2", "FCGR3A", "GZMB", "TRGV9"],
       'MAIT' : ["KLRB1", "LTB", "GZMK", "IL7R"],
       'MAIT_17': ["CCR6", "IL27RA", "RORC", "STAT3"]
    },
    
    'ILC' : {
       'NK_CD16high' : ["FCGR3A", "CX3CR1", "S100A4"],
       'NK_CD56dimCD16' : ["PRF1", "GZMB"],
       'NK_CD56high' : ["NCAM1", "XCL1", "XCL2", "IL18R1"],
       'NK_IFN1response' : ["IFIT2", "IFIT3", "OASL", "ISG15"],
       'NK_adaptive' : ["GZMK", "IL32", "CD3E", "IL2RB"],
       'NK_lowRibocontent' : ["MT-ND1", "MT-CYB", "RPS2"],
    },

    'Progenitors' : {
        'Progenitors': ["CD34", "SOX4",  "KIT"],
        'HSC_LMP' : ["CSF3R", "MDK", "BEX1", "SELL"],
        'HSC_MEMP' : ["PRG2", "CLC", "HDC", "MS4A2"],
        'HSC_MMP' : ["GATA1", "KLF1", "CNRIP1", "TFR2"],
        'B_Progenitors' : ["TNFRSF18", "JCHAIN", "HCST"], 
        'T_Progenitors' : ["CD3G", "IL7R", "BCL11B", "IL32"]
    },
    
    'Cycling_cells': {
        'Proliferation': ['MKI67', 'TOP2A'], 
        'Plasma': ['MZB1', 'JCHAIN', 'XBP1', 'IRF4', 'PRDM1', 'FKBP11'],
         'DC': ['CD1C', 'FCER1A', 'CLEC9A', 'SERPINA1'],
        'T': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B'],
         'ILC': ['NCAM1', 'FCGR3A']
    },
    
    'Platelets': ["PPBP"],  # To remove
    'RBC': ["HBA1", "HBB"]  # To remove
}

##### SELECTED CELL TYPES #################

shap_cell_types = ["Mono", "T_CD4_Naive", "T_CD4_NonNaive", "T_CD8_Naive", "T_CD8_NonNaive", "B", "Plasma", "UTC", "ILC", "pDC", "DC"]
