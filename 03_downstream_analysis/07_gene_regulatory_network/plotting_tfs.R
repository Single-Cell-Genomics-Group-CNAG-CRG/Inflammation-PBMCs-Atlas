require(dplyr)
require(data.table)
require(RColorBrewer)
require(ggplot2)
require(colorRamp2)

# Inputs
inputpath <- "data/tfs/toplot_tfs.rds"
color_palette_path <- "config/colors_palette.R"

# Outputs
outputpath <- ""

# Arguments
spectra_of_interest <- c("SPECTRA_110", "SPECTRA_13", "SPECTRA_130", "SPECTRA_24",
                         "SPECTRA_34", "SPECTRA_46", "SPECTRA_56", "SPECTRA_71",
                         "SPECTRA_86", "SPECTRA_98")

# Load data
toplot <- readRDS(inputpath)
source(color_palette_path)

xx <-
    toplot %>%
    tidyr::pivot_wider(names_from = "annot", values_from = "t") %>%
    tibble::column_to_rownames("target") %>%
    select(contains(spectra_of_interest)) %>%
    filter(!if_all(everything(), is.na))

col_annot <- data.frame(
    tf = as.factor(sapply(strsplit(names(xx), ":"), "[", 2)),
    disease = as.factor(sapply(strsplit(names(xx), ":"), "[", 3)),
#        SPECTRA = as.factor(sapply(strsplit(names(xx), ":"), "[", 1)),
    cell_type = as.factor(sapply(strsplit(names(xx), ":"), "[", 4)),
    row.names = names(xx)
) %>%
arrange(cell_type, tf, disease)

# Colors
my_colour_annot <- list(
    tf = tfs_colors[names(tfs_colors) %in% levels(col_annot$tf)],
    disease = disease_colors[names(disease_colors) %in% levels(col_annot$disease)],
    cell_type = cell_colors[names(cell_colors) %in% levels(col_annot$cell_type)]
)

my_breaks <- c(seq(-10, 0, length.out=ceiling(palette_length/2) + 1),
            seq(0.05, 15, length.out=floor(palette_length/2)))

pdf(outputpath, width = 15, height = 12)
pheatmap::pheatmap(
    xx[, rownames(col_annot)],
    border_color = "gray",
    na_col = "white",
    color = my_color,
    breaks = my_breaks,
    annotation_colors = my_colour_annot,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    cellheight = 10,
    cellwidth = 10,
    annotation_col = col_annot,
    cex = 1,
    fontsize = 9,
    labels_col = col_annot$disease,
    main = "IFN_response",
    gaps_col = c(6, 10, 20, 25, 29, 31, 33, 38)
    )
dev.off()
