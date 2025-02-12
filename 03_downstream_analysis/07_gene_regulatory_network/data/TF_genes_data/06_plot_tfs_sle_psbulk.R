require(dplyr)
require(ggpubr)
require(data.table)

# Inputs
inputpath_perez2022 <- "tmp/STAT1_SP1_perez2022.csv"
inputpath_scgt00 <- "tmp/STAT1_SP1_scgt00.csv"

# Load
perez2022 <- read.csv(inputpath_perez2022, header = TRUE)
scgt00 <- read.csv(inputpath_scgt00, header = TRUE)

# Plotting flare ---------------------
toplot_flare <-
    perez2022 %>%
    mutate(
        Flare = case_when(Flare == "not_F" ~ "notF",
        TRUE ~ as.character(Flare))
    ) %>%
    filter(
        Flare %in% c("notF", "F")
    )

flare_palette <- c(
  "F" = "#e76f51",
  "notF" = "#e9c46a"
)

plot_flare_sp1 <- ggboxplot(toplot_flare, x = "Level2", y = "SP1", fill = "Flare") +
    facet_wrap(~Level1, scales = "free") +
    scale_fill_manual(values = flare_palette) +
    theme_bw() +
    ggtitle(label = "Flare SP1") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

ggsave(
    plot_flare_sp1,
    filename = "flare_sp1.pdf",
    device = "pdf",
    width = 20,
    height = 15,
    path = "figures/lupus/")


plot_flare_stat1 <- ggboxplot(toplot_flare, x = "Level2", y = "STAT1", fill = "Flare") +
    facet_wrap(~Level1, scales = "free") +
    scale_fill_manual(values = flare_palette) +
    theme_bw() +
    ggtitle(label = "Flare STAT1") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

ggsave(
    plot_flare_stat1,
    filename = "flare_stat1.pdf",
    device = "pdf",
    width = 20,
    height = 15,
    path = "figures/lupus/")


# Plotting sledai ---------------------
SLEDAI_palette <- c(
    "High" = "#103783",
    "Low" = "#9bafd9"
)

toplot_sledai_perez <-
    perez2022 %>%
    filter(
        SLEDAI_category %in% c("High", "Low")
    )

toplot_sledai_scgt00 <-
    scgt00 %>%
    filter(
        SLEDAI_category %in% c("High", "Low")
    )

toplot_sledai <- rbind(toplot_sledai_perez, toplot_sledai_scgt00)

plot_sledai_stat1 <- ggboxplot(toplot_sledai, x = "Level2", y = "STAT1", fill = "SLEDAI_category") +
    facet_grid(studyID ~ Level1, scales = "free_x") +
    scale_fill_manual(values = SLEDAI_palette) +
    theme_bw() +
    ggtitle(label = "SLEDAI STAT1") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

ggsave(
    plot_sledai_stat1,
    filename = "sledai_stat1.pdf",
    device = "pdf",
    width = 40,
    height = 25,
    path = "figures/lupus/")


plot_sledai_sp1 <- ggboxplot(toplot_sledai, x = "Level2", y = "SP1", fill = "SLEDAI_category") +
    facet_grid(studyID ~ Level1, scales = "free_x") +
    scale_fill_manual(values = SLEDAI_palette) +
    theme_bw() +
    ggtitle(label = "SLEDAI SP1") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

ggsave(
    plot_sledai_sp1,
    filename = "sledai_sp1.pdf",
    device = "pdf",
    width = 40,
    height = 25,
    path = "figures/lupus/")




# Plotting Response ------------------
response_palette <- c(
    "R" = "#70e000",
    "NR" = "#007200"
)

toplot_response <-
    scgt00 %>%
    filter(
        Responder %in% c("R", "NR")
    )

plot_response_stat1 <- ggboxplot(toplot_response, x = "Level2", y = "STAT1", fill = "Responder") +
    facet_wrap(~ Level1, scales = "free") +
    scale_fill_manual(values = response_palette) +
    theme_bw() +
    ggtitle(label = "Response STAT1") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(
    plot_response_stat1,
    filename = "response_stat1.pdf",
    device = "pdf",
    width = 20,
    height = 15,
    path = "figures/lupus/")

plot_response_sp1 <- ggboxplot(toplot_response, x = "Level2", y = "SP1", fill = "Responder") +
    facet_wrap(~ Level1, scales = "free") +
    scale_fill_manual(values = response_palette) +
    theme_bw() +
    ggtitle(label = "Response SP1") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(
    plot_response_sp1,
    filename = "response_sp1.pdf",
    device = "pdf",
    width = 20,
    height = 15,
    path = "figures/lupus/")
