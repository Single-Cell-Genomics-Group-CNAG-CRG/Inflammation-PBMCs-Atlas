require(dplyr)
require(ggpubr)
require(data.table)

# Inputs
inputpath_perez2022 <- "tmp/STAT1_SP1_perez2022_level1.csv"
inputpath_scgt00 <- "tmp/STAT1_SP1_scgt00_level1.csv"

inputpath_perez2022_level2 <- "tmp/STAT1_SP1_perez2022_level2.csv"
inputpath_scgt00_level2 <- "tmp/STAT1_SP1_scgt00_level2.csv"


# Load
perez2022_level1 <- read.csv(inputpath_perez2022, header = TRUE)
scgt00_level1 <- read.csv(inputpath_scgt00, header = TRUE)

perez2022_level2 <- read.csv(inputpath_perez2022_level2, header = TRUE)
scgt00_level2 <- read.csv(inputpath_scgt00_level2, header = TRUE)

# Flare ---------------------
flare_palette <- c(
  "F" = "#e76f51",
  "notF" = "#e9c46a"
)

data_flare_level1 <-
    perez2022_level1 %>%
    mutate(
        Flare = case_when(Flare == "not_F" ~ "notF",
        TRUE ~ as.character(Flare))
    ) %>%
    filter(
        Flare %in% c("notF", "F")
    ) %>%
    tidyr::pivot_longer(
        cols = c(STAT1, SP1),
        names_to = "tfs",
        values_to = "activity"
        ) %>%
    filter(
        !is.na(activity)
    )

pvals_flare <- data_flare_level1 %>%
  group_by(Level1, tfs) %>%
  filter(n_distinct(Flare) == 2) %>% # Ensure that Flare has two levels
  filter(all(table(Flare) >= 4)) %>% # Ensure that each level has at least two observations
  summarise(
    t_test = list(broom::tidy(t.test(activity ~ Flare))),
    .groups = "drop"
  ) %>%
  tidyr::unnest(t_test) %>%
  select(Level1, tfs, p.value) %>%
  mutate(
    pval_adj = p.adjust(p.value, method = "BH")
  ) %>%
  filter(
    pval_adj < 0.05
  )

data_flare_level2 <-
    perez2022_level2 %>%
    mutate(
        Flare = case_when(Flare == "not_F" ~ "notF",
        TRUE ~ as.character(Flare))
    ) %>%
    filter(
        Flare %in% c("notF", "F")
    ) %>%
    tidyr::pivot_longer(
        cols = c(STAT1, SP1),
        names_to = "tfs",
        values_to = "activity"
        ) %>%
    filter(
        !is.na(activity),
        Level1 %in% pvals_flare$Level1
    )

plot_flare <- ggboxplot(data_flare_level2, x = "Level2", y = "activity", fill = "Flare") +
    facet_wrap(~Level1 + tfs, scales = "free") +
    scale_fill_manual(values = flare_palette) +
    theme_bw() +
    ggtitle(label = "Flare") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )

ggsave(
    plot_flare,
    filename = "flare.pdf",
    device = "pdf",
    width = 20,
    height = 15,
    path = "figures/lupus/")



# Plotting sledai ---------------------
SLEDAI_palette <- c(
    "High" = "#103783",
    "Low" = "#9bafd9"
)

data_sledai_perez_level1 <-
    perez2022_level1 %>%
    filter(
        SLEDAI_category %in% c("High", "Low")
    ) %>%
    tidyr::pivot_longer(
        cols = c(STAT1, SP1),
        names_to = "tfs",
        values_to = "activity"
        ) %>%
    filter(
        !is.na(activity)
    )

data_sledai_scgt_level1 <-
    scgt00_level1 %>%
    filter(
        SLEDAI_category %in% c("High", "Low")
    ) %>%
    tidyr::pivot_longer(
        cols = c(STAT1, SP1),
        names_to = "tfs",
        values_to = "activity"
        ) %>%
    filter(
        !is.na(activity)
    )

pvals_sledai_perez <- data_sledai_perez_level1 %>%
  group_by(Level1, tfs) %>%
  filter(n_distinct(SLEDAI_category) == 2) %>% # Ensure that Flare has two levels
  filter(all(table(SLEDAI_category) >= 4)) %>% # Ensure that each level has at least two observations
  summarise(
    t_test = list(broom::tidy(t.test(activity ~ SLEDAI_category))),
    .groups = "drop"
  ) %>%
  tidyr::unnest(t_test) %>%
  select(Level1, tfs, p.value) %>%
  mutate(
    pval_adj = p.adjust(p.value, method = "BH")
  ) %>%
  filter(
    pval_adj < 1
  )

pvals_sledai_scgt <- data_sledai_scgt_level1 %>%
  group_by(Level1, tfs) %>%
  filter(n_distinct(SLEDAI_category) == 2) %>% # Ensure that Flare has two levels
  filter(all(table(SLEDAI_category) >= 4)) %>% # Ensure that each level has at least two observations
  summarise(
    t_test = list(broom::tidy(t.test(activity ~ SLEDAI_category))),
    .groups = "drop"
  ) %>%
  tidyr::unnest(t_test) %>%
  select(Level1, tfs, p.value) %>%
  mutate(
    pval_adj = p.adjust(p.value, method = "BH")
  ) %>%
  filter(
    pval_adj < 1
  )

############# PLOT SLEDAI
# INCLUDE ggboxplot() and ggsave()



# Plotting Response ------------------
response_palette <- c(
    "R" = "#70e000",
    "NR" = "#007200"
)

data_response_level1 <-
    scgt00_level1 %>%
    filter(
        Responder %in% c("R", "NR")
    ) %>%
    tidyr::pivot_longer(
        cols = c(STAT1, SP1),
        names_to = "tfs",
        values_to = "activity"
        ) %>%
    filter(
        !is.na(activity)
    )

data_response_level2 <-
    scgt00_level2 %>%
    filter(
        Responder %in% c("R", "NR")
    ) %>%
    tidyr::pivot_longer(
        cols = c(STAT1, SP1),
        names_to = "tfs",
        values_to = "activity"
        ) %>%
    filter(
        !is.na(activity)
    )

pvals_response <- data_response_level1 %>%
  group_by(Level1, tfs) %>%
  filter(n_distinct(Responder) == 2) %>% # Ensure that Flare has two levels
  filter(all(table(Responder) >= 3)) %>% # Ensure that each level has at least two observations
  summarise(
    t_test = list(broom::tidy(t.test(activity ~ Responder))),
    .groups = "drop"
  ) %>%
  tidyr::unnest(t_test) %>%
  select(Level1, tfs, p.value) %>%
  mutate(
    pval_adj = p.adjust(p.value, method = "fdr")
  ) %>%
  filter(
    p.value < 0.1
  ) %>%
  arrange(p.value)

pvals_response


toplot_response <-
    data_response_level2 %>%
    semi_join(pvals_response, by = c("Level1", "tfs"))

plot_response <- ggboxplot(toplot_response, x = "Level2", y = "activity", fill = "Responder") +
    facet_wrap(~ Level1 + tfs, scales = "free") +
    scale_fill_manual(values = response_palette) +
    theme_bw() +
    ggtitle(label = "Response") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(
    plot_response,
    filename = "response.pdf",
    device = "pdf",
    width = 20,
    height = 15,
    path = "figures/lupus/")
