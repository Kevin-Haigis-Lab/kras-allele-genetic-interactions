# Build Figure 1

FIGNUM <- 1
VERSION <- 1
FIGFILENAME <- glue("figure_{FIGNUM}_{VERSION}.svg")
FIG_DIMENSIONS <- get_figure_dimensions(2, "medium")


library(patchwork)

#### ---- A. KRAS mutation lollipop ---- ####

# Panel A.
# The distribution of mutations along the KRAS amino acid sequence.
# original script: "src/90_05_kras-allele-distribution.R"

panel_A <- read_fig_proto("lollipop-kras_2", FIGNUM)


#### ---- B. KRAS allele frequency ---- ####

# Panel B.
# Barplots of the KRAS allele frequency across the cancers.
# original script: "src/90_05_kras-allele-distribution.R"

panel_B <- read_fig_proto("allele_dist_barplot_stackplot", FIGNUM)


#### ---- C. Mutational signatures probability of causing KRAS allele ---- ####

# Panel C.
# The probability that each allele was caused by each detectable mutational
# signature.
# original script: "src/50_30_mutsignatures_prob-causing-allele.R"

panel_C <- read_fig_proto("probability-mutsig-caused-allele", FIGNUM)


#### ---- D. Probability of select signatures causing KRAS allele ---- ####

# Panel D.
# The distribution of the probabilities that alleles were cause by a selection
# of the mutational signatures.
# original script: "src/50_30_mutsignatures_prob-causing-allele.R"

panel_D <- read_fig_proto("contribution-of-select-signatures", FIGNUM)


#### ---- E. Predicting the alleles from the mutational signatures ---- ####

# Panel E.
# The observed vs. predicted KRAS allele frequencies.
# original script: "src/50_10_observed-predicted-kras-alleles.R"

panel_E_proto_list <- read_fig_proto("obs_pred_plot_stats", FIGNUM)
panel_E <- wrap_plots(panel_E_proto_list, nrow = 1, guides = "collect")




#### ---- F. Predicting the alleles from the mutational signatures ---- ####

# Panel F.
# The observed vs. predicted KRAS allele frequencies.
# original script: "src/50_10_observed-predicted-kras-alleles.R"

panel_F_proto_list <- read_fig_proto("obs_pred_plot_g12_stats", FIGNUM)
panel_F <- wrap_plots(panel_F_proto_list, nrow = 1, guides = "collect")




#### ---- Figure assembly ---- ####

# ROW 1
panel_B_mod <- wrap_plots(panel_B, ncol = 2)
row_1 <- panel_A + panel_B_mod + plot_layout(ncol = 2, widths = c(1, 2, 2))

# ROW 2
row_2 <- (panel_C - panel_D) + plot_layout(widths = c(2, 1))

# COMPLETE FIGURE
full_figure <- row_1 / row_2 / panel_E / panel_F

save_figure(
    full_figure,
    figure_num = 1,
    dim = FIG_DIMENSIONS
)
