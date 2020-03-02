# A script to make simple plot of WDR26 expression vs. gene effect in COAD.

GRAPHS_DIR <- "10_57_coad_depmap_wdr26"
reset_graph_directory(GRAPHS_DIR)


GENE <- "WDR26"

dat <- model1_tib %>%
    filter(cancer == "COAD" & hugo_symbol == !!GENE) %>%
    select(hugo_symbol, cancer, data) %>%
    unnest(data)

fit <- lm(gene_effect ~ rna_expression, data = dat)
p_val <- glance(fit)$p.value
m <- coef(fit)["rna_expression"]
b <- coef(fit)["(Intercept)"]
R2 <- glance(fit)$r.squared

equation_dat <- tibble(
    x = c(5.6),
    y = c(-0.26, -0.29, -0.32),
    label = c(glue("*y* = {round(m, 3)}*x* + {round(b, 3)}"),
              glue("p-value: {round(p_val, 2)}"),
              glue("*R*<sup>2</sup> = {round(R2, 2)}"))
)

p <- dat %>%
    mutate(allele = factor_alleles(allele),
           copy_number_label = factor(copy_number_label,
                                      levels = c("norm", "amp"))) %>%
    ggplot(aes(x = rna_expression, y = gene_effect)) +
    geom_point(aes(color = allele, size = copy_number_label), alpha = 0.7) +
    geom_smooth(method = "lm", color = "grey20", linetype = 2, size = 0.6) +
    geom_richtext(aes(x = x, y = y, label = label),
                  data = equation_dat,
                  size = 2.6, family = "arial", hjust = 1.0,
                  fill = NA, label.color = NA,
                  label.padding = unit(rep(0, 4), "pt")) +
    scale_color_manual(values = short_allele_pal) +
    scale_size_manual(values = c("norm" = 1.3, "amp" = 2)) +
    theme_bw(base_size = 7, base_family = "arial") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_markdown()
    ) +
    labs(
        x = "mRNA expression (*log* TPM)",
        y = "dependency score",
        title = GENE,
        size = "copy number",
        color = "KRAS allele"
    )
ggsave_wrapper(
    p,
    plot_path(GRAPHS_DIR, glue("{GENE}-rna-dep.svg")),
    "small"
)

saveRDS(
    p,
    get_fig_proto_path("coad_depmap_wdr26-rna-v-dep.svg", 15, supp = TRUE)
)
