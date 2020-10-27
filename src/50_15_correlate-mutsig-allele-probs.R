# Check for correlation between signature contribution and KRAS allele prob.

GRAPHS_DIR <- "50_15_correlate-mutsig-allele-probs"
reset_graph_directory(GRAPHS_DIR)

theme_set(theme_bw(base_size = 7, base_family = "Arial"))

library(Hmisc)

mutsig_df_compact <- mutsig_noartifact_df %.% {
  select(cancer, tumor_sample_barcode, description, contribution)
  group_by(cancer, tumor_sample_barcode, description)
  summarise(contribution = sum(contribution))
  group_by(cancer, tumor_sample_barcode)
  nest()
  ungroup()
  rename(mutsig_data = data)
}


allele_signature_corr_data <- ranked_allele_predictions %.%
  {
    filter(is_tested)
    group_by(cancer, tumor_sample_barcode)
    nest()
    ungroup()
    rename(alleleprob_data = data)
    left_join(mutsig_df_compact, by = c("cancer", "tumor_sample_barcode"))
    unnest(alleleprob_data)
    unnest(mutsig_data)
    group_by(cancer, kras_allele, description)
    nest()
    ungroup()
  }

d <- allele_signature_corr_data %.%
  {
    filter(kras_allele == "G12C" & description == "4" & cancer == "LUAD")
    select(data)
    unnest(data)
    mutate(
      real_kras_allele = factor_alleles(real_kras_allele),
    )
  }


fit_me <- lme4::lmer(
  allele_prob ~ 1 + (1 + contribution | real_kras_allele),
  data = d
)

fixef_intercept <- lme4::fixef(fit_me)[["(Intercept)"]]
me_ranef <- lme4::ranef(fit_me)$real_kras_allele %.% {
  rownames_to_column(var = "real_kras_allele")
  as_tibble()
  janitor::clean_names()
  mutate(adj_intercept = intercept + fixef_intercept)
}

eg_p <- d %>%
  ggplot(aes(contribution, allele_prob)) +
  geom_point(
    aes(color = real_kras_allele),
    size = 1,
    alpha = 0.6
  ) +
  geom_abline(
    aes(
      slope = contribution,
      intercept = adj_intercept,
      color = real_kras_allele
    ),
    data = me_ranef,
    alpha = 0.7,
    size = 0.5
  ) +
  scale_color_manual(
    values = short_allele_pal,
    drop = TRUE
  ) +
  theme(
    legend.position = "right",
    legend.background = element_blank()
  ) +
  labs(
    x = "mutational signature 4 contribution",
    y = "probability of G12C",
    color = "observed\nKRAS allele",
    title = "Correlation between signature 4 and KRAS G12C in LUAD",
    subtitle = "Lines represent a multilevel linear model fit with a varying intercept and slope\nfor each observed KRAS allele"
  )
ggsave_wrapper(
  eg_p,
  plot_path(GRAPHS_DIR, "luad_g12c-prob_sig4_me-model.svg"),
  "small"
)



#### ---- Probability spline models ---- ####


logisitc_spline_signature_allele <- function(cancer,
                                             kras_allele,
                                             description,
                                             data,
                                             n_knots = 3) {
  model_data <- data %>%
    transmute(
      is_allele = as.numeric(real_kras_allele == !!kras_allele),
      contribution
    )

  spline_fit <- try(
    spline_fit <- rcspline.plot(
      x = model_data$contribution,
      y = model_data$is_allele,
      model = "logistic",
      nk = n_knots,
      show = "prob",
      noprint = TRUE,
      statloc = "none"
    ),
    silent = TRUE
  )

  if (inherits(spline_fit, "try-error")) {
    message(glue("Skipping {kras_allele} & sig. {description} in {cancer}."))
    return(NULL)
  }

  spline_fit_data <- tibble(
    x = spline_fit$x,
    y = spline_fit$xbeta[, 1],
    y_lower = spline_fit$lower[, 1],
    y_upper = spline_fit$upper[, 1]
  )

  spline_knots <- tibble(knots = spline_fit$knots)

  p <- spline_fit_data %>%
    ggplot(aes(x = x, y = y)) +
    geom_ribbon(
      aes(ymin = y_lower, ymax = y_upper),
      fill = "grey50",
      alpha = 0.3
    ) +
    geom_line(group = "a") +
    geom_point(
      aes(x = contribution, y = is_allele),
      data = model_data,
      shape = 3,
      size = 1,
      alpha = 0.5,
      color = "black"
    ) +
    geom_vline(
      aes(xintercept = knots),
      data = spline_knots,
      lty = 2,
      alpha = 0.75,
      color = "blue"
    ) +
    labs(
      x = glue("signature {description} contribution"),
      y = glue("probability of {kras_allele}"),
      title = glue("{kras_allele} & sig {description} in {cancer}")
    )
  fn <- as.character(glue(
    "prob-splie_nk{n_knots}_{cancer}_{kras_allele}_sig{description}.svg"
  ))
  ggsave_wrapper(
    p,
    plot_path(GRAPHS_DIR, fn),
    "small"
  )
}

allele_signature_corr_data %>%
  pwalk(logisitc_spline_signature_allele)


