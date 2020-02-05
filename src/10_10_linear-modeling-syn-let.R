
set.seed(0)

#### ---- Linear Model 1 ---- ####
# The first model is a simple two-step process. First, the depletion effect is
# regressed on the RNA expression of the target gene. If this fails to explain
# the gene effect (this model is not significant), then an ANOVA and pair-wise
# test is performed over the alleles.


info(logger, "Beginning linear model 1.")


# put the alleles in the correct order of turning into factors
#   if there is WT, put that first
#   all of the other alleles are in alphabetical order
get_allele_factor_levels <- function(als) {
    als <- unique(als)
    if (any(als == "WT")) {
        all_but_wt <- sort(als[als != "WT"])
        return(c("WT", all_but_wt))
    } else {
        return(sort(als))
    }
}


# prepare the data for the first linear model
#   put alleles as factors
#   dummy variable for if the gene is mutated
#   scale the RNA expression to mean of 0 and std. dev. of 1
lm1_prepare_data <- function(tib, hugo_symbol) {
    allele_levels <- get_allele_factor_levels(tib$allele)
    mod_tib <- tib %>%
        ungroup() %>%
        mutate(
            allele = factor(allele, levels = allele_levels),
            is_altered = as.numeric(is_altered),
            rna_scaled = scale(rna_expression)[, 1],
            gene_effect_scaled = scale(gene_effect)[, 1]
        )

    if (any(table(mod_tib$dep_map_id) > 1)) {
        stop(glue("Multiple data points for a cell line in {hugo_symbol}."))
    }

    if (all(is.na(mod_tib$rna_scaled))) { mod_tib$rna_scaled <- 0 }
    return(mod_tib)
}


# linear model on RNA expression
lm_on_rna <- function(data, ...) {
    fit <- lm(gene_effect ~ rna_scaled, data = data)
    return(fit)
}


# does the linear model on RNA have a significant p-value? (returns Boolean)
rna_pvalue_is_significant <- function(pval, cutoff = 0.01) {
    if (is.na(pval)) { return(FALSE) }
    if (pval < cutoff) { return(TRUE) }
    return(FALSE)
}


# ANOVA on the alleles
#   only done if the RNA expression of the gene does not explain its effect
anova_wrapper <- function(data, rna_pvalue, scale_gene_effect = FALSE, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    if (scale_gene_effect) {
        res <- aov(gene_effect_scaled ~ allele, data = data)
    } else {
        res <- aov(gene_effect ~ allele, data = data)
    }
    return(res)
}


# Kruskal-Wallis rank sum test on the alleles
#   only done if the RNA expression of the gene does not explain its effect
kruskal_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    res <- kruskal.test(gene_effect ~ allele, data = data)
    return(res)
}


# conduct a pair-wise comparison on each pair of alleles
#   only done if the RNA expression of the gene does not explain its effect
pairwise_wrapper <- function(data, rna_pvalue, scale_gene_effect = FALSE, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    if (scale_gene_effect) {
        res <- pairwise.t.test(data$gene_effect_scaled,
                               data$allele,
                               p.adjust.method = "BH")
    } else {
        res <- pairwise.t.test(data$gene_effect,
                               data$allele,
                               p.adjust.method = "BH")
    }
    return(res)
}

# conduct the first modeling attempt
model1_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    nest() %>%
    mutate(data = purrr::map2(data, hugo_symbol, lm1_prepare_data)) %>%
    ungroup() %>%
    mutate(
        rna_lm_fit = map(data, lm_on_rna),
        rna_pvalue = map_dbl(rna_lm_fit, ~ tidy(.x)$p.value[2]),
        allele_aov = map2(data, rna_pvalue, anova_wrapper),
        allele_pairwise = map2(data, rna_pvalue, pairwise_wrapper)
    )

info(logger, "Caching results of model 1.")
cache("model1_tib")



#### ---- RNAi Screen modeling ---- ####


# conduct the first modeling attempt
rnai_model1_tib <- rnai_model_data %>%
    group_by(cancer, allele, hugo_symbol) %>%
    filter(n_distinct(dep_map_id) >= 3) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(n_distinct(allele) > 1) %>%
    nest() %>%
    mutate(data = purrr::map2(data, hugo_symbol, lm1_prepare_data)) %>%
    ungroup() %>%
    mutate(
        rna_lm_fit = map(data, lm_on_rna),
        rna_pvalue = map_dbl(rna_lm_fit, ~ tidy(.x)$p.value[2]),
        allele_aov = map2(data, rna_pvalue, anova_wrapper),
        allele_pairwise = map2(data, rna_pvalue, pairwise_wrapper)
    )

info(logger, "Caching results of RNAi model 1.")
cache("rnai_model1_tib")
