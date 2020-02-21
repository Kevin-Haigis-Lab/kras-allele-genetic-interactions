# Survival analysis

library(survminer)
library(survival)


GRAPHS_DIR <- "70_10_survival-analysis"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)



#### ---- Prepare TCGA survival data ---- ####

# KRAS information for all cancer samples.
tumor_sample_data <- cancer_full_muts_df %>%
    group_by(cancer, tumor_sample_barcode) %>%
    slice(1) %>%
    ungroup() %>%
    select(cancer, tumor_sample_barcode, ras_allele) %>%
    unique() %>%
    arrange(cancer) %>%
    mutate(
        kras_allele = str_remove(ras_allele, "KRAS_"),
        kras_mut = as.numeric(ras_allele != "WT")
    ) %>%
    select(cancer, tumor_sample_barcode, kras_mut, kras_allele)


# TCGA survival data with KRAS annotation.
tcga_survival_data_anno <- inner_join(
    tcga_survival_data,
    tumor_sample_data,
    by = c("cancer" = "cancer", "patient_id" = "tumor_sample_barcode")
) %>%
    filter(!is.na(os_months) & !is.na(os_status)) %>%
    mutate(
        time = os_months,
        status = as.numeric(os_status == "DECEASED")
    )


#### ---- Prepare MMRF survival data ---- ####

mmrf_survival_data_anno <- inner_join(
    mmrf_survival_data,
    tumor_sample_data,
    by = c("cancer" = "cancer", "public_id" = "tumor_sample_barcode")
) %>%
    select(cancer, public_id, kras_mut, kras_allele,
           date_of_death, time_to_os, censor_flag_overall_survival,
           overall_survival_censored_date, time_to_os_event_censored) %>%
    mutate(
        time = overall_survival_censored_date,
        status = as.numeric(censor_flag_overall_survival)
    )


#### ---- Combined survival data ---- ####

# Merged survival data.
# Note that a lot of features will be all `NA` for MM.
# Some may be in the original data and can be retrieved during munging.
cancer_survival_df <- bind_rows(
    {
        tcga_survival_data_anno %>%
            dplyr::rename(tumor_sample_barcode = patient_id) %>%
            select(cancer, tumor_sample_barcode,
                   age, sex, weight, ethnicity,
                   path_t_stage, path_m_stage,
                   kras_mut, kras_allele,
                   time, status)
    },
    {
        mmrf_survival_data_anno %>%
            dplyr::rename(tumor_sample_barcode = public_id) %>%
            select(cancer, tumor_sample_barcode,
                   kras_mut, kras_allele,
                   time, status)
    }
)


#### ---- Survival Curve Helpers ---- ####

# Get p-value from results of `survdiff()`.
survdiff_pval <- function(fit) {
    pchisq(fit$chisq, length(fit$n) - 1, lower.tail = FALSE)
}


# Alter a palette's names with the covariate and "=" sign.
alter_pal_for_ggsurvplot <- function(pal, covariate_name) {
    names(pal) <- paste(covariate_name, names(pal), sep = "=")
    return(pal)
}


# Use 'patchwork' to put the curve and table together.
patch_ggsurvplot <- function(ggsurv_obj, layout_heights = c(3, 1)) {
    survival_curve <- ggsurv_obj$plot
    survival_tbl <- ggsurv_obj$table
    p <- (survival_curve / survival_tbl) +
        plot_layout(heights = layout_heights, guides = "collect")
    return(p)
}


# Style the survival curve and the table from `ggsurvplot()`.
style_ggsurvminer_plot <- function(ggsurv_obj,
                                   table_y = "number at risk") {
    x_expand <- c(0.02, 0.01)
    ggsurv_obj$plot <- ggsurv_obj$plot +
        scale_x_continuous(expand = expand_scale(mult = x_expand)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.01))) +
        theme_classic(base_size = 7, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(),
            axis.title.x = element_blank()
        )
    ggsurv_obj$table <- ggsurv_obj$table +
        scale_x_continuous(expand = expand_scale(mult = x_expand)) +
        theme_classic(base_size = 7, base_family = "arial") +
        theme(
            legend.position = "none",
            plot.title = element_blank(),
            axis.ticks = element_blank()
        ) +
        labs(y = table_y)
    return(ggsurv_obj)
}


write_summary <- function(fit_summary, fpath) {
    sink(fpath, append = TRUE)
    print(fit_summary)
    sink()
}


write_survfit_summary <- function(fit, title, fname, clear_file = FALSE) {
    fpath <- table_path(GRAPHS_DIR, fname)

    if (clear_file & file.exists(fpath)) file.remove(fpath)

    cat(title, "\n", file = fpath, append = TRUE)
    write_summary(summary(fit)$table, fpath)
    cat(str_rep("=", 80), "\n\n", file = fpath, append = TRUE)
}


write_survdiff_summary <- function(fit, title, fname, clear_file = FALSE) {
    fpath <- table_path(GRAPHS_DIR, fname)

    if (clear_file & file.exists(fpath)) file.remove(fpath)

    cat(title, "\n", file = fpath, append = TRUE)
    write_summary(fit, fpath)
    cat(str_rep("=", 80), "\n\n", file = fpath, append = TRUE)
}


write_coxph_summary <- function(fit, title, fname, clear_file = FALSE) {
    fpath <- table_path(GRAPHS_DIR, fname)

    if (clear_file & file.exists(fpath)) file.remove(fpath)

    cat(title, "\n", file = fpath, append = TRUE)
    write_summary(summary(fit), fpath)
    cat(str_rep("=", 80), "\n\n", file = fpath, append = TRUE)
}


#### ---- Overall survival curves ---- ####

# A model of all cancers compared to each other.
all_cancer_survival_model <- function(data) {
    fit <- survfit(Surv(time = time, event = status) ~ cancer, data = data)
    write_survfit_summary(fit, "All cancers survival model",
                          "all_cancer_model.txt", clear_file = TRUE)
    p <- ggsurvplot(
        fit = fit,
        data = data,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        surv.median.line = "hv",
        fontsize = 3,
        font.family = "arial",
        palette = unname(cancer_palette)
    )
    p <- style_ggsurvminer_plot(p)
    surv_plot <- patch_ggsurvplot(p)
    ggsave_wrapper(
        surv_plot,
        plot_path(GRAPHS_DIR, glue("cancer-survival-curve_ALL.svg")),
        "wide"
    )
}


# Survival curve for a cancer.
cancer_survival_model <- function(cancer, data) {
    fit <- survfit(Surv(time = time, event = status) ~ 1, data = data)
    write_survfit_summary(fit, glue("{cancer} survival model"),
                          "cancer_models.txt")
    p <- ggsurvplot(
        fit = fit,
        data = data,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        surv.median.line = "hv",
        fontsize = 3,
        font.family = "arial",
        palette = cancer_palette[[cancer]]
    ) +
        ggtitle(cancer)
    p <- style_ggsurvminer_plot(p)
    surv_plot <- patch_ggsurvplot(p)
    ggsave_wrapper(
        surv_plot,
        plot_path(GRAPHS_DIR, glue("cancer-survival-curve_{cancer}.svg")),
        "wide"
    )
}


cancer_survival_df %T>%
    all_cancer_survival_model() %>%
    group_by(cancer) %>%
    nest() %>%
    pwalk(cancer_survival_model)


#### ---- KRAS survival curves ---- ####

# Survival curve for a cancer's data separated by KRAS mutated or WT.
kras_mutated_survival_analysis <- function(cancer, data) {
    fit <- survfit(Surv(time = time, event = status) ~ kras_mut,
                   data = data)
    fit_diff <- survdiff(Surv(time = time, event = status) ~ kras_mut,
                         data = data)

    title <- glue("{cancer} KRAS WT vs. mutant")
    fname <- "wt-vs-mutant.txt"
    write_survfit_summary(fit, title, fname)
    write_survdiff_summary(fit_diff, title, fname)

    p <- ggsurvplot(
        fit = fit,
        data = data,
        pval = TRUE,
        conf.int = TRUE,
        linetype = "strata",
        risk.table = TRUE,
        risk.table.col = "strata",
        surv.median.line = "hv",
        fontsize = 3,
        font.family = "arial",
        palette = c("grey50", "black")
    ) +
        ggtitle(cancer)
    p <- style_ggsurvminer_plot(p)
    surv_plot <- patch_ggsurvplot(p)
    ggsave_wrapper(
        surv_plot,
        plot_path(GRAPHS_DIR, glue("kras-mutated-survival-curve_{cancer}.svg")),
        "wide"
    )

    return(tibble(
        cancer = cancer,
        data = list(data),
        kras_mut_model = list(fit)
    ))
}


cancer_survival_df %>%
    group_by(cancer) %>%
    nest() %>%
    pmap_df(kras_mutated_survival_analysis)



#### ---- KRAS allele survival curves ---- ####


# Get the alleles with more than `other_min` samples.
get_alleles_to_keep <- function(df, other_min) {
    df %>%
        group_by(kras_allele) %>%
        summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
        ungroup() %>%
        filter(num_allele_samples > !!other_min) %>%
        pull(kras_allele)
}


# Group the KRAS alleles into "Other" if fewer than `other_min` samples.
group_kras_alleles <- function(df, other_min) {
    alleles_to_keep <- get_alleles_to_keep(df, other_min)
    df %>%
        mutate(kras_allele_grp = case_when(
            kras_allele == "WT" ~ "WT",
            kras_allele %in% !!alleles_to_keep ~ kras_allele,
            TRUE ~ "Other"
        ))
}


# Wrapper around standard survival curve analysis workflow.
allele_group_survival_analysis <- function(cancer,
                                           data,
                                           do_survdiff = TRUE,
                                           do_coxph = TRUE,
                                           p_val_sig = NULL,
                                           model_output_file_template = NULL,
                                           model_output_title_template = NULL,
                                           plot_name_template = NULL,
                                           curve_palette = NULL,
                                           plot_file_size = "wide") {
    title <- glue(model_output_title_template)
    fname <- glue(model_output_file_template)

    write_out <- (
        !is.null(model_output_file_template) &
        !is.null(model_output_title_template)
    )

    fit <- survfit(Surv(time = time, event = status) ~ kras_allele_grp,
               data = data)
    if (write_out) write_survfit_summary(fit, title, fname)

    if (do_survdiff) {
        fit_diff <- survdiff(
            Surv(time = time, event = status) ~ kras_allele_grp,
            data = data
        )
        if (write_out) write_survdiff_summary(fit_diff, title, fname)
    }

    if (do_coxph) {
        fit_coxph <- coxph(
            Surv(time = time, event = status) ~ kras_allele_grp,
            data = data
        )
        if (write_out) write_coxph_summary(fit_coxph, title, fname)
    }

    # If no name for saving plots, return early.
    if (is.null(plot_name_template)) return()

    # If using a cut-off, if the p-value is not low enough, return early.
    if (!is.null(p_val_sig) & do_survdiff) {
        if (survdiff_pval(fit_diff) >= p_val_sig) return()
    }

    p <- ggsurvplot(
        fit = fit,
        data = data,
        pval = TRUE,
        conf.int = FALSE,
        risk.table = TRUE,
        risk.table.col = "strata",
        surv.median.line = "hv",
        fontsize = 3,
        font.family = "arial",
        palette = curve_palette
    ) +
        ggtitle(cancer)
    p <- style_ggsurvminer_plot(p)
    surv_plot <- patch_ggsurvplot(p)
    ggsave_wrapper(
        surv_plot,
        plot_path(GRAPHS_DIR, glue(plot_name_template)),
        plot_file_size
    )
}


# Survival curve for a cancer's data separated by KRAS alleles.
kras_alleles_survival_analysis <- function(cancer, data, other_min = 10) {

    data <- group_kras_alleles(data, other_min)

    allele_group_survival_analysis(
        cancer = cancer,
        data = data,
        do_survdiff = TRUE,
        do_coxph = TRUE,
        model_output_file_template = "kras-alleles.txt",
        model_output_title_template = glue("{cancer} KRAS alleles"),
        plot_name_template = glue("kras-allele-survival-curve_{cancer}.svg"),
        curve_palette = alter_pal_for_ggsurvplot(short_allele_pal,
                                                 "kras_allele_grp")
    )
}


kras_allele_vs_rest_survival_analysis <- function(cancer, data,
                                                  other_min = 10) {
    data <- group_kras_alleles(data, other_min)
    all_alleles <- unique(data$kras_allele_grp)

    pal <- alter_pal_for_ggsurvplot(
        c(short_allele_pal, "rest" = "grey40"),
        "kras_allele_grp"
    )

    for (allele in all_alleles[!all_alleles %in% c("WT", "Other")]) {
        mod_data <- data %>%
            mutate(kras_allele_grp = ifelse(kras_allele_grp == allele,
                                            allele, "rest"))

        filename <- glue("{cancer}_{allele}-vs-rest.txt")
        title <- glue("{cancer} - {allele} vs. the rest")
        plot_title <- glue("allele-vs-rest_{cancer}_{allele}.svg")

        allele_group_survival_analysis(
            cancer = cancer,
            data = mod_data,
            do_survdiff = TRUE,
            do_coxph = FALSE,
            model_output_file_template = filename,
            model_output_title_template = title,
            plot_name_template = plot_title,
            curve_palette = pal
        )
    }
}


allele_combination_tibble <- function(x, m = 2) {
    combn(x, m, stringsAsFactors = FALSE) %>%
        t() %>%
        as.data.frame() %>%
        as_tibble() %>%
        dplyr::rename(allele_1 = V1, allele_2 = V2)
}


kras_allele_vs_each_allele_survival_analysis <- function(cancer, data,
                                                         other_min = 10,
                                                         p_val_sig = 0.1) {
    data <- group_kras_alleles(data, other_min)
    all_alleles <- unique(data$kras_allele_grp)

    pal <- alter_pal_for_ggsurvplot(short_allele_pal, "kras_allele_grp")

    allele_combs <- allele_combination_tibble(all_alleles) %>%
        filter(allele_1 != "Other" & allele_2 != "Other")

    for (i in seq(1, nrow(allele_combs))) {
        allele_1 <- allele_combs$allele_1[[i]]
        allele_2 <- allele_combs$allele_2[[i]]
        mod_data <- data %>%
            filter(kras_allele_grp %in% c(allele_1, allele_2))

        filename <- glue("{cancer}_{allele_1}-vs-{allele_2}.txt")
        title <- glue("{cancer} - {allele_1} vs. {allele_2}")
        plot_title <- glue(
            "allele-vs-allele_{cancer}_{allele_1}-vs-{allele_2}.svg"
        )

        allele_group_survival_analysis(
            cancer = cancer,
            data = mod_data,
            do_survdiff = TRUE,
            do_coxph = FALSE,
            p_val_sig = p_val_sig,
            model_output_file_template = filename,
            model_output_title_template = title,
            plot_name_template = plot_title,
            curve_palette = pal
        )
    }
}


cancer_survival_df %>%
    group_by(cancer) %>%
    nest() %T>%
    pwalk(kras_alleles_survival_analysis) %T>%
    pwalk(kras_allele_vs_rest_survival_analysis) %T>%
    pwalk(kras_allele_vs_each_allele_survival_analysis)
