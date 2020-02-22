# General functions used in survival analysis.


# Get p-value from results of `survdiff()`.
survdiff_pval <- function(fit) {
    pchisq(fit$chisq, length(fit$n) - 1, lower.tail = FALSE)
}


# Get p-value from the Log rank test resulting from using `coxph()`.
coxph_logtest_pval <- function(fit) {
    summary(fit)$logtest[["pvalue"]]
}


# Get p-value from the Wald test resulting from using `coxph()`.
coxph_logtest_pval <- function(fit) {
    summary(fit)$waldtest[["pvalue"]]
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


# Style the survival curve from `ggsurvplot()` or ``ggadjustedcurves()`,
style_ggsurvminer_surv_curve <- function(p, x_expand = c(0.02, 0.01)) {
    p_new <- p +
        scale_x_continuous(expand = expand_scale(mult = x_expand)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.01))) +
        theme_classic(base_size = 7, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            strip.background = element_blank()
        )
    return(p_new)
}


# Style the table from `ggsurvplot()`.
style_ggsurvminer_table <- function(p, table_y, x_expand = c(0.02, 0.01)) {
    p_new <- p +
        scale_x_continuous(expand = expand_scale(mult = x_expand)) +
        theme_classic(base_size = 7, base_family = "arial") +
        theme(
            legend.position = "none",
            plot.title = element_blank(),
            axis.ticks = element_blank()
        ) +
        labs(y = table_y)
    return(p_new)
}


# Style the survival curve and the table from `ggsurvplot()`.
style_ggsurvminer_plot <- function(ggsurv_obj,
                                   table_y = "number at risk") {
    x_expand <- c(0.02, 0.01)
    ggsurv_obj$plot <- style_ggsurvminer_surv_curve(ggsurv_obj$plot, x_expand)
    ggsurv_obj$table <- style_ggsurvminer_table(ggsurv_obj$table, x_expand)
    return(ggsurv_obj)
}


# Write a model's summary to file.
write_summary <- function(fit_summary, fpath) {
    sink(fpath, append = TRUE)
    print(fit_summary)
    sink()
}


# Write summary for a `survfit` model.
write_survfit_summary <- function(fit, title, fname, clear_file = FALSE) {
    fpath <- table_path(GRAPHS_DIR, fname)

    if (clear_file & file.exists(fpath)) file.remove(fpath)

    cat(title, "\n", file = fpath, append = TRUE)
    write_summary(summary(fit)$table, fpath)
    cat(str_rep("=", 80), "\n\n", file = fpath, append = TRUE)
}


# Write summary for a `survdiff` model.
write_survdiff_summary <- function(fit, title, fname, clear_file = FALSE) {
    fpath <- table_path(GRAPHS_DIR, fname)

    if (clear_file & file.exists(fpath)) file.remove(fpath)

    cat(title, "\n", file = fpath, append = TRUE)
    write_summary(fit, fpath)
    cat(str_rep("=", 80), "\n\n", file = fpath, append = TRUE)
}


# Write summary for a `coxph` model.
write_coxph_summary <- function(fit, title, fname, clear_file = FALSE) {
    fpath <- table_path(GRAPHS_DIR, fname)

    if (clear_file & file.exists(fpath)) file.remove(fpath)

    cat(title, "\n", file = fpath, append = TRUE)
    write_summary(summary(fit), fpath)
    cat(str_rep("=", 80), "\n\n", file = fpath, append = TRUE)
}
