# Observed vs. predicted KRAS allele frequency.
# This is another analysis with the same goal. I just want to make sure
#   it was done correctly. Hopefully, I get the same results as the first time.

GRAPHS_DIR <- "50_11_observed-predicted-kras-alleles_v2"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)

artifact_signatures <- c("sigAR1", "sigAR2")


#### ---- Observed vs. Predicted KRAS allele frequency          ---- ####
#### ---- using the mutational signature contributions directly ---- ####


#### ---- Matrix A: KRAS allele x mutational signature ---- ####

kras_by_mutsig_matrix <- kras_trinucleotide_contexts %>%
    left_join(mutational_signature_spectra, by = "tricontext") %>%
    select(kras_allele, signature, composition) %>%
    mutate(signature = paste0("sig", signature)) %>%
    filter(!signature %in% artifact_signatures) %>%
    group_by(kras_allele, signature) %>%
    summarise(composition = sum(composition)) %>%
    ungroup() %>%
    pivot_wider(kras_allele,
                names_from = signature,
                values_from = composition) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    column_to_rownames("kras_allele") %>%
    as.matrix()

kras_by_mutsig_matrix[, 1:5]


residual_signature <- matrix(1 - colSums(kras_by_mutsig_matrix), nrow = 1)
rownames(residual_signature) <- "remainder"


kras_by_mutsig_matrix <- rbind(kras_by_mutsig_matrix, residual_signature)

kras_by_mutsig_matrix[, 1:5]


#### ---- Matrix B: mutational signature x tumor sample barcode ---- ####

MIN_NUM_MUTATION_PER_SAMPLE <- 10

REMOVE_TSB <- trinucleotide_mutations_df %>%
    count(tumor_sample_barcode) %>%
    filter(n < MIN_NUM_MUTATION_PER_SAMPLE) %>%
    pull(tumor_sample_barcode)

length(REMOVE_TSB)


mutsig_by_tsb_matrix <- mutational_signatures_df %>%
    filter(!tumor_sample_barcode %in% REMOVE_TSB) %>%
    select(tumor_sample_barcode, signature, contribution) %>%
    mutate(signature = paste0("sig", signature)) %>%
    filter(!signature %in% artifact_signatures) %>%
    pivot_wider(tumor_sample_barcode,
                names_from = signature,
                values_from = contribution) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    column_to_rownames("tumor_sample_barcode") %>%
    as.matrix() %>%
    t()

mutsig_by_tsb_matrix[, 1:5]


#### ---- Organize matrices A and B ---- ####

# Set column names of A equal to row names of B.
all_sigs <- unique(paste0("sig", mutational_signature_spectra$signature))
all_sigs <- setdiff(all_sigs, artifact_signatures)
all(all_sigs %in% colnames(kras_by_mutsig_matrix))
all(all_sigs %in% rownames(mutsig_by_tsb_matrix))

kras_by_mutsig_matrix <- kras_by_mutsig_matrix[, all_sigs]
mutsig_by_tsb_matrix <- mutsig_by_tsb_matrix[all_sigs, ]

identical(colnames(kras_by_mutsig_matrix), rownames(mutsig_by_tsb_matrix))

kras_by_mutsig_matrix[1:5, 1:5]
mutsig_by_tsb_matrix[1:5, 1:5]


#### ---- Matrix C: KRAS allele x tumor sample barcode ---- ####

kras_by_tsb_matrix <- kras_by_mutsig_matrix %*% mutsig_by_tsb_matrix

kras_by_tsb_matrix[, 1:5]
colSums(kras_by_tsb_matrix[, 1:5])


#### ---- Analysis of results ---- ####

# Tumor sample information to add back to results.
sample_information <- trinucleotide_mutations_df %>%
    distinct(cancer, tumor_sample_barcode, target, kras_allele) %>%
    dplyr::rename(real_kras_allele = kras_allele)


# Matrix to long tibble.
kras_by_tsb_df <- kras_by_tsb_matrix %>%
    as.data.frame() %>%
    rownames_to_column("kras_allele") %>%
    as_tibble() %>%
    pivot_longer(-kras_allele,
                 names_to = "tumor_sample_barcode",
                 values_to = "probability") %>%
    left_join(sample_information, by = "tumor_sample_barcode")


#### ---- Allele frequency as most likely per sample ---- ####

# With "remainder".
kras_by_tsb_df %>%
    group_by(tumor_sample_barcode, cancer, real_kras_allele, target) %>%
    summarise(pred_kras_allele = kras_allele[which.max(probability)]) %>%
    ungroup() %>%
    filter(pred_kras_allele != "remainder" & real_kras_allele != "WT") %T>%
    print() %>%
    filter(pred_kras_allele == real_kras_allele)


allele_freq_tibble <- function(realk, predk) {
    all_alleles <- sort(unique(c(realk, predk)))
    freqs <- bind_rows(table(realk), table(predk)) %>%
        t() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        rownames_to_column("kras_allele") %>%
        as_tibble() %>%
        set_names(c("kras_allele", "real", "pred")) %>%
        mutate(real = ifelse(is.na(real), 0, real),
               pred = ifelse(is.na(pred), 0, pred))
}


# without remainder
kras_freq_per_sample <- kras_by_tsb_df %>%
    filter(kras_allele != "remainder") %>%
    group_by(tumor_sample_barcode, cancer, real_kras_allele, target) %>%
    summarise(pred_kras_allele = kras_allele[which.max(probability)]) %>%
    ungroup() %>%
    filter(real_kras_allele != "WT") %>%
    group_by(cancer) %>%
    summarise(freqs = list(allele_freq_tibble(real_kras_allele,
                                              pred_kras_allele))) %>%
    unnest(freqs) %>%
    group_by(cancer) %>%
    mutate(real_freq  = real / sum(real),
           pred_freq  = pred / sum(pred)) %>%
    ungroup()

kras_freq_per_sample_scatter <- kras_freq_per_sample %>%
    ggplot(aes(x = real_freq, y = pred_freq)) +
    facet_wrap(. ~ cancer, scales = "free") +
    geom_abline(lty = 2, color = "grey80") +
    geom_point(size = 0.7, color = "grey30") +
    ggrepel::geom_text_repel(
        aes(label = kras_allele),
        size = 2, family = "Arial",
        segment.size = 0.1, segment.alpha = 0.5, segment.color = "grey50",
        seed = 0
    ) +
    theme_bw(base_size = 7, base_family = "arial") +
    theme(
        strip.background = element_blank()
    ) +
    labs(x = "observed", y = "predicted")
ggsave_wrapper(kras_freq_per_sample_scatter,
               plot_path(GRAPHS_DIR, "kras_freq_per_sample_scatter.svg"),
               "medium")


#### ---- Allele frequency averaged across samples ---- ####

# The minimum number of samples with a KRAS allele for the KRAS allele to be
# included in the calculation (per cancer).
MIN_NUMBER_OF_KRAS_ALLELE_SAMPLES <- 10

real_kras_allele_freq <- trinucleotide_mutations_df %>%
    filter(kras_allele != "WT") %>%
    distinct(cancer, tumor_sample_barcode, kras_allele) %>%
    count(cancer, kras_allele, name = "real_allele_count") %>%
    group_by(cancer) %>%
    mutate(real_allele_freq = real_allele_count / sum(real_allele_count)) %>%
    ungroup()

kras_freq_sample_avg <- kras_by_tsb_df %>%
    filter(kras_allele != "remainder") %>%
    left_join(real_kras_allele_freq,
              by = c("cancer", "kras_allele")) %>%
    mutate(real_allele_freq = ifelse(is.na(real_allele_freq),
                                     0, real_allele_freq)) %>%
    filter(real_allele_count >= MIN_NUMBER_OF_KRAS_ALLELE_SAMPLES) %>%
    group_by(cancer, kras_allele, real_allele_freq) %>%
    summarise(avg_probability = mean(probability),
              sd_probability = sd(probability)) %>%
    ungroup() %>%
    group_by(cancer) %>%
    mutate(avg_probability = avg_probability / sum(avg_probability),
           sd_probability = sd_probability / sum(avg_probability)) %>%
    ungroup()


plot_obs_pred_pointrange_scatter <- function(cancer, data) {
    max_val <- max(c(data$real_allele_freq, data$avg_probability))

    ggplot(data, aes(x = real_allele_freq, y = avg_probability)) +
        geom_abline(lty = 2, color = "grey80") +
        geom_pointrange(
            aes(ymin = ymin, ymax = ymax),
            fatten = 1.4, size = 0.6, color = "grey30"
        ) +
        ggrepel::geom_text_repel(
            aes(label = kras_allele),
            size = 4, family = "Arial",
            segment.size = 0.7, segment.alpha = 0.5, segment.color = "grey50",
            point.padding = unit(2, "mm"),
            min.segment.length = unit(5, "mm"),
            seed = 0
        ) +
        theme_bw(base_size = 10, base_family = "arial") +
        theme(
            strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_x_continuous(limits = c(0, max_val),
                           expand = expansion(mult = c(0, 0.02))) +
        scale_y_continuous(limits = c(0, max_val),
                           expand = expansion(mult = c(0, 0.02))) +
        coord_fixed() +
        labs(x = "observed", y = "predicted", title = cancer)
}

save_kras_freq_pointrange_scatter <- function(cancer, plt) {
    ggsave_wrapper(
        plt,
        plot_path(GRAPHS_DIR,
                  glue("{cancer}_kras_freq_sample_avg_scatter.svg")),
        "small"
    )
    invisible(NULL)
}

kras_freq_sample_avg %>%
    mutate(ymin = avg_probability - sd_probability,
           ymax = avg_probability + sd_probability) %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(plt = map2(cancer, data, plot_obs_pred_pointrange_scatter),
           tmp = map2(cancer, plt, save_kras_freq_pointrange_scatter))


#### ---- Observed vs. Predicted KRAS allele frequency       ---- ####
#### ---- using the trinucleotide context mutation frequency ---- ####

# TODO: calculate the estimated KRAS allele freq. using the frequency of finding
# the same trinucleotide mutation.


trinucleotide_mutations_df %>%
    filter(kras_allele != "WT") %>%
    distinct(tumor_sample_barcode, kras_allele) %>%
    count(kras_allele, name = "real_allele_count") %>%
    mutate(real_allele_freq = real_allele_count / sum(real_allele_count))