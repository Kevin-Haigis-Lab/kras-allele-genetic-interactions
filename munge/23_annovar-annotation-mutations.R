# Annotate the coding mutation using ANNOVAR

#### ---- Functions for preparing ANNOVAR input ---- ####

# Split he `genomic_position` column into the input for ANNOVAR.
split_genomic_position_col <- function(v) {
  # Split the data and add column names.
  df <- as_tibble(as.data.frame(str_split_fixed(v, ":|,", 4)))
  colnames(df) <- c("Chr", "Start", "Ref", "Alt")

  # Fix Start and End position.
  df$Start <- as.numeric(df$Start)
  df$End <- df$Start + (str_length(df$Ref) - 1)

  # Put in the correct order for ANNOVAR.
  df %<>% select(Chr, Start, End, Ref, Alt)

  if (any(is.na(df))) stop("There are missing values in the data frame.")

  return(df)
}


# Make the column with the comments for ANNOVAR input.
make_comments_col <- function(df) {
  tibble(comments_col = paste("comments:", seq(1, nrow(df))))
}

transform_to_annovar_input <- function(df) {
  bind_cols(
    split_genomic_position_col(df$genomic_position),
    make_comments_col(df)
  )
}


#### ---- Functions for preparing ANNOVAR output ---- ####

# Is the mutation clinically significant.
is_clinisig_pathogenic <- function(cs) {
  if (is.na(cs)) {
    return(FALSE)
  }
  cs_split <- str_to_lower(unlist(str_split(cs, ",|\\|")))
  any(cs_split %in% c("pathogenic", "drug response", "likely pathogenic"))
}


# Add logical columns for each prediction for if they predict harmful mutations.
access_predictions <- function(avdf) {
  avdf %>%
    mutate(
      is_sift_pred = sift_pred == "D",
      is_polyphen2_hdiv_pred = polyphen2_hdiv_pred %in% c("P", "D"),
      is_polyphen2_hvar_pred = polyphen2_hvar_pred %in% c("P", "D"),
      is_lrt_pred = lrt_pred == "D",
      is_mutation_taster_pred = mutation_taster_pred %in% c("A", "D"),
      is_mutation_assessor_pred = mutation_assessor_pred %in% c("L", "N"),
      is_fathmm_pred = fathmm_pred == "D",
      is_provean_pred = provean_pred == "D",
      is_fathmm_mkl_coding_pred = fathmm_mkl_coding_pred == "D",
      is_meta_svm_pred = meta_svm_pred == "D",
      is_meta_lr_pred = meta_lr_pred == "D",
      is_gerp = as.numeric(gerp_rs) > quantile(as.numeric(gerp_rs), 0.70, na.rm = TRUE),
      is_phast_cons20way_mammalian = phast_cons20way_mammalian > 0.70,
      si_phy_29way_log_odds = as.numeric(si_phy_29way_log_odds),
      is_si_phy_29way_log_odds =
        si_phy_29way_log_odds > quantile(si_phy_29way_log_odds, 0.70, na.rm = TRUE),
      is_icgc = !is.na(icgc_occurrence),
      is_cosmic = !is.na(cosmic68wgs),
      is_clinsig = purrr::map_lgl(clinsig, is_clinisig_pathogenic)
    )
}


# Replace the "." from ANNOVAR
replace_annovar_nas <- function(df) {
  df[df == "."] <- NA
  return(df)
}


#### ---- Run ANNOVAR annotation ---- ####

ProjectTemplate::cache("cancer_coding_av_muts_df",
  depends = "cancer_coding_muts_df",
  {
    # Make input for ANNOVAR
    annovar_in_file <- file.path("data", "annovar", "annovar_input.tsv")

    transform_to_annovar_input(cancer_coding_muts_df) %>%
      write_tsv(annovar_in_file,
        append = FALSE,
        col_names = FALSE
      )


    # Run ANNOVAR script and wait for output.

    annovar_out_file <- file.path("data", "annovar", "annovar_out")
    annovar_out_file_final <- paste0(annovar_out_file, ".hg19_multianno.txt")

    system(glue("sbatch munge/24_annovar-annotation-mutations.sh {annovar_in_file} {annovar_out_file}"))

    annovar_log_file <- file.path(
      "logs", "annovar", "24_annovar-annotation-mutations.log"
    )

    ANNOVAR_DONE <- FALSE
    while (!ANNOVAR_DONE) {
      # If file exists and the specific finish tag exists, then continue with
      # this R script.
      if (file.exists(annovar_log_file)) {
        annovar_log <- unlist(tail(readLines(annovar_log_file)))
        if (any(str_detect(annovar_log, "ANNOVAR FINISHED"))) {
          if (file.exists(annovar_out_file_final)) ANNOVAR_DONE <- TRUE
        }
      }
      Sys.sleep(10)
    }
    # Provide a short hold to make sure ANNOVAR is done and some O2 lag.
    Sys.sleep(120)

    annovar_df <- read_tsv(annovar_out_file_final, progress = FALSE) %>%
      janitor::clean_names() %>%
      mutate(matched_row = as.numeric(str_extract(otherinfo, "[:digit:]+"))) %>%
      replace_annovar_nas() %>%
      access_predictions()

    cancer_coding_av_muts_df <- cancer_coding_muts_df %>%
      mutate(matched_row = seq(1, n())) %>%
      left_join(annovar_df, by = "matched_row")

    return(cancer_coding_av_muts_df)
  }
)
