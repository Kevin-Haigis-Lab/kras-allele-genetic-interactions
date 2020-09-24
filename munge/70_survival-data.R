# Preparation of the survival data from TCGA and MMRF.


DATA_DIR <- file.path("data", "cbioportal")


#### ---- TCGA ---- ####

# Extract the cancer name from the directory.
extract_cancer_from_tcga <- function(x) {
  str_to_upper(str_remove(basename(x), "_tcga"))
}

# Read in and prepare the TCGA clinical data from cBioPortal.
get_tcga_clinical <- function(tcga_dir, skip_rows = 4, cancer = NULL) {
  if (is.null(cancer)) cancer <- extract_cancer_from_tcga(tcga_dir)

  path <- file.path(DATA_DIR, tcga_dir, "data_clinical_patient.txt")
  read_tsv(path, skip = skip_rows) %>%
    janitor::clean_names() %>%
    mutate(cancer = !!cancer) %>%
    select(
      cancer, patient_id, age, days_to_birth, sex, weight, ethnicity,
      days_last_followup, days_to_initial_pathologic_diagnosis,
      path_t_stage, path_m_stage,
      os_status, os_months, dss_status, dss_months,
      dfs_status, dfs_months, pfs_status, pfs_months
    )
}

ProjectTemplate::cache("tcga_survival_data", {
  dirs <- c("coad_tcga", "luad_tcga", "paad_tcga")
  tcga_survival_data <- bind_rows(
    purrr::map(dirs, get_tcga_clinical)
  )
  return(tcga_survival_data)
})


#### ---- MMRF ---- ####

# Column names for MMRF survival data.
mmrf_survival_columns <- list(
  "public_id" = "Public ID",
  "deathdy" = "Date of Death",
  "lstalive" = "Last known alive",
  "lvisitdy" = "Last Visit Date from D_VISIT",
  "pddy" = "Date of disease progression",
  "lpddy" = "Last PD date",
  "pdflag" = "Flag if patient had a PD",
  "lastdy" = "Last date on study",
  "trtsdy" = "First date of treatment",
  "trtedy" = "Last date of treatment",
  "ttfpd" = "Time to first PD",
  "ttfpdw" = "Time to first PD (wk)",
  "pfsdy" = "PFS date",
  "ttpfs" = "Time to PFS",
  "censpfs" = "Censor flag: progression-free survival",
  "pfscdy" = "PFS censored date",
  "ttcpfs" = "Time to PFS event (censored)",
  "ttpfsw" = "Time to PFS (wk)",
  "ttcpfsw" = "Time to PFS event (censored) (wk)",
  "ttos" = "Time to OS",
  "censos" = "Censor flag: overall survival",
  "oscdy" = "Overall survival censored date",
  "ttcos" = "Time to OS event (censored)",
  "ttosw" = "Time to OS (wk)",
  "ttcosw" = "Time to OS event (censored) (wk)",
  "ttfrsp" = "Time to first response",
  "censfrsp" = "Censor flag: first response",
  "frspcdy" = "First response date (censored)",
  "ttcfrsp" = "Time to first response (censored)",
  "ttfrspw" = "Time to first response (wk)",
  "ttcfrspw" = "Time to first response (censored) (wk)",
  "mmstatus" = "MM status (derived)"
)


mmrf_perpatient_columns <- list(
  "PUBLIC_ID" = "Public ID",
  "D_PT_gender" = "Gender",
  "D_PT_age" = "Age",
  "D_PT_iss" = "ISS Disease Stage"
)

# Read in and prepare the MMRF clinical data.
get_mmrf_clinical <- function(mmrf_dir, skip_rows = 0, cancer = "MM") {
  cancer <- str_to_upper(cancer)
  path <- file.path(
    DATA_DIR, mmrf_dir,
    "CoMMpass_IA14_FlatFiles",
    "MMRF_CoMMpass_IA14_STAND_ALONE_SURVIVAL.csv"
  )
  df <- read_csv(path, skip = skip_rows) %>%
    mutate(cancer = !!cancer) %>%
    janitor::clean_names() %>%
    select(cancer, !!!names(mmrf_survival_columns)) %>%
    unique()
  colnames(df) <- c("cancer", unname(mmrf_survival_columns))
  df %>%
    janitor::clean_names()
}


get_mmrf_patient <- function(mmrf_dir, skip_rows = 0, cancer = "MM") {
  cancer <- str_to_upper(cancer)
  path <- file.path(
    DATA_DIR, mmrf_dir,
    "CoMMpass_IA14_FlatFiles",
    "MMRF_CoMMpass_IA14_PER_PATIENT.csv"
  )
  df <- read_csv(path, skip = skip_rows) %>%
    mutate(cancer = !!cancer) %>%
    select(cancer, !!!names(mmrf_perpatient_columns)) %>%
    unique()
  colnames(df) <- c("cancer", unname(mmrf_perpatient_columns))
  df %>%
    janitor::clean_names() %>%
    mutate(sex = ifelse(gender == 1, "M", "F")) %>%
    select(-gender)
}


ProjectTemplate::cache("mmrf_survival_data", {
  mmrf_survival_data <- get_mmrf_clinical("mm_mmrf") %>%
    left_join(get_mmrf_patient("mm_mmrf"),
      by = c("cancer", "public_id")
    )
  return(mmrf_survival_data)
})
