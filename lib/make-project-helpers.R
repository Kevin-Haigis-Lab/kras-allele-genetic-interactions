
ssrc <- function(n1, n2) {
  all_files <- list.files("src", pattern = "R$", full.names = TRUE)
  all_files <- unlist(all_files)

  n1_pad <- str_pad(as.character(n1), width = 2, side = "left", pad = "0")
  n2_pad <- str_pad(as.character(n2), width = 2, side = "left", pad = "0")
  prefix <- as.character(glue("^{n1_pad}_{n2_pad}"))

  f <- all_files[str_detect(basename(all_files), prefix)]
  if (length(f) == 1) {
    message(glue("Running '{f}'"))
    source(f)
  } else if (length(f) == 0) {
    message(glue("No 'src' file found for '{n1_pad}_{n2_pad}...'"))
  } else {
    message("Multiple files found:")
    message(paste(f, collapse = "\n"))
  }
}


# Write all installed packages to a file.
write_package_list <- function(out_file) {
  x <- as.data.frame(installed.packages())
  x <- janitor::clean_names(x)
  x <- dplyr::distinct(x, package, version)
  readr::write_tsv(x, out_file)
}


# Install all packages from a file.
# The file should be a TSV with two columns: 'package' and 'version'.
install_pacakges <- function(package_df_file) {

  get_all_installed_packages <- function() {
    as.data.frame(installed.packages())$Package
  }

  installed_pkgs <- get_all_installed_packages()

  if (!"devtools" %in% installed_pkgs) {
    install.packages("devtools")
  }

  github_installs <- function(loc_pkg) {
    if (!(basename(loc_pkg) %in% installed_pkgs)) {
      devtools::install_github(loc_pkg)
    }
  }

  github_installs("jhrcook/wext")
  github_installs("jhrcook/jhcutils")
  github_installs("easystats/easystats")
  github_installs("moodymudskipper/nakedpipe")

  installed_pkgs <- get_all_installed_packages()
  err_packages <- c()

  read_tsv(package_df_file) %>%
    filter(!(package %in% installed_pkgs)) %>%
    pmap2(package, version, function(pkg, v) {
      message(paste0("Installing '", pkg, "'"))
      a1 <- try(devtools::install_version(package = pkg, version = v))
      if (inherits(a1, "try-error")) {
        a2 <- try(BiocManager::install(pkg))
      }
      if (inherits(a1, "try-error") & inherits(a2, "try-error")) {
        err_packages <- c(err_packages, pkg)
      }
      installed_pkgs <- get_all_installed_packages()
    })

    if (length(err_packages) > 0) {
      message("The following packages could not be installed.")
      cat(err_packages, sep = "\n")
    }
}
