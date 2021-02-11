# Add any project specific configuration here.
ProjectTemplate::add.config(
  apply.override = FALSE
)

# Add project specific configuration that can be overridden from load.project()
ProjectTemplate::add.config(
  apply.override = TRUE
)

default_libraries <- c(
  "stats", "org.Hs.eg.db", "glue", "conflicted", "assertr", "testthat",
  "glmnet", "parallel", "caret", "ggfortify", "tidygraph", "jhcutils",
  "clisymbols", "magrittr", "nakedpipe", "ggpubr", "ggraph", "ggtext",
  "patchwork", "ggplot2", "broom", "tidyverse"
)

for (lib in default_libraries) {
  a <- require(lib, character.only = TRUE, quietly = TRUE)
  if (!a) {
    renv::install(lib)
    library(lib, character.only = TRUE)
  }
}


#### ---- Conflicts ---- ####
# Declare which namespaces to use for conflicting functions.

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("cache", "ProjectTemplate")
conflict_prefer("rename", "dplyr")
conflict_prefer("parLapply", "parallel")
conflict_prefer("which", "Matrix")


#### ---- Options ---- ####

# dplyr
options(dplyr.summarise.inform = FALSE)

## ggrepel
options(ggrepel.max.overlaps = Inf)
