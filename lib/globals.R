# Add any project specific configuration here.
ProjectTemplate::add.config(
  apply.override = FALSE
)

# Add project specific configuration that can be overridden from load.project()
ProjectTemplate::add.config(
  apply.override = TRUE
)

#### ---- Packages used frequently ---- ####

library(stats)
library(org.Hs.eg.db)
library(glue)
library(conflicted)
library(assertr)
library(testthat)
library(glmnet)
library(parallel)
library(caret)
library(ggfortify)
library(tidygraph)
library(jhcutils)
library(clisymbols)
library(magrittr)
library(nakedpipe)
library(ggpubr)
library(ggraph)
library(ggtext)
library(patchwork)
library(ggplot2)
library(broom)
library(tidyverse)


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
