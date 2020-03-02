# Script that runs first from "lib/".

#### ---- Conflicts ---- ####

# Declare which namespaces to use for conflicting functions.

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")