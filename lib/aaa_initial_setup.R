# Script that runs first from "lib/".

# Declare the `select()` function from 'dplyr' to override that
# from 'AnnotationDbi'.
conflicted::conflict_prefer("select", "dplyr")
