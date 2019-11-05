
library(tidyverse)


replace_pheatmap <- function(line) {
    old_img_regex <- "COAD_pheatmap\\.svg|LUAD_pheatmap\\.svg|PAAD_pheatmap\\.svg"
    if (str_detect(line, old_img_regex)) {
        old_img_name <- str_extract(line, old_img_regex)
        cancer <- str_extract(old_img_name, "COAD|LUAD|PAAD")
        new_img_name <- paste0(cancer, "_CRISPR_pheatmap.svg")
        new_line <- str_replace(line, old_img_name, new_img_name)
        return(new_line)
    } else {
        return(line)
    }
}


posts <- list.files("reports/content/post", pattern = "md$", full.names = TRUE)

for (post in posts) {
    post_lines <- readLines(post)
    fixed_post_lines <- purrr::map(post_lines, replace_pheatmap)
    writeLines(unlist(fixed_post_lines), post)
}
