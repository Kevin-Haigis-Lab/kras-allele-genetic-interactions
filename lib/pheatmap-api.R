# Functions for interfacing with the 'pheatmap' library.

#### ---- Saving a 'pheatmap' heatmap ---- ####

#' Save a pheatmap object
save_pheatmap_svg <- function(x, save_path, size = NA, width = NA, height = NA) {
  size <- decide_size(size = size[[1]], width = width, height = height)
  svg(save_path, width = size[[1]], height = size[[2]])
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, save_path, size = NA, width = NA, height = NA) {
  size <- decide_size(size = size[[1]], width = width, height = height)
  pdf(save_path, width = size[[1]], height = size[[2]])
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#### ---- Modify a pheatmap plot ---- ####

italicize_pheatmap_rownames <- function(ph) {
  mod_gpar <- ph$gtable$grobs[[5]]$gp
  mod_gpar$fontface <- as.integer(3)
  mod_gpar$font <- as.integer(3)
  ph$gtable$grobs[[5]]$gp <- mod_gpar
  return(ph)
}

remove_bold_in_pheatmap <- function(ph) {
  one <- integer(1)
  for (i in c(7, 9)) {
    mod_gpar <- ph$gtable$grobs[[i]]$gp
    mod_gpar$fontface <- one
    mod_gpar$font <- one
    ph$gtable$grobs[[i]]$gp <- mod_gpar
  }
  return(ph)
}


reduce_height_of_column_anno <- function(ph,
                                         new_height = grid::unit(1.5, "bigpts"),
                                         new_y = grid::unit(7.75, "bigpts")) {
  mod_gpar <- ph$gtable$grobs[[6]]
  mod_gpar$height <- new_height
  mod_gpar$y <- rep(new_y, length(mod_gpar$y))
  ph$gtable$grobs[[6]] <- mod_gpar
  return(ph)
}
