geom_subplot <- function(mapping, width = rel(0.95), height = rel(0.95),
                         data = waiver(), x_scale = identity, y_scale = identity,
                         position = "identity", reference = NULL,
                         ply.aes = TRUE, .ref = FALSE) {
  
  missing <- c(is.null(mapping$x), is.null(mapping$y), is.null(mapping$group),
               is.null(mapping$subplot))
  if (any(missing)) {
    stop(paste("Missing required aesthetics in geom_subplot:",
               paste(c("x", "y", "group", "subplot")[missing], collapse = ", ")))
  }
  
  if (position != "identity" & position != "merge") {
    stop("geom_subplot only supports position = 'identity' or 'merge'",
         call. = FALSE)
  }
  
  if (position == "merge") {
    merge.overlaps <- TRUE
  } else {
    merge.overlaps <- FALSE
  }
  
  layer <- extract_layer(mapping$subplot, parent.frame())
  layer$data <- data
  mapping$subplot <- NULL
  layer$embed <- list(width = width, height = height, x_scale = x_scale,
                      y_scale = y_scale, merge.overlaps = merge.overlaps, major.aes = mapping)
  layer$assign_subplots <- assign_subplots
  layer$combine_subplots <- combine_subplots
  if (.ref) layer$combine_subplots <- combine_refs
  
  if (is.null(reference)) {
    if (ply.aes) {
      ply_aes(sp_layer(layer))
    } else {
      sp_layer(layer)
    }
  } else {
    ref.layer <- reference(layer, "subplot", mapping, width, height, position)
    if (ply.aes) {
      list(ref.layer, ply_aes(sp_layer(layer)))
    } else {
      list(ref.layer, sp_layer(layer))
    }
  }
}