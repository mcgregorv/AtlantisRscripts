extract_layer <- function(subplot_aes, env) {
  layer <- eval(subplot_aes, env)
  if ("sp_plot" %in% class(layer) | "ggsubplot" %in% class(layer)) {
    stop("Cannot place subplots inside subplots", call. = FALSE)
  }
  if ("ggplot" %in% class(layer)) {
    # propogate data and aesthetics
    layer$data <- NULL
    mapping <- layer$mapping
    layer <- propogate_aes(layer, mapping)
  }
  if (!("proto" %in% class(layer))) {
    stop("subplot aes should be a ggplot2 layer object", call. = FALSE)
  }
  layer
}
