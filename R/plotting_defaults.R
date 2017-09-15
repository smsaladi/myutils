#' Plotting defaults

#' @export
geom_abline <- partial(ggplot2::geom_abline, size = 0.25, linetype = "dashed")

#' @export
geom_hline <- partial(ggplot2::geom_hline, size = 0.25, linetype = "dashed")

#' @export
geom_vline <- partial(ggplot2::geom_vline, size = 0.25, linetype = "dashed")

#' @export
geom_path <- partial(ggplot2::geom_path, lineend = "round")

#' @export
geom_text <- partial(ggplot2::geom_text, size = 6*5/14)

#' @export
annotate_text <- partial(ggplot2::annotate, geom = "text", size = 6*5/14)

#' @export
scale_x_continuous <- partial(ggplot2::scale_x_continuous, expand = c(0, 0))

#' @export
scale_x_discrete <- partial(ggplot2::scale_x_discrete, expand = c(0, 0))

#' @export
scale_y_continuous <- partial(ggplot2::scale_y_continuous, expand = c(0, 0))

#' @export
draw_plot_label <- partial(cowplot::draw_plot_label, size = 8)
