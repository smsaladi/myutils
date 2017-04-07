#' @import grid
#' @import magrittr
#' @import ggplot2

#' @export
#' @rdname geom_tile2
#' geom rect with alpha for tranparency only
geom_rect2 <- function(mapping = NULL, data = NULL,
                      stat = "identity", position = "identity",
                      ...,
                      na.rm = FALSE,
                      show.legend = NA,
                      inherit.aes = TRUE) {
    layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomRect2,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            na.rm = na.rm,
            ...
        )
    )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomRect2 <- ggproto("GeomRect2", Geom,
                    default_aes = aes(colour = NA, fill = "grey35", size = 0.5, linetype = 1,
                                      alpha = NA),
                    
                    required_aes = c("xmin", "xmax", "ymin", "ymax"),
                    
                    draw_panel = function(self, data, panel_scales, coord) {
                        if (!coord$is_linear()) {
                            aesthetics <- setdiff(
                                names(data), c("x", "y", "xmin", "xmax", "ymin", "ymax")
                            )
                            
                            polys <- plyr::alply(data, 1, function(row) {
                                poly <- ggplot2:::rect_to_poly(row$xmin, row$xmax, row$ymin, row$ymax)
                                aes <- as.data.frame(row[aesthetics],
                                                     stringsAsFactors = FALSE)[rep(1,5), ]
                                
                                GeomPolygon$draw_panel(cbind(poly, aes), panel_scales, coord)
                            })
                            
                            ggname("bar", do.call("grobTree", polys))
                        } else {
                            coords <- coord$transform(data, panel_scales)
                            ggname("geom_rect", rectGrob(
                                coords$xmin, coords$ymax,
                                width = coords$xmax - coords$xmin,
                                height = coords$ymax - coords$ymin,
                                default.units = "native",
                                just = c("left", "top"),
                                gp = gpar(
                                    col = coords$colour,
                                    fill = alpha(coords$fill, coords$alpha),
                                    lwd = coords$size * .pt,
                                    lty = coords$linetype,
                                    lineend = "butt"
                                )
                            ))
                        }
                    },
                    
                    draw_key = draw_key_polygon
)
