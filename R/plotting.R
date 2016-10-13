# A set of functions written to make data munging and preparation easier
#' @import magrittr
#' @import ggplot2
#' @importFrom grid convertUnit
#' @importFrom grid gpar
#' @importFrom foreach foreach

# necessary for PPV (e.g. Figure 2B)
#' @export
prep_for_polygon <-
    function(df, x = "thresholds", y = "ppv", min_y = "min_ppv") {
    # order by thresholds
    df <- df[order(df[,x], decreasing = TRUE),]

    # min y value
    min_y_val <- df[, min_y][1]

    # top of polygon
    df$y_greater_min <- df[,y] > min_y_val
    df <- df[seq_len(nrow(df)) %>% rep(each = 2) %>% tail(-1) %>% head(-1),]
    rownames(df) <- seq(length = nrow(df))
    df$id <-  rep(seq(nrow(df)/2), each = 2)
    df$type <- y

    # bottom of polygon
    df <- df[rep(seq_len(nrow(df)), each = 2),]
    df[grepl(".1", rownames(df), fixed = TRUE),]$type <- min_y
    df[df$type == min_y, y] <- min_y_val

    # order for geom_polygon (needs to be closed)
    foreach(thisid = unique(df$id), .combine = rbind) %do% {
        thisdata <- filter(df, id == thisid)
        thisdata <-
            thisdata[order(thisdata[, x], thisdata[,y], decreasing = TRUE),]
        rbind(thisdata[1,], thisdata[2,], thisdata[4,], thisdata[3,])
    }
}


# http://stackoverflow.com/a/34859307/2320823
GeomStepHist <-
    ggproto("GeomStepHist", GeomPath, required_aes = c("x"),
            draw_panel = function(data, panel_scales, coord, direction) {
                data <- as.data.frame(data)[order(data$x), ]

                n <- nrow(data)
                i <- rep(1:n, each = 2)
                newdata <- rbind(
                    transform(data[1, ], x = x - width/2, y = 0),
                    transform(data[i, ],
                              x = c(rbind(data$x - data$width/2,
                                          data$x + data$width/2))),
                    transform(data[n, ], x = x + width/2, y = 0)
                )
                rownames(newdata) <- NULL

                GeomPath$draw_panel(newdata, panel_scales, coord)
            }
)

#' @export
geom_step_hist <- function(mapping = NULL, data = NULL, stat = "bin",
                           direction = "hv", position = "stack", na.rm = FALSE,
                           show.legend = NA, inherit.aes = TRUE, ...) {
    layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomStepHist,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            direction = direction,
            na.rm = na.rm,
            ...
        )
    )
}


# GeomPolygon2, GeomViolin2, geom_violin2
# Used to create a violin plot where alpha pertain only to fill (not outline)
# http://stackoverflow.com/questions/34754357/fill-transparency-with-geom-violin/34762360#34762360
GeomPolygon2 <-
    ggproto("GeomPolygon2", Geom,
            draw_panel =
                function(data, panel_scales, coord) {
                    n <- nrow(data)
                    if (n == 1)
                        return(zeroGrob())
                    munched <- coord_munch(coord, data, panel_scales)
                    munched <- munched[order(munched$group), ]
                    first_idx <- !duplicated(munched$group)
                    first_rows <- munched[first_idx, ]
                    ggplot2:::ggname("geom_polygon",
                                     polygonGrob(
                                         munched$x,
                                         munched$y,
                                         default.units = "native",
                                         id = munched$group,
                                         gp = gpar(
                                             col = first_rows$colour,
                                             fill = alpha(first_rows$fill,
                                                          first_rows$alpha),
                                             lwd = first_rows$size * .pt,
                                             lty = first_rows$linetype)
                                         )
                                     )
                    },
            default_aes = aes(colour = "NA", fill = "grey20",
                              size = 0.5, linetype = 1,
                              alpha = NA),
            handle_na = function(data, params) { data },
            required_aes = c("x", "y"),
            draw_key = draw_key_polygon)

GeomViolin2 <-
    ggproto("GeomViolin", Geom,
            setup_data =
                function(data, params) {
                    data$width <- data$width %||%
                        params$width %||% (resolution(data$x, FALSE) * 0.9)
                    plyr::ddply(data, "group", transform,
                                xmin = x - width / 2,
                                xmax = x + width / 2
                                )
                    },
            draw_group = function(self, data, ..., draw_quantiles = NULL) {
                data <- transform(data,
                                  xminv = x - violinwidth * (x - xmin),
                                  xmaxv = x + violinwidth * (xmax - x))
                newdata <- rbind(
                    plyr::arrange(transform(data, x = xminv), y),
                    plyr::arrange(transform(data, x = xmaxv), -y))
                newdata <- rbind(newdata, newdata[1,])
                if (length(draw_quantiles) > 0) {
                    stopifnot(all(draw_quantiles >= 0),
                              all(draw_quantiles <= 1))
                    quantiles <-
                        create_quantile_segment_frame(data, draw_quantiles)
                    aesthetics <- data[
                        rep(1, nrow(quantiles)),
                        setdiff(names(data), c("x", "y")),
                        drop = FALSE]
                    both <- cbind(quantiles, aesthetics)
                    quantile_grob <- GeomPath$draw_panel(both, ...)
                    ggplot2:::ggname("geom_violin",
                                     grobTree(
                                         GeomPolygon2$draw_panel(newdata, ...),
                                         quantile_grob)
                                     )
                    } else {
                        ggplot2:::ggname("geom_violin",
                                         GeomPolygon2$draw_panel(newdata, ...))
                        }
                },
            draw_key = draw_key_polygon,
            default_aes = aes(weight = 1, colour = "grey20", fill = "white",
                              size = 0.5, alpha = NA, linetype = "solid"),
            required_aes = c("x", "y")
)

#' @export
geom_violin2 <- function(mapping = NULL, data = NULL, stat = "ydensity",
                         draw_quantiles = NULL, position = "dodge",
                         trim = TRUE, scale = "area",
                         na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                         ...) {
    layer(data = data,
          mapping = mapping,
          stat = stat,
          geom = GeomViolin2,
          position = position,
          show.legend = show.legend,
          inherit.aes = inherit.aes,
          params = list(trim = trim,
                        scale = scale,
                        draw_quantiles = draw_quantiles,
                        na.rm = na.rm,
                        ...))
}

#' @export
find_cell <- function(table, row, col, name="core-fg") {
    l <- table$layout
    which(l$t == row & l$l == col & l$name == name)
}

#' @export
change_cell <- function(table, row, col, name="core",
                        tnew=NA, bnew=NA, lnew=NA, rnew=NA, znew=TRUE) {
    l <- table$layout

    ind_fg = which(l$t == row & l$l == col & l$name == paste0(name, '-fg'))
    ind_bg = which(l$t == row & l$l == col & l$name == paste0(name, '-bg'))

    if (!is.na(tnew)) {
        table$layout[ind_fg, 't'] <- tnew
        table$layout[ind_bg, 't'] <- tnew
    }
    if (!is.na(bnew)) {
        table$layout[ind_fg, 'b'] <- bnew
        table$layout[ind_bg, 'b'] <- bnew
    }
    if (!is.na(lnew)) {
        table$layout[ind_fg, 'l'] <- lnew
        table$layout[ind_bg, 'l'] <- lnew
    }
    if (!is.na(rnew)) {
        table$layout[ind_fg, 'r'] <- rnew
        table$layout[ind_bg, 'r'] <- rnew
    }
    if (znew) {
        zpos <- max(table$layout$z)
        table$layout[ind_fg, 'z'] <- zpos + 1
        table$layout[ind_bg, 'z'] <- zpos
    }
    table
}

## Add secondary y-axis
## http://stackoverflow.com/a/36761846/2320823
#' @export
secondary_y_axis <- function(g_base, g_secondary) {
    require(grid)
    require(gtable)
    g1 <- ggplotGrob(g_base)
    g2 <- ggplotGrob(g_secondary)

    ## Get the position of the plot panel in g1
    pp <- c(subset(g1$layout, name == "panel", se = t:r))

    # Title grobs have margins.
    # The margins need to be swapped.
    # Function to swap margins -
    # taken from the cowplot package:
    # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R
    vinvert_title_grob <- function(grob) {
        heights <- grob$heights
        grob$heights[1] <- heights[3]
        grob$heights[3] <- heights[1]
        grob$vp[[1]]$layout$heights[1] <- heights[3]
        grob$vp[[1]]$layout$heights[3] <- heights[1]

        grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust
        grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust
        grob$children[[1]]$y <- unit(1, "npc") - grob$children[[1]]$y
        grob
    }

    # Copy xlab from g2 and swap margins
    index <- which(g2$layout$name == "xlab")
    xlab <- g2$grobs[[index]]
    xlab <- vinvert_title_grob(xlab)

    # Put xlab at the top of g1
    g1 <- gtable_add_rows(g1, g2$heights[g2$layout[index, ]$t], pp$t - 1)
    g1 <- gtable_add_grob(g1, xlab, pp$t, pp$l, pp$t, pp$r,
                          clip = "off", name = "topxlab")

    # Get "feet" axis (axis line, tick marks and tick mark labels) from g2
    index <- which(g2$layout$name == "axis-b")
    xaxis <- g2$grobs[[index]]

    # Move the axis line to the bottom (Not needed in your example)
    xaxis$children[[1]]$y <- unit.c(unit(0, "npc"), unit(0, "npc"))

    # Swap axis ticks and tick mark labels
    ticks <- xaxis$children[[2]]
    ticks$heights <- rev(ticks$heights)
    ticks$grobs <- rev(ticks$grobs)

    # Move tick marks
    ticks$grobs[[2]]$y <- ticks$grobs[[2]]$y - unit(1, "npc") + unit(3, "pt")

    # Sswap tick mark labels' margins
    ticks$grobs[[1]] <- vinvert_title_grob(ticks$grobs[[1]])

    # Put ticks and tick mark labels back into xaxis
    xaxis$children[[2]] <- ticks

    # Add axis to top of g1
    g1 <- gtable_add_rows(g1, g2$heights[g2$layout[index, ]$t], pp$t)
    gtable_add_grob(g1, xaxis, pp$t + 1, pp$l, pp$t + 1, pp$r,
                    clip = "off", name = "axis-t")
}


#' @export
find_cell <- function(table, row, col, name="core-fg") {
    l <- table$layout
    which(l$t == row & l$l == col & l$name == name)
}

#' @export
change_cell <- function(table, row, col, name="core",
                        tnew=NA, bnew=NA, lnew=NA, rnew=NA, znew=TRUE) {
    l <- table$layout

    ind_fg = which(l$t == row & l$l == col & l$name == paste0(name, '-fg'))
    ind_bg = which(l$t == row & l$l == col & l$name == paste0(name, '-bg'))

    if (!is.na(tnew)) {
        table$layout[ind_fg, 't'] <- tnew
        table$layout[ind_bg, 't'] <- tnew
    }
    if (!is.na(bnew)) {
        table$layout[ind_fg, 'b'] <- bnew
        table$layout[ind_bg, 'b'] <- bnew
    }
    if (!is.na(lnew)) {
        table$layout[ind_fg, 'l'] <- lnew
        table$layout[ind_bg, 'l'] <- lnew
    }
    if (!is.na(rnew)) {
        table$layout[ind_fg, 'r'] <- rnew
        table$layout[ind_bg, 'r'] <- rnew
    }
    if (znew) {
        zpos <- max(table$layout$z)
        table$layout[ind_fg, 'z'] <- zpos + 1
        table$layout[ind_bg, 'z'] <- zpos
    }
    table
}

# adjustment for logticks on the outside and 
# no tickmarks outside the plotted area
# https://groups.google.com/forum/#!topic/ggplot2/OlGA8Gm9O7w
GeomLogticks2 <- ggproto("GeomLogticks", Geom,
  extra_params = "",
  handle_na = function(data, params) {
    data
  },

  draw_panel = function(data, panel_scales, coord, base = 10, sides = "bl",
    scaled = TRUE, short = unit(0.1, "cm"), mid = unit(0.2, "cm"),
    long = unit(0.3, "cm"), lineend = "butt")
  {
    ticks <- list()

    # Convert these units to numbers so that they can be put in data frames
    short <- convertUnit(short, "cm", valueOnly = TRUE)
    mid   <- convertUnit(mid,   "cm", valueOnly = TRUE)
    long  <- convertUnit(long,  "cm", valueOnly = TRUE)

    if (grepl("[b|t]", sides)) {

      # Get positions of x tick marks
      xticks <- ggplot2:::calc_logticks(
        base = base,
        minpow = floor(panel_scales$x.range[1]),
        maxpow = ceiling(panel_scales$x.range[2]),
        start = 0,
        shortend = short,
        midend = mid,
        longend = long
      )

      if (scaled)
        xticks$value <- log(xticks$value, base)

      # adapted from 
      # https://groups.google.com/forum/#!topic/ggplot2/OlGA8Gm9O7w
      xticks <- with(xticks, xticks[value > panel_scales$x.range[[1]] &
                                        value < panel_scales$x.range[[2]],])

      names(xticks)[names(xticks) == "value"] <- "x"   # Rename to 'x' for coordinates$transform
      xticks <- coord$transform(xticks, panel_scales)
      
      # Make the grobs
      if (grepl("b", sides)) {
        ticks$x_b <- with(data, segmentsGrob(
          x0 = unit(xticks$x, "native"), x1 = unit(xticks$x, "native"),
          y0 = unit(xticks$start, "cm"), y1 = unit(xticks$end, "cm"),
          gp = gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt,
                    lineend = lineend)
        ))
      }
      if (grepl("t", sides)) {
        ticks$x_t <- with(data, segmentsGrob(
          x0 = unit(xticks$x, "native"), x1 = unit(xticks$x, "native"),
          y0 = unit(1, "npc") - unit(xticks$start, "cm"), y1 = unit(1, "npc") - unit(xticks$end, "cm"),
          gp = gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt,
                    lineend = lineend)
        ))
      }
    }


    if (grepl("[l|r]", sides)) {
      yticks <- ggplot2:::calc_logticks(
        base = base,
        minpow = floor(panel_scales$y.range[1]),
        maxpow = ceiling(panel_scales$y.range[2]),
        start = 0,
        shortend = short,
        midend = mid,
        longend = long
      )

      if (scaled)
        yticks$value <- log(yticks$value, base)

      names(yticks)[names(yticks) == "value"] <- "y"   # Rename to 'y' for coordinates$transform
      yticks <- coord$transform(yticks, panel_scales)

      # Make the grobs
      if (grepl("l", sides)) {
        ticks$y_l <- with(data, segmentsGrob(
          y0 = unit(yticks$y, "native"), y1 = unit(yticks$y, "native"),
          x0 = unit(yticks$start, "cm"), x1 = unit(yticks$end, "cm"),
          gp = gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt,
                    lineend = lineend)
        ))
      }
      if (grepl("r", sides)) {
        ticks$y_r <- with(data, segmentsGrob(
          y0 = unit(yticks$y, "native"), y1 = unit(yticks$y, "native"),
          x0 = unit(1, "npc") - unit(yticks$start, "cm"), x1 = unit(1, "npc") - unit(yticks$end, "cm"),
          gp = gpar(col = alpha(colour, alpha), lty = linetype, lwd = size * .pt,
                    lineend = lineend)
        ))
      }
    }

    gTree(children = do.call("gList", ticks))
  },

  default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = 1)
)

#' @export
annotation_logticks2 <-
    function(base = 10, sides = "bl", scaled = TRUE, short = unit(0.1, "cm"),
            mid = unit(0.2, "cm"), long = unit(0.3, "cm"), colour = "black",
            size = 0.5, linetype = 1, alpha = 1, color = NULL, lineend = "butt",
            ...) {
    if (!is.null(color)) 
        colour <- color
    layer(data = data.frame(x = NA), mapping = NULL, stat = StatIdentity,
          geom = GeomLogticks2, position = PositionIdentity, show.legend = FALSE,
          inherit.aes = FALSE,
          params = list(base = base, sides = sides,
                        scaled = scaled, short = short, mid = mid, long = long,
                        colour = colour, size = size, linetype = linetype,
                        alpha = alpha, lineend = lineend, ...))
}
