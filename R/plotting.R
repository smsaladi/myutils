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

#' modified from GGally package
#' @export
ggcorr <- function(data, method = c("pairwise", "pearson"), cor_matrix = NULL, 
          nbreaks = NULL, digits = 2, name = "", low = "#3B9AB2", mid = "#EEEEEE", 
          high = "#F21A00", midpoint = 0, palette = NULL, geom = "tile", 
          min_size = 2, max_size = 6, label = FALSE, label_alpha = FALSE, 
          label_color = "black", label_round = 1, label_size = 4,
          limits = c(-1,  1), drop = is.null(limits) || identical(limits, FALSE),
          layout.exp = 0, legend.position = "right", title.position = "top", legend.size = 9,
          diagColors = NULL, ...) 
{
    if (is.numeric(limits)) {
        if (length(limits) != 2) {
            stop("'limits' must be of length 2 if numeric")
        }
    }
    if (is.logical(limits)) {
        if (limits) {
            limits <- c(-1, 1)
        }
        else {
            limits <- NULL
        }
    }
    if (length(geom) > 1 || !geom %in% c("blank", "circle", "text", "tile")) {
        stop("incorrect geom value")
    }
    if (length(method) == 1) {
        method = c(method, "pearson")
    }
    if (!is.null(data)) {
        if (!is.data.frame(data)) {
            data = as.data.frame(data)
        }
        x = which(!sapply(data, is.numeric))
        if (length(x) > 0) {
            warning(paste("data in column(s)", paste0(paste0("'", 
                                                             names(data)[x], "'"), collapse = ", "), "are not numeric and were ignored"))
            data = data[, -x]
        }
    }
    if (is.null(cor_matrix)) {
        cor_matrix = cor(data, use = method[1], method = method[2])
    }
    m = cor_matrix
    m = reshape2::melt(m * lower.tri(m))
    names(m) = c("x", "y", "coefficient")
    m$coefficient[m$coefficient == 0] = NA
    if (!is.null(nbreaks)) {
        x = seq(-1, 1, length.out = nbreaks + 1)
        if (!nbreaks %% 2) {
            x = sort(c(x, 0))
        }
        m$breaks = cut(m$coefficient, breaks = unique(x), include.lowest = TRUE, 
                       dig.lab = digits)
    }
    if (is.null(midpoint)) {
        midpoint = median(m$coefficient, na.rm = TRUE)
        message(paste("Color gradient midpoint set at median correlation to", 
                      round(midpoint, 2)))
    }
    m$label = round(m$coefficient, label_round)
    p = ggplot(na.omit(m), aes(x, y))
    if (geom == "tile") {
        if (is.null(nbreaks)) {
            p = p + geom_tile(aes(
                fill = coefficient#, height = abs(coefficient),
               # width = abs(coefficient)
                ), color = "white")
        }
        else {
            p = p + geom_tile(aes(fill = breaks), color = "white")
        }
        if (is.null(nbreaks) && !is.null(limits)) {
            p = p + scale_fill_gradient2(name, low = low, mid = mid, 
                                         high = high, midpoint = midpoint, limits = limits)
        }
        else if (is.null(nbreaks)) {
            p = p + scale_fill_gradient2(name, low = low, mid = mid, 
                                         high = high, midpoint = midpoint)
        }
        else if (is.null(palette)) {
            x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))
            p = p + scale_fill_manual(name, values = x, drop = drop)
        }
        else {
            p = p + scale_fill_brewer(name, palette = palette, 
                                      drop = drop)
        }
    }
    else if (geom == "circle") {
        p = p + geom_point(aes(size = abs(coefficient) * 1.25), 
                           color = "grey50")
        if (is.null(nbreaks)) {
            p = p + geom_point(aes(size = abs(coefficient), color = coefficient))
        }
        else {
            p = p + geom_point(aes(size = abs(coefficient), color = breaks))
        }
        p = p + scale_size_continuous(range = c(min_size, max_size)) + 
            guides(size = FALSE)
        r = list(size = (min_size + max_size)/2)
        if (is.null(nbreaks) && !is.null(limits)) {
            p = p + scale_color_gradient2(name, low = low, mid = mid, 
                                          high = high, midpoint = midpoint, limits = limits)
        }
        else if (is.null(nbreaks)) {
            p = p + scale_color_gradient2(name, low = low, mid = mid, 
                                          high = high, midpoint = midpoint)
        }
        else if (is.null(palette)) {
            x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))
            p = p + scale_color_manual(name, values = x, drop = drop) + 
                guides(color = guide_legend(override.aes = r, title.position = title.position))
        }
        else {
            p = p + scale_color_brewer(name, palette = palette, 
                                       drop = drop) + guides(color = guide_legend(override.aes = r, title.position = title.position))
        }
    }
    else if (geom == "text") {
        if (is.null(nbreaks)) {
            p = p + geom_text(aes(label = label, color = coefficient), 
                              size = label_size)
        }
        else {
            p = p + geom_text(aes(label = label, color = breaks), 
                              size = label_size)
        }
        if (is.null(nbreaks) && !is.null(limits)) {
            p = p + scale_color_gradient2(name, low = low, mid = mid, 
                                          high = high, midpoint = midpoint, limits = limits)
        }
        else if (is.null(nbreaks)) {
            p = p + scale_color_gradient2(name, low = low, mid = mid, 
                                          high = high, midpoint = midpoint)
        }
        else if (is.null(palette)) {
            x = colorRampPalette(c(low, mid, high))(length(levels(m$breaks)))
            p = p + scale_color_manual(name, values = x, drop = drop)
        }
        else {
            p = p + scale_color_brewer(name, palette = palette, 
                                       drop = drop)
        }
    }
    if (label) {
        if (isTRUE(label_alpha)) {
            p = p + geom_text(aes(x, y, label = label, alpha = abs(coefficient)), 
                              color = label_color, size = label_size, show.legend = FALSE)
        }
        else if (label_alpha > 0) {
            p = p + geom_text(aes(x, y, label = label, show_guide = FALSE), 
                              alpha = label_alpha, color = label_color, size = label_size)
        }
        else {
            p = p + geom_text(aes(x, y, label = label), color = label_color, 
                              size = label_size)
        }
    }
    textData <- m[m$x == m$y & is.na(m$coefficient), ]
    xLimits <- levels(textData$y)
    textData$diagLabel <- textData$x
    if (!is.null(diagColors)) {
        textData <- left_join(textData, diagColors, by = c("x" = "feature"))
    }
    print(textData)
    if (!is.numeric(layout.exp) || layout.exp < 0) {
        stop("incorrect layout.exp value")
    }
    else if (layout.exp > 0) {
        layout.exp <- as.integer(layout.exp)
        textData <- rbind(textData[1:layout.exp, ], textData)
        spacer <- paste(".ggally_ggcorr_spacer_value", 1:layout.exp, 
                        sep = "")
        textData$x[1:layout.exp] <- spacer
        textData$diagLabel[1:layout.exp] <- NA
        xLimits <- c(spacer, levels(m$y))
    }
    if (is.null(diagColors)) {
        p = p + geom_text(data = textData, aes_string(label = "diagLabel"), ..., 
                      na.rm = TRUE)
    } else {
        p = p + geom_text(data = textData,
                          aes_string(label = "diagLabel", color = "colorgrp"),
                          ..., na.rm = TRUE)
    }
    p = p +
        scale_x_discrete(breaks = NULL, limits = xLimits) +
        scale_y_discrete(breaks = NULL, limits = levels(m$y)) + 
        labs(x = NULL, y = NULL) + coord_equal() +
        theme(panel.background = element_blank(), 
              legend.key = element_blank(), legend.position = legend.position,
              legend.title = element_text(size = legend.size),
              legend.text = element_text(size = legend.size))
    return(p)
}
