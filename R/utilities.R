# A set of functions written to make data munging and preparation easier
library(magrittr)
library(tidyverse)

library(caret)
library(pROC)

library(foreach)


## Used for model training
#' @export
write_dataset <-
    function(df_to_write, filename, outcome_col = "outcome", query_col = "qid",
             csv = FALSE, feature_select = NULL, features_for_training = FALSE)
    {
        if (!is.null(feature_select) && features_for_training) {
            columns <- data.frame(names = colnames(df_to_write))
            columns <- columns[columns$names %in% as.character(feature_select),]
            df_to_write <- subset(df_to_write, select = as.character(columns))
        }
        
        if (!csv) {
            save_cols <- df_to_write %>% select_(outcome_col, query_col)
            df_to_write <- df_to_write %>%
                select(-one_of(outcome_col, query_col)) %>%
                as.matrix %>%
                svmlight.file %>%
                data.frame %>%
                cbind(save_cols, .)
        }
        
        if (!is.null(feature_select) && !features_for_training) {
            columns <- data.frame(names = colnames(df_to_write))
            columns <- columns[columns$names %in% as.character(feature_select),]
            df_to_write <- subset(df_to_write, select = as.character(columns))
        }
        
        if (csv) {
            data_fn <- paste0(filename, '.csv')
            write.csv(df_to_write, file = data_fn, row.names = FALSE, quote = FALSE)
            data_fn
        } else {
            data_fn <- paste0(filename, '.svmlight')
            write.table(df_to_write, file = data_fn, row.names = FALSE,
                        col.names = FALSE, quote = FALSE)
            system(paste("perl -p -i -e 's/ \\d+:NA//g'", data_fn))
            system(paste("perl -p -i -e 's/ \\d+:NaN//g'", data_fn))
            data_fn
        }
    }

## Used for model training
#' @export
calc_performance <- function(df, score_col, outcome_col,
                             subgrp_col, group_col) {
    require(foreach)
    require(reshape2)
    require(dplyr)
    require(magrittr)
    
    cor_df <- foreach(thisgrouping = unique(df[, subgrp_col] %>% unlist),
                      .combine = bind_rows, .multicombine = TRUE) %dopar% {
                          thisgroup <- filter(df, grouping == thisgrouping)
                          pairs <- ConDis.matrix2(thisgroup[, score_col] %>% unlist,
                                                  thisgroup[, outcome_col] %>% unlist)
                          data.frame(subgrp_col = thisgrouping,
                                     group_col = unique(thisgroup[, group_col] %>% unlist),
                                     table(pairs))
                      } %>%
        dcast(group_col + subgrp_col ~ pairs, value.var = "Freq")
    colnames(cor_df) <- c(group_col, subgrp_col, "disconcor", "invalid", "concor")
    cor_df %>%
        mutate(total = disconcor + concor) %>%
        group_by_(group_col) %>%
        summarize(concor = sum(concor, na.rm = TRUE),
                  disconcor = sum(disconcor, na.rm = TRUE),
                  total = sum(total, na.rm = TRUE),
                  kendall = (concor - disconcor)/total)
}

# http://stackoverflow.com/a/12135122/2320823
#' @export
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)

#' @export
substrRight <- function(x, n){
    substr(x, nchar(x) - n + 1, nchar(x))
}

#' @export
perc.rank <- function(x) trunc(rank(x))/length(x)

#' @export
perc_rank_single <- function(x, xo)  length(x[x <= xo])/length(x)*100

#' @export
rem_extrema <- function(x, max = TRUE, min = TRUE) {
    x <- data.frame(V1 = x, V2 = x)
    if (max & min) {
        x[x$V1 != min(x$V1) & x$V1 != max(x$V1), ]$V1
    }
    else if (max) {
        x[x$V1 != max(x$V1), ]$V1
    }
    else if (min) {
        x[x$V1 != min(x$V1), ]$V1
    }
    else {
        x
    }
}

#' @export
getname <- function(vector, token = 1) {
    vector <- as.character(vector)
    unlist(lapply(vector, function(x) strsplit(x, "[[:space:]]|(?=[|])", perl = TRUE)[[1]][token]))
}

#' @export
fd_bw <- function(x) {
    2 * IQR(x) * length(x) ^ (-1/3)
}

#' @export
level_ranks <- function(x, extra_starting_level = FALSE) {
    x <- as.numeric(factor(x))
    x_levels <- unique(x)

    if (extra_starting_level) {
        x_levels <- c(min(x_levels) - 1, x_levels)
    }
    map_df <- data.frame(old = x_levels,
                         new = percent_rank(x_levels))
    map_df$new[match(x, map_df$old)]
}

# modified from from asbio::ConDis.matrix
# http://www.inside-r.org/packages/cran/asbio/docs/ConDis.matrix
ConDis.matrix2 <- function(Y1, Y2)
{
    n <- length(Y1)
    foreach(i = 1:n, .combine = c) %do% {
        # changed bounds from 1:n to i:n
        foreach(j = i:n, .combine = c) %do% {
            ifelse((Y1[i] > Y1[j] & Y2[i] > Y2[j]) | (Y1[i] < Y1[j] & Y2[i] < Y2[j]), 1,
                   ifelse((Y1[i] < Y1[j] & Y2[i] > Y2[j]) | (Y1[i] > Y1[j] & Y2[i] < Y2[j]), -1, 0))
        }
    }
}

# analytic method (didn't end up using this)
#' @export
spearman_95ci <- function(spearman, nobs) {
    # http://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation
    # http://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_Spearmans_Rank_Correlation.pdf
    base <- atanh(spearman)
    delta <- 1.96/sqrt(nobs - 3)
    c(tanh(base - delta), tanh(base + delta))
}

#' @export
estimate_overlap <- function(A, B, method = 'locfit', n.grid = 2 ^ 13) {
    require(gss)
    # if both values are only positive
    # take the log
    #     if(all(A >= 0) & all(A >= 0)) {
    #         A <- log(A)
    #         B <- log(B)
    #     }

    # if the length of the sets are non-zero, the calculate a KDE
    if (length(A) == 0 | length(B) == 0) {
        return(NA)
    }

    # rangeA <- range(A)
    # rangeB <- range(B)
    # range_all <- c(min(rangeA[[1]], rangeB[[1]]), max(rangeA[[2]], rangeB[[2]]))

    # Calculate Extent
    range_all <- c(min(A, B), max(A, B))
    # Set up mesh
    mesh <- seq(range_all[[1]], range_all[[2]], length = n.grid)

    if (method == 'gss') {
        # fit each density independently, set extent with `domain`
        fitA <- gss::ssden(~A, domain = data.frame(A = range_all))
        fitB <- gss::ssden(~B, domain = data.frame(B = range_all))

        # calculate density at each grid point
        densityA <- gss::dssden(fitA, mesh)
        densityB <- gss::dssden(fitB, mesh)
    } else if (method == 'locfit') {
        # calculate density at each grid point
        densityA <- locfit::density.lf(A, ev = mesh, maxit = 10000)$y
        densityB <- locfit::density.lf(B, ev = mesh, maxit = 10000)$y
    } else if (method == 'logspline') {
        # Doesn't work for edge cases, e.g. freqens
        #       tryCatch({
        fitA <- logspline::logspline(A)
        fitB <- logspline::logspline(B)

        densityA <- logspline::dlogspline(mesh, fitA)
        densityB <- logspline::dlogspline(mesh, fitB)
        #      },  error = function(e) err <- TRUE)
    } else {
        stop("Unrecognized method")
    }
    # normalize density to 1
    densityA <- densityA / sum(densityA)
    densityB <- densityB / sum(densityB)
    # calculate overall overlap
    sum(pmin(densityA, densityB))
}

#' @export
calc_auc_roc <- function(data, grouping_col, outcome_col, score_col = "ml21", method = "auc") {
    require(foreach)

    foreach(this_group = data[, grouping_col] %>% unique, .combine = rbind) %do% {
        this_subset <- data[data[, grouping_col] == this_group, ]

        # within each grouping and go through each activity
        foreach(this_threshold = this_subset[, outcome_col] %>% unique %>% rem_extrema(max = FALSE), .combine = rbind) %dopar% {
            # take all activities greater than the current one as true
            this_response <- this_subset[, outcome_col] >= this_threshold
            this_perc_rank <- perc_rank_single(this_subset[, outcome_col], this_threshold)

            results_df <- data.frame(grouping_col = grouping_col,
                                     group = this_group,
                                     outcome_threshold = this_threshold,
                                     outcome_percentile = this_perc_rank,
                                     stringsAsFactors = FALSE)

            if (method == "auc") {
                ci_obj <- ci.auc(response = this_response,
                                 predictor = this_subset[, score_col],
                                 direction = "<",
                                 method = "delong")
                data.frame(
                    results_df,
                    count = nrow(this_subset),
                    pos_count = sum(this_response),
                    lower95 = ci_obj[[1]],
                    auc = ci_obj[[2]],
                    upper95 = ci_obj[[3]],
                    stringsAsFactors = FALSE
                )
            } else if (method == "roc") {
                data.frame(
                    results_df,
                    my_roc(response = this_response,
                           predictor = this_subset[, score_col],
                           direction = "<"),
                    stringsAsFactors = FALSE
                )
            }
        }
    }
}

#' @export
my_roc <- function(...) {
    require(pROC)
    require(dplyr)
    pROC_outcome <- roc(...)
    pROC_outcome$ppv <- foreach(thisthreshold = pROC_outcome$thresholds,
                                .combine = c, .multicombine = TRUE) %dopar% {
        truepos <- sum(pROC_outcome$cases >= thisthreshold)
        predpos <- sum(c(pROC_outcome$cases, pROC_outcome$controls) >= thisthreshold)
        truepos/predpos
    }
    pROC_outcome$npv <- foreach(thisthreshold = pROC_outcome$thresholds,
                                .combine = c, .multicombine = TRUE) %dopar% {
        trueneg <- sum(pROC_outcome$controls <= thisthreshold)
        predneg <- sum(c(pROC_outcome$cases, pROC_outcome$controls) <= thisthreshold)
        trueneg/predneg
    }
    pROC_outcome$threshold_percentiles <- percent_rank(pROC_outcome$thresholds)
    pROC_outcome$min_ppv <- length(pROC_outcome$cases) / (length(pROC_outcome$cases) + length(pROC_outcome$controls))
    pROC_outcome$min_npv <- length(pROC_outcome$controls) / (length(pROC_outcome$cases) + length(pROC_outcome$controls))
    #    pROC_outcome$pos_count <- length(pROC_outcome$cases)
    pROC_outcome[c("thresholds", "threshold_percentiles", "sensitivities","specificities","ppv", "min_ppv", "npv", "min_npv")]
}

# necessary for PPV (e.g. Figure 2B)
#' @export
prep_for_polygon <- function(df, x = "thresholds", y = "ppv", min_y = "min_ppv") {
    # prep_for_polygon(data.frame(thresholds = c(1, 2, -3), ppv = c(1, 5, 10)), min_ppv = 1)
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
        thisdata <- thisdata[order(thisdata[, x], thisdata[,y], decreasing = TRUE),]
        rbind(thisdata[1,], thisdata[2,], thisdata[4,], thisdata[3,])
    }
}


# http://stackoverflow.com/a/34859307/2320823
#' @export
GeomStepHist <- ggproto("GeomStepHist", GeomPath,
                        required_aes = c("x"),

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
#' @export
GeomPolygon2 <- ggproto("GeomPolygon2", Geom,
                        draw_panel = function(data, panel_scales, coord) {
                            n <- nrow(data)
                            if (n == 1) return(zeroGrob())
                            munched <- coord_munch(coord, data, panel_scales)
                            munched <- munched[order(munched$group), ]
                            first_idx <- !duplicated(munched$group)
                            first_rows <- munched[first_idx, ]
                            ggplot2:::ggname("geom_polygon",
                                             polygonGrob(munched$x, munched$y, default.units = "native",
                                                         id = munched$group,
                                                         gp = gpar(
                                                             col = first_rows$colour,
                                                             fill = alpha(first_rows$fill, first_rows$alpha),
                                                             lwd = first_rows$size * .pt,
                                                             lty = first_rows$linetype
                                                         )
                                             )
                            )
                        },
                        default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1,
                                          alpha = NA),
                        handle_na = function(data, params) {
                            data
                        },
                        required_aes = c("x", "y"),
                        draw_key = draw_key_polygon
)

#' @export
GeomViolin2 <- ggproto("GeomViolin", Geom,
                       setup_data = function(data, params) {
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
                                             xmaxv = x + violinwidth * (xmax - x)
                           )
                           newdata <- rbind(
                               plyr::arrange(transform(data, x = xminv), y),
                               plyr::arrange(transform(data, x = xmaxv), -y)
                           )
                           newdata <- rbind(newdata, newdata[1,])
                           if (length(draw_quantiles) > 0) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[
                                   rep(1, nrow(quantiles)),
                                   setdiff(names(data), c("x", "y")),
                                   drop = FALSE
                                   ]
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_violin", grobTree(
                                   GeomPolygon2$draw_panel(newdata, ...),
                                   quantile_grob)
                               )
                           } else {
                               ggplot2:::ggname("geom_violin", GeomPolygon2$draw_panel(newdata, ...))
                           }
                       },
                       draw_key = draw_key_polygon,
                       default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                                         alpha = NA, linetype = "solid"),
                       required_aes = c("x", "y")
)

#' @export
geom_violin2 <- function(mapping = NULL, data = NULL, stat = "ydensity",
                         draw_quantiles = NULL, position = "dodge",
                         trim = TRUE, scale = "area",
                         na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                         ...) {
    layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomViolin2,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            trim = trim,
            scale = scale,
            draw_quantiles = draw_quantiles,
            na.rm = na.rm,
            ...
        )
    )
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
preprocess_vars <- function(test_data, xtrans) {
    prediction_features <- names(xtrans$mean)

    # features in the test dataset
    test_data_features <- colnames(test_data)

    # features requested for the ML but not in the input dataset
    # keep in mind setdiff gives the asymmetric difference
    # (http://stat.ethz.ch/R-manual/R-patched/library/base/html/sets.html)
    missing_features <- setdiff(prediction_features, colnames(test_data))

    for (x in missing_features) {
        test_data[[x]] <- NA
    }
    rm(x)

    # keep only the columns we want for the ML (get rid of the extras)
    test_data <- test_data[, prediction_features]

    # do the preprocessing (scaling and centering)
    transformed <- predict(xtrans, test_data)

    transformed[, !(colnames(transformed) %in% missing_features)]
}

#' @export
cor.mtest <- function(mat, conf.level = 0.95){
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
            p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
            lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
            uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
        }
    }
    return(list(p.mat, lowCI.mat, uppCI.mat))
}

#' @export
prediction_fn <- function(test_data, xtrans, weights) {
    require(data.table, quietly = TRUE)
    # In case there were features that were not calculated,
    # e.g. RNAss, add columns with NA for them
    # Contribution will not be taken into account for score calculation
    
    # features we use for the ML
    prediction_features <- names(weights)
    
    # features requested for the ML but not in the input dataset
    # keep in mind setdiff gives the asymmetric difference
    # (http://stat.ethz.ch/R-manual/R-patched/library/base/html/sets.html)
    missing_features <- setdiff(prediction_features, colnames(test_data))
    
    for (x in missing_features) {
        test_data[[x]] <- NA
    }
    
    # keep only the columns we want for the ML (get rid of the extras)
    test_data <- test_data[, prediction_features]
    
    # do the preprocessing (scaling and centering)
    if (!is.null(xtrans)) {
        library(caret)
        test_data <- predict(xtrans, test_data)
    }
    
    # remove features without a value
    test_data[is.na(test_data)] <- 0
    
    # reorder vector
    weights <- weights[colnames(test_data)]
    for (i in seq_along(test_data))
        set(test_data, j = i, value = test_data[[i]] * weights[[i]])
    
    # prediction
    test_data %>% rowSums(na.rm = TRUE)
}

#' @export
read_svmlight_model <- function(filename, feat_names) {
    weights <- read_delim(filename, skip = 11, delim = " ", col_names = FALSE,
                          col_types = cols(.default = col_character(),
                                           X1 = col_integer())) %>%
        select(-X1) %>%
        gather(id, weight) %>%
        filter(weight != "#") %>%
        separate(weight, into = c("idx", "weight"), sep = ":") %>%
        mutate(weight = as.numeric(weight)) %>%
        select(-id, -idx) %>% unlist
    names(weights) <- feat_names
    weights
}


# Adapted from klaR
# http://www.inside-r.org/packages/cran/klaR/docs/svmlight
#' @export
svmlight.file <- function(x, train = FALSE, ...)
{
    if (is.vector(x)) x <- t(x)
    erg <- x
    sn <- 1:nrow(x)
    if (!train) erg[sn, 1] <- paste("1:", x[sn, 1], sep = "")
    if (ncol(x) > 1) {
        j <- 2:ncol(x)
        erg[ , -1] <- matrix(paste(j - train, t(x[,j]), sep = ":"),
                             ncol = ncol(x) - 1, byrow = TRUE)
    }
    return(erg)
}
