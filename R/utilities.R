# A set of functions written to make data munging and preparation easier
#' @import magrittr
#' @import tidyverse
#' @import caret
#' @importFrom pROC roc ci.auc
#' @importFrom foreach foreach
#' @importFrom grid polygonGrob

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
            write.csv(df_to_write, file = data_fn, row.names = FALSE,
                      quote = FALSE)
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
                          pairs <- ConDis.matrix2(thisgroup[, score_col] %>%
                                                      unlist,
                                                  thisgroup[, outcome_col] %>%
                                                      unlist)
                          data.frame(subgrp_col = thisgrouping,
                                     group_col = unique(thisgroup[, group_col] %>%
                                                            unlist),
                                     table(pairs))
                      } %>%
        dcast(group_col + subgrp_col ~ pairs, value.var = "Freq")
    colnames(cor_df) <- c(group_col, subgrp_col,
                          "disconcor", "invalid", "concor")
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
specify_decimal <- function(x, k = 1L) {
    format(round(x, k), nsmall = k)
}

#' @export
substrRight <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
}

#' @export
perc.rank <- function(x) trunc(rank(x))/length(x)

#' @export
rel_diff <- function(mut, wt, func = range) {
    2*(mut - wt)/abs((wt + mut))
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
            ifelse((Y1[i] > Y1[j] & Y2[i] > Y2[j]) |
                       (Y1[i] < Y1[j] & Y2[i] < Y2[j]), 1,
                   ifelse((Y1[i] < Y1[j] & Y2[i] > Y2[j]) |
                              (Y1[i] > Y1[j] & Y2[i] < Y2[j]), -1, 0))
        }
    }
}

# analytic method (didn't end up using this)
# http://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation
# http://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_Spearmans_Rank_Correlation.pdf
#' @export
spearman_95ci <- function(spearman, nobs) {
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

rem_extrema <- function(x, max = TRUE, min = TRUE) {
    tmp <- x
    if (max & min) {
        x[tmp != min(tmp) & tmp != max(tmp)]
    } else if (max) {
        x[tmp != max(tmp)]
    } else if (min) {
        x[tmp != min(tmp)]
    } else {
        x
    }
}

perc_rank_single <- function(x, xo) {
    length(x[x <= xo])/length(x)*100
}

#' @export
calc_auc_roc <- function(data, grouping_col, outcome_col = "outcome",
                         score_col = "score", method = "auc") {
    foreach(this_group = select_(data, grouping_col) %>% distinct %>% unlist,
            .combine = bind_rows, .multicombine = TRUE) %do% {
        this_subset <- data[data[, grouping_col] %>% unlist == this_group, ]

        # within each grouping and go through each activity
        foreach(this_threshold = this_subset[, outcome_col] %>%
                    unlist %>% unique %>% rem_extrema(max = FALSE),
                .combine = bind_rows, .multicombine = TRUE) %dopar% {
            # take all activities greater than the current one as true
            this_response <- 
                this_subset[, outcome_col] %>% unlist >= this_threshold
            this_perc_rank <- perc_rank_single(
                this_subset[, outcome_col] %>% unlist, this_threshold)

            results_df <- tibble(grouping_col = grouping_col,
                                 group = this_group,
                                 outcome_threshold = this_threshold,
                                 outcome_percentile = this_perc_rank)

            if (method == "auc") {
                ci_obj <- ci.auc(
                    response = this_response,
                    predictor = this_subset[, score_col] %>% unlist,
                    direction = "<",
                    method = "delong")
                results_df %>% 
                    mutate(count = nrow(this_subset),
                           pos_count = sum(this_response),
                           lower95 = ci_obj[[1]],
                           auc = ci_obj[[2]],
                           upper95 = ci_obj[[3]])
            } else if (method == "roc") {
                results_df %>%
                    data.frame(my_roc(response = this_response,
                               predictor = this_subset[, score_col] %>% unlist,
                               direction = "<"))
            }
        }
    }
}

#' @export
find_closest <- function(col, vals) {
    keep <- rep(FALSE, length(col))
    for (i in vals)
        keep <- keep | col == col[[which.min(abs(col - i))]]
    keep
}

#' @export
my_roc <- function(...) {
    require(pROC)
    require(dplyr)
    pROC_outcome <- roc(...)
    pROC_outcome$ppv <- foreach(thisthreshold = pROC_outcome$thresholds,
                                .combine = c, .multicombine = TRUE) %dopar% {
        truepos <- sum(pROC_outcome$cases >= thisthreshold)
        predpos <-
            sum(c(pROC_outcome$cases, pROC_outcome$controls) >= thisthreshold)
        truepos/predpos
    }
    pROC_outcome$npv <- foreach(thisthreshold = pROC_outcome$thresholds,
                                .combine = c, .multicombine = TRUE) %dopar% {
        trueneg <- sum(pROC_outcome$controls <= thisthreshold)
        predneg <-
            sum(c(pROC_outcome$cases, pROC_outcome$controls) <= thisthreshold)
        trueneg/predneg
    }
    pROC_outcome$threshold_percentiles <- percent_rank(pROC_outcome$thresholds)
    pROC_outcome$min_ppv <- length(pROC_outcome$cases) /
        (length(pROC_outcome$cases) + length(pROC_outcome$controls))
    pROC_outcome$min_npv <- length(pROC_outcome$controls) /
        (length(pROC_outcome$cases) + length(pROC_outcome$controls))
    #    pROC_outcome$pos_count <- length(pROC_outcome$cases)
    pROC_outcome[c("thresholds", "threshold_percentiles", "sensitivities",
                   "specificities", "ppv", "min_ppv", "npv", "min_npv")]
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
