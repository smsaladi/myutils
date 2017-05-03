#' Used for model training
#'
#' groups should be specified within the dataframe provided
#' using dplyr (i.e., group_by(df))
#' @export
write_dataset <- function(df, filename, outcome_col = "outcome") {
    data_fn <- paste0(filename, '.svmlight')

    df$qid <- group_indices(df)

    group_cols <- groups(df) %>%
        sapply(deparse)

    df %<>%
        ungroup %>%
        arrange(qid) %>%
        mutate(qid = paste0("qid:", qid)) %>%
        select(-one_of(group_cols))

    save_cols <- df %>%
        select(one_of(outcome_col, "qid"))

    df %>%
        ungroup %>%
        # remove outcome and grouping cols
        select(-one_of(outcome_col, "qid")) %>%
        # cast all into numeric
        as.matrix %>%
        # transform into SVMlight format
        svmlight.file %>%
        as_tibble %>%
        bind_cols(save_cols, .) %>%
        write.table(file = data_fn, row.names = FALSE,
                    col.names = FALSE, quote = FALSE)

    # remove missing feature NA's
    c("perl -p -i -e 's/ \\d+:NA//g'",
      "perl -p -i -e 's/ \\d+:NaN//g'") %>%
        paste(data_fn) %>%
        lapply(system)

    data_fn
}

#' @export
count_uniq <- function(x)
    length(unique(x))

#' @export
replace_all_na <- function(df, repl = 0) {
    df[is.na(df)] <- repl
    df
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

# http://stats.stackexchange.com/a/143447/62183
nclass.FD <- function(x) {
    h <- stats::IQR(x)
    if (h == 0)
        h <- stats::mad(x, constant = 2)
    if (h > 0)
        ceiling(diff(range(x))/(2 * h * length(x) ^ (-1/3)))
    else
        1L
}

ppv_from_sens <- function(pos, neg, sens, spec) {
    A = pos * sens
    D = neg * spec
    B = neg - D
    A / (A + B)
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

# function to split groups in a grouped_df
#' @export
split_groups <- function(df) {
        do(df, {
            df_sub <- .
            groups(df) %>%
                lapply(function(x) distinct_(df_sub, x)) %>%
                bind_cols %>%
                bind_rows(df_sub, .)
        })
}

#' returns concordant - disconcordant pairs
#' number of pairs
#' @export
ConDis_fast <- function(Y1, Y2) {
    Y1 <- unlist(Y1, use.names = FALSE)
    Y2 <- unlist(Y2, use.names = FALSE)
    # pairwise complete obs
    keep <- !(is.na(Y1) | is.na(Y2))
    Y1 <- Y1[keep]
    Y2 <- Y2[keep]
    # This is how SVMrank would have calculated kendall's tau
    pairs <- DescTools::ConDisPairs(table(Y1, Y2))
    count <- pairs$C + pairs$D
    data_frame(valid_count = count,
               con = pairs$C,
               dis = pairs$D)
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
        require(gss)
        # fit each density independently, set extent with `domain`
        fitA <- gss::ssden(~A, domain = data.frame(A = range_all))
        fitB <- gss::ssden(~B, domain = data.frame(B = range_all))

        # calculate density at each grid point
        densityA <- gss::dssden(fitA, mesh)
        densityB <- gss::dssden(fitB, mesh)
    } else if (method == 'locfit') {
        require(locfit)
        # calculate density at each grid point
        densityA <- locfit::density.lf(A, ev = mesh, maxit = 10000)$y
        densityB <- locfit::density.lf(B, ev = mesh, maxit = 10000)$y
    } else if (method == 'logspline') {
        require(logspline)
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
calc_auc_roc <- function(df, method = "auc") {
    df %>%
        mutate(count = n()) %>%
        distinct(outcome, count) %>%
        mutate(outcome_percentile = percent_rank(outcome)*100) %>%
        filter(outcome != min(outcome)) %>%
        rowwise() %>%
        # go through activities
        do({
            thresh <- .
            # take all activities greater than the current one as true
            this_response <- df$outcome >= thresh$outcome
            this_score <- df$score

            if (method == "auc") {
                thresh %<>%
                    as_tibble %>%
                    mutate(pos_count = sum(this_response),
                           neg_count = sum(!this_response))

                # no positive cases
                if (thresh$pos_count == 0 || thresh$neg_count == 0)
                    return(thresh)

                roc_obj <- roc(response = this_response,
                               predictor = this_score,
                               direction = "<")

                bestCI <- c("threshold", "specificity", "sensitivity", "ppv") %>%
                    coords(roc_obj, "best", ret = .)
                power <- power.roc.test(roc_obj)$power

                thresh %<>%
                    mutate(auc = as.numeric(roc_obj$auc),
                           best_thresh = tryCatch({bestCI[["threshold"]]},
                                                  error = function(x) NA),
                           best_spec = tryCatch({bestCI[["specificity"]]},
                                                error = function(x) NA),
                           best_sens = tryCatch({bestCI[["sensitivity"]]},
                                                error = function(x) NA),
                           power = power)

                # only one case
                if (thresh$pos_count == 1 || thresh$neg_count == 1)
                    return(thresh)

                ci_obj <- ci.auc(roc_obj, method = "delong")

                thresh %>%
                    mutate(lower95 = ci_obj[[1]],
                           upper95 = ci_obj[[3]])
            } else if (method == "roc") {
                my_roc(response = this_response,
                       predictor = this_score,
                       direction = "<") %>%
                    mutate(outcome = thresh$outcome,
                           outcome_percentile = thresh$outcome_percentile)
            }
        }) %>%
        ungroup
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
    pROC_out <- roc(...)

    pROC_out %$%
        data_frame(sensitivities, specificities, thresholds) %>%
        rowwise() %>%
        mutate(truepos = sum(pROC_out$cases >= thresholds),
               predpos = sum(c(pROC_out$cases, pROC_out$controls) >= thresholds),
               trueneg = sum(pROC_out$controls <= thresholds),
               predneg = sum(c(pROC_out$cases, pROC_out$controls) <= thresholds)) %>%
        ungroup() %>%
        arrange(thresholds) %>%
        mutate(threshold_percentiles = percent_rank(thresholds),
               ppv = truepos/predpos,
               min_ppv = ppv[[1]],
               npv = trueneg/predneg,
               min_npv = npv[[1]]) %>%
        select(-truepos, -predpos, -trueneg, -predneg)
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

#' @export
great_eq <- function(cond, data) {
    if (is.na(data)) {
        is.na(cond)
    } else {
        data >= cond
    }
}

#' @export
less_eq <- function(cond, data) {
    if (is.na(data)) {
        is.na(cond)
    } else {
        data <= cond
    }
}

#' @export
eq_fun <- function(cond, data) {
    if (is.na(cond)) {
        is.na(data)
    } else {
        data == cond
    }
}

#' @export
ne_fun <- function(cond, data) {
    if (any(is.na(cond)) | any(is.na(data))) {
        stop("NA handling not implemented")
    }
    data != cond
}
