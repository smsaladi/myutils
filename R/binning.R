#' Cut Numeric Values Into Evenly Distributed Groups (Bins)
#'
#' \code{bins} - Cuts points in vector x into evenly distributed groups (bins).
#' \code{bins} takes 3 separate approaches to generating the cuts, picks the one
#' resulting in the least mean square deviation from the ideal cut -
#' \code{length(x) / target.bins} points in each bin - and then  merges small bins
#' unless excat.groups is \code{TRUE}
#'
#' The gains are computed using incremental analytical expresions derived for moving
#' a value from one bin to the next, splitting a bin into two or merging two bins.
#' @param x Vector of numbers
#' @param target.bins Number of groups desired; this is also the max number of groups.
#' @param max.breaks   Used for initial cut. If \code{exact.groups} is \code{FALSE}, bins are merged
#'                     until there's no bins with fewer than \code{length(x) / max.breaks} points.
#'                     In \code{bins}, one of \code{max.breaks} and \code{minpts} must be supplied.
#' @param exact.groups if TRUE, the result will have exactly the number of target.bins;
#'                     if FALSE, the result may contain fewer than target.bins bins
#' @param verbose      Indicates verbose output.
#' @param errthresh    If the error is below the provided value, stops after the first rough estimate of the bins.
#' @param minpts       Minimum number of points in a bin.
#'                     In \code{bins}, one of \code{max.breaks} and \code{minpts} must be supplied.
#' @return A list containing the following items (not all of them may be present):
#' \itemize{
#'    \item{gain}{ - Error gain obtained as the result of the function call.}
#' }
#'
#' @examples
#' \dontrun{
#'    # Seriously skewed x:
#'    x <- floor(exp(rnorm(200000 * 1.3)))
#'    cuts <- bins(x, target.bins = 10, minpts = 2000)
#'    cuts$breaks <- bins.getvals(cuts)
#'    cuts$binct
#'    #   [0, 0]    [1, 1]    [2, 2]    [3, 3]    [4, 4]    [5, 5]    [6, 7]   [8, 10]
#'    # 129868     66611     28039     13757      7595      4550      4623      2791
#'    #   [11, 199]
#'    # 2166
#'
#'    # Centered x:
#'    x <- rep(c(1:10,20,31:40), c(rep(1, 10), 100, rep(1,10)))
#'    cuts <- bins(x, target.bins = 3, minpts = 10)
#'    cuts$binct
#'    # [1, 10] [20, 20] [31, 40]
#'    #      10      100       10
#' }
#' @seealso \code{\link{binr}}
#' @export
#' @rdname bins
bins <- function(x, ...) {
    UseMethod("bins")
}

#' @export
bins.default <- function(x, target.bins) {
    x_narm <- x[!is.na(x)]
    target.bins <- min(target.bins, length(unique(x_narm)))
    probs <- seq_len(target.bins) / target.bins
    structure(quantile(x, probs = probs[1:(target.bins - 1)], na.rm = TRUE),
              hasna = any(is.na(x)),
              class = "binr")
}

#' @export
bins.data.frame <- function(df, ...) {
    structure(
        lapply(df, function(x, ...) bins(x, ...), ...),
        class = "binr"
    )
}

predict_data_frame <- function(obj, data, labels = TRUE, makewide = TRUE) {
    if (all(names(data) %in% names(obj))) {
        return(bind_cols(lapply(
            names(data),
            function(x)
                as.data.frame(
                    predict_vector(
                        obj[[x]], data[[x]], name = x,
                        labels = labels, makewide = makewide
                        )
                    )
            )))
    } else {
        stop("data.frame has columns without bin cuts")
    }
}

predict_vector <- function(obj, data, labels = TRUE, makewide = TRUE, name = NA) {
    out <- cut(data, c(-Inf, obj, Inf),
               labels = if (labels) NULL else FALSE)
    # Add a level for NAs (if they exist and are in the training data)
    # Otherwise they are left as NA
    if (NA %in% out && attr(obj, 'hasna')) {
        if (labels) {
            levels(out) <- c(levels(out), "[NA]")
            out[is.na(out)] <- "[NA]"
        } else {
            out[is.na(out)] <- max(out, na.rm = TRUE) + 1
        }
    }

    if (makewide) {
        if (is.na(name))
            stop("must provide name if makewide == TRUE")
        # set name of column
        out <- setNames(list(out), name)
        df_tmp <- mutate_all(as.data.frame(out), as.factor)
        out <- predict(dummyVars(~ ., data = df_tmp), newdata = df_tmp)
        # if any NA's remain, it's because obj$hasna is FALSE
        out[is.na(out)] <- 0
    }

    return(out)
}


#' @export
predict.binr <- function(obj, data, labels = TRUE, makewide = TRUE) {
    if (is.data.frame(data)) {
        predict_data_frame(obj, data, labels = labels, makewide = makewide)
    } else {
        predict_vector(obj, data, labels = labels, makewide = makewide)
    }
}
