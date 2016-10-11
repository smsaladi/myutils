# http://www.math.umt.edu/gideon/Kendelemslp.pdf

# computing the two unique vectors with ties present: the function is tauuniq
#' @export
tauuniq <- function(x, y) {
    n <- length(x)
    e <- 1:n
    xrr <- n + 1 - rank(x)
    xtp <- x[order(y, x)]
    xtn <- x[order(y, xrr)]
    rkyp <- order(xtp, e)
    rkyn <- order(xtn, n:1)
    cbind(rkyp, rkyn)
}

# calculation of Kendall's tau on unique max-min vectors
# the function is rtau
rtau <- function(x, y) {
    ot <- tauuniq(x, y)
    rkyp <- ot[, 1]
    rkyn <- ot[, 2]
    dyp <- 0
    dyn <- 0
    n <- length(x)
    n2 <- (n*(n - 1))/2
    n1 <- n - 1
    for (i in 1:n1) {
        j <- i + 1
        tempp <- rkyp[i] - rkyp[j:n]
        tempn <- rkyn[i] - rkyn[j:n]
        dyp <- dyp + sum(tempp < 0)
        dyn <- dyn + sum(tempn < 0)
    }
    (dyp + dyn)/n2 - 1
}
# output is slope and intercept, function name tauslp in positions 1 and 2
tauslp <- function(x, y) {
    rat <- c(outer(y, y, "-")/outer(x, x, "-"))
    ratv <- rat[!is.na(rat)]
    slp <- median(ratv)
    res <- y - slp * x
    aint <- median(res)
    res <- res - aint
    ck <- rtau(x,res)
    ck1 <- sum(res)
    ck2 <- median(res)
    c(slp, aint, ck, ck1, ck2)
}
