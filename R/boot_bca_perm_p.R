# A set of functions written to make data munging and preparation easier
#' @import magrittr
#' @import tidyverse
#' @import boot
#' @import parallel

# summary metric of interest
cor_spear <- purrr::partial(cor, method = "spearman",
                            use = "pairwise.complete.obs")

# single bootstrap for confidence intervals
boot_cor <- function(data, indices, x_col, y_col, SUMMARY_FUN) {
    SUMMARY_FUN(data[indices, x_col] %>% unlist,
                data[indices, y_col] %>% unlist)
}

# single permutation test for p-value significance
shuff_cor_pval <- function(index, data, x_col, y_col, SUMMARY_FUN) {
    base::sample(data[, y_col] %>% unlist,
                 size = nrow(data), replace = FALSE) %>%
        SUMMARY_FUN(data[, x_col] %>% unlist)
}

# summary statistic after all shuffles (p-value)
calc_pval <- function(orig_df, shuf, x_col, y_col, SUMMARY_FUN) {
    stat <- SUMMARY_FUN(orig_df[, x_col] %>% unlist,
                        orig_df[, y_col] %>% unlist)

    # adjustment to account for the original data that does meet the criteria
    perm_p <- (sum(shuf >= stat) + 1) / (length(shuf) + 1)
    names(perm_p) <- "perm_p"
    perm_p
}

#' @export
boot_perm <- function(data, x_col = "score", y_col = "outcome",
                      fn_prefix, count, cores, SUMMARY_FUN = cor_spear,
                      force = FALSE) {
    # check if outcome is saved
    boot_perm_summary_fn <- paste0(fn_prefix, "_cor_summary.rds")
    if (file.exists(boot_perm_summary_fn) & !(force))
        return(readRDS(boot_perm_summary_fn))

    # bootstrapping
    boot_data <- boot(data = data[, c(x_col, y_col)],
                      statistic = boot_cor, R = count,
                      x_col = x_col, y_col = y_col,
                      SUMMARY_FUN = SUMMARY_FUN,
                      parallel = "multicore", ncpus = cores)
    saveRDS(boot_data, file = paste0(fn_prefix, "_boot.rds"))
    boot_ci <- boot.ci(boot_data, type = "bca")
    rm(boot_data)
    gc()

    # permutation testing
    shuf_data_stat <- mclapply(seq(count), shuff_cor_pval,
                               data = data[, c(x_col, y_col)],
                               x_col = x_col, y_col = y_col,
                               SUMMARY_FUN = SUMMARY_FUN,
                               mc.cores = cores) %>% unlist
    saveRDS(shuf_data_stat, file = paste0(fn_prefix, "_perm.rds"))
    perm_p <- calc_pval(data, shuf_data_stat,
                        x_col = x_col, y_col = y_col, SUMMARY_FUN = SUMMARY_FUN)
    rm(shuf_data_stat)
    gc()
    result <- c(boot_ci, perm_p) %$%
        tibble(R = R,
               t0 = t0,
               bca_level = bca[[1]],
               bca_order1 = bca[[2]],
               bca_order2 = bca[[3]],
               bca_low = bca[[4]],
               bca_high = bca[[5]],
               perm_p = perm_p)

    saveRDS(result, boot_perm_summary_fn)
    result
}
