#' @export
read_ranklib <- function(model_dat_fn, xtrans) {
  structure(model_dat_fn,
            xtrans = xtrans,
            class = "ranklib")
}

#' @export
predict.ranklib <- function(obj, newdata, exec_jar = "~/apps/ranklib/RankLib.jar") {
  require(caret, quietly = TRUE)
  
  # set up files for svmlight
  data_fn <- tempfile("ranklib_", fileext = ".dat")
  out_fn <- paste0(data_fn, ".predictions")
  
  # In case there were features that were not calculated,
  # e.g. RNAss, add columns with NA for them
  # Contribution will not be taken into account for score calculation
  
  # features we use for the ML
  prediction_features <- names(attr(obj, "xtrans")$mean)
  
  # features in the test dataset
  newdata_features <- colnames(newdata)
  
  # features requested for the ML but not in the input dataset
  # keep in mind setdiff gives the asymmetric difference
  missing_features <- setdiff(prediction_features, colnames(newdata))
  
  for (x in missing_features) {
    newdata[[x]] <- NA
  }
  rm(x)
  
  data_fn <- newdata %>%
    # keep only the columns we want for the ML (get rid of the extras)
    select(one_of(prediction_features)) %>%
    # do the preprocessing (scaling and centering)
    predict(attr(obj, "xtrans"), .) %>%
    # add necessary "metadata" columns 
    mutate(outcome = 1,
           grp = "grp") %>%
    group_by(grp) %>%
    # write the file in svmlight format
    write_dataset(data_fn)
  
  # remove features with no value (NA)
  system(paste("perl -p -i -e 's/ \\d+:NA//g'", data_fn))
  system(paste("perl -p -i -e 's/ \\d+:NaN//g'", data_fn))
  
  system2("java",
          args = c("-jar", exec_jar,
                   "-load", obj,
                   "-rank", data_fn,
                   "-score", out_fn,
                   "-thread", "1"),
          stdout = FALSE, stderr = "")
  
  if (file.exists(out_fn)) {
    preds <- read_tsv(out_fn, col_names = c("qid", "idx", "score"),
                      col_types = cols(.default = col_integer(),
                                       score = col_double()))
    
    # remove input and predictions file
    unlink(c(data_fn, out_fn))
    
    preds$score
  } else {
    NA
  }
}
