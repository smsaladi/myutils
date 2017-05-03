# *************************************************
#                     Setup
# *************************************************

#' @importFrom cowplot theme_cowplot

.onAttach <- function(libname, pkgname) {
    # switch the default theme to my defaults
    theme_set(
        theme_cowplot(font_size = 7, line_size = 0.25) +
            theme(line = element_line(lineend = "round"),
                  legend.title = element_blank(),
                  legend.key = element_blank(),
                  legend.position = "None",
                  plot.margin = unit(c(0, 0, 0, 0), "in")
                  )
        )
}
