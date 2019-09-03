#' @import dplyr

# Function to plot a treemap based on  graphics -------------------------------------
#' @title Draw a simple treemap from pr2.
#' @description
#' Can also be used for any data frame with same taxo structure
#'
#' @param pr2 data frame - can be the whole database or a filtered version of it.
#' It needs at three columns, 2 taxonomy columns and one for example with sequences
#' @param taxo_rank integer between 1 and 7 - defines the first grouping level, level i-1 is used also (1=kingdom etc...)
#'
#' @examples
#' # Will plot at the order and family levels
#' pr2_treemap (pr2_filtered, 5)
#' @export
pr2_treemap <- function(pr2, taxo_rank) {

  # Define the levels
  level1 = pr2.env$taxo_levels[taxo_rank]
  level2 = pr2.env$taxo_levels[taxo_rank+1]
  # Group
  pr2_class <- pr2 %>%
                group_by_(level1,level2) %>%
                summarise(sequence_number= n())

  # Do a simple treemap
  treemap::treemap(pr2_class, index=c(level1,level2),
        vSize="sequence_number",
        title="",asp=1, lowerbound.cex.labels= 0.2, fontsize.labels = 12,
        palette="Blues",
        format.legend = list(scientific = FALSE, big.mark = " "))
}



