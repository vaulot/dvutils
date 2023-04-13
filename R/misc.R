#' @import stringr


# Function to convert between Unix and Dos text format --------------------
#' @title Convert between Unix and Dos text format
#' @description
#' Read the Unix text file and writes back to DOS format
#' @return
#' TRUE if successful
#' @examples
#' file_unix2dos("my_dos.txt")
#' @export
#' @md
file_unix2dos <- function (filename) {

  big_text <- readr::read_file(filename)
  big_text <- stringr::str_replace_all(big_text, "\n$", "\r\n")
  readr::write_file(big_text, filename)
  remove(big_text)  # cleaning up
  return(TRUE)
}


# filename_append -----------------------------------------------------

#' @title Append a string at the end of a file name before the exension
#' @description
#' The extension is unchanged
#' @return
#' New file name
#' @examples
#' # Will result in "mfile.pr2.txt",
#' new_name <- filename_append("mfile.txt",".pr2")
#' @export
#' @md

filename_append <- function (filename, new_end){
  file_name_out <- str_c(fs::path_ext_remove(filename),
                                  new_end,
                                  '.',
                                  fs::path_ext(filename))
  }

# filename_change_ext -----------------------------------------------------

#' @title Change the extension of the file name
#' @description
#' @return
#' New file name
#' @examples
#' # Will result in "mfile.pr2",
#' new_name <- filename_change_ext("mfile.txt","pr2")
#' @export
#' @md

filename_change_ext <- function (filename, new_ext){
  file_name_out <- str_c(fs::path_ext_remove(filename),
                                  '.',
                                  new_ext)
  }


# Fix a bib file  --------------------
#' @title Fix bib library file created by Mendeley
#' @description
#' Does the following :
#'
#' * Transforms <i>...</i> tag to textit
#' * Replace "<="  (found in the title of some papers) by Tex symbol (leq)
#' * Put Vaulot in bold -  if dv_bold=TRUE
#'
#' Note: the bib file must be exported using the "Escape Latex character" setting OFF
#' @return
#' TRUE if successful
#' @examples
#' latex_fix_bibfile(dv_bold=TRUE)
#' @md
#'
latex_fix_bibfile <- function (filename="C:/Daniel/Paper reprints Bibtex/library.bib",dv_bold=FALSE) {

  big_text <- readr::read_file(filename)

  # Transform italics

#The following lines are for the case when the bib file has been exported with escaping latex character
  #    big_text <- str_replace_all(big_text, "\\{\\\\textless\\}i\\{\\\\textgreater\\}", "\\\\textit\\{")
  #    big_text <- str_replace_all(big_text, "\\{\\\\textless\\}/i\\{\\\\textgreater\\}", "\\}")

# When export without escaping latex
  big_text <- stringr::str_replace_all(big_text, "<i>([^<>]+)</i>", "\\\\textit\\{\\1\\}")

  # Look for <= sign
    big_text <- stringr::str_replace_all(big_text, "\u2264", "$\\\\leq$")
  # Look for Vaulot, D and make it bold
    if (dv_bold)
      {big_text <- stringr::str_replace_all(big_text, "Vaulot, D\\.|Vaulot, Daniel", "\\\\textbf\\{Vaulot, D\\.\\}")}

  readr::write_file(big_text, filename)
  remove(big_text)  # cleaning up

  return(TRUE)

}

# Unfill a table --------------------
#' @title Remove repeated values in vector
#' @description
#' Repeated values are replaced by NA. This is useful for formatting tables for papers.
#' @return
#' Vector with repeated values replaced by NA
#' @examples
#' x <- c("A","A","A","B","B","C","C")
#' unfill(x)
#' \dontrun{
#' mutate(df, across(c("division", "class"), ~ unfill(.x)))
#' }
#' @export
#' @md

unfill <- function(x) {
  same <- x == dplyr::lag(x)
  ifelse(!is.na(same) & same, NA, x)
}
