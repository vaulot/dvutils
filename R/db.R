#' @import RMySQL
#' @import RSQLite
#' @import DBI
#' @import dbplyr
#' @import stringr

#' @title Put a string between single quotes
#' @description
#' Replaces the dbplyr escape function that returns vectors
#' @param sql_value String or value to be put between single quotes
#' @return
#' A string between single quotes
#' @examples
#' db_sql_escape("Eukaryota")
#' @export
#' @md

db_sql_escape <- function(sql_value) {
  str_c("'", sql_value,"'")
}


# Build a string to update a field in a database  -----------------------------------

#' @title Build a string to update a field in a database
#' @description
#' @param field_name Character
#' @param field_value Numeric or string
#' @param append_comma TRUE  append a ", " if there more than one field value
#'
#' @return
#' A string in the format "field_name = 'field_value'"
#' @examples
#' query <- update_field_string("pr2_accession", "ZZ220100", append_comma=FALSE)
#' @export
#' @md

db_update_field_string <- function(field_name, field_value, append_comma=TRUE) {

  single_quote <- "'"

  # Convert first to character
  field_value <- as.character(field_value)

  # Remove any ' occuring in any of the fields which will generate an error
  field_value <- str_replace_all(field_value, "'","")

  # Convert everything to ASCII
  field_value <- iconv(field_value, "latin1", "ASCII", sub="?")

  # Check whether the value is empty and then return empty string
  field_string <- case_when(is.na(field_value) ~ "",
                            TRUE ~ str_c(field_name," = ",
                                         single_quote,field_value,single_quote,
                                         if_else(append_comma,", ", "")))

  # If use escape() it convert NA to NULL so no need to check
  # However, we cannot use escape because it collapses the vector into a list...
  field_string <- str_replace_all(field_string, "'NA'","NULL")

  # if (append_comma) {field_string <- str_c(field_string, ", ")}

  return(field_string)
}

# Establish a connection to a SQLite database --------------------------------------------------
#' @title Establish a connection to a SQLite database
#' @description
#' @param file_name Character, name of the file to store the database
#' @return
#' A database connection
#' @examples
#' pr2_con <- db_connect("PR2.sqlite")
#' @export
#' @md

db_connect_sqlite <- function(file_name)  {

db <- dbConnect(SQLite(),  file_name)
    return(db)
}

# Establish a connection to a MySQL database --------------------------------------------------
#' @title Establish a connection database
#' @description
#' All connections parameters are stored into a my.conf file
#' @param db_info list(default.file= ,groups= ,dbname= )
#' @return
#' A database connection
#' @examples
#' pr2_con <- db_connect(db_info("pr2_local"))
#' test_con <- db_connect(list(default.file="my.cnf", groups="root" , dbname= "test"))
#' @export
#' @md

db_connect <- function(db_info)  {

db <- dbConnect(MySQL(),  default.file=db_info$default.file,
                          groups=db_info$groups,
                          dbname=db_info$dbname)

    return(db)
}

# Get the connection info for a specific database --------------------------------------------------
#' @title Get the connection info for a specific database
#' @description
#' To be used with db_connect (see the example).
#' For the database on Scrol server, one can go though SSH (e.g. when using EduRoam or at NTU).
#' However at present this connection does not work.
#' It is first necessary to set up a tunnel through Putty (see https://blog.devolutions.net/2017/4/how-to-configure-an-ssh-tunnel-on-putty).
#' Same thing for SSL (Google database).
#' The connection info are now stored in a file (my.cnf). See https://github.com/r-dbi/RMySQL/issues/131
#' @param database_name One of the following choices "rcc", "pr2_google", "metapr2_google", "rcc_local", "pr2_local", "metapr2_local"
#' @param file_cnf by default "C:/daniel.vaulot@gmail.com/Databases/MySQL/my.cnf" but can be moved to another location (e.g. when using dvutils on a server)
#' @return
#' A list with the database connection info "user","password","dbname","host"
#' @examples
#' db_connect(db_info("rcc"))
#' @export
#' @md

db_info <- function(database_name,
                    file_cnf = "C:/daniel.vaulot@gmail.com/Databases/MySQL/my.cnf")  {



# To do, try https://theautomatic.net/2019/06/25/how-to-hide-a-password-in-r-with-the-keyring-package/

  db_rcc       <- list( dbname='rcc_database',
                        groups="scrol",
                        default.file=file_cnf)
  db_pr2_google <- list(dbname='pr2',
                        default.file=file_cnf,
                        groups="google")
  db_metapr2_google <- list(dbname='metapr2',
                        default.file=file_cnf,
                        groups="google")
  db_rcc_local <- list(dbname='rcc_local',
                       groups="local",
                       default.file=file_cnf)
  db_pr2_local <- list(dbname='pr2_local',
                       groups="local",
                       default.file=file_cnf)
  db_metapr2_local <- list(dbname='metapr2_local',
                       groups="local",
                       default.file=file_cnf)

  dbinfo <- NULL
  if (database_name == "rcc") {db_info=db_rcc}
  if (database_name == "pr2_google") {db_info=db_pr2_google}
  if (database_name == "metapr2_google") {db_info=db_metapr2_google}
  if (database_name == "rcc_local") {db_info=db_rcc_local}
  if (database_name == "pr2_local") {db_info=db_pr2_local}
  if (database_name == "metapr2_local") {db_info=db_metapr2_local}

    return(db_info)
}

# Disconnect a connection to a database --------------------------------------------------
#' @title Disconnect an existing connection
#' @description
#' @param db_con Database connection
#' @examples
#' db_disconnect(pr2_con)
#' @export
#' @md

db_disconnect <- function(db_con)  {
    dbDisconnect(db_con)
    return(TRUE)
}



# Write new records to database --------------------------------------------------

#' Write new records to a  table.
#'
#' The data frame must have the same fields than the table (previous name = db_write)
#' @inheritParams db_connect
#' @param table_name Character - Name of the table
#' @param table_df Dataframe - Must have the same columns than the database table
#' @return TRUE if successful
#' @export
#'
#' @examples
#' db_append_records(db_con, "taxo", taxo_new_records)
#' @md

db_append_records <- function(database_info, table_name, table_df)  {
    db <- dbConnect(MySQL(),  user=database_info$user,
                              password=database_info$password,
                              dbname=database_info$dbname,
                              host=database_info$host)
    success <- dbWriteTable(db,table_name, table_df, row.names=FALSE, append=TRUE, nrows=(nrow(table_df)) )
    dbDisconnect(db)
    return(success)
}
# Read query from database --------------------------------------------------
#' @title Read a query into a dataframe.
#'
#' @inheritParams db_connect
#' @param query Character - SQL Query
#' @return
#' The value returned by the query
#'
#' @examples
#' taxo <- db_get_query(pr2_db, taxo_query)
#' @md
#' @export

db_get_query <- function(database_info, query)  {
    db <- dbConnect(MySQL(),  user=database_info$user,
                              password=database_info$password,
                              dbname=database_info$dbname,
                              host=database_info$host)
    table <- dbGetQuery(db,query)
    dbDisconnect(db)
    return(table)
}


# Execute query vector from database --------------------------------------------------
#' Execute a vector of query.
#'
#' @inheritParams db_connect
#' @param query Character vector - Queries
#' @return
#' A vector with the number of rows affected by each query
#' @export
#'
#' @examples
#' @md
#' db_execute_query_vector(pr2_db, taxo_query_vector)

db_execute_query_vector <- function(database_info, query_vector)  {
    db <- RMySQL::dbConnect(MySQL(),  user=database_info$user,
                              password=database_info$password,
                              dbname=database_info$dbname,
                              host=database_info$host)
	query_result <- vector(mode="numeric", length=0)
    for (i in 1:length(query_vector)) {
		# print(str_c("i = ", i, "query = ",query_vector[i]) )
		query_result[i] <- dbExecute(db,query_vector[i] )
	}
	return(query_result)
    dbDisconnect(db)
}

