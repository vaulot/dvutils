#' @import RMySQL
#' @import RSQLite
#' @import duckdb
#' @import DBI
#' @import dbplyr
#' @import stringr
#' @import yaml


# Establish a connection to a SQLite database --------------------------------------------------
#' @title Establish a connection to a SQLite database
#' @description Establish a connection to a SQLite database
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

# Establish a connection to a database --------------------------------------------------
#' @title Establish a connection database
#' @description
#' The **db_config** must be read before using read_db_config
#' @param db_name name of the database ("pr2", "metapr2", "rcc")
#' @return
#' A database connection
#' @examples
#' db_config <- read_db_config("C:/daniel.vaulot/Databases/google_pc.yaml") # Must always use db_config
#' pr2_con <- db_connect("pr2")
#' @export
#' @md

db_connect <- function(db_name)  {

  if (db_config$db_type == "mysql") {db <- dbConnect(MySQL(),
                                                   default.file=db_config$db_config_file,
                                                   groups=db_config$db_groups,
                                                   dbname=db_name) }

  # if (db_config$db_type == "sqlite") {db <- dbConnect(SQLite(),  db_config$db_file) }

  if (db_config$db_type == "duckdb") {
      if (db_name == "pr2") dbdir = db_config$db_file_pr2
      if (db_name == "metapr2") dbdir = db_config$db_file_metapr2
      db <- dbConnect(duckdb(), dbdir = dbdir, read_only = FALSE)
      }

  return(db)

}

# Get the connection info for a specific database --------------------------------------------------
#' @title Get the connection info from a yaml file
#' @description
#' To be used with db_connect (see the example).
#' The MySQL connection info are stored in a file (my.cnf). See https://github.com/r-dbi/RMySQL/issues/131
#' @param db_config_yaml  The yaml file containing the info
#' @return
#' A list with the database connection info
#'  - db_type: mysql, duckdb ...
#'  - db_config$db_config_file: for MySQL, the location of my.cnf
#'  - db_config$db_groups: for MySQL
#'  - db_file_pr2: for duckdb
#'  - db_file_metapr2: for duckdb
#' @examples
#' db_config <- read_db_config("C:/daniel.vaulot/Databases/google_pc.yaml") # Must always use db_config
#' db_config <- read_db_config(db_config_yaml="C:/daniel.vaulot/Databases/duckdb_pc.yaml")
#' db_connect("pr2")
#' @export
#' @md

read_db_config <- function(db_config_yaml = "C:/daniel.vaulot@gmail.com/Databases/google_pc.yaml") {

  db_config <- yaml::read_yaml(db_config_yaml)

  return(db_config)

}

# Disconnect a connection to a database --------------------------------------------------
#' @title Disconnect an existing connection
#' @description Disconnect an existing connection
#' @param db_con Database connection
#' @examples
#' db_disconnect(pr2_con)
#' @export
#' @md

db_disconnect <- function(db_con)  {
    dbDisconnect(db_con)
    return(TRUE)
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

db_get_query <- function(database, query)  {
    db <- db_connect(database)
    table <- dbGetQuery(db,query)
    dbDisconnect(db)
    return(table)
}


