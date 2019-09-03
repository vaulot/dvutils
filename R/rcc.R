#' @import dplyr
#' @import stringr
#' @import readxl
#' @import tidyr


# GenBank processing ------------------------------------------------------

#' @title Prepare list of GenBank entries to include in rcc database
#'
#' @description
#' Input is an Excel file exported from Geneious with three additional columns
#'   * rcc_id
#'   * gene_location
#'   * gene_name
#'
#' Entries already in the rcc database are discarded
#' @param genbank_file Name of excel file exported from Geneious
#'
#' @return
#' Data frame ready for import into rcc database. .
#'
#' @examples
#' df <- rcc_genbank("RCC from Genbank 2017 12 26.xlsx")
#' @md
#' @export

rcc_genbank <- function(genbank_file)  {

      # genbank_file <- "C:/Daniel/Cultures/sequences/GenBank/from GenBank/RCC from Genbank 2017 12 26.xlsx"
      genbank_df <- read_excel(genbank_file)

      # replace all spaces in column names by underscores
      names(genbank_df)<-names(genbank_df)%>% stringr::str_replace_all("\\s","_")

      genbank_df <- genbank_df %>% transmute(genbank_accession = Accession,
                                             genbank_date = Created_Date,
                                             rcc_id=rcc_id ,
                                             gene_name = gene_name,
                                             gene_location = gene_location,
                                             genbank_description = Description,
                                             genbank_taxonomy = Taxonomy,
                                             genbank_organism = Organism,
                                             sequence_length = Sequence_Length)
      genbank_df$date_added <- format(Sys.Date(),'%Y-%m-%d %H:%M:%S')
      genbank_df <- genbank_df %>% filter(!(genbank_accession %in% read_field_db(db_info("rcc"), "sequences", "genbank_accession")))

      return(genbank_df)
}


# Orders ------------------------------------------------------------------

#' @title Format "orders.csv" for import into rcc database
#'
#' @description
#' This function reads the file "orders.csv" exported from the web site
#' and import it into the rcc database
#' @param orders_file Name of file "orders.csv" exported from the web site
#'
#' @return
#' Data frame ready for import into rcc database
#'
#' @examples
#' df <- rcc_orders("orders.csv")
#' @export

rcc_orders <- function(orders_file)  {
    orders<-read.csv(orders_file, sep=";", quote='"', stringsAsFactors=FALSE)
    # rename the columns
    orders <- orders %>% rename(customer_web_id = web_user_id,
                                custom_registration_number = cr_number,
                                order_web_id = order_id)
    # replace carriage returns by -
    orders$culture_use <- str_replace_all(orders$culture_use, "[\r\n]", " - ")
    # convert date
    date <- strptime(orders$order_date, "%m/%d/%Y - %H:%M")
    orders$order_date <- format(date,'%Y-%m-%d %H:%M:%S')
    # filters the orders to have only the new orders
    orders <- orders %>% filter(!(order_web_id %in% read_field_db(db_info("rcc"),"orders", "order_web_id")))
    return(orders)
}

# Order details ------------------------------------------------------------------

#' @title Format "order_details.csv" for import into rcc database
#'
#' @description
#' This function reads the file "order_details.csv" exported from the web site
#' and import it into the rcc database
#' @param order_details_file Name of file "order_details.csv" exported from the web site
#'
#' @return
#' Data frame ready for import into rcc database
#'
#' @examples
#' df <- rcc_order_details("order_details.csv")
#' @export

rcc_order_details <- function(order_details_file)  {
    # read the file
      order_details<-read.csv(order_details_file, sep=";", quote='"', stringsAsFactors=FALSE)
    # extract rcc number
      buffer <- as.data.frame(str_split_fixed(order_details$line_item_label, "_", 2))
      order_details$rcc_id <- as.numeric(str_replace(buffer$V1, "rcc", ""))
    # extract conditioning and charge_category
      buffer <- as.data.frame(str_split_fixed(order_details$conditioning, " - ", 3))
      order_details$conditioning <- buffer$V2
      order_details$order_type <- buffer$V3
    # rename the columns,  filter the file, filters the orders to have only the new items, and reorder by order_id and rcc_id
      order_details <-  order_details %>%
                        rename(order_detail_web_id = order_detail_id, order_web_id=order_id, culture_number=quantity) %>%
                        select(-created_date, -product_id, -line_item_label) %>%
                        filter(!(order_web_id %in% read_field_db(db_info("rcc"), "order_details", "order_web_id"))) %>%
                        filter(order_web_id > 705) %>%
                        arrange(order_web_id, rcc_id, conditioning)
      return(order_details)
}


# Customers ------------------------------------------------------------------

#' @title Format "customers.csv" for import into rcc database
#'
#' @description
#' This function reads the file "customers.csv" exported from the web site
#' and import it into the rcc database
#' @param customers_file Name of file "customers.csv" exported from the web site
#'
#' @return
#' Data frame ready for import into rcc database
#'
#' @examples
#' df <- rcc_customers("customers.csv")
#' @export

rcc_customers <- function(customers_file)  {
  # read the file
    customers<-read.csv(customers_file, sep=";", quote='"', stringsAsFactors=FALSE)
  # rename the columns and filter
    customers <- customers %>%  rename(customer_web_id = web_id, user_name = name) %>%
                                filter(!(customer_web_id %in% read_field_db(db_info("rcc"), "customers", "customer_web_id")))
      return(customers)
}

# Export tables from MySQL database ---------------------------------------

#
#' @title Export tables from MySQL database
#' @description
#' This function exports the data from the rcc database to
#' 3 files to be imported into the RCC web site version 1.0
#' * web_sequences.csv
#' * web_strains.csv
#' * web_strains_products.csv
#'
#' @param export_directory By defaut "C:/Databases/web_export/"
#'
#' @return
#' TRUE if succesful
#'
#' @examples
#' rcc_export()
#' @md
#' @export

rcc_export <- function(export_directory = "C:/Databases/web_export/") {

   db_con_rcc <- db_connect(db_info("rcc"))

   cultures <- tbl(db_con_rcc, "cultures") %>% collect()
   taxonomy <- tbl(db_con_rcc, "taxonomy") %>% collect()
   sequences <- tbl(db_con_rcc, "sequences") %>% collect()
   projects <- tbl(db_con_rcc, "projects") %>% collect()
   products <- tbl(db_con_rcc, "products") %>% collect()

   db_disconnect(db_con_rcc)

  # I do not do how to do this with dplyr....
   product_query = "SELECT
                    cultures.rcc_id,
                    products.`code`,
                    products.description,
                    products.price,
                    products.available_for,
                    cultures.species,
                    cultures.distributed,
                    cultures.lost
                    FROM
                    cultures ,
                    products
                    ORDER BY
                    cultures.rcc_id ASC,
                    products.`code` ASC"
   products <- db_get_query(db_info("rcc"), product_query)


   # For cultures which are not distributed only the collaborator conditioning is available and
   # for lost cultures no conditionning at all
   # (status=0)
   # products <- products %>% filter(!((distributed==0)&(code!="01")))
   products <- products %>% transmute (
                   RCC = rcc_id,
                   SKU=paste("rcc",rcc_id,"_",code, sep=""),
                   Title = paste("RCC",rcc_id," ",species," - ", description, sep = ""),
                   Status = ifelse((((distributed==0)&(code!="01"))|(lost!=0)),0,1),
                   Price = price,
                   Culture_form = code
                )
   products_file = paste(export_directory, "web_strains_products.csv", sep="")
   write.table(products, file=products_file, sep="\t",
               quote=TRUE, row.names = FALSE, na="", fileEncoding = "UTF-8")


   projects <- projects %>% group_by(rcc_id) %>% summarise(project_names=paste(project_name, collapse="; "))

   products <- products %>% group_by(RCC) %>% summarise(SKU=paste(SKU, collapse=",")) %>% rename(rcc_id=RCC)


   cultures <- left_join(cultures, taxonomy)
   cultures <- left_join(cultures, projects)
   cultures <- left_join(cultures, products)

   cultures_export <- cultures %>% transmute(
          Numero_Roscoff = rcc_id ,
          Publish_to_web = distributed ,
          Date_entered_catalog = date_entered_catalog,
          Category=common_name ,
          Domain = domain ,
          Division = division ,
          Class = class ,
          Order = order ,
          Family = family ,
          Genus = genus ,
          Species = str_split_fixed(species ,"_",n=2)[,2],
          Ecotype = clade ,
          Strain = strain_name ,
          Other_Name = strain_name_synonyms ,
          Clonal = clonal ,
          Heterotrophic = heterotrophic ,
          Axenic = axenic ,
          Toxic = toxic ,
          Symbiotic = symbiotic,
          Frozen = active_transfer_stopped ,
          Cryo_strain_substitute = active_transfer_stopped_substitute ,
          Lost = lost ,
          Date_Lost = lost_date,
          Size = length ,
          Cell_shape = cell_shape ,
          Cell_assemblage = cell_assemblage ,
          Pigment_phenotype = phenotype_pigment ,
          Ocean_origin = sampling_ocean ,
          Region_origin = sampling_regional_sea ,
          Country_origin = sampling_country ,
          Cruise_Isolation = sampling_cruise ,
          Station_Isolation = sampling_station ,
          Depth_Isolation = sampling_depth ,
          Date_Isolation = sampling_date ,
          Lat_Deg_Isolation = sampling_lat_deg ,
          Lat_Mn_Isolation = sampling_lat_mn ,
          Lat_NS_Isolation = sampling_lat_ns ,
          Long_Deg_Isolation = sampling_long_deg ,
          Long_Mn_Isolation = sampling_long_mn ,
          Long_EW_Isolation = sampling_long_ew ,
          Isolator = isolation_by ,
          Temperature_Isolation = isolation_growth_temperature,
          isolated_from_strain = isolated_from_strain,
          Medium_SBR = rcc_medium ,
          Temperature_SBR = rcc_temperature ,
          Light_SBR = rcc_light ,
          Remark_Web = remark_web,
          Genome_sequenced = sequenced_genome,
          Projects = project_names ,
          SKU = SKU,
          RCC_Image = sprintf("sites/default/files/images/RCC%04d.jpg", rcc_id),
          Latitude  = mapply(lat_long_dec,sampling_lat_deg,sampling_lat_mn,sampling_lat_ns,"latitude"),
          Longitude = mapply(lat_long_dec,sampling_long_deg,sampling_long_mn,sampling_long_ew,"longitude")
          )
   # Compute latitude and longitude in decimal - Replace by mapply
   # cultures_export$Latitude_2  =  ifelse(cultures$sampling_lat_ns=="N",1,-1)*
   #                              (cultures$sampling_lat_deg+cultures$sampling_lat_mn/60)
   # cultures_export$Longitude_2 =  ifelse(cultures$sampling_long_ew=="E",1, -1)*
   #                              (cultures$sampling_long_deg+cultures$sampling_long_mn/60)

   # Correct species for any group that contains X_sp
   cultures_export$Species[str_detect(cultures_export$Species,"X_sp")] <- "sp"

   cultures_file = paste(export_directory, "web_strains.csv", sep="")
   write.table(cultures_export, file=cultures_file, sep="\t",
               quote=TRUE, row.names = FALSE, na="",  fileEncoding = "UTF-8")


    cultures_taxo <- select(cultures, rcc_id, class, order, genus, species)
    sequences <- left_join(sequences, cultures_taxo)
    sequences <- sequences %>% transmute(
            id_sequence=sequence_id ,
            RCC=rcc_id ,
            DateAdded=format(strptime(date_added, "%Y-%m-%d"), "%d/%m/%Y - %H:%M:%S") ,
            Accession=genbank_accession ,
            Genome_link=genome_link ,
            Length=sequence_length ,
            Gene=gene_name ,
            Gene_location=gene_location ,
            Genbank_date=genbank_date ,
            Genbank_organism=genbank_organism ,
            Genbank_taxonomy=genbank_taxonomy ,
            Description=genbank_description,
            Class = class ,
            Order = order ,
            Genus = genus ,
            Species = str_split_fixed(species ,"_",n=2)[,2]
            )
   sequences_file = paste(export_directory, "web_sequences.csv", sep="")
   write.table(sequences, file=sequences_file, sep="\t",
               quote=TRUE, row.names = FALSE, na="",  fileEncoding = "UTF-8")

}
