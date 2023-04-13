
# Libraries tidyr ---------------------------------------------------------

  library("stringr")

# Helper function to add remarks

  str_append<- function(x, y) ifelse(is.na(x), y, str_c(x, y, sep= "; "))

# Define PR2 global variables ----------------------------

  pr2.env <- list()
  pr2.env$pr2_directory = str_c(here::here("tests", "testthat", "output", "pr2"), "/")
  # pr2.env$genbank_directory = str_c("GenBank/")
  # pr2.env$pr2database_directory <- str_c(pr2.env$pr2_directory,"pr2-database/")
  # pr2.env$R_directory <- str_c(here::here("R"), "/")

  pr2.env$version = "5.0.0"
  pr2.env$version_directory = str_c(pr2.env$pr2_directory, "versions/",pr2.env$version,"/")


  dir.create(pr2.env$version_directory)

  pr2.env$date = format(Sys.time(), "%Y-%m-%d")

  # do not forget to add in the annotator table in the MySQL database
  pr2.env$editor = "D. Vaulot"

  pr2.env$ambiguities_max = 20
  pr2.env$sequence_length_min = 500
  pr2.env$sequence_length_max = 15000 # only for imported GenBank sequences
  pr2.env$sequence_N_repeat = "NN"

# Useful strings
  pr2.env$taxo_levels <- list()
  pr2.env$taxo_levels[[8]] = c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
  pr2.env$taxo_levels[[9]] = c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")

  pr2.env$taxo_levels_number = 9

  pr2.env$ambig_regex ="[^ATGCU]"

# create a regex containing the names of the levels (kingdom|supergroup|....)
  pr2.env$taxo_levels_regex <- str_c(pr2.env$taxo_levels[[pr2.env$taxo_levels_number]], collapse="|")
  pr2.env$taxo_levels_regex <- str_c("(",pr2.env$taxo_levels_regex, ")")

# The Eukref must be searched in this direction (first sp. and cf. and then the Genus species)
  pr2.env$genus_sp_regex <- "(^[A-Z]{1}[a-z]+)[ ]sp[.]"
  pr2.env$genus_cf_regex <- "(^[A-Z]{1}[a-z]+)[ ]cf[.][ ]([a-z]+)"
  pr2.env$genus_species_regex <- "(^[A-Z]{1}[a-z]+)[ ]([a-z]+)[ ]?"

# pr2_traits used
  pr2.env$traits_used <- c("mixoplankton")

# Files for export
  pr2.env$file_head <- str_c(pr2.env$version_directory, "pr2_version_", pr2.env$version)

  pr2.env$file_merged_excel <- str_c(pr2.env$file_head,"_merged", ".xlsx")

  pr2.env$file_chimera_excel <- str_c(pr2.env$file_head,"_chimera", ".xlsx")
  pr2.env$file_taxonomy_excel <- str_c(pr2.env$file_head,"_taxonomy", ".xlsx")
  pr2.env$file_unassigned_excel <- str_c(pr2.env$file_head,"_unassigned", ".xlsx")
  pr2.env$file_updated_excel <- str_c(pr2.env$file_head,"_updated", ".xlsx")

# Files for SSU export
  pr2.env$file_SSU_fasta_UTAX = str_c(pr2.env$file_head,"_SSU_UTAX", ".fasta.gz")
  pr2.env$file_SSU_fasta_taxo_long = str_c(pr2.env$file_head,"_SSU_taxo_long", ".fasta.gz")
  pr2.env$file_SSU_fasta_mothur =  str_c(pr2.env$file_head,"_SSU_mothur", ".fasta.gz")
  pr2.env$file_SSU_taxo_mothur <- str_c(pr2.env$file_head,"_SSU_mothur",".tax.gz")
  pr2.env$file_SSU_fasta_dada2 =  str_c(pr2.env$file_head,"_SSU_dada2", ".fasta.gz")


  pr2.env$file_sqlite <- str_c(pr2.env$file_head,".sqlite")
