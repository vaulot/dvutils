# dvutils, a R package

A set of R utilities developped by Daniel Vaulot (vaulot@gmail.com)

* Graphics
* Databases
* Bioinfo
    * interface with Phyloseq
    * interface with dada2
    * fasta file read and write
    * Genbank parsing
* Interface and management to pecific databases
    * PR2 (18S rRNA sequences)
    * MetaPR2 (metabarcodes)
    * RCC (Roscoff Culture Collection)

# Installation

``` r
install.packages(devtools)
devtools::install_github("vaulot/dvutils")
```

# Functions


Function | Aim
--- | ----
blast_18S_reformat|Process a tabular output from blastn (BLAST+)
blast_summary|Write a summary for a tabular output from blastn (BLAST+)
dada2_assign|Assign sequences using dada2 wang assigner
db_append_records|Write new records to a table.
db_connect|Establish a connection database
db_connect_sqlite|Establish a connection to a SQLite database
db_disconnect|Disconnect an existing connection
db_execute_query_vector|Execute a vector of query.
db_get_query|Read a query into a dataframe.
db_info|Get the connection info for a specific database
db_sql_escape|Put a string between single quotes
db_update_field_string|Build a string to update a field in a database
fasta_filter|Filter a fasta file
fasta_read|Read a fasta file into a data frame
fasta_write|Write a fasta file with the taxonomy
fastq_subsample|Subsample fastq files
filename_append|Append a string at the end of a file name before the exension
filename_change_ext|Change the extension of the file name
file_unix2dos|Convert between Unix and Dos text format
genbank_download|Download sequences from GenBank
genbank_download_parse|Download and parse sequences from GenBank
genbank_features|Read features of sequences from GenBank
genbank_field|Read a single field from a set of existing GenBank files
genbank_locus|Read locus from a set of existing GenBank files
genbank_taxonomy|Read NCBI taxo
get_primer_position|Get primer position on sequence set
gg_bar_discrete|Do a simple barplot
gg_boxplot|Plot a box plot
gg_density|Plot simple density plot to compare factors
gg_hist|Plot simple histogram
gg_violin|Plot a violin and a box plot
kmer|Get all kmers in a sequence
kmer_set|Get all kmers in a set of sequence
latex_fix_bibfile|Fix bib library file created by Mendeley
lat_long_dec|Convert lat and long to decimal
map_distribution|Draw distribution map from data frame
map_get_world|Get the world map for the rworldmap package
map_leaflet|Map with leaflet package (OpenStreetMap)
map_US|Background US map using the maps package to extract a specific country
map_world|Background map using the maps package
map_world_google|Background map using the ggmap package
map_zoom_europe|Zoom a map on Europe
metapr2_export_asv|Exports the metapr2 database
pcr_sequences|In silico amplification
phyloseq_filter_abundant_taxa|Filter a phyloseq table keeping only abundant taxa
phyloseq_filter_autotrophic_taxa|Filter a phyloseq table keeping only autotrophic taxa
phyloseq_import_mothur|Create a phyloseq file from mothur database file
pr2_buid_taxons|Build taxon table (long format)
pr2_build_taxonomy|Buikd taxonomy table (wide format) from taxon table (long format)
pr2_export|Export the PR2 database (one file)
pr2_export_all|Export the PR2 database (all files)
pr2_export_sqlite|Export the PR2 database to a SQLite file
pr2_extract_one_taxo_level|Extract one level of taxonomy
pr2_read|Reads the whole PR2 database into a data frame
pr2_sequence_label|Create a simple label for a sequence
pr2_sequence_reassign|Reassign pr2 sequences
pr2_taxo_check|Check taxonomy
pr2_taxo_list|Build a list of taxa
pr2_taxo_read|Reads the PR2 taxonomic database into a data frame
pr2_taxo_update|Update pr2 taxonomy
pr2_taxo_X|Add _X for taxon that are repeated on same line
pr2_treemap|Draw a simple treemap from pr2.
rcc_customers|Format "customers.csv" for import into rcc database
rcc_export|Export tables from MySQL database
rcc_genbank|Prepare list of GenBank entries to include in rcc database
rcc_orders|Format "orders.csv" for import into rcc database
rcc_order_details|Format "order_details.csv" for import into rcc database
theme_dviz_grid|Theme - grid lines along major axis ticks, no axes
theme_dviz_hgrid|Theme - horizontal grid lines only
theme_dviz_vgrid|Theme - vertical grid lines only
treemap_dv|Do a simple treemap
XStringSet_to_df|Transforms a DNA or AA String set into a data frame
