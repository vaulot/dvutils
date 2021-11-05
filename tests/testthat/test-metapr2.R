# MetaPR2 ---------------------------------------------------------------------

context("metapr2")


test_that("metapr2 asv export ", {
  # Export all the asv in a single fasta
  # metapr2_export_asv()

  # Export as specific data set as a phyloseq file

  # metapr2_export_asv(dataset_id_selected = 23, export_phyloseq = TRUE)

  # Export a specific genus as a fasta file and an excel file

  # asv_set <- metapr2_export_asv(taxo_level = genus, taxo_name="Pseudohaptolina",
  #                               taxonomy_full= FALSE, boot_min = 100, export_xls = TRUE)

  # Export without any filter

  # metapr2_export_asv (taxo_level = kingdom, taxo_name="Eukaryota",
  #                      boot_level = class_boot, boot_min = 0,
  #                      directory = "output/metapr2/47/",
  #                      dataset_id_selected = 47,
  #                      filter_metadata = NULL,
  #                      export_long_xls=FALSE,
  #                      export_wide_xls=TRUE,
  #                      export_sample_xls=TRUE,
  #                      export_phyloseq = FALSE,
  #                      export_fasta=TRUE,
  #                      taxonomy_full = TRUE,
  #                      use_hash = FALSE,
  #                      sum_reads_min = 0)

  # Export with filter for water only samples

  # metapr2_export_asv (taxo_level = kingdom, taxo_name="Eukaryota",
  #                     boot_level = class_boot, boot_min = 0,
  #                     directory = "output/metapr2/47_filter/",
  #                     dataset_id_selected = 47,
  #                     filter_metadata = "((substrate == 'water') & is.na(substrate_description )) | (substrate == 'sediment trap material')",
  #                     export_long_xls=TRUE,
  #                     export_wide_xls=TRUE,
  #                     export_sample_xls=TRUE,
  #                     export_phyloseq = TRUE,
  #                     export_fasta=TRUE,
  #                     taxonomy_full = TRUE,
  #                     use_hash = FALSE,
  #                     sum_reads_min = 0)

  # Export the basic dataset for shiny

   metapr2_export_qs (set_type = "basic",
                      directory = "output/metapr2/")
})

