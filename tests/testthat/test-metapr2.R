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
  #                      directory = "output/metapr2/6/",
  #                      dataset_id_selected = 6,
  #                      filter_metadata = NULL,
  #                      export_long_xls=FALSE,
  #                      export_wide_xls=TRUE,
  #                      export_sample_xls=TRUE,
  #                      export_phyloseq = FALSE,
  #                      export_fasta=TRUE,
  #                      taxonomy_full = TRUE,
  #                      use_hash = TRUE,
  #                      sum_reads_min = 0,
  #                      sample_reads_min = 35000
  #                     )

  # Export with filter for water only samples

#   metapr2_export_asv (taxo_level = domain, taxo_name="Eukaryota",
#                       boot_level = class_boot, boot_min = 0,
#                       directory = "output/metapr2/47_filter/",
#                       dataset_id_selected = 47,
#                       filter_metadata = "((substrate == 'water') & is.na(substrate_description )) | (substrate == 'sediment trap material')",
#                       export_long_xls=TRUE,
#                       export_wide_xls=TRUE,
#                       export_sample_xls=TRUE,
#                       export_phyloseq = TRUE,
#                       export_fasta=TRUE,
#                       taxonomy_full = TRUE,
#                       use_hash = FALSE,
#                       sum_reads_min = 0)

metapr2_export_asv (gene = "SSU",
                    taxo_level = class, taxo_name="Coscinodiscophyceae",
                    boot_level = order_boot, boot_min = 0,
                    directory = "output/metapr2/PacBio/",
                    dataset_id_selected = c(393),
                    filter_metadata = NULL,
                    export_long_xls=FALSE,
                    export_wide_xls=FALSE,
                    export_sample_xls=FALSE,
                    export_phyloseq = FALSE,
                    export_fasta=TRUE,
                    export_fasta_sum_reads=TRUE,
                    taxonomy_full = TRUE,
                    use_hash = FALSE,
                    sum_reads_min = 0)


# test_that("metapr2 asv qs ", {
#   metapr2_export_qs (set_type = "public 1.0",
#                      directory = "output/metapr2/",
#                      do_cluster = TRUE)
#   asv_set <- qs::qread("tests/testthat/output/metapr2/asv_set_cluster.qs")
#   global <- qs::qread("tests/testthat/output/metapr2/global.qs")
# })



})
