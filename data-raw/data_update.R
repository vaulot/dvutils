library (dplyr)

colors_euk_division <- readxl::read_excel("data-raw/colors.xlsx", sheet = "euk division")

colors_euk_division <-structure(colors_euk_division$color_name,.Names=colors_euk_division$division)


colors_euk_class <- readxl::read_excel("data-raw/colors.xlsx", sheet = "euk classes")

colors_euk_class <-structure(colors_euk_class$color_name,.Names=colors_euk_class$class)

colors_euk_genus <- readxl::read_excel("data-raw/colors.xlsx", sheet = "euk genus")

colors_euk_genus <-structure(colors_euk_genus$color_name,.Names=colors_euk_genus$genus)

usethis::use_data(colors_euk_division, colors_euk_class, colors_euk_genus, overwrite = TRUE, internal=FALSE)
