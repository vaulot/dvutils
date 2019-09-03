#' @import ggplot2


# Plot parameters  --------------------------------------------------------

  depth_color<-c("1.5"="#c6dbef","5"="#9ecae1", "10"="#6baed6", "20"="#4292c6",  "30"="#2171b5", "40"="#08519c", "60"="#08306b")
  depth_linetype<-c("1.5"=1,"5"=1, "10"=1, "20"=2, "30"=2, "40"=2, "60"=2)
  depth_shape <-c("1.5"=15,"5"=15, "10"=16, "20"=16, "30"= 17, "40"=17, "60"=18)

  scaling_factor=15
  cell_label <- expression (paste("cell.",mL^-1))
  pico_label <- expression (paste("Pico cell.",mL^-1))
  nano_label <- expression (paste("Nano cell.",mL^-1))
  cell_breaks=c(10, 100,1000,10000,100000)
  cell_limits_phyto = c(10, 100000)
  cell_limits_bact = c(100000, 1500000)
  julian_day_limits = c(1, 365)
  depth_limits = c(60, 0)  # Need to invert for inverted scale....
  date_limits = c(as.Date("2014-01-01"), as.Date("2015-02-01"))


# Default BW theme ----------------------------------------
#' @export
#' @md
gg_theme_bw_dv1 <- function() {

        theme_bw(scaling_factor) +
        theme(panel.border = element_rect(colour = "black"),
              axis.line = element_line(colour = "black") ,
              legend.title=element_text(size=scaling_factor),
              legend.key=element_blank(),
              axis.title = element_text(size=scaling_factor),
              legend.text=element_text(size=scaling_factor*0.5),
              legend.key.height = unit(1, "cm"),
              legend.key.width = unit(3, "cm"),
              axis.text = element_text(size=0.8*scaling_factor),
              panel.background = element_rect(fill="white"),
              legend.position = "top",
              legend.box = "horizontal"
              )
  }

# Histogram ----------------------------------------
#' @title Plot simple histogram
#' @description
#' @param df Data frame - data to be plotted
#' @param x Character - Which column use for x axis (quoted variable), should correspond to a **numeric column**
#' @param bin Numeric - Size of the bins (default = 0.2)
#' @return a ggplot object
#' @examples
#' @export
#' @md
gg_hist <- function(df,x, bin = 0.2) {

  p <- ggplot(df, aes_string(x=x)) +
        geom_histogram(aes(y=..density..), alpha = 0.5, position="identity", binwidth = bin) +
        geom_density(alpha=0.6 )
  return(p)
  }

# Density -----------------------------------------
#' @title Plot simple density plot to compare factors
#' @description
#' @param df Data frame - data to be plotted
#' @param x Character - Column corresponding x axis (quoted variable), should correspond to a **numeric column**
#' @param factor Character - Factor for the different densities
#' @return a ggplot object
#' @examples
#' @export
#' @md
gg_density <- function(df,x, factor) {

  p <- ggplot(df, aes_string(x = x, fill = factor)) +
          geom_density(alpha=0.6 )
  return(p)
  }

# Boxplot -----------------------------------------
#' @title Plot a box plot
#' @description
#' @param df Data frame - data to be plotted
#' @param y Character - Which column use for y axis (quoted variable), should correspond to a **numeric column**
#' @param factor Character - Which column use for x axis (quoted variable), should correspond to a **character column**
#' @return a ggplot object
#' @examples
#' @export
#' @md
gg_boxplot <- function(df,y, factor) {

  p <- ggplot(df, aes_string( y= y, x = factor, fill = factor)) +
          geom_boxplot(alpha=0.6)
  return(p)
  }

# Violin -----------------------------------------
#' @title Plot a violin and a box plot
#' @description
#' @param df Data frame - data to be plotted
#' @param y Character - Which column use for y axis (quoted variable), should correspond to a **numeric column**
#' @param factor Character - Which column use for x axis (quoted variable), should correspond to a **character column**
#' @return a ggplot object
#' @examples
#' @export
#' @md
gg_violin <- function(df,y, factor) {

  p <- ggplot(df, aes_string(y=y, x = factor, fill = factor)) +
    geom_violin() +
    geom_boxplot(width=0.1)
  return(p)
  }

# Barplot -----------------------------------------
#' @title Do a simple barplot
#' @description
#' @param df Data frame
#' @param x Character - Which column use for x axis (quoted variable), should correspond to a **character column**
#' @param y Character - Which column use for y axis (quoted variable), should correspond to a **numeric column**
#' @param fill Character - Which column use for the slices of the bar plot (quoted variable), should correspond to a **character column**
#' @param bar_labels TRUE or FALSE - to add the label inside the graph itself to the bar
#' @param flip TRUE or FALSE - to flip to horizontal mode
#' @param legend TRUE or FALSE - if TRUE the  two following items must be passed to the function, if FALSE then automatic legends
#' @param legend_colors List - Correspondance between fill categories and colors
#' @param legend_label Character - Title of the legend
#' @return
#' A ggplot2 object
#' @examples
#' p <- dv_bar_discrete (cultures_division,x="division", y="n_strains", fill="domain",
#'                      ymin=0, ymax=400, xlab="Division", ylab="Number of strains",
#'                      bar_labels=TRUE, flip=TRUE,
#'                      legend=TRUE,
#'                      legend_colors = c("deposit" = "yellow", "collaborator" = "blue",
#'                      "public" = "red", "education" = "green", "private"="black"),
#'                      legend_label = "order type")
#' @export
#' @md

gg_bar_discrete <- function (df,x, y, fill, xlimits, ymin, ymax, xlab, ylab,
                               bar_labels=FALSE, flip=FALSE,
                               legend=FALSE, legend_colors=c("none"="black", legend_label="")) {

  p <- ggplot(df, aes_string(x=x, y=y, fill=fill) ) +
      geom_bar(stat = "identity") +
      scale_y_continuous(limits = c(ymin,ymax)) +
      xlab(xlab) + ylab(ylab) +
      scale_x_discrete(limits = xlimits)   # This line used to order the x scale, but does not work

  if (bar_labels) { p <- p +   geom_text(aes_string(label=y), vjust=0.2, hjust=-0.2)}
  if (flip) { p <- p +   coord_flip()  }
  if (legend) {p <- p +   scale_fill_manual(legend_label, values = legend_colors)}

  return(p)
}

# Function to plot a treemap ----------------------------------------------
#' @title Do a simple treemap
#' @param df Data frame
#' @param index Character vector - what variables are used for categories
#' @param size Character - what variable is use for size of boxes
#' @param title Character - Title
#' @return
#' A treemap object and plots the treemap
#' @examples
#' treemap_dv(dada2_species, c("Division", "Class"),"n_seq","Dada2" )
#' @export
#' @md

treemap_dv <- function(df, index, size, title){
  treemap::treemap(df, index=index,
        vSize=size,
        title=title,asp=1, lowerbound.cex.labels= 0.2, fontsize.labels = 12,
        palette="Blues",
        format.legend = list(scientific = FALSE, big.mark = " "))
}



# Themes from dataviz book ------------------------------------------------

#' @title Theme - horizontal grid lines only
#' @description
#' From data viz book : http://serialmentor.com/dataviz/
#' @param font_size font size
#' @param font_family font family
#'
#' @return
#' Theme
#' @export

theme_dviz_hgrid <- function(font_size = 14, font_family = "") {
  color = "grey90"
  line_size = 0.5

  # Starts with theme_cowplot and then modify some parts
  # Programing note ; %+replace% replaces the entire element
  cowplot::theme_cowplot(font_size = font_size, font_family = font_family) %+replace%
    theme(
      # make horizontal grid lines
      panel.grid.major   = element_line(colour = color,
                                        size = line_size),
      panel.grid.major.x = element_blank(),

      # adjust axis tickmarks
      axis.ticks        = element_line(colour = color, size = line_size),

      # adjust x axis
      axis.line.x       = element_line(colour = color, size = line_size),
      # no y axis line
      axis.line.y       = element_blank()
      )
}


#' @title Theme - vertical grid lines only
#' @description
#' From data viz book : http://serialmentor.com/dataviz/
#' @param font_size font size
#' @param font_family font family
#'
#' @return
#' Theme
#' @export

theme_dviz_vgrid <- function(font_size = 14, font_family = "") {
  color = "grey90"
  line_size = 0.5

  # Starts with theme_cowplot and then modify some parts
  cowplot::theme_cowplot(font_size = font_size, font_family = font_family) %+replace%
    theme(
      # make vertical grid lines
      panel.grid.major   = element_line(colour = color,
                                        size = line_size),
      panel.grid.major.y = element_blank(),

      # adjust axis tickmarks
      axis.ticks        = element_line(colour = color, size = line_size),

      # adjust y axis
      axis.line.y       = element_line(colour = color, size = line_size),
      # no x axis line
      axis.line.x       = element_blank()
    )
}


#' @title Theme - grid lines along major axis ticks, no axes
#' @description
#' From data viz book : http://serialmentor.com/dataviz/
#' @param font_size font size
#' @param font_family font family
#'
#' @return
#' Theme
#' @export

theme_dviz_grid <- function(font_size = 14, font_family = "") {
  color = "grey90"
  line_size = 0.5

  # Starts with theme_cowplot and then modify some parts
  cowplot::theme_cowplot(font_size = font_size, font_family = font_family) %+replace%
    theme(
      # make horizontal grid lines
      panel.grid.major   = element_line(colour = color,
                                        size = line_size),

      # adjust axis tickmarks
      axis.ticks        = element_line(colour = color, size = line_size),

      # no x or y axis lines
      axis.line.x       = element_blank(),
      axis.line.y       = element_blank()
    )
}
