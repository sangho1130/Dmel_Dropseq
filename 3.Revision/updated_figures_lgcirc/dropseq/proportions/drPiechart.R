  
drPiechart <- function(columnNames, Values, PlotTitle, outputPdf){
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  
  data <- data.frame(
    group = columnNames,
    value = Values
  )
  data$group <- factor(data$group, columnNames)
  
  colourCount <- length(Values)
  getPalette <- colorRampPalette(brewer.pal(9, "Reds")) #or Spectral

  pie <- ggplot(data, aes(x="", y=value, fill=factor(group))) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + 
    scale_fill_manual(values = getPalette(colourCount)) + 
    labs(title = PlotTitle) + 
    coord_polar(theta = "y", direction = -1) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()) 

  ggsave(outputPdf, units = 'cm', height = 8, width = 16)
}

drPiechart(c('CC 1', 'CC 2'), 
           c(163, 113), 
           '', 'cell_count_CC_LG.pdf')
drPiechart(c('CC 1', 'CC 2'), 
           c(91, 60), 
           '', 'cell_count_CC_Circ.pdf')

drPiechart(c('PM 1', 'PM 2', 'PM 3', 'PM 4'), 
           c(6599, 1812, 457, 578), 
           '', 'cell_count_PM_LG.pdf')
drPiechart(c('PM 1', 'PM 2', 'PM 3', 'PM 4'), 
           c(2096, 0, 0, 0), 
           '', 'cell_count_PM_Circ.pdf')

drPiechart(c('PH 1', 'PH 2', 'PH 3', 'PH 4', 'PH 5', 'PH 6'), 
           c(64, 72, 336, 4267, 1321, 448), 
           '', 'cell_count_PH_LG.pdf')
drPiechart(c('PH 1', 'PH 2', 'PH 3', 'PH 4', 'PH 5', 'PH 6'), 
           c(65, 0, 0, 54, 0, 0), 
           '', 'cell_count_PH_Circ.pdf')

