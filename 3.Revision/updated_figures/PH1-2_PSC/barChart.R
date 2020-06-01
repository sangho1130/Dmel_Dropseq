library(ggplot2)

barChart <- function(Categories, Values, GOcolors, PlotTitle, outputPdf){
  data <- data.frame(
    group = Categories,
    value = -log10(Values),
    colors = GOcolors
  )
  data$group <- factor(data$group, rev(Categories))
  data
  bar <- ggplot(data, aes(x=group, y=value)) +
    geom_bar(width = 0.75, stat = "identity", col = 'black', fill = rev(as.character(data$colors))) + 
    coord_flip() +
    labs(title = PlotTitle, x = '', y = '-log10 FDR') + 
    theme_bw(base_size = 7) + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black', size = 7),
          axis.ticks = element_line(colour = 'black')) 
  
  ggsave(outputPdf, units = 'cm', height = 3, width = 8)
}

barChart(c('Notch signaling pathway',
           'MAPK signaling pathway - fly'), #Keywords
         c(0.00112056,
           0.025338896), 
         c('#e9137c', 
           '#e9137c'), 
         'KEGG pathway', 'gProfiler_Dl+N+.pdf')
