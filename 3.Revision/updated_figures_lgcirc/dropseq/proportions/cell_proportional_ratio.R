library(ggplot2)
library(reshape2)

countdf <- read.delim('cell_proportional_ratio.txt', sep = '\t', check.names = F)
colnames(countdf) <- c('celltype', 'L', 'C')
countdf.m <- melt(countdf)
colnames(countdf.m) <- c('celltype', 'Origin', 'Proportion')
countdf.m$celltype <- factor(countdf.m$celltype, levels = rev(c("PSC", "PH", "PM", "LM", "CC", "GST-rich", "Adipohemocyte")))
countdf.m$Origin <- factor(countdf.m$Origin, levels = c('C', 'L'))

ggplot(countdf.m, aes(celltype, Proportion, fill = Origin)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = c('#7ac5cd', '#ffa500')) +
  scale_y_continuous(breaks = seq(0, 10, 5)/10) +
  coord_flip() + 
  labs(title = '', x = '', y = 'Proportional ratio') + 
  theme_bw(base_size = 7) +
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'None')
#ggsave('cell_proportional_ratio.pdf', units = 'cm', width = 3.5, height = 3)

