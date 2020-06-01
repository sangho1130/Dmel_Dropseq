library(ggplot2)
library(reshape2)

countdf <- read.delim('cell_count.txt')
colnames(countdf) <- c('celltype', 'L', 'C')
countdf.m <- melt(countdf)
colnames(countdf.m) <- c('celltype', 'Origin', 'Proportion')
countdf.m$celltype <- factor(countdf.m$celltype, levels = rev(c('PSC', 'PH', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte')))

ggplot(countdf.m, aes(Origin, Proportion, fill = celltype)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = rev(c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'))) +
  labs(title = '', x = '', y = 'Prop. (%)') + 
  theme_bw(base_size = 7) +
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'None')
#ggsave('cell_count.pdf', units = 'cm', width = 2, height = 3)

