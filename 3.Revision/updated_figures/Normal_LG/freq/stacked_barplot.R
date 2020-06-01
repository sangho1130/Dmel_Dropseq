library(reshape2)
library(scales)
library(ggplot2)
library(plyr)
library(RColorBrewer)
label <- readRDS('../rdata/label.Rds')
head(label)

### PH
ph_label <- subset(label, anno_simple == 'PH')
head(ph_label)
ph_label_summ <- ph_label[, c(8, 9)]
ph_label_summ <- droplevels(ph_label_summ)

result <- data.frame(matrix(nrow = 11, ncol = 7))
colnames(result) <- c('Subclustering', levels(ph_label_summ$new_subclustering))
result$Subclustering <- levels(ph_label_summ$Subclustering)
for (new in levels(ph_label_summ$new_subclustering)) {
  tmp <- subset(ph_label_summ, new_subclustering == new)
  result[, new] <- summary(tmp$Subclustering)
}

for (new in levels(ph_label_summ$new_subclustering)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
result_m$Subclustering <- factor(result_m$Subclustering, levels = levels(ph_label_summ$Subclustering))

getPalette <- colorRampPalette(brewer.pal(9, 'GnBu'))


p <- ggplot(result_m ,aes(x = variable, y = value, fill = Subclustering)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c(getPalette(12)[2:13])) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
#ggsave('ph.freqs.pdf', units = 'cm', width = 10, height = 10)


ph_label <- subset(label, anno_simple == 'PH')
head(ph_label)
ph_label_summ <- ph_label[, c(4, 9)]
ph_label_summ <- droplevels(ph_label_summ)
ph_label_summ$Library <- factor(ph_label_summ$Library, levels = unique(ph_label_summ$Library))

result <- data.frame(matrix(nrow = 14, ncol = 7))
colnames(result) <- c('Library', levels(ph_label_summ$new_subclustering))
result$Library <- unique(ph_label_summ$Library)

for (new in unique(ph_label_summ$new_subclustering)) {
  tmp <- subset(ph_label_summ, new_subclustering == new)
  result[, new] <- summary(tmp$Library)
}

for (new in unique(ph_label_summ$new_subclustering)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
result_m$Library <- factor(result_m$new_subclustering, levels = unique(ph_label_summ$new_subclustering))
head(result_m)

p <- ggplot(result_m ,aes(x = variable, y = value, fill = Library)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('#eee5c0', '#ebd2b0', '#eeba86', '#eea160', '#ec8c3c', 
                               '#e4d3df', '#e9c0c8', '#e191b2', '#df5598', '#e60c7d', 
                               '#c2e2dd', '#8ed2d3', '#63afbc', '#519692')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
#ggsave('ph.freqs.libs.pdf', units = 'cm', width = 10, height = 10)


### PM
pm_label <- subset(label, anno_simple == 'PM')
head(pm_label)
pm_label_summ <- pm_label[, c(8, 9)]
pm_label_summ <- droplevels(pm_label_summ)

result <- data.frame(matrix(nrow = 10, ncol = 5))
colnames(result) <- c('Subclustering', levels(pm_label_summ$new_subclustering))
result$Subclustering <- levels(pm_label_summ$Subclustering)
for (new in levels(pm_label_summ$new_subclustering)) {
  tmp <- subset(pm_label_summ, new_subclustering == new)
  result[, new] <- summary(tmp$Subclustering)
}

for (new in levels(pm_label_summ$new_subclustering)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
result_m$Subclustering <- factor(result_m$Subclustering, levels = levels(pm_label_summ$Subclustering))

getPalette <- colorRampPalette(brewer.pal(9, 'OrRd'))

p <- ggplot(result_m ,aes(x = variable, y = value, fill = Subclustering)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c(getPalette(11)[2:12])) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('pm.freqs.pdf', units = 'cm', width = 8, height = 10)


pm_label <- subset(label, anno_simple == 'PM')
head(pm_label)
pm_label_summ <- pm_label[, c(4, 9)]
pm_label_summ <- droplevels(pm_label_summ)
pm_label_summ$Library <- factor(pm_label_summ$Library, levels = unique(pm_label_summ$Library))

result <- data.frame(matrix(nrow = 14, ncol = 7))
colnames(result) <- c('Library', levels(pm_label_summ$new_subclustering))
result$Library <- levels(pm_label_summ$Library)
for (new in levels(pm_label_summ$new_subclustering)) {
  tmp <- subset(pm_label_summ, new_subclustering == new)
  result[, new] <- summary(tmp$Library)
}

for (new in levels(pm_label_summ$new_subclustering)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
result_m$Library <- factor(result_m$Library, levels = levels(pm_label_summ$Library))

p <- ggplot(result_m ,aes(x = variable, y = value, fill = Library)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('#eee5c0', '#ebd2b0', '#eeba86', '#eea160', '#ec8c3c', 
                               '#e4d3df', '#e9c0c8', '#e191b2', '#df5598', '#e60c7d', 
                               '#c2e2dd', '#8ed2d3', '#63afbc', '#519692')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('pm.freqs.libs.pdf', units = 'cm', width = 8, height = 10)



### PSC

psc_label <- subset(label, anno_simple == 'PSC')
head(ph_label)
psc_label_summ <- psc_label[, c(4, 9)]
psc_label_summ <- droplevels(psc_label_summ)
psc_label_summ$Library <- factor(psc_label_summ$Library, levels = unique(psc_label_summ$Library))
head(psc_label_summ)

result <- data.frame(matrix(nrow = 13, ncol = 2))
colnames(result) <- c('Library', levels(psc_label_summ$new_subclustering))
result$Library <- unique(psc_label_summ$Library)

for (new in unique(psc_label_summ$new_subclustering)) {
  tmp <- subset(psc_label_summ, new_subclustering == new)
  result[, new] <- summary(tmp$Library)
}

for (new in unique(psc_label_summ$new_subclustering)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
#result_m$Library <- factor(result_m$new_subclustering, levels = unique(psc_label_summ$new_subclustering))
head(result_m)

p <- ggplot(result_m ,aes(x = variable, y = value, fill = Library)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('#eee5c0', '#ebd2b0', '#eea160', '#ec8c3c', 
                               '#e4d3df', '#e9c0c8', '#e191b2', '#df5598', '#e60c7d', 
                               '#c2e2dd', '#8ed2d3', '#63afbc', '#519692')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
#ggsave('psc.freqs.libs.pdf', units = 'cm', width = 5, height = 10)



### GST
gst_label <- subset(label, anno_simple == 'GST-rich')
head(gst_label)
gst_label_summ <- gst_label[, c(4, 9)]
gst_label_summ <- droplevels(gst_label_summ)
gst_label_summ$Library <- factor(gst_label_summ$Library, levels = unique(gst_label_summ$Library))
head(gst_label_summ)

result <- data.frame(matrix(nrow = 13, ncol = 2))
colnames(result) <- c('Library', levels(gst_label_summ$new_subclustering))
result$Library <- unique(gst_label_summ$Library)

for (new in unique(gst_label_summ$new_subclustering)) {
  tmp <- subset(gst_label_summ, new_subclustering == new)
  result[, new] <- summary(tmp$Library)
}

for (new in unique(gst_label_summ$new_subclustering)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
head(result_m)

p <- ggplot(result_m ,aes(x = variable, y = value, fill = Library)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('#eee5c0', '#ebd2b0', '#eeba86', '#eea160', 
                               '#e4d3df', '#e9c0c8', '#e191b2', '#df5598', '#e60c7d', 
                               '#c2e2dd', '#8ed2d3', '#63afbc', '#519692')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
#ggsave('gst.freqs.libs.pdf', units = 'cm', width = 5, height = 10)


