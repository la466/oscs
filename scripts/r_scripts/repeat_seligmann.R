
# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/other/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)

theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2, size=rel(0.8)),
      axis.title.x = element_text(vjust = -0.2, size=rel(0.8)),
      axis.text.x = element_text(size=rel(0.8)),
      axis.text.y = element_text(size=rel(0.8)),
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      # legend.margin = unit(0, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
    ))
}

file <- read.csv('outputs/seligmann_pollock_repeat/correlation.csv', head=TRUE)
file$padj_mean_use <- p.adjust(file$p_mean_use, method="fdr")
file$sig_factor <- ifelse(file$padj_mean_use < 0.05, 'P < 0.05', 'P \u2265 0.05')

colours <- c('blue', 'black')
ylab <- expression(italic(rho))

library(ggplot2)
library(Cairo)
plot <- ggplot(file) +
  geom_point(aes(x=gc, y=cor_mean_use, col=sig_factor, shape=sig_factor), size=1.2) +
  scale_shape_manual(values=c(4,16)) +
  scale_colour_manual(values=colours) +
  scale_y_continuous(limits = c(-1,1)) +
  labs(x="GC", y=ylab) +
  theme_Publication() +
  theme(
    legend.position = c(0.9,0.95),
    legend.direction = "vertical",
    legend.title=element_blank(),
    legend.key.size = unit(1, 'lines')
  )
ggsave('graphs/other/codon_contribution_to_stop_correlations_seligmann.pdf', plot=plot, device=cairo_pdf, width=7, height=7)



## Analysis

sink('outputs/r_outputs/repeat_seligmann.csv')

file$padj_mean_use <- p.adjust(file$p_mean_use, method="fdr")

pos_mean <- file[file$cor_mean_use > 0,]
pos_sig_mean <- pos_mean[pos_mean$p_mean_use < 0.05,]
pos_sig_adj_mean <- pos_mean[pos_mean$padj_mean_use < 0.05,]
cat(sprintf('Number of genomes with significant excesses: %s / %s\n', nrow(pos_mean), nrow(file)))
cat(sprintf('Number of genomes with significant excesses after correction: %s / %s\n', nrow(pos_sig_adj_mean), nrow(file)))


neg_mean <- file[file$cor_mean_use < 0,]
neg_sig_mean <- neg_mean[neg_mean$p_mean_use < 0.05,]
neg_sig_adj_mean <- neg_mean[neg_mean$padj_mean_use < 0.05,]
cat(sprintf('Number of genomes with significant reduction: %s\n', nrow(neg_sig_mean)))
cat(sprintf('Number of genomes with significant reduction after control: %s\n', nrow(neg_sig_adj_mean)))


sink()
