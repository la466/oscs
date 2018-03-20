# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/trna/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)

theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text.x = element_text(size=rel(1)),
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



anticodon_repetoire_plot <- function(){
  file <- read.csv('outputs/trna/anticodon_repetoire.csv', head=T)
  pdf('graphs/trna/anticodon_repetoire.pdf')
  plot(file$gc3, file$unique_anticodons, pch=16, cex=0.6, xlab="GC3", ylab="Anticodon repetoire")
  dev.off()
}

anticodon_sparing_plot <- function(){
  file <- read.csv('outputs/trna/anticodon_sparing.csv', head=T)
  codons <- c("TTT","TTG","TTC","TTA","TGT","TGG","TGC","TGA","TCT","TCG","TCC","TCA","TAT","TAG","TAC","TAA","GTT","GTG","GTC","GTA","GGT","GGG","GGC","GGA","GCT","GCG","GCC","GCA","GAT","GAG","GAC","GAA","CTT","CTG","CTC","CTA","CGT","CGG","CGC","CGA","CCT","CCG","CCC","CCA","CAT","CAG","CAC","CAA","ATT","ATG","ATC","ATA","AGT","AGG","AGC","AGA","ACT","ACG","ACC","ACA","AAT","AAG","AAC","AAA")
  pdf('graphs/trna/anticodon_sparing.pdf', width=10)
  plot(file$gc3, file$index, pch=16, cex=0.4, yaxt="n", xlab="GC3", ylab="Anticodon")
  axis(side = 2, at = seq(1:64), labels=codons, cex.axis=0.3, las=2)
  dev.off()
}

median_shift_probability_plot <- function(frame) {
  library(ggplot2)
  file <- read.csv('outputs/trna/median_shift_probability.csv', head=T)
  col <- paste('median_cost_', frame, sep='')

  fit <- lm(file[[col]]~file$gc3)
  plot <- ggplot(file) +
    geom_point(aes(x=file$gc3, y=file[[col]]), size=0.8) +
    labs(x="Genome GC3", y='CDS median frameshift probability') +
    # geom_hline(yintercept=0, lty=2) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], col="red") +
    theme_Publication() +
    theme(legend.position = "none")

  return(plot)
}

cost_osc_density_correlation_plot <- function(frame) {
  library(ggplot2)
  file <- read.csv('outputs/trna/cds_cost_osc_density_correlations.csv', head=T)
  col <- paste(frame, '_cor', sep='')

  fit <- lm(file[[col]]~file$gc3)
  plot <- ggplot(file) +
    geom_point(aes(x=file$gc3, y=file[[col]]), size=1.5) +
    labs(x="Genome GC3", y=expression(paste(rho))) +
    geom_hline(yintercept=0, lty=2) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], col="red") +
    theme_Publication() +
    theme(axis.title = element_text(face = "bold",size = rel(1)),legend.position = "none")

  return(plot)
}



z_median_probability_plot <- function(frame, model) {

  file <- read.csv('outputs/trna/median_shift_probability_z.csv', head=T)
  zcol <- paste('z_', frame, '_', model, sep="")
  prob_col <- paste('prob_', frame, sep="")

  fit <- lm(file[[prob_col]]~file[[zcol]])

  library(ggplot2)
  plot <- ggplot(file) +
    geom_point(aes(x=file[[zcol]], y=file[[prob_col]]), size=0.5) +
    labs(x="Genome Z", y='CDS median frameshift probability') +
    geom_vline(xintercept=0, lty=2) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], col="red") +
    theme_Publication() +
    theme(legend.position = "none")



  gt <- ggplot_gtable(ggplot_build(plot))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  # grid.draw(gt)
  return(gt)
}




compress_tiff <- function(filepath){
  library(tiff)
  tiff <- readTIFF(filepath)
  writeTIFF(tiff, filepath, compression="LZW")
}



## Run

# anticodon_repetoire_plot()
# anticodon_sparing_plot()

pdf('graphs/trna/gc3_osc_density_median_cost_correlation.pdf')
cost_osc_density_correlation_plot('p1')
dev.off()


## Combined plot
library(gridExtra)
library(cowplot)
plots <- plot_grid(
  median_shift_probability_plot("p1"),
  z_median_probability_plot('p1', 'ss'),
  ncol = 1, nrow = 2, labels = c("A", "B")
)

filepath <- "graphs/trna/shiftability.tiff"
ggsave(filepath, plot=plots, width=6, height=10)
compress_tiff(filepath)

library(gridExtra)
library(cowplot)
plots <- plot_grid(
  cost_osc_density_correlation_plot("p1"),
  z_median_probability_plot('p1', 'cs'),
  ncol = 1, nrow = 2, labels = c("A", "B")
)


# ## Save as tiff to reduce file size
ggsave(file='graphs/trna/z_cds_median_frameshift_probability.tiff', plots, width=6, height=10)
compress_tiff('graphs/trna/z_cds_median_frameshift_probability.tiff')


sink('outputs/r_outputs/trna.csv')

# Correlation between gc3 and the correlation between cds cost and osc density
cat('Correlation between GC3 and the correlation between cds cost and osc density\n')
file <- read.csv('outputs/trna/cds_cost_osc_density_correlations.csv', head=T)
print(cor.test(file$gc3, file$p1_cor, method="spearman"))

# Correlation between gc3 and shift probability
cat('Correlation between GC3 and the probability of frameshifting\n')
file <- read.csv('outputs/trna/median_shift_probability.csv', head=T)
print(cor.test(file$gc3, file$median_cost_p1, method="spearman"))

## Partial correlation
library(ppcor)
file <- read.csv('outputs/trna/median_shift_probability_z.csv', head=T)
file <- file[complete.cases(file),]

cat('\nPartial correlation between z score and frameshift prob with gc3 for synonymous site model\n')
print(pcor.test(file$z_p1_ss, file$prob_p1, file$gc3, method=c("spearman")))
cat('\nPartial correlation between z score and frameshift prob with gc3 for codon shuffle model\n')
print(pcor.test(file$z_p1_cs, file$prob_p1, file$gc3, method=c("spearman")))

sink()
