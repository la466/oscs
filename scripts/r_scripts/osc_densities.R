
# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/osc_densities/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)

theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold", hjust = 0.5),
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


combined_density_plot <- function(frame) {

  t11 <- read.csv('outputs/osc_densities/combined_osc_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/combined_osc_densities_t4.csv', head=TRUE)
  t11 <- t11[t11$trans_table != 4,]

  t11$trans_table <- 'Table 11'
  t4$trans_table <- 'Table 4'


  if(frame == 'both'){
    density_metric <- paste(frame, '_density', sep='')
  } else {
    density_metric <- paste('X', frame, '_density', sep='')
  }

  data1 <- data.frame(t11$gc, t11[[density_metric]], t11$trans_table)
  data2 <- data.frame(t4$gc, t4[[density_metric]], t4$trans_table)

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)
  data.lo <- loess(data$density~data$gc, parametric = F)

  ylabel = paste('OSC density per 100 codons', sep='')

  colours <- c('blue', 'black')
  j <- order(data$gc)

  library(ggplot2)
  ggplot(data) +
    geom_point(aes(x=gc, y=density, color=factor(table), shape=factor(table)), size=1.5) +
    geom_line(aes(x=data$gc[j], y=data.lo$fitted[j]), col="red")+
    scale_y_continuous(limits = c(0,30)) +
    scale_shape_manual(values=c(16,17)) +
    scale_color_manual(values=c("black", "blue"))+
    labs(x="GC", y=ylabel) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme_Publication(base_size = 12) +
    theme(
      legend.position=c(0.85, 0.9),
      legend.title=element_blank(),
      legend.direction="vertical",
      legend.text=element_text(size=12),
      legend.key.size = unit(1, 'lines')
    )
}

combined_density_vioplot <- function(frame){

  t11 <- read.csv('outputs/osc_densities/combined_osc_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/combined_osc_densities_t4.csv', head=TRUE)
  t11 <- t11[t11$trans_table != 4,]

  if(frame == 'both'){
    density_metric <- paste(frame, '_density', sep='')
  } else {
    density_metric <- paste('X', frame, '_density', sep='')
  }

  data1 <- data.frame(t11$gc, t11[[density_metric]], t11$trans_table)
  data2 <- data.frame(t4$gc, t4[[density_metric]], t4$trans_table)

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)
  data.lo <- loess(data$density~data$gc, parametric = F)

  resids <- data.lo$residuals
  resids4 <- data.frame(resids[data$table==4], 4)
  colnames(resids4) <- c('resids', 'table')
  resids11 <- data.frame(resids[data$table==11], 11)
  colnames(resids11) <- c('resids', 'table')
  all_resids <- rbind(resids4, resids11)

  library(ggplot2)
  ggplot(all_resids, aes(factor(table), resids)) +
    theme_Publication() +
    geom_violin(aes(fill = factor(table))) +
    geom_boxplot(width=.1, outlier.colour=NA) +
    scale_fill_manual(values=c('blue', 'grey')) +
    geom_hline(yintercept=0, lty=2) +
    labs(y="Loess residual", x='') +
    scale_x_discrete(labels=c('Table 4', 'Table 11')) +
    theme(legend.position = "none")
}

combined_density_plots <- function() {
  library(gridExtra)
  plots <- grid.arrange(
    arrangeGrob(combined_density_plot('both'), combined_density_vioplot('both'), top=textGrob("Both", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    arrangeGrob(combined_density_plot(1), combined_density_vioplot(1), top=textGrob("+1", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    arrangeGrob(combined_density_plot(2), combined_density_vioplot(2), top=textGrob("+2", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    ncol=1, nrow=3
  )
  ggsave(file='graphs/osc_densities/combined_osc_densities.pdf', plots, width=8, height=12)
}


codon_density_plot <- function(codon, frame) {

  t11 <- read.csv('outputs/osc_densities/off_frame_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/off_frame_densities_t4.csv', head=TRUE)

  # Remove table 4 genomes from t11 set
  t11 <- t11[t11$trans_table != 4,]

  density_metric <- paste(codon, '_', frame, sep='')

  data1 <- data.frame(t11$gc, t11[[density_metric]], paste('Table ', t11$trans_table, sep=""))
  data2 <- data.frame(t4$gc, t4[[density_metric]], paste('Table ', t4$trans_table, sep=""))

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)

  data.lo <- loess(data$density~data$gc, parametric = F)

  ylabel = paste('OSC density per 100 codons', sep='')

  colours <- c('black', 'blue')
  j <- order(data$gc)

  library(ggplot2)
  ggplot(data) +
    geom_point(aes(x=gc, y=density, color=factor(table), shape=factor(table)), size=1.5) +
    geom_line(aes(x=data$gc[j], y=data.lo$fitted[j]), col="red")+
    scale_y_continuous(limits = c(0,7.5)) +
    scale_color_manual(values=colours) +
    scale_shape_manual(values=c(16, 17))+
    labs(x="GC", y=ylabel) +
    theme_Publication() +
    theme(
      legend.position=c(0.85, 0.9),
      legend.title=element_blank(),
      legend.direction="vertical",
      legend.text=element_text(size=12),
      legend.key.size = unit(1, 'lines')
    )

}

codon_density_vioplot <- function(codon, frame){

  t11 <- read.csv('outputs/osc_densities/off_frame_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/off_frame_densities_t4.csv', head=TRUE)

  # Remove table 4 genomes from t11 set
  t11 <- t11[t11$trans_table != 4,]

  density_metric <- paste(codon, '_', frame, sep='')

  data1 <- data.frame(t11$gc, t11[[density_metric]], t11$trans_table)
  data2 <- data.frame(t4$gc, t4[[density_metric]], t4$trans_table)

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)

  data.lo <- loess(data$density~data$gc, parametric = F)
  resids <- data.lo$residuals
  resids4 <- data.frame(resids[data$table==4], 4)
  colnames(resids4) <- c('resids', 'table')
  resids11 <- data.frame(resids[data$table==11], 11)
  colnames(resids11) <- c('resids', 'table')
  all_resids <- rbind(resids4, resids11)

  library(ggplot2)
  ggplot(all_resids, aes(factor(table), resids)) +
    theme_Publication() +
    geom_violin(aes(fill = factor(table))) +
    geom_boxplot(width=.1, outlier.colour=NA) +
    scale_fill_manual(values=c('blue', 'grey')) +
    geom_hline(yintercept=0, lty=2) +
    labs(y="Loess residual", x='') +
    scale_x_discrete(labels=c('Table 4', 'Table 11')) +
    theme(legend.position = "none")
}

codon_density_plots <- function() {
  library(gridExtra)
  plots <- grid.arrange(
    arrangeGrob(codon_density_plot('TAA', 1), codon_density_vioplot('TAA', 1), top=textGrob("+1 TAA", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    arrangeGrob(codon_density_plot('TAG', 1), codon_density_vioplot('TAG', 1), top=textGrob("+1 TAG", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    arrangeGrob(codon_density_plot('TGA', 1), codon_density_vioplot('TGA', 1), top=textGrob("+1 TGA", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    ncol=1, nrow=3
  )
  ggsave(file='graphs/osc_densities/+1_osc_densities.pdf', plots, width=8, height=12)
}

combined_density_test <- function(frame){
  t11 <- read.csv('outputs/osc_densities/combined_osc_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/combined_osc_densities_t4.csv', head=TRUE)
  t11 <- t11[t11$trans_table != 4,]


  if(frame == 'both'){
    density_metric <- paste(frame, '_density', sep='')
  } else {
    density_metric <- paste('X', frame, '_density', sep='')
  }

  data1 <- data.frame(t11$gc, t11[[density_metric]], t11$trans_table)
  data2 <- data.frame(t4$gc, t4[[density_metric]], t4$trans_table)

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)
  data.lo <- loess(data$density~data$gc, parametric = F)

  resids <- data.lo$residuals

  krus <- kruskal.test(resids~data$table)

  cat('\n==========================================\n')
  cat(sprintf('Combined OSC for frame %s\n', frame))
  print(krus)
  cat(sprintf('Mean residual t4: %s\n', mean(resids[data$table==4])))
  cat(sprintf('Mean residual t11: %s\n', mean(resids[data$table==11])))
}

codon_density_test <- function(codon, frame){

  t11 <- read.csv('outputs/osc_densities/off_frame_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/off_frame_densities_t4.csv', head=TRUE)

  # Remove table 4 genomes from t11 set
  t11 <- t11[t11$trans_table != 4,]

  density_metric <- paste(codon, '_', frame, sep='')

  data1 <- data.frame(t11$gc, t11[[density_metric]], t11$trans_table)
  data2 <- data.frame(t4$gc, t4[[density_metric]], t4$trans_table)

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)

  data.lo <- loess(data$density~data$gc, parametric = F)

  resids <- data.lo$residuals

  krus <- kruskal.test(resids~data$table)


  cat('\n==========================================\n')
  cat(sprintf('%s frame %s\n', codon, frame))
  print(krus)
  cat(sprintf('Mean residual t4: %s\n', mean(resids[data$table==4])))
  cat(sprintf('Mean residual t11: %s\n', mean(resids[data$table==11])))
}

codon_density_test_mean_t11_resid <- function(codon, frame){

  t11 <- read.csv('outputs/osc_densities/off_frame_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/off_frame_densities_t4.csv', head=TRUE)

  # Remove table 4 genomes from t11 set
  t11 <- t11[t11$trans_table != 4,]

  density_metric <- paste(codon, '_', frame, sep='')

  data1 <- data.frame(t11$gc, t11[[density_metric]], t11$trans_table)
  data2 <- data.frame(t4$gc, t4[[density_metric]], t4$trans_table)

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)

  data.lo <- loess(data$density~data$gc, parametric = F)

  resids <- data.lo$residuals
  return(mean(resids[data$table==11]))
}

codon_density_test_mean_t4_resid <- function(codon, frame){

  t11 <- read.csv('outputs/osc_densities/off_frame_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/off_frame_densities_t4.csv', head=TRUE)

  # Remove table 4 genomes from t11 set
  t11 <- t11[t11$trans_table != 4,]

  density_metric <- paste(codon, '_', frame, sep='')

  data1 <- data.frame(t11$gc, t11[[density_metric]], t11$trans_table)
  data2 <- data.frame(t4$gc, t4[[density_metric]], t4$trans_table)

  colnames(data1) <- c('gc', 'density', 'table')
  colnames(data2) <- c('gc', 'density', 'table')

  data <- rbind(data1, data2)

  data.lo <- loess(data$density~data$gc, parametric = F)

  resids <- data.lo$residuals
  return(mean(resids[data$table==4]))
}

# Run
combined_density_plots()
codon_density_plots()

# Density tests
sink('outputs/r_outputs/osc_densities.csv')

combined_density_test('both')
combined_density_test(1)
combined_density_test(2)

codon_density_test('TGA', 1)
codon_density_test('TGA', 2)
codon_density_test('TAA', 1)
codon_density_test('TAG', 1)

codon_density_test('TAC', 1)
codon_density_test('TAT', 1)
codon_density_test('TGC', 1)
codon_density_test('TGG', 1)
codon_density_test('TGT', 1)

sink()
