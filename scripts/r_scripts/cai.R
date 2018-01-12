

# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/cai/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)


library(reshape2)
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

draw_plot <- function(all, xcol, ycol, group, xlabel, ylabel, colours) {

  # Order the x variable
  j<-order(xcol)

  # Loess regression
  data.lo <- lm(ycol~xcol)

  # Plot x against y, with the loess curve, coloured by group
  library(ggplot2)
  ggplot(all, aes(x=xcol, y=ycol, color=factor(group), shape=factor(group))) +
    geom_point(size=1.5) +
    geom_smooth(method="lm", fill=NA, size=0.5) +
    scale_color_manual(values=colours)+
    scale_shape_manual(values = c(16,17)) +
    labs(x=xlabel, y=ylabel) +
    # guides(colour = guide_legend(override.aes = list(size=3)), size="size") +
    guides(fill = guide_legend(override.aes = list(size=3)), colour = guide_legend(override.aes = list(linetype=c(0,0)), shape=c(16,17))) +
    theme_Publication(base_size = 12) +
    theme(
      legend.position=c(0.85, 0.9),
      legend.title=element_blank(),
      legend.direction="vertical",
      legend.key.size = unit(0.8, 'lines'),
      legend.text=element_text(size=12)
    )
}

get_data1 <- function(){
  file <- read.csv('outputs/cai/cai_osc_densities.csv', head=T)
  high <- data.frame(file$gc, file$high_osc_density_p1, 'high')
  low <- data.frame(file$gc, file$low_osc_density_p1, 'low')
  colnames(high) <- c('gc', 'density', 'group')
  colnames(low) <- c('gc', 'density', 'group')
  all <- rbind(high, low)
  return(all)
}

all1 <- get_data1()

plot <- draw_plot(all1, all1$gc, all1$density, all1$group, 'GC', 'OSC density per 100 codons', c('black', 'blue'))
ggsave('graphs/cai/osc_density_cai_p1.pdf', plot=plot, width=7, height=7)

get_data_p1 <- function(){
  file <- read.csv('outputs/cai/cai_osc_densities.csv', head=T)
  high <- data.frame(file$gc, file$high_osc_density_p1, 'high')
  low <- data.frame(file$gc, file$low_osc_density_p1, 'low')
  colnames(high) <- c('gc', 'density', 'group')
  colnames(low) <- c('gc', 'density', 'group')
  all <- rbind(high, low)
  return(all)
}
get_data_m1 <- function(){
  file <- read.csv('outputs/cai/cai_osc_densities.csv', head=T)
  high <- data.frame(file$gc, file$high_osc_density_m1, 'high')
  low <- data.frame(file$gc, file$low_osc_density_m1, 'low')
  colnames(high) <- c('gc', 'density', 'group')
  colnames(low) <- c('gc', 'density', 'group')
  all <- rbind(high, low)
  return(all)
}

data_p1 <- get_data_p1()

loess_test <- function(data, xcol, ycol) {

  # Order the x variable
  j<-order(data[[xcol]])

  # Loess regression
  data.lo <- loess(data[[ycol]]~data[[xcol]], parametric = F)

  # Get residuals
  resids <- data.lo$residuals

  krus <- kruskal.test(resids~data$group)
  cat('Loess test between OSC densities for high CAI genes and low CAI genes\n')
  print(krus)
  cat(sprintf('Mean residual high: %s\n', mean(resids[data$group=='high'])))
  cat(sprintf('Mean residual low: %s\n', mean(resids[data$group=='low'])))
}

sink('outputs/r_outputs/cai.csv')
loess_test(data_p1, "gc", "density")
sink()
