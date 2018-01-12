

# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/misc/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)


theme_Publication <- function(base_size=16) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(
      plot.title = element_text(face = "bold" , size = rel(1), hjust =0.5, vjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -3),
      axis.text = element_text(),
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = c(0.9,0.9),
      legend.background = element_rect(fill=NA),
      legend.title=element_blank(),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0")
    ))
}


file <- read.csv('outputs/simulation_tests/means.csv', head=T)

sampling <- function(frame) {
  samples <- c(5,10,25,50,75,100,150,200,300,400,500)
  data <- data.frame(mean=numeric(), sd=numeric(), sample=numeric())

  count <- 0
  for (sample in samples) {
    count <- count + 1
    col <- paste('frame_',frame, '_', sample, sep='')
    sample_data <- data.frame(mean(file[[col]]), sd(file[[col]]))
    sample_data$new <- count
    colnames(sample_data) <- c('mean', 'sd', 'sample')
    data <- rbind(data, sample_data)



  }

  if(frame == "both"){
    frame_label <- 'Both'
  } else {
    frame_label <- paste('+', frame, sep='')
  }

  library(ggplot2)
  plot <- ggplot(data, aes(x=sample, y=mean)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, size=0.7) +
    geom_line() +
    geom_point(size=2) +
    scale_x_continuous(breaks=seq(1,11,1), limits=c(1,11), labels=samples) +
    labs(x="Number of randomisations", y="Mean OSC density per 100 codons", title=paste(frame_label)) +
    theme_Publication(base_size=12)

  return(plot)

}

plot1 <- sampling(1)
plot2 <- sampling(2)
plot3 <- sampling('both')

library(gridExtra)
ggsave('graphs/misc/repeat_determination_codon_shuffle.pdf', plot=grid.arrange(plot3, plot1, plot2, ncol=1), width=180, height=340, units='mm')
