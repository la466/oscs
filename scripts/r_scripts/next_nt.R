

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
      # legend.position = "bottom",
      legend.direction = "vertical",
      legend.key.size= unit(0.4, "cm"),
      # legend.margin = unit(0, "cm"),
      legend.title = element_blank(),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
    ))
}

library(reshape2)
file <- read.csv('outputs/other/nt_after_osc.csv', head=T)
file <- subset(file, select = -c(gc) )
file <- melt(file, id.vars=c("acc","gc2"))
colnames(file) <- c('acc','gc2','nt','value')
file$nt <- toupper(file$nt)


plot <- ggplot(file) +
  geom_point(aes(x=gc2, y=log(value), colour=factor(nt), shape=factor(nt)), size=1.2) +
  scale_color_manual(values = c("#F5AB35", "red", "blue", "black")) +
  scale_shape_manual(values = c(15,16,17,18)) +
  geom_hline(yintercept=0, lty=2) +
  labs(x="GC2", y="Log ratio\nproportion of nt +1 to OSC:proportion of codon second positions with nt") +
  theme_Publication() +
  theme(legend.position = c(0.1,0.2))

ggsave('graphs/other/next_nt.pdf', plot=plot, width=7, height=7)
