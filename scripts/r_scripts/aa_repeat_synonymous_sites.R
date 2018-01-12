

# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/aa_repeat_synonymous_sites/"), showWarnings = FALSE)
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


use_plot <- function(frame, codon) {

  file <- read.csv('outputs/aa_repeat_synonymous_sites/aa_repeat_osc_synonymous_sites.csv', head=T)

  if(frame == 1) {
    query_nt = 'a'
    ylabel <- 'Log(A3:A6) proportions'
  } else {
    query_nt = 't'
    ylabel <- 'Log(T3:T6) proportions'
  }

  col1 <- paste(codon, '_', frame, '_', query_nt, '3', sep='')
  col2 <- paste(codon, '_', frame, '_', query_nt, '6', sep='')

  file$ratio <- log((file[[col1]]/file[[col2]]))

  file$ratio[!is.finite(file$ratio)] <- NA

  title <- paste('+', frame, ' ', codon, sep='')

  fit <- lm(file$ratio~file$gc3)

  if(frame == 1){
    ylims = c(-4,3)
  } else {
    ylims = c(-4,2)
  }

  library(ggplot2)
  ggplot(file) +
    geom_point(aes(x=gc3, y=ratio), size=1) +
    labs(x="GC3", y=ylabel, title=title) +
    scale_y_continuous(limits = ylims) +
    geom_hline(yintercept=0, lty=2) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color="red") +
    theme_Publication(base_size = 12)

}

use_plot_at <- function(frame, codon) {

  file <- read.csv('outputs/aa_repeat_synonymous_sites/aa_repeat_osc_synonymous_sites_only_at.csv', head=T)

  query_nt = 'a'
  ylabel <- 'Log(A3:A6) proportions'

  col1 <- paste(codon, '_', frame, '_', query_nt, '3', sep='')
  col2 <- paste(codon, '_', frame, '_', query_nt, '6', sep='')


  file$ratio <- log((file[[col1]]/file[[col2]]))

  file$ratio[!is.finite(file$ratio)] <- NA

  title <- paste('+', frame, ' ', codon, ' (A/T restricted)', sep='')

  fit <- lm(file$ratio~file$gc3)

  library(ggplot2)
  ggplot(file) +
    geom_point(aes(x=gc3, y=ratio), size=1) +
    labs(x="GC3", y=ylabel, title=title) +
    geom_hline(yintercept=0, lty=2) +
    theme_Publication(base_size = 12) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color="red")

}

aa_repeat_all_nts <- function() {
  library(gridExtra)
  plots <- arrangeGrob(
    use_plot(1, 'TAA'),
    use_plot(1, 'TAG'),
    use_plot_at(1, 'TAA'),
    use_plot(2, 'TAA'),
    use_plot(2, 'TAG'),
    use_plot(2, 'TGA'),
    ncol=3
  )
  ggsave(file='graphs/aa_repeat_synonymous_sites/synonymous_site_use.pdf', plots, width=12, height=8)
}


use_stats <- function(frame, codon) {

  file <- read.csv('outputs/aa_repeat_synonymous_sites/aa_repeat_osc_synonymous_sites.csv', head=T)

  if(frame == 1) {
    query_nt = 'a'
  } else {
    query_nt = 't'
  }

  col1 <- paste(codon, '_', frame, '_', query_nt, '3', sep='')
  col2 <- paste(codon, '_', frame, '_', query_nt, '6', sep='')

  file$ratio <- log((file[[col1]]/file[[col2]]))

  file$ratio[!is.finite(file$ratio)] <- NA

  cat('\n==================================================\n')
  cat(sprintf('+%s %s\n', frame, codon))
  print(wilcox.test(file[[col1]], file[[col2]], paired=T))
  cat(sprintf('Mean site 3: %s\n', mean(file[[col1]])))
  cat(sprintf('Mean site 6: %s\n', mean(file[[col2]])))
  print(cor.test(file$gc3, file$ratio, method="spearman"))
  cat(sprintf('Number where site 3 > site 6: %s\n', nrow(file[file$ratio > 0,])))
  cat(sprintf('Prop where site 3 > site 6: %s\n', nrow(file[file$ratio > 0,]) / nrow(file)))

}


use_stats_at <- function(frame, codon) {

  file <- read.csv('outputs/aa_repeat_synonymous_sites/aa_repeat_osc_synonymous_sites_only_at.csv', head=T)


  if(frame == 1){
    query_nt = 'a'
  } else {
    query_nt = 't'
  }


  col1 <- paste(codon, '_', frame, '_', query_nt, '3', sep='')
  col2 <- paste(codon, '_', frame, '_', query_nt, '6', sep='')

  file$ratio <- log((file[[col1]]/file[[col2]]))

  file$ratio[!is.finite(file$ratio)] <- NA

  cat('\n==================================================\n')
  cat(sprintf('+%s %s\n', frame, codon))
  print(wilcox.test(file[[col1]], file[[col2]], paired=T))
  cat(sprintf('Mean site 3: %s\n', mean(file[[col1]])))
  cat(sprintf('Mean site 6: %s\n', mean(file[[col2]])))
  print(cor.test(file$gc3, file$ratio, method="spearman"))
  cat(sprintf('Number where site 3 > site 6: %s\n', nrow(file[file$ratio > 0,])))
  cat(sprintf('Prop where site 3 > site 6: %s\n', nrow(file[file$ratio > 0,]) / nrow(file)))
}

one_tail_test <- function(frame, codon) {

  file <- read.csv('outputs/aa_repeat_synonymous_sites/aa_repeat_osc_synonymous_sites.csv', head=T)

  if(frame == 1) {
    query_nt = 'a'
  } else {
    query_nt = 't'
  }

  col1 <- paste(codon, '_', frame, '_', query_nt, '3', sep='')
  col2 <- paste(codon, '_', frame, '_', query_nt, '6', sep='')

  test <- wilcox.test(file[[col1]], file[[col2]], paired=T, alternative = "g")

  return(test$p.value)
}




## Run

aa_repeat_all_nts()

sink('outputs/r_outputs/aa_repeat_synonymous_sites.csv')
cat(sprintf('Tests between A use at sites 3/6 of codon repeats that can encode an OSC\n'))
cat(sprintf('col1 = site 3 use, col2 = site 6 use\n\n'))


use_stats(1, 'TAA')
use_stats(1, 'TAG')
use_stats(2, 'TAA')
use_stats(2, 'TAG')
use_stats(2, 'TGA')

cat('\n==================================================\n')
cat('\nOnly allowing for A/T\n')

use_stats_at(1, 'TAA')
use_stats_at(1, 'TAG')
# Cant do for +2

TAA1 <- one_tail_test(1, 'TAA')
TAG1 <- one_tail_test(1, 'TAG')
TAA2 <- one_tail_test(2, 'TAA')
TAG2 <- one_tail_test(2, 'TAG')
TGA2 <- one_tail_test(2, 'TGA')

pvals <- c(TAA1, TAG1, TAA2, TAG2, TGA2)
p.adj <- p.adjust(pvals, method="fdr")

cat('\n==================================================\n')
cat('Combined P value test\n')
cat(pchisq(-2*sum(log(p.adj)), 10, lower.tail=FALSE))

sink()
