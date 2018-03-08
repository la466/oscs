

# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/simulation_synonymous_site/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)


processFile = function(filepath) {
  t4_genomes <- c()
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    t4_genomes <- c(t4_genomes, line)
  }
  close(con)
  return(t4_genomes)
}

# Get table 4 genomes
t4 <- processFile("outputs/gene_filtering/table_4_genomes_in_list.txt")

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



# Excesses plot for when the stop codons are combined
excess_plot <- function(frame) {
  file <- read.csv('outputs/simulation_synonymous_site_analysis/combined_stops.csv', head=T)

  col <- paste('osc_', frame, '_z', sep='')
  pcol <- paste('osc_', frame, '_pval', sep='')

  file$padj <- p.adjust(file[[pcol]], method="fdr")
  file$sig <- ifelse(file$padj < 0.05 & file[[col]] > 0, 0, ifelse(file$padj < 0.05 & file[[col]] < 0, 1, 2))

  library(ggplot2)
  ggplot(file) +
    geom_point(aes(x=gc, y=file[[col]], color=factor(sig))) +
    scale_y_continuous(limits = c(-200, 100)) +
    scale_color_manual(name="", values=c('blue', 'red', 'black'), labels = c("Z > 0, p < 0.05", "Z < 0, p < 0.05", "p > 0.05")) +
    labs(x="GC", y='Z') +
    geom_hline(yintercept=0, lty=2) +
    theme_Publication() +
    theme(legend.justification=c(0,1),
          legend.position=c(0,0.26),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.direction='vertical',
          legend.key.size = unit(0.6, 'lines'),
          legend.text=element_text(size=8)
    )

}

# Violin plot for genome GC content grouped by significant excesses or not
excess_vioplot <- function(frame){
  file <- read.csv('outputs/simulation_synonymous_site_analysis/combined_stops.csv', head=T)

  col <- paste('osc_', frame, '_z', sep='')
  pcol <- paste('osc_', frame, '_pval', sep='')

  file$padj <- p.adjust(file[[pcol]], method="fdr")

  file$sig <- ifelse(file$padj < 0.05 & file[[col]] > 0, 1, ifelse(file$padj < 0.05 & file[[col]] < 0, 2, 0))
  file <- file[file$sig > 0,]

  library(ggplot2)
  ggplot(file, aes(factor(sig), gc)) +
    theme_Publication() +
    geom_violin(aes(fill = factor(sig))) +
    geom_boxplot(width=.1, outlier.colour=NA) +
    scale_fill_manual(values=c('blue', 'red')) +
    labs(y="GC", x='Z score (p < 0.05)') +
    scale_x_discrete(labels=c('Z > 0', 'Z < 0')) +
    theme(legend.position = "none")
}

# Create a plot combining the excess plots and violin plots
excess_plots <- function() {
  library(ggplot2)
  library(gridExtra)
  plots <- grid.arrange(
    arrangeGrob(excess_plot('both'), excess_vioplot('both'), top=textGrob("Both", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    arrangeGrob(excess_plot(1), excess_vioplot(1), top=textGrob("+1", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    arrangeGrob(excess_plot(2), excess_vioplot(2), top=textGrob("+2", gp=gpar(fontsize=14,fontface="bold")), ncol=2),
    ncol=1, nrow=3
  )
  ggsave(file='graphs/simulation_synonymous_site/osc_excesses.pdf', plots, width=8, height=12, device=cairo_pdf)
}



# Excess plot for codon in given frame
excess_codon_plot <- function(frame, codon, min=NULL, max=NULL) {
  file <- read.csv('outputs/simulation_synonymous_site_analysis/all_codons.csv', head=T)

  col <- paste(codon, '_', frame, '_z', sep='')
  pcol <- paste(codon, '_', frame, '_pval', sep='')

  file$padj <- p.adjust(file[[pcol]], method="fdr")

  file$sig <- ifelse(file$padj < 0.05 & file[[col]] > 0, 0, ifelse(file$padj < 0.05 & file[[col]] < 0, 1, 2))

  if(frame == "both") {
    title <- paste('Both', codon)
  } else {
    title <- paste('+', frame, ' ', codon, sep='')
  }

  if(is.null(min)){
    min = min(file[[col]])
  }
  if(is.null(max)){
    max = max(file[[col]])
  }

  ggplot(file) +
    geom_point(aes(x=gc, y=file[[col]], color=factor(sig)), size=0.8) +
    scale_y_continuous(limits=c(min,max), breaks=scales::pretty_breaks(n = 10)) +
    scale_color_manual(name="", values=c('blue', 'red', 'black'), labels = c("Z > 0, p < 0.05", "Z < 0, p < 0.05", "p > 0.05")) +
    labs(x="GC", y='Z', title=title) +
    geom_hline(yintercept=0, lty=2) +
    theme_Publication() +
    theme(legend.justification=c(0,1),
          legend.position=c(0.7,1.12),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.direction='vertical',
          legend.key.size = unit(0.6, 'lines'),
          legend.text=element_text(size=8)
    )

}

# Get the minimum Z value for the stop codons
get_min_z <- function(){
  file <- read.csv('outputs/simulation_synonymous_site_analysis/all_codons.csv', head=T)
  min = 0
  for (frame in c('both', 1, 2)) {
    for (codon in c('TAA', 'TAG', 'TGA')) {
      col <- paste(codon, '_', frame, '_z', sep='')
      if(min(file[[col]])< min){
        min = min(file[[col]])
      }
    }
  }
  return(min)
}

# Get the maximum Z value for the stop codons
get_max_z <- function(){
  file <- read.csv('outputs/simulation_synonymous_site_analysis/all_codons.csv', head=T)

  max = 0
  for (frame in c('both', 1, 2)) {
    for (codon in c('TAA', 'TAG', 'TGA')) {
      col <- paste(codon, '_', frame, '_z', sep='')
      if(max(file[[col]]) > max){
        max = max(file[[col]])
      }
    }
  }
  return(max)
}

# Generate plot for each of the stop codons in each of the reading frames
excess_individual_osc_plots <- function(min, max) {
  plots <- arrangeGrob(
    # excess_codon_plot('both', 'TAA', min, max),
    # excess_codon_plot('both', 'TAG', min, max),
    # excess_codon_plot('both', 'TGA', min, max),
    excess_codon_plot(1, 'TAA', min, max),
    excess_codon_plot(1, 'TAG', min, max),
    excess_codon_plot(1, 'TGA', min, max),
    excess_codon_plot(2, 'TAA', min, max),
    excess_codon_plot(2, 'TAG', min, max),
    excess_codon_plot(2, 'TGA', min, max),
    ncol=3
  )
  ggsave(file='graphs/simulation_synonymous_site/osc_individual_excesses_+1_+2.pdf', plots, width=12, height=8)
}


# Vioplot of genomes GC content with significnat or not signincant excesses
individual_codon_vioplot <- function(frame, codon, grouped=FALSE, count) {

  file <- read.csv('outputs/simulation_synonymous_site_analysis/all_codons.csv', head=T)

  col <- paste(codon, '_', frame, '_z', sep='')
  pcol <- paste(codon, '_', frame, '_pval', sep='')

  file$padj <- p.adjust(file[[pcol]], method="fdr")

  file$sig <- ifelse(file$padj < 0.05 & file[[col]] > 0, 1, ifelse(file$padj < 0.05 & file[[col]] < 0, 2, 0))
  file <- file[file$sig > 0,]

  if(frame=="both"){
    title <- paste('Both ', codon, sep='')
  } else {
    title <- paste('+', frame, ' ', codon, sep='')
  }

  library(ggplot2)
  ggplot(file, aes(factor(sig), gc)) +
    theme_Publication() +
    geom_violin(aes(fill = factor(sig))) +
    geom_boxplot(width=.1, outlier.colour=NA) +
    scale_fill_manual(values=c('blue', 'red')) +
    labs(y="GC", x='Z score (p < 0.05)', title=title) +
    scale_x_discrete(labels=c('Z > 0', 'Z < 0')) +
    theme(legend.position = "none")
}

# Plot the vioplots
osc_vioplots <- function() {
  plots <- arrangeGrob(
    # individual_codon_vioplot('both', 'TAA'),
    # individual_codon_vioplot('both', 'TAG'),
    # individual_codon_vioplot('both', 'TGA'),
    individual_codon_vioplot(1, 'TAA'),
    individual_codon_vioplot(1, 'TAG'),
    individual_codon_vioplot(1, 'TGA'),
    individual_codon_vioplot(2, 'TAA'),
    individual_codon_vioplot(2, 'TAG'),
    individual_codon_vioplot(2, 'TGA'),
    ncol=3
  )

  ggsave(file='graphs/simulation_synonymous_site/osc_individual_gc_vioplots.pdf', plots, width=12, height=8, device=cairo_pdf)
}

# Calculate the number of genomes with an excess of combined OSCs for the given frame
combined_stops_excess <- function(frame) {
  file <- read.csv('outputs/simulation_synonymous_site_analysis/combined_stops.csv', head=T)

  col <- paste('osc_', frame, '_z', sep='')
  pcol <- paste('osc_', frame, '_pval', sep='')

  file$padj <- p.adjust(file[[pcol]], method="fdr")

  cat("\n")
  print(frame)
  print(cor.test(file$gc, file[[col]], method="spearman"))
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,]))
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,])/nrow(file)*100)

}



# Calculate the number of genomes with excesses for frame and codon
all_codons_excess <- function(frame, codon) {
  file <- read.csv('outputs/simulation_synonymous_site_analysis/all_codons.csv', head=T)

  col <- paste(codon, '_', frame, '_z', sep='')
  pcol <- paste(codon, '_', frame, '_pval', sep='')

  head(file)

  file$padj <- p.adjust(file[[pcol]], method="fdr")

  cat("\n")
  print(frame)
  print(codon)
  print(cor.test(file$gc, file[[col]], method="spearman"))
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,]))
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,])/nrow(file)*100)
}






## Run

# Codons together
excess_plots()

# Individual codons
min = get_min_z()
max = get_max_z()
excess_individual_osc_plots(min, max)
osc_vioplots()

sink('outputs/r_outputs/synonymous_site_model.txt')
combined_stops_excess('both')
combined_stops_excess(1)
combined_stops_excess(2)
all_codons_excess("both", 'TAA')
all_codons_excess("both", 'TAC')
all_codons_excess("both", 'TAG')
all_codons_excess("both", 'TAT')
all_codons_excess("both", 'TGA')
all_codons_excess("both", 'TGC')
all_codons_excess("both", 'TGG')
all_codons_excess("both", 'TGT')
all_codons_excess(1, 'TAA')
all_codons_excess(1, 'TAC')
all_codons_excess(1, 'TAG')
all_codons_excess(1, 'TAT')
all_codons_excess(1, 'TGA')
all_codons_excess(1, 'TGC')
all_codons_excess(1, 'TGG')
all_codons_excess(1, 'TGT')
all_codons_excess(2, 'TAA')
all_codons_excess(2, 'TAC')
all_codons_excess(2, 'TAG')
all_codons_excess(2, 'TAT')
all_codons_excess(2, 'TGA')
all_codons_excess(2, 'TGC')
all_codons_excess(2, 'TGG')
all_codons_excess(2, 'TGT')
sink()
