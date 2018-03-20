

# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/markov_models/"), showWarnings = FALSE)
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
            axis.title = element_text(size = rel(1)),
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
excess_plot <- function(model, frame, main=NULL) {
  filePath <- paste('outputs/mmodels_analysis/', model, '/combined_stops.csv', sep="")
  file <- read.csv(filePath, head=T)

  col <- paste('osc_', frame, '_z', sep='')
  pcol <- paste('osc_', frame, '_pval', sep='')

  file$padj <- p.adjust(file[[pcol]], method="fdr")
  file$sig <- ifelse(file$padj < 0.05 & file[[col]] > 0, 0, ifelse(file$padj < 0.05 & file[[col]] < 0, 1, 2))

  library(ggplot2)
  plot <- ggplot(file) +
    geom_point(aes(x=gc, y=file[[col]] )) +
    scale_y_continuous(limits = c(-0.5, 1.5)) +
    labs(x="GC", y='Z') +
    geom_hline(yintercept=0, lty=2) +
    theme_Publication() +
    theme(legend.justification=c(0,1),
          legend.position=c(0,1),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.direction='vertical',
          legend.key.size = unit(0.6, 'lines'),
          legend.text=element_text(size=8)
    )

  if(!is.null(main)){
    plot <- plot + ggtitle(main)
  }
  plot
}



# Create a plot combining the excess plots and violin plots
excess_plots <- function(model1, model2) {
  library(ggplot2)
  library(gridExtra)
  plots <- grid.arrange(
    arrangeGrob(excess_plot(model1, 'both', '23MM'), excess_plot(model2, 'both', '53MM'), ncol=2, top="Both"),
    arrangeGrob(excess_plot(model1, 1), excess_plot(model2, 1), ncol=2, top="+1"),
    arrangeGrob(excess_plot(model1, 2), excess_plot(model2, 2), ncol=2, top="+2"),
    ncol=1, nrow=3
  )
  filePath <- paste('graphs/markov_models/osc_excesses.pdf', sep='')
  ggsave(file=filePath, plots, width=8, height=12, device=cairo_pdf)
}


# Excess plot for codon in given frame
excess_codon_plot <- function(model, frame, codon, min=NULL, max=NULL) {

  filePath <- paste('outputs/mmodels_analysis/', model, '/all_codons.csv', sep="")
  file <- read.csv(filePath, head=T)


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

  library(ggplot2)
  ggplot(file) +
    geom_point(aes(x=gc, y=file[[col]], color=factor(sig)), size=0.8) +
    scale_y_continuous(limits=c(min,max), breaks=scales::pretty_breaks(n = 10)) +
    scale_color_manual(name="", values=c("0" = 'blue', "1"= 'red', "2" = 'black'), labels = c("0" = "Z > 0, p < 0.05", "1" = "Z < 0, p < 0.05", "2" = "p > 0.05")) +
    labs(x="GC", y='Z', title=title) +
    geom_hline(yintercept=0, lty=2) +
    theme_Publication() +
    theme(legend.justification=c(0,1),
          legend.position=c(0,1),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.direction='vertical',
          legend.key.size = unit(0.6, 'lines'),
          legend.text=element_text(size=8)
    )
}



# Get the minimum Z value for the stop codons
get_min_z <- function(model){
  filepath <- paste('outputs/mmodels_analysis/', model, '/all_codons.csv', sep="")
  file <- read.csv(filepath, head=T)
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
get_max_z <- function(model){
  filepath <- paste('outputs/mmodels_analysis/', model, '/all_codons.csv', sep="")
  file <- read.csv(filepath, head=T)

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
excess_individual_osc_plots <- function(model, min, max) {
  plots <- arrangeGrob(
    excess_codon_plot(model, 1, 'TAA', min, max),
    excess_codon_plot(model, 1, 'TAG', min, max),
    excess_codon_plot(model, 1, 'TGA', min, max),
    excess_codon_plot(model, 2, 'TAA', min, max),
    excess_codon_plot(model, 2, 'TAG', min, max),
    excess_codon_plot(model, 2, 'TGA', min, max),
    ncol=3
  )
  ggsave(file=paste('graphs/markov_models/', model, '_osc_individual_excesses_+1_+2.pdf', sep=''), plots, width=12, height=8)
}

# Calculate the number of genomes with an excess of combined OSCs for the given frame
combined_stops_excess <- function(model, frame) {
  filePath <- paste('outputs/mmodels_analysis/', model, '/combined_stops.csv', sep="")
  file <- read.csv(filePath, head=T)


  col <- paste('osc_', frame, '_z', sep='')
  pcol <- paste('osc_', frame, '_pval', sep='')

  file$padj <- p.adjust(file[[pcol]], method="fdr")


  cat('\n')
  print(frame)
  print(cor.test(file$gc, file[[col]], method="spearman"))
  print('Excess')
  print(nrow(file[file[[col]] > 0,]))
  print(nrow(file[file[[col]] > 0 ,])/nrow(file)*100)
  print('Significant excess')
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,]))
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,])/nrow(file)*100)
}




# Calculate the number of genomes with excesses for frame and codon
all_codons_excess <- function(model, frame, codon) {
  filePath <- paste('outputs/mmodels_analysis/', model, '/all_codons.csv', sep="")
  file <- read.csv(filePath, head=T)

  col <- paste(codon, '_', frame, '_z', sep='')
  pcol <- paste(codon, '_', frame, '_pval', sep='')

  head(file)

  file$padj <- p.adjust(file[[pcol]], method="fdr")

  cat('\n')
  print(frame)
  print(codon)
  print(cor.test(file$gc, file[[col]], method="spearman"))
  print('Excess')
  print(nrow(file[file[[col]] > 0,]))
  print(nrow(file[file[[col]] > 0 ,])/nrow(file)*100)
  print('Significant excess')
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,]))
  print(nrow(file[file[[col]] > 0 & file$padj < 0.05,])/nrow(file)*100)
}


## Run

excess_plots("23mm", '53mm')
excess_individual_osc_plots('23mm', get_min_z('23mm'), get_max_z('23mm'))
excess_individual_osc_plots('53mm', get_min_z('53mm'), get_max_z('53mm'))

sink('outputs/r_outputs/23mm_model.txt')
combined_stops_excess('23mm', 'both')
combined_stops_excess('23mm', 1)
combined_stops_excess('23mm', 2)
all_codons_excess("23mm", "both", 'TAA')
all_codons_excess("23mm", "both", 'TAC')
all_codons_excess("23mm", "both", 'TAG')
all_codons_excess("23mm", "both", 'TAT')
all_codons_excess("23mm", "both", 'TGA')
all_codons_excess("23mm", "both", 'TGC')
all_codons_excess("23mm", "both", 'TGG')
all_codons_excess("23mm", "both", 'TGT')
all_codons_excess("23mm", 1, 'TAA')
all_codons_excess("23mm", 1, 'TAC')
all_codons_excess("23mm", 1, 'TAG')
all_codons_excess("23mm", 1, 'TAT')
all_codons_excess("23mm", 1, 'TGA')
all_codons_excess("23mm", 1, 'TGC')
all_codons_excess("23mm", 1, 'TGG')
all_codons_excess("23mm", 1, 'TGT')
all_codons_excess("23mm", 2, 'TAA')
all_codons_excess("23mm", 2, 'TAC')
all_codons_excess("23mm", 2, 'TAG')
all_codons_excess("23mm", 2, 'TAT')
all_codons_excess("23mm", 2, 'TGA')
all_codons_excess("23mm", 2, 'TGC')
all_codons_excess("23mm", 2, 'TGG')
all_codons_excess("23mm", 2, 'TGT')
sink()

sink('outputs/r_outputs/53mm_model.txt')
combined_stops_excess('53mm', 'both')
combined_stops_excess('53mm', 1)
combined_stops_excess('53mm', 2)
all_codons_excess("53mm", "both", 'TAA')
all_codons_excess("53mm", "both", 'TAC')
all_codons_excess("53mm", "both", 'TAG')
all_codons_excess("53mm", "both", 'TAT')
all_codons_excess("53mm", "both", 'TGA')
all_codons_excess("53mm", "both", 'TGC')
all_codons_excess("53mm", "both", 'TGG')
all_codons_excess("53mm", "both", 'TGT')
all_codons_excess("53mm", 1, 'TAA')
all_codons_excess("53mm", 1, 'TAC')
all_codons_excess("53mm", 1, 'TAG')
all_codons_excess("53mm", 1, 'TAT')
all_codons_excess("53mm", 1, 'TGA')
all_codons_excess("53mm", 1, 'TGC')
all_codons_excess("53mm", 1, 'TGG')
all_codons_excess("53mm", 1, 'TGT')
all_codons_excess("53mm", 2, 'TAA')
all_codons_excess("53mm", 2, 'TAC')
all_codons_excess("53mm", 2, 'TAG')
all_codons_excess("53mm", 2, 'TAT')
all_codons_excess("53mm", 2, 'TGA')
all_codons_excess("53mm", 2, 'TGC')
all_codons_excess("53mm", 2, 'TGG')
all_codons_excess("53mm", 2, 'TGT')
sink()
