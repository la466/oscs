
# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/osc_densities/"), showWarnings = FALSE)
dir.create(file.path("outputs/", "r_outputs/"), showWarnings = FALSE)

combined_density_test <- function(frame){
  
  t11 <- read.csv('outputs/osc_densities/combined_osc_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/combined_osc_densities_t4.csv', head=TRUE)
  
  file <- read.csv('outputs/simulation_synonymous_site_analysis_t4/all_codons.csv', head=T)
  t4_list <- read.csv('outputs/genome_sorting/genome_list_t4.csv', head=T)
  t4_list <- t4_list[t4_list$genus == "Mycoplasma",]
  mycos_list <- file[file$acc %in% t4_list$acc,]
  mycos_list$padj <- p.adjust(mycos_list$TGA_1_pval, method="fdr")
  mycos_list <- mycos_list[order(-mycos_list$TGA_1_z),]
  mycos_list <- mycos_list$acc[1:9]

  
  # Remove table 4 genomes from t11 set
  t11 <- t11[t11$trans_table != 4,]
  
  t4_not_myco <- t4[t4$gen != "Mycoplasma",]
  t4_myco <- t4[t4$acc %in% mycos_list,]
  
  t4 <- rbind(t4_not_myco, t4_myco)

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
combined_density_test(1)

codon_density_test <- function(codon, frame){

  t11 <- read.csv('outputs/osc_densities/off_frame_densities.csv', head=TRUE)
  t4 <- read.csv('outputs/osc_densities/off_frame_densities_t4.csv', head=TRUE)
  
  file <- read.csv('outputs/simulation_synonymous_site_analysis_t4/all_codons.csv', head=T)
  t4_list <- read.csv('outputs/genome_sorting/genome_list_t4.csv', head=T)
  t4_list <- t4_list[t4_list$genus == "Mycoplasma",]
  mycos_list <- file[file$acc %in% t4_list$acc,]
  mycos_list$padj <- p.adjust(mycos_list$TGA_1_pval, method="fdr")
  mycos_list <- mycos_list[order(-mycos_list$TGA_1_z),]
  mycos_list <- mycos_list$acc[1:9]

  density_metric <- paste(codon, '_', frame, sep='')

  # Remove table 4 genomes from t11 set
  t11 <- t11[t11$trans_table != 4,]

  t4_not_myco <- t4[t4$gen != "Mycoplasma",]
  t4_myco <- t4[t4$acc %in% mycos_list,]

  t4 <- rbind(t4_not_myco, t4_myco)

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

codon_density_test('TGA', 1)




# Density tests
sink('outputs/r_outputs/osc_densities_restricted_mycoplasma.csv')

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

