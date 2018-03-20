
# Remove all variables
closeAllConnections()
rm(list=ls(all=TRUE))

# Create directory if not already created
dir.create(file.path("graphs/osc_densities/"), showWarnings = FALSE)
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

combined_density_test <- function(frame){
  t11 <- read.csv('outputs/osc_densities/combined_osc_densities.csv', head=TRUE)
  # t4 <- read.csv('outputs/osc_densities/combined_osc_densities_t4.csv', head=TRUE)
  t4 <- t11[t11$trans_table == 4,]
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
  # t4 <- read.csv('outputs/osc_densities/combined_osc_densities_t4.csv', head=TRUE)
  t4 <- t11[t11$trans_table == 4,]
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

# Run

# Density tests
sink('outputs/r_outputs/osc_densities_original_t4.csv')
combined_density_test('both')
combined_density_test(1)
combined_density_test(2)
codon_density_test('TGA', 1)
codon_density_test('TGA', 2)
sink()
