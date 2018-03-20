##Scripts

All scripts to be run using using python 3, except **1/get_genomes** and **9/extract_trna** which need to run using python 2.


#####1_get_genomes

- **get_genomes**: Download genomes from EMBL using source file 'misc/bacteria_accessions.txt'

#####2_sort_genomes

- **sort_genomes**: Sort the genomes, limiting to 1 per genus greater than 500,000bp. Take all table 4 genomes.

#####3_embl_to_fasta

- **embl_to_fasta**: Filter CDSs for each of the sorted genomes and write to fasta format.

#####4_osc_stats

- **osc_stats**: Stats on the OSC densities for genome sample.

- **osc_stats_t4**: Stats on the OSC densities for the extended table 4 genome sample.

#####5_repeat_seligmann

- **codon_contribution_use**: Repeat the analysis used by Seligmann and Pollock (2004).

#####6_mmodels

- **23mm**: Second order Markov model based on Morgens et al (2013). Outputs OSC densities for 200 simulations for each genome.

- **23mm_analysis**: Analysis of results from the second order Markov model.

- **23mm_excessess**: Output excess stats of the second order Markov model.

- **53mm**: Fifth order Markov model based on Morgens et al (2013). Outputs OSC densities for 200 simulations for each genome.

- **53mm_analysis**: Analysis of results from the fifth order Markov model.

- **53mm_excessess**: Output excess stats of the fifth order Markov model.



#####7_models

- **codon_shuffle/codon_shuffle**: Randomise CDSs by shuffling the codons within each CDS. Outputs OSC densities for 200 simulations for each genome.

- **codon_shuffle/codon_shuffle_analysis**: Analysis of results from codon shuffle.

- **codon_shuffle/excesses**: Output excess stats of codon shuffle model.

- **synonymous_site_simulation/synonymous_site_simulation**: Randomise synonymous sites within the genome, preserving amino acid identities and genome codon usage biases. Outputs OSC densities for 200 simulations for each genome. Codon count excludes one-fold degenerates.

- **synonymous_site_simulation/synonymous_site_simulation_analysis**: Analysis of results from synonymous site simulations.

- **synonymous_site_simulation/excesses**: Output excess stats of synonymous site model.

- **synonymous_site_simulation/synonymous_site_simulation_t4**: For table 4 genomes. Randomise synonymous sites within the genome, preserving amino acid identities and genome codon usage biases. Outputs OSC densities for 200 simulations for each genome. Codon count excludes one-fold degenerates.

- **synonymous_site_simulation/synonymous_site_simulation_analysis_t4**: For table 4 genomes. Analysis of results from synonymous site simulations.

- **synonymous_codon_simulation/synonymous_codon_simulation**: Randomise synonymous sites within the genome, preserving amino acid identities and genome codon usage biases. Outputs OSC densities for 200 simulations for each genome. Codon count excludes one-fold degenerates.

- **synonymous_codon_simulation/synonymous_codon_simulation_analysis**: Analysis of results from synonymous site simulations.

- **synonymous_codon_simulation/excesses**: Output excess stats of synonymous codon model.

- **ecoli_codon_shuffle_500_repeats**: Generate 500 repeats for E coli sing codon shuffle model.

- **test_simulation_number**: Test the number of simulations to run. Conducted using codon shuffle model.

#####8_aa_repeat_synonymous_site_use

- **aa_repeat_synonymous_site_use**: Compare the use of nucleotides at synonymous sites for amino acids that when repeated have the oppertunity to encode an OSC.

#####9_trna

- **extract_trna**: Extract trna from file to fasta format. Input file must be located at raw_files/trna/trna_seqeunce.fasta.

- **costs**: Calculate the correlation between the costs of frameshift in CDS and OSC density for each genome.**

- **median_shift_probability**: Calculate the median shift probability for each CDS.

- **median_shift_probability_z**: Group the shift proabilities with the genome Z score from the synonymous site simulation model.

- *test_anticodon_repetoire: Sanity check comparing with Warnceke (2010). Get the number of anticodons, anticodon gene copy number for each genome.*

- *test_anticodon_sparing: Sanity check comparing with Warnceke (2010). Get the anticodons that are spared in each genone.*


#####other

- **cai/get_genomes**: Get genomes that contain the correct genes.

- **cai/write_genome_files**: Write genome files for use with CodonW.

- **cai/run_codonw**: Run CodonW.

- **cai/cai_analysis**: Analyse CodonW outputs.

- **next_nt**: Get the use of nucleotides that come after an OSC.

- **compare_z**: Compare Z scores for the codon shuffle and synoymous site models.


#####r_scripts

- **osc_densities**: Comapre OSC densities between table 11 and table 4 genomes.

- **analysis_simulation_codon_shuffle**: Analysis and plots for the codon shuffle model.

- **analysis_simulation_synonymous_sites**: Analysis and plots for the synonymous site simulaion model.

- **analysis_simulation_synonymous_codons**: Analysis and plots for the codon shuffle model (permitting changes between coding blocks).

- **analysis_simulation_synonymous_codons**: Analysis and plots for the Markov models.

- **repeat_determination**: Plot the variation in Z scores for the codon shuffle model.

- **aa_repeat_synonymous_sites**: Analysis of the use of nucleotides from amino acids repeats whos codons have the oppertunity to encode an OSC.

- **trna** Analysis for the tRNA data.

- **cai** Analysis and plot scripts for CAI data.

- **next_nt** Analysis and plot scripts for the data on the nucleotide following an OSC.

- **repeat_seligmann** Analysis and plot scripts for the data repeating the analysis used by Seligmann (2004).
