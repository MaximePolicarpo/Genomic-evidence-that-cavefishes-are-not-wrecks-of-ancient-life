#### Libraries  ---------------------------------

library(plotrix)
library(data.table)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(ggpubr)
library(gtools)
library(forcats)
library(FSA)
library("patchwork")
library("quantable")
library("reticulate")
library("sjPlot")
library(gridExtra)
library(ape)
library(phytools)
library(seqinr)
library(ggthemes)

############### SECTION VISION ############### 

#### Import datas  ---------------------------------

setwd("~/Desktop/New_Cavefishes_Analysis/Analysis_Scripts/Vision")
Vision_Sequences <- read.table("Vision_Sequences_Table.tsv", header=TRUE, sep="\t")
Vision_LoFs <- read.table("LoF_Table_Vision.tsv", header=TRUE, sep="\t")


Typhlichthys_subterraneus_table <- Vision_Sequences %>% 
  filter(Species == "Typhlichthys_subterraneus") %>%
  filter(Gene != "Not_Found")

Percopsis_transmontana_table <- Vision_Sequences %>% 
  filter(Species == "Percopsis_transmontana") %>%
  filter(Gene != "Not_Found")

Lamprologus_lethops_table <- Vision_Sequences %>% 
  filter(Species == "Lamprologus_lethops") %>%
  filter(Gene != "Not_Found")




#### Perform alignments  ---------------------------------

setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Vision")

#trimNonHomologousFragments
system("for i in *.fa ; do java -jar /Users/maxime/Downloads/macse_v2.03.jar -prog trimNonHomologousFragments -seq $i -min_homology_to_keep_seq 0 ; done")
system("rm *_AA.fa ; rm *_detail_NT.fa ; rm *.csv")

#alignSequences
system("for i in *_NT.fa ; do java -jar /Users/maxime/Downloads/macse_v2.03.jar -prog alignSequences -seq $i -out_NT $i.NT_ALN.fasta -out_AA $i.AA_ALN.fasta ; done")
system("rm *.AA_ALN.fasta")

#exportAlignment and rename sequences 
system("for i in *.NT_ALN.fasta ; do java -jar /Users/maxime/Downloads/macse_v2.03.jar -prog exportAlignment -align $i -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS - -out_NT $i.align_noFS_NT.fasta -out_AA $i.align_noFS_AA.fasta ; done")
system("for i in *_noFS_AA.fasta ; do sed -i '' 's/Danio rerio.*/Danio_rerio/g' $i ; sed -i '' 's/Typhlichthys_subterraneus.*/Typhlichthys_subterraneus/g' $i ; sed -i '' 's/Gadus_morhua.*/Gadus_morhua/g' $i ; sed -i '' 's/Percopsis_transmontana.*/Percopsis_transmontana/g' $i ; sed -i '' 's/Lamprologus_lethops.*/Lamprologus_lethops/g' $i ; sed -i '' 's/Lamprologus_tigripic.*/Lamprologus_tigripic/g' $i ; sed -i '' 's/Neolamprologus_brichardi.*/Neolamprologus_brichardi/g' $i ; sed -i '' 's/Astyanax_mexicanus_Surface.*/Astyanax_mexicanus_Surface/g' $i ; sed -i '' 's/Astyanax_mexicanus_Cave.*/Astyanax_mexicanus_Cave/g' $i ; sed -i '' 's/Pygocentrus_nattereri.*/Pygocentrus_nattereri/g' $i ; sed -i '' 's/Lucifuga_dentata.*/Lucifuga_dentata/g' $i ; sed -i '' 's/Lucifuga_holguinensis.*/Lucifuga_holguinensis/g' $i ; sed -i '' 's/Brotula_barbata.*/Brotula_barbata/g' $i ; sed -i '' 's/Carapus_acus.*/Carapus_acus/g' $i ; sed -i '' 's/Lamprogrammus_exutus.*/Lamprogrammus_exutus/g' $i ; done")
system("for i in *_noFS_NT.fasta ; do sed -i '' 's/Danio rerio.*/Danio_rerio/g' $i ; sed -i '' 's/Typhlichthys_subterraneus.*/Typhlichthys_subterraneus/g' $i ; sed -i '' 's/Gadus_morhua.*/Gadus_morhua/g' $i ; sed -i '' 's/Percopsis_transmontana.*/Percopsis_transmontana/g' $i ; sed -i '' 's/Lamprologus_lethops.*/Lamprologus_lethops/g' $i ; sed -i '' 's/Lamprologus_tigripic.*/Lamprologus_tigripic/g' $i ; sed -i '' 's/Neolamprologus_brichardi.*/Neolamprologus_brichardi/g' $i ; sed -i '' 's/Astyanax_mexicanus_Surface.*/Astyanax_mexicanus_Surface/g' $i ; sed -i '' 's/Astyanax_mexicanus_Cave.*/Astyanax_mexicanus_Cave/g' $i ; sed -i '' 's/Pygocentrus_nattereri.*/Pygocentrus_nattereri/g' $i ; sed -i '' 's/Lucifuga_dentata.*/Lucifuga_dentata/g' $i ; sed -i '' 's/Lucifuga_holguinensis.*/Lucifuga_holguinensis/g' $i ; sed -i '' 's/Brotula_barbata.*/Brotula_barbata/g' $i ; sed -i '' 's/Carapus_acus.*/Carapus_acus/g' $i ; sed -i '' 's/Lamprogrammus_exutus.*/Lamprogrammus_exutus/g' $i ; done")
system("for i in *_noFS_NT.fasta ; do sed -i '' 's/Danio_rerio.*/Danio_rerio/g' $i ; done")
system("for i in *_noFS_AA.fasta ; do sed -i '' 's/Danio_rerio.*/Danio_rerio/g' $i ; done")
system("for i in *_noFS_NT.fasta ; do sed -i '' 's/!/-/g' $i ; done")
system("for i in *_noFS_AA.fasta ; do sed -i '' 's/!/-/g' $i ; done")





#Concatenate every alignments using AMAS

#system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i *_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Vision_Genes_Concatenated_Alignment.fa")





#### Run aaml on alignments and mutpred2  ---------------------------------

setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Vision")


#Create proper nwk file for every alignments
my_species_tree <- read.tree("Species_tree_unrooted.nwk")

alignment_files <- list.files()

alignment_files <- alignment_files[grepl("noFS_AA.fasta", alignment_files)]


for(alignments in alignment_files){
  file_name_align_nwk <- paste0(alignments, ".nwk")
  
  species_present <- read.alignment(alignments, format = "fasta")$nam
  
  current_species_tree <- keep.tip(my_species_tree, species_present)
  
  write.tree(current_species_tree, file = file_name_align_nwk, append = FALSE,
             digits = 10, tree.names = FALSE)
}


#Create ctl files for every alignments

system('for i in *_noFS_AA.fasta ; do sed "s/exorh.fas/$i/g" ../aaml_control_file.ctl | sed "s/exorh.tree/$i.nwk/g" | sed "s/exorh.amml/$i.aaml/g" > $i.CTL_AMML ; done')

#Run codeml and save rst file

system("for i in *.CTL_AMML ; do yes | codeml $i ; mv rst $i.RST_FILE ; done")



#For each species and genes, extract the ancestral aa sequence and the species sequence

system("./From_RST_to_Mp2.bash")
system("mv *_input ../../Vision_Genes_Mutpred/")


#Run homemade python script to generate mutpred2 input files for each species


setwd("~/Desktop/These_Supplemental_Analysis/Vision_Genes_Mutpred")

system('for file in *_input ; do if grep -q "subterraneus" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Typhlichthys_subterraneus_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "rerio" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Danio_rerio_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "morhua" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Gadus_morhua_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "transmontana" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Percopsis_transmontana_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "lethops" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Lamprologus_lethops_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "tigripic" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Lamprologus_tigripic_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "brichardi" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Neolamprologus_brichardi_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "Lamprogrammus" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Lamprogrammus_exutus_Vision_Mutpred2.input ; done')
system('for file in *_input ; do if grep -q "Carapus" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Carapus_acus_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "barbata" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Brotula_barbata_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "holguinensis" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Lucifuga_holguinensis_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "dentata" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Lucifuga_dentata_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "mexicanus_Surface" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Astyanax_SF_Vision_Mutpred2.input ; done') 
system('for file in *_input ; do if grep -q "mexicanus_Cave" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Astyanax_CF_Vision_Mutpred2.input ; done')
system('for file in *_input ; do if grep -q "Pygocentrus" $file ; then python ../From_rst_to_mutpred2.py -f $file ; fi >> Pygocentrus_nattereri_Vision_Mutpred2.input ; done')
system("rm *.Typhlichthys_subterraneus_input ; rm *.Danio_rerio_input ; rm *.Gadus_morhua_input ; rm *.Percopsis_transmontana_input ; rm *.Lamprologus_lethops_input ; rm *.Lamprologus_tigripic_input  ; rm *.Neolamprologus_brichardi_input ; rm *.Lamprogrammus_exutus_input ; rm *.Carapus_acus_input ; rm *.Brotula_barbata_input ; rm *.Lucifuga_holguinensis_input ; rm *.Lucifuga_dentata_input ; rm *.Astyanax_mexicanus_Surface_input ; rm *.Astyanax_mexicanus_Cave_input ; rm *.Pygocentrus_nattereri_input")


#Run mutpred2


system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Astyanax_CF_Vision_Mutpred2.input -o Astyanax_CF_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Gadus_morhua_Vision_Mutpred2.input -o AGadus_morhua_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Lucifuga_holguinensis_Vision_Mutpred2.input -o Lucifuga_holguinensis_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Astyanax_SF_Vision_Mutpred2.input -o Astyanax_SF_Vision.rslt &")	
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Lamprogrammus_exutus_Vision_Mutpred2.input -o Lamprogrammus_exutus_Vision.rslt &")	
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Neolamprologus_brichardi_Vision_Mutpred2.input -o Neolamprologus_brichardi_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Brotula_barbata_Vision_Mutpred2.input -o Brotula_barbata_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Lamprologus_lethops_Vision_Mutpred2.input -o Lamprologus_lethops_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Percopsis_transmontana_Vision_Mutpred2.input -o Percopsis_transmontana_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Carapus_acus_Vision_Mutpred2.input -o Carapus_acus_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Lamprologus_tigripic_Vision_Mutpred2.input -o Lamprologus_tigripic_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Pygocentrus_nattereri_Vision_Mutpred2.input -o Pygocentrus_nattereri_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Danio_rerio_Vision_Mutpred2.input -o Danio_rerio_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Lucifuga_dentata_Vision_Mutpred2.input -o Lucifuga_dentata_Vision.rslt &")
system("nohup nice -n 10 /home/policarpo/mutpred2.0/run_mutpred2.sh -i Typhlichthys_subterraneus_Vision_Mutpred2.input -o Typhlichthys_subterraneus_Vision.rslt &")





#### Pseudogene number  ---------------------------------

Vision_Sequences %>% filter(Gene != "Not_Found") %>%
  group_by(Species, Gene.Type) %>%
  summarise(n())


#### LoF stats  ---------------------------------


#LoF distribution

labels_lof <- Vision_LoFs %>% 
  filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops") %>%
  group_by(LoF) %>%
  summarise(value = n()) %>%
  unite(n_string, c(value, LoF), sep=" ") %>%
  pull(n_string)


Vision_LoFs %>% group_by(LoF) %>%
  filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops") %>%
  summarise(value = n()) %>%
  ggplot(., aes(x="", y=value, fill=LoF)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)  +
  theme_minimal() +
  theme_void() +
  scale_fill_manual(name=paste("LoF Type", "(",nrow(Vision_LoFs %>% filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops")) ,"mutations )"),
                    values=c("gray", "black", "purple", "orange", "red"),
                    labels=labels_lof) + 
  ggtitle(label = "LoF mutations distribution")


#Frameshift distribution

Vision_LoFs %>% filter(LoF == "Deletion" | LoF == "Insertion") %>%
  group_by(LoF, Nb_Bases) %>%
  summarise(value=n()) %>%
  ggplot(., aes(x=Nb_Bases, y=value, fill=LoF)) +
  geom_bar(stat='identity', position='dodge') + theme_minimal() +
  xlab("Indel size") +
  ylab("Number") +
  ggtitle("Distribution of indel size") +
  scale_fill_discrete(name="Indel type", labels=c("Deletions", "Insertions")) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=18,face="bold"))



#CDS repartition

Vision_LoFs %>%
  filter(LoF != "Splice" & LoF != "Start_loss" & LoF != "Stop_loss") %>%
  mutate(pos_cds_perc = (Pos_CDS/Danio_rerio_Cds_length)*100) %>%
  ggplot(., aes(y=1, x=pos_cds_perc)) +
  geom_point(aes(color=LoF)) +
  theme_void() +
  geom_hline(yintercept = 1, color = "gray") 







#### Estimation LoF Mutation rate in Typhlichthys subterraneus (vision genes) ---------------------------------


######### PROBA STOP


#transition/transversion ratio esitmated on paml on concatenated Vision genes with species of this study
#paml file : Concatenated_Alignment/This_study_species/Vision_Genes_Concatenated_Alignment.codeml


R_ratio = 1.68141 

#To avoid taking sequences full of frameshifts or stop codons, we take Percopsis sequences as reference


Concat_vision_seq <- Percopsis_transmontana_table %>%
  pull(Sequence) %>%
  paste(., collapse = '') %>%
  as.character(.)


n=0
x=0

for(codon_stop in seq(3, nchar(Concat_vision_seq), 3)){
  codon_start <- codon_stop-2
  codon_seq <- substr(Concat_vision_seq, codon_start, codon_stop)
  
  if(codon_seq == "CAA" | codon_seq == "CGA" | codon_seq == "CAG"){
    x <- x + R_ratio/((3*R_ratio)+6)
    n <- n + 3
  } else if(codon_seq == "GAA" | codon_seq == "AAA" | codon_seq == "AGA" | codon_seq == "TGT" | codon_seq == "TGC" | codon_seq == "AAG" | codon_seq == "GAG" | codon_seq == "TTG" | codon_seq == "TCG" | codon_seq == "GGA"){
    x <- x + 1.0/((3*R_ratio)+6)
    n <- n + 3
    
  } else if(codon_seq == "TTA" | codon_seq == "TCA" | codon_seq == "TAC" | codon_seq == "TAT" | codon_seq == "TTA" | codon_seq == "TAC"){
    x <- x + 2.0/((3*R_ratio)+6)
    n <- n + 3
    
  } else if(codon_seq == "TGG"){
    x <= x + (2.0*R_ratio)/((3*R_ratio)+6)
    n <- n + 3
    
  } else {
    
    x <- x + 0.0
    n <- n + 3
    
  }
  
  
  
}

proba_stop_gain <- x/(n/3)




#### PROBA Frameshift


nb_stop <- nrow(Vision_LoFs %>% 
                  filter(Species == "Typhlichthys_subterraneus") %>%
                  filter(LoF == "Stop_gain"))

nb_fs <- nrow(Vision_LoFs %>% 
                filter(Species == "Typhlichthys_subterraneus") %>%
                filter(LoF == "Deletion" | LoF == "Insertion"))

proba_frameshift <- proba_stop_gain * (nb_fs/nb_stop)



###### PROBA Splice

intron_number <- Typhlichthys_subterraneus_table %>%
  mutate(Intron_count = Exon.Count - 1) %>%
  pull(Intron_count) %>%
  sum(.)


base_number <- Typhlichthys_subterraneus_table %>%
  pull(CDS.length) %>%
  sum(.)


proba_splice <- 4 * (intron_number / base_number)


###### PROBA Start and Stop loss

proba_start_loss <- 3 * (nrow(Typhlichthys_subterraneus_table) / base_number)
proba_stop_loss <- 0.852 * proba_start_loss

###### PROBA LoF


LoF_proba <- proba_stop_loss + proba_start_loss + proba_splice + proba_frameshift + proba_stop_gain





#### Simulation Number neutral genes Typhlichthys_subterraneus  ---------------------------------


Total_Number_Genes <- nrow(Vision_Sequences %>% 
                             filter(Species == "Typhlichthys_subterraneus") %>%
                             filter(Gene != "Not_Found")) 

Total_Number_Pseudo <- nrow(Vision_Sequences %>% 
                              filter(Species == "Typhlichthys_subterraneus") %>%
                              filter(Gene != "Not_Found") %>% 
                              filter(Gene.Type == "Pseudogene"))




Number_Lof <- nrow(Vision_LoFs %>% filter(Species == "Typhlichthys_subterraneus"))

Nb_lof_intron <- nrow(Vision_LoFs %>% filter(Species == "Typhlichthys_subterraneus") %>% filter(LoF == "Splice"))
Nb_lof_cds <- Number_Lof - Nb_lof_intron



## Simulated numbers of LoF ##


#c(0.4354839, 0.5, 0.6, 0.7, 0.8, 0.9, 1) * 62


Table_rslt_simu_final <- as.data.frame(NULL)
nb_simu <- 10000
for(NB_neutral_gene in c(62 ,56, 50, 43, 37, 31, 27)){
  
  nb_genes_0_lof <- c()
  nb_genes_1_lof <- c()
  nb_genes_2_lof <- c()
  nb_genes_3_lof <- c()
  nb_genes_4_lof <- c()
  nb_genes_5_lof <- c()
  
  for(i in rep(1, nb_simu)){
    current_table <- sample_n(Typhlichthys_subterraneus_table, NB_neutral_gene) #pick x random genes evolving as neutral
    
    simulation_rslt <- bind_rows(
      sample_n(current_table, size=Nb_lof_cds, weight=CDS.length, replace = TRUE),
      sample_n(current_table, size=Nb_lof_intron, weight=Exon.Count, replace = TRUE))
    
    
    table_Lof_nb <- simulation_rslt %>% group_by(Gene) %>% summarise(LoF_nb = n()) %>%
      group_by(LoF_nb) %>% summarise(Count = n())
    
    
    nb_genes_0_lof <- c(nb_genes_0_lof, Total_Number_Genes - length(unique(simulation_rslt %>% pull(Gene))))
    nb_genes_1_lof <- c(nb_genes_1_lof, table_Lof_nb %>% filter(LoF_nb == 1) %>% pull(Count))
    nb_genes_2_lof <- c(nb_genes_2_lof, table_Lof_nb %>% filter(LoF_nb == 2) %>% pull(Count))
    nb_genes_3_lof <- c(nb_genes_3_lof, table_Lof_nb %>% filter(LoF_nb == 3) %>% pull(Count))
    nb_genes_4_lof <- c(nb_genes_4_lof, table_Lof_nb %>% filter(LoF_nb == 4) %>% pull(Count))
    nb_genes_5_lof <- c(nb_genes_5_lof, sum(table_Lof_nb %>% filter(LoF_nb >= 5) %>% pull(Count)))
    
  }
  
  mean_0_lof <- (sum(nb_genes_0_lof)/nb_simu) * 100 / Total_Number_Genes
  mean_1_lof <- (sum(nb_genes_1_lof)/nb_simu) * 100 / Total_Number_Genes
  mean_2_lof <- (sum(nb_genes_2_lof)/nb_simu) * 100 / Total_Number_Genes
  mean_3_lof <- (sum(nb_genes_3_lof)/nb_simu) * 100 / Total_Number_Genes
  mean_4_lof <- (sum(nb_genes_4_lof)/nb_simu) * 100 / Total_Number_Genes
  mean_5_lof <- (sum(nb_genes_5_lof)/nb_simu) * 100 / Total_Number_Genes
  
  
  LoF_count_vector <- c(mean_0_lof, mean_1_lof, mean_2_lof, mean_3_lof, mean_4_lof, mean_5_lof)
  LoF_nb_vector <- c(0, 1, 2, 3, 4, 5)
  Nb_neutral_gene_vector <- rep(NB_neutral_gene, 6)
  
  Table_rslt_simu_ongoing <- as.data.frame(cbind(Nb_neutral_gene_vector, LoF_nb_vector, LoF_count_vector))
  
  Table_rslt_simu_final <- rbind(Table_rslt_simu_final, Table_rslt_simu_ongoing)
  
}



## Observed numbers of LoF ##

obs_table <- Vision_LoFs %>% filter(Species == "Typhlichthys_subterraneus") %>%
  group_by(Gene) %>% summarise(LoF_nb=n()) %>%
  group_by(LoF_nb) %>%
  summarise(Count = n())

obs_1_lof <- (obs_table %>% filter(LoF_nb == 1) %>% pull(Count)) * 100 / Total_Number_Genes
obs_2_lof <- (obs_table %>% filter(LoF_nb == 2) %>% pull(Count)) * 100 / Total_Number_Genes
obs_3_lof <- (obs_table %>% filter(LoF_nb == 3) %>% pull(Count)) * 100 / Total_Number_Genes
obs_4_lof <- (obs_table %>% filter(LoF_nb == 4) %>% pull(Count)) * 100 / Total_Number_Genes
obs_5_lof <- (sum(obs_table %>% filter(LoF_nb >= 5) %>% pull(Count))) * 100 / Total_Number_Genes
obs_0_lof <- (Total_Number_Genes - Total_Number_Pseudo) * 100 / Total_Number_Genes




## Graph ##

graph_LoF_sitrib <- Table_rslt_simu_final %>% 
  ggplot(., aes(x=LoF_nb_vector, y=LoF_count_vector, group=as.factor(Nb_neutral_gene_vector))) +
  geom_bar(stat="identity", position="dodge", aes(fill= as.factor(Nb_neutral_gene_vector))) +
  theme_minimal() +
  scale_fill_discrete(name = "Neutral gene number", labels = c("62" ,"55", "50", "45", "40", "35", "27")) +
  geom_segment(aes(x = -0.45, y = obs_0_lof, xend = 0.5, yend = obs_0_lof, color = "red")) +
  geom_segment(aes(x = 0.55, y = obs_1_lof, xend = 1.4, yend = obs_1_lof, color = "red")) +
  geom_segment(aes(x = 1.55, y = obs_2_lof, xend = 2.4, yend = obs_2_lof, color = "red")) +
  geom_segment(aes(x = 2.55, y = obs_3_lof, xend = 3.4, yend = obs_3_lof, color = "red")) +
  geom_segment(aes(x = 3.55, y = obs_4_lof, xend = 4.4, yend = obs_4_lof, color = "red")) +
  geom_segment(aes(x = 4.55, y = obs_5_lof, xend = 5.4, yend = obs_5_lof, color = "red")) +
  ggtitle("Simulations with 62 genes in Typhlichthys subterraneus") +
  xlab("Number of LoF") +
  ylab("Percentage of genes") 







#### Simulation Number neutral genes Lucifuga dentata  ---------------------------------

Total_Number_Genes_dent <- nrow(Dentata_Table)

Total_Number_Pseudo_dent <- nrow(Dentata_Table %>% 
                              filter(V6 == "Pseudogene"))




Number_Lof_dent <- 22
Nb_lof_intron_dent <- 6
  
  
Nb_lof_cds_dent <- Number_Lof_dent - Nb_lof_intron_dent



## Simulated numbers of LoF ##



#c(0.4354839, 0.5, 0.6, 0.7, 0.8, 0.9, 1) * 76

Table_rslt_simu_final_dent <- as.data.frame(NULL)
nb_simu <- 10000
for(NB_neutral_gene in c(76 ,68, 61, 53, 46, 38, 33)){
  
  nb_genes_0_lof <- c()
  nb_genes_1_lof <- c()
  nb_genes_2_lof <- c()
  nb_genes_3_lof <- c()
  nb_genes_4_lof <- c()
  nb_genes_5_lof <- c()
  
  for(i in rep(1, nb_simu)){
    current_table_dent <- sample_n(Dentata_Table, NB_neutral_gene) #pick x random genes evolving as neutral
    
    simulation_rslt_dent <- bind_rows(
      sample_n(current_table_dent, size=Nb_lof_cds_dent, weight=V5, replace = TRUE),
      sample_n(current_table_dent, size=Nb_lof_intron_dent, weight=V8, replace = TRUE))
    
    
    table_Lof_nb_dent <- simulation_rslt_dent %>% group_by(V2) %>% summarise(LoF_nb = n()) %>%
      group_by(LoF_nb) %>% summarise(Count = n())
    
    
    nb_genes_0_lof <- c(nb_genes_0_lof, Total_Number_Genes_dent - length(unique(simulation_rslt_dent %>% pull(V2))))
    nb_genes_1_lof <- c(nb_genes_1_lof, table_Lof_nb_dent %>% filter(LoF_nb == 1) %>% pull(Count))
    nb_genes_2_lof <- c(nb_genes_2_lof, table_Lof_nb_dent %>% filter(LoF_nb == 2) %>% pull(Count))
    nb_genes_3_lof <- c(nb_genes_3_lof, table_Lof_nb_dent %>% filter(LoF_nb == 3) %>% pull(Count))
    nb_genes_4_lof <- c(nb_genes_4_lof, table_Lof_nb_dent %>% filter(LoF_nb == 4) %>% pull(Count))
    nb_genes_5_lof <- c(nb_genes_5_lof, sum(table_Lof_nb_dent %>% filter(LoF_nb >= 5) %>% pull(Count)))
    
  }
  
  mean_0_lof <- (sum(nb_genes_0_lof)/nb_simu) * 100 / Total_Number_Genes_dent
  mean_1_lof <- (sum(nb_genes_1_lof)/nb_simu) * 100 / Total_Number_Genes_dent
  mean_2_lof <- (sum(nb_genes_2_lof)/nb_simu) * 100 / Total_Number_Genes_dent
  mean_3_lof <- (sum(nb_genes_3_lof)/nb_simu) * 100 / Total_Number_Genes_dent
  mean_4_lof <- (sum(nb_genes_4_lof)/nb_simu) * 100 / Total_Number_Genes_dent
  mean_5_lof <- (sum(nb_genes_5_lof)/nb_simu) * 100 / Total_Number_Genes_dent
  
  
  LoF_count_vector_dent <- c(mean_0_lof, mean_1_lof, mean_2_lof, mean_3_lof, mean_4_lof, mean_5_lof)
  LoF_nb_vector_dent <- c(0, 1, 2, 3, 4, 5)
  Nb_neutral_gene_vector_dent <- rep(NB_neutral_gene, 6)
  
  Table_rslt_simu_ongoing_dent <- as.data.frame(cbind(Nb_neutral_gene_vector_dent, LoF_nb_vector_dent, LoF_count_vector_dent))
  
  Table_rslt_simu_final_dent <- rbind(Table_rslt_simu_final_dent, Table_rslt_simu_ongoing_dent)
  
}



## Observed numbers of LoF ##


obs_table_dent <- cbind(c(1, 2, 3, 4, 5), c(16, 3, 0, 0, 0))
colnames(obs_table_dent) <- c("LoF_nb", "Count")
obs_table_dent <- as.data.frame(obs_table_dent)

obs_1_lof_dent <- (obs_table_dent %>% filter(LoF_nb == 1) %>% pull(Count)) * 100 / Total_Number_Genes_dent
obs_2_lof_dent <- (obs_table_dent %>% filter(LoF_nb == 2) %>% pull(Count)) * 100 / Total_Number_Genes_dent
obs_3_lof_dent <- (obs_table_dent %>% filter(LoF_nb == 3) %>% pull(Count)) * 100 / Total_Number_Genes_dent
obs_4_lof_dent <- (obs_table_dent %>% filter(LoF_nb == 4) %>% pull(Count)) * 100 / Total_Number_Genes_dent
obs_5_lof_dent <- (sum(obs_table_dent %>% filter(LoF_nb >= 5) %>% pull(Count))) * 100 / Total_Number_Genes_dent
obs_0_lof_dent <- (Total_Number_Genes_dent - Total_Number_Pseudo_dent) * 100 / Total_Number_Genes_dent



# Graph with the two species #




Table_rslt_simu_final_dent <- Table_rslt_simu_final_dent %>% mutate(LoF_count_minus = -LoF_count_vector_dent)




Table_rslt_simu_final_dent %>% 
  ggplot(., aes(x=LoF_nb_vector_dent, y=LoF_count_minus, group=as.factor(Nb_neutral_gene_vector_dent))) +
  geom_bar(stat="identity", position="dodge", aes(fill= as.factor(Nb_neutral_gene_vector_dent)), show.legend = FALSE) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_blank(), panel.background=element_blank()) +
  #theme_classic() +
  geom_bar(data=Table_rslt_simu_final, aes(x=LoF_nb_vector, y=LoF_count_vector, group=as.factor(Nb_neutral_gene_vector), fill= as.factor(Nb_neutral_gene_vector)), stat="identity", position="dodge",  show.legend = FALSE) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#009E73","#0072B2", "#F0E442", "#D55E00", "#0072B2", "#CC79A7", "#D55E00", "#CC79A7")) +
  #coord_flip() +
  scale_x_discrete(name ="Number of LoF mutations per gene", limits=c("0","1","2", "3", "4", "5")) +
  geom_segment(aes(x = -0.45, y = obs_0_lof, xend = 0.5, yend = obs_0_lof, colour = "red")) +
  geom_segment(aes(x = 0.55, y = obs_1_lof, xend = 1.4, yend = obs_1_lof, colour = "red")) +
  geom_segment(aes(x = 1.55, y = obs_2_lof, xend = 2.4, yend = obs_2_lof, colour = "red")) +
  geom_segment(aes(x = 2.55, y = obs_3_lof, xend = 3.4, yend = obs_3_lof, colour = "red")) +
  geom_segment(aes(x = 3.55, y = obs_4_lof, xend = 4.4, yend = obs_4_lof, colour = "red")) +
  geom_segment(aes(x = 4.55, y = obs_5_lof, xend = 5.4, yend = obs_5_lof, colour = "red")) +
  geom_segment(aes(x = -0.45, y = -obs_0_lof_dent, xend = 0.5, yend = -obs_0_lof_dent)) +
  geom_segment(aes(x = 0.55, y = -obs_1_lof_dent, xend = 1.4, yend = -obs_1_lof_dent)) +
  geom_segment(aes(x = 1.55, y = -obs_2_lof_dent, xend = 2.4, yend = -obs_2_lof_dent)) +
  geom_segment(aes(x = 2.55, y = -obs_3_lof_dent, xend = 3.4, yend = -obs_3_lof_dent)) +
  geom_segment(aes(x = 3.55, y = -obs_4_lof_dent, xend = 4.4, yend = -obs_4_lof_dent)) +
  geom_segment(aes(x = 4.55, y = -obs_5_lof_dent, xend = 5.4, yend = -obs_5_lof_dent)) +
  coord_flip()
  




#### Intersection Dentata - Subterraneus  ---------------------------------

# Calcul analytique #

Dentata_Table <- read.table("Dentata_Table.tsv", header=FALSE, sep="\t")

dentata_genes <- Dentata_Table %>% pull(V2)
typhlichthys_genes <- Typhlichthys_subterraneus_table %>% pull(Gene)

nb_genes_common <- length(intersect(dentata_genes, typhlichthys_genes))
nb_genes_dentata <- length(dentata_genes)
nb_genes_typhlichthys <- length(typhlichthys_genes)
nb_pseudo_common <- length(intersect(Dentata_Table %>% filter(V6 == "Pseudogene") %>% pull(V2),
                                     Typhlichthys_subterraneus_table %>% filter(Gene.Type == "Pseudogene") %>% pull(Gene)))

nb_pseudo_dent <- 19
nb_pseudo_typhlichthys <- 27


#Expected number of pseudogenes in common:
(nb_pseudo_typhlichthys/nb_genes_typhlichthys) * (nb_pseudo_dent/nb_genes_dentata) * nb_genes_common




# Simulation #

pseudo_in_common <- intersect(Dentata_Table %>% filter(V6 == "Pseudogene") %>% pull(V2), Typhlichthys_subterraneus_table %>% filter(Gene.Type == "Pseudogene") %>% pull(Gene))
genes_in_common <- intersect(dentata_genes, typhlichthys_genes)





subset_table <- Typhlichthys_subterraneus_table %>% filter(Gene %in% genes_in_common)
nb_subset_pseudo_subt <- nrow(Typhlichthys_subterraneus_table %>% filter(Gene %in% genes_in_common) %>% filter(Gene.Type == "Pseudogene"))
nb_subset_pseudo_dent <- nrow(Dentata_Table %>% filter(V2 %in% genes_in_common) %>% filter(V6 == "Pseudogene"))


rslt_simu_inter <- c()
neutal_vector_sim <- c()
for(neutral_portion in seq(nb_subset_pseudo_subt, nrow(subset_table), 1)){
  for(i in rep(1, 10000)){
    
    subset_table_picked <- sample_n(subset_table, neutral_portion)
    
    nb_intersect <- length(intersect(sample_n(subset_table_picked, nb_subset_pseudo_dent) %>% pull(Gene),
                                     sample_n(subset_table_picked, nb_subset_pseudo_subt) %>% pull(Gene)))
    
    rslt_simu_inter <- c(rslt_simu_inter, nb_intersect)
    neutal_vector_sim <- c(neutal_vector_sim, neutral_portion)
    
  }
}



Simu_df_nolength <- as.data.frame(cbind(rslt_simu_inter, neutal_vector_sim))


test_no_gene_length <- Simu_df_nolength %>% group_by(neutal_vector_sim) %>%
  summarise(n_10_obs = (sum(rslt_simu_inter == nb_pseudo_common)),
            n_simu = n(),
            p_value = n_10_obs / n_simu) 

test_no_gene_length %>% 
  ggplot(., aes(x=neutal_vector_sim, y=p_value)) +
  geom_line() +
  ggtitle("Maximum likelyhood of observing 10 common pseudogenes depending on the number of neutral genes among 58 genes") +
  xlab("Number of neutral genes") +
  ylab("Probability") +
  theme_minimal() +
  geom_hline(yintercept = 0.05, color="red") +
  annotate("text", x=20, y=0.07, label="0.05", color="red")




#### SAME SIMULATION BUT TAKING INTO ACCOUNT GENE LENGTH 


#dent_subset <- Dentata_Table %>% filter(V2 %in% genes_in_common) %>%
#  dplyr::select(V2, V5)
#colnames(dent_subset) <- c("Gene", "Dent_length")
#
#subset_table_V2 <- merge(subset_table, dent_subset, by="Gene") %>%
#  mutate(Total_length = CDS.length + Dent_length)


#Compute mean length of pseudos and non-pseudo taking Danio rerio as ref


pseudo_in_common <- intersect(Dentata_Table %>% filter(V6 == "Pseudogene") %>% pull(V2), Typhlichthys_subterraneus_table %>% filter(Gene.Type == "Pseudogene") %>% pull(Gene))
genes_in_common <- intersect(dentata_genes, typhlichthys_genes)


#some genes dont have the same annotation, so change them individually
genes_in_common[6] <- "opn6a"
genes_in_common[15] <- "rh1.1"
genes_in_common[16] <- "rh2.1"

Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
  filter(Gene %in% pseudo_in_common) %>% pull(CDS.length) %>% mean()
  

Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
  filter(Gene %in% genes_in_common) %>% pull(CDS.length) %>% mean()



subset_table_V2 <- Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
  filter(Gene %in% genes_in_common)
  
  
  

#simu

rslt_simu_inter <- c()
neutal_vector_sim <- c()
for(neutral_portion in seq(nb_subset_pseudo_subt, nrow(subset_table_V2), 1)){
  for(i in rep(1, 10000)){
    
    subset_table_picked <- sample_n(subset_table_V2, neutral_portion)
    
    nb_intersect <- length(intersect(sample_n(subset_table_picked, nb_subset_pseudo_dent, weight=CDS.length) %>% pull(Gene),
                                     sample_n(subset_table_picked, nb_subset_pseudo_subt, weight=CDS.length) %>% pull(Gene)))
    
    rslt_simu_inter <- c(rslt_simu_inter, nb_intersect)
    neutal_vector_sim <- c(neutal_vector_sim, neutral_portion)
    
  }
}



Simu_df <- as.data.frame(cbind(rslt_simu_inter, neutal_vector_sim))


test <- Simu_df %>% group_by(neutal_vector_sim) %>%
  summarise(n_10_obs = (sum(rslt_simu_inter == nb_pseudo_common)),
            n_simu = n(),
            p_value = n_10_obs / n_simu) 

graph_inter <- test %>% 
  ggplot(., aes(x=neutal_vector_sim, y=p_value)) +
  geom_line() +
  ggtitle("Maximum likelyhood of observing 10 common pseudogenes depending on the number of neutral genes among 58 genes") +
  xlab("Number of neutral genes") +
  ylab("Probability") +
  theme_minimal() +
  geom_hline(yintercept = 0.05, color="red") +
  annotate("text", y=0.06, x=25, label="0.05", color="red")


(graph_LoF_sitrib / graph_inter)


#+
#geom_hline(yintercept = 0.05, color="red") +
#annotate("text", x=20, y=0.07, label="0.05", color="red")

#compute 95 perc confience ineterval




test %>% arrange(desc(n_10_obs))
test %>% arrange(neutal_vector_sim)


# Probability of finind n common pseudogenes



Simu_df %>% group_by(neutal_vector_sim) %>%
  summarise(moyenne = mean(rslt_simu_inter),
            sd = sd(rslt_simu_inter),se = sd/sqrt(n),
            n=n()) %>%
  ggplot(., aes(x=neutal_vector_sim, y=moyenne)) +
  geom_point() +
  geom_errorbar(aes(ymin = moyenne - se,
                    ymax = moyenne + se)) +
  theme_minimal() +
  ylab("Mean number of common pseudogenes found") +
  xlab("Number of neutral genes") +
  geom_hline(yintercept = 10, color="red", linetype="dashed")
  
  



####

Simu_df %>% filter(neutal_vector_sim == 58) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  mutate(nb_obs_prop = nb_obs/10000) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs_prop)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  ylab("Number of observations") +
  geom_hline(yintercept = 0.05, col="red", linetype="dashed")#+
  #ggtitle("Number of common pseudogenes in 10,000 simulations with 58 neutral evolving genes")
  


Simu_df_nolength %>% filter(neutal_vector_sim == 58) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  mutate(nb_obs_prop = nb_obs/10000) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs_prop)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  ylab("Number of observations") +
  geom_hline(yintercept = 0.05, col="red", linetype="dashed") #+
#ggtitle("Number of common pseudogenes in 10,000 simulations with 58 neutral evolving genes")




#### SAME SIMULATION BUT TAKING INTO ACCOUNT GENE LENGTH + NUMBER OF INTRON



subset_table_V2 <- subset_table_V2 %>%
  mutate(Intron_nb = Exon.Count - 1)

#Ponderate cds length and intron length by the LoF proba rate in exon or in intron 


subset_table_V2 <- subset_table_V2 %>%
  mutate(exon_weight = (proba_stop_gain + proba_frameshift + proba_stop_loss + proba_start_loss) * CDS.length) %>%
  mutate(intron_weight = proba_splice * Intron_nb) %>%
  mutate(Proba_mut = intron_weight+exon_weight)



rslt_simu_inter <- c()
neutal_vector_sim <- c()
for(neutral_portion in seq(nb_subset_pseudo_subt, nrow(subset_table_V2), 1)){
  for(i in rep(1, 10000)){
    
    subset_table_picked <- sample_n(subset_table_V2, neutral_portion)
    
    nb_intersect <- length(intersect(sample_n(subset_table_picked, nb_subset_pseudo_dent, weight=Proba_mut) %>% pull(Gene),
                                     sample_n(subset_table_picked, nb_subset_pseudo_subt, weight=Proba_mut) %>% pull(Gene)))
    
    rslt_simu_inter <- c(rslt_simu_inter, nb_intersect)
    neutal_vector_sim <- c(neutal_vector_sim, neutral_portion)
    
  }
}

Simu_df_exon_intron <- as.data.frame(cbind(rslt_simu_inter, neutal_vector_sim))


test_exon_intron <- Simu_df_exon_intron %>% group_by(neutal_vector_sim) %>%
  summarise(n_10_obs = (sum(rslt_simu_inter == nb_pseudo_common)),
            n_simu = n(),
            p_value = n_10_obs / n_simu) 





Simu_df_exon_intron %>% filter(neutal_vector_sim == 58) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)) +
  ylab("Number of observations") 


Simu_df_exon_intron %>% filter(neutal_vector_sim == 58) %>%
  group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  mutate(nb_obs_prop = nb_obs/10000) %>%
  ggplot(., aes(x=rslt_simu_inter, y=nb_obs_prop)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  ylab("Number of observations") +
  geom_hline(yintercept = 0.05, col="red", linetype="dashed")




Simu_df_nolength %>% filter(neutal_vector_sim == 58) %>% pull(rslt_simu_inter) %>% mean(.)
Simu_df %>% filter(neutal_vector_sim == 58) %>% pull(rslt_simu_inter) %>% mean(.)
Simu_df_exon_intron %>% filter(neutal_vector_sim == 58) %>% pull(rslt_simu_inter) %>% mean(.)

#+
#ggtitle("Number of common pseudogenes in 10,000 simulations with 58 neutral evolving genes")


## Three on the same graph ##

d1_inter <- as.data.frame(Simu_df_nolength %>% filter(neutal_vector_sim == 58) %>% group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  mutate(nb_obs_prop = nb_obs/10000) %>% mutate(Simu_name = "Same_proba"))

d2_inter <- as.data.frame(Simu_df %>% filter(neutal_vector_sim == 58) %>% group_by(rslt_simu_inter) %>%
  summarise(nb_obs = n()) %>%
  mutate(nb_obs_prop = nb_obs/10000) %>% mutate(Simu_name = "Gene_length"))

d3_inter <- as.data.frame(Simu_df_exon_intron %>% filter(neutal_vector_sim == 58) %>% group_by(rslt_simu_inter) %>%
                            summarise(nb_obs = n()) %>%
                            mutate(nb_obs_prop = nb_obs/10000) %>% mutate(Simu_name = "Gene_length_plus_intron_nb"))


d1_2_3 <- bind_rows(d1_inter, d2_inter, d3_inter)

d1_2_3 %>% ggplot(., aes(x=rslt_simu_inter, y=nb_obs_prop, fill=Simu_name)) +
  geom_bar(stat="identity", position="dodge", show.legend=FALSE) +
  theme_minimal() +
  scale_x_continuous("Number of common pseudogenes", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  ylab("Proportion of observations") +
  geom_hline(yintercept = 0.05, col="red", linetype="dashed")
  






##


#### Simulation for datation Typhlichthys subterraneus ---------------------------------



setwd("~/Desktop/These_Supplemental_Analysis/Datations/Simulations_Subterraneus/")

#Simulations taking account gene lengths
#system("/Users/maxime/miniconda3/bin/python Script_simu_subteraneus.py -g 1000000 -r 10000 -m 1e-8 -s 0.1278851 -N 500.0")


#Analytical models assuming small populations sizes and equal gene lengths

#system("python2 Simulation_Analytique_Subterraneus.py")
#system("/Users/maxime/miniconda3/bin/python Analatycal_gen_subteraneus.py")


#Analytical models assuming small populations sizes and only genes found as pseudogene in 1 species

Pseudogene_list <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn4m3", "opn6b", "opn7a", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "tmt3b", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1", "gnb3b", "gc2", "gc3", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")
Pseudogene_list_wo_sub <- c("opn4m1", "opn4m2", "opn4m3", "opn6b", "opn7a", "opn7b", "opn7c", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt2a", "tmt3a", "tmt3b","cryaa", "crybb1", "crybgx", "crygm5", "crygn2", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh1", "gnb3b", "gc2", "gc3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")



nb_gene_neutre_subt <- nrow(Typhlichthys_subterraneus_table %>% filter(Gene %in% Pseudogene_list_wo_sub))
moy_length_neutre_subt <- Typhlichthys_subterraneus_table %>% filter(Gene %in% Pseudogene_list_wo_sub) %>% pull(CDS.length) %>% mean(.)
Nb_pseudo_subset_subt <- nrow(Typhlichthys_subterraneus_table %>% filter(Gene %in% Pseudogene_list_wo_sub) %>% filter(Gene.Type == "Pseudogene"))

system("python2 Simulation_Analytique_Subterraneus_subset.py")
system("/Users/maxime/miniconda3/bin/python Analatycal_gen_subteraneus_subset.py")




#Graphs

Binomial_all_neutral_subt <- read.table(file = "Matrice_loi_binomiale.csv", sep=",", header=FALSE)
Binomial_subset_neutral_subt <- read.table(file = "Matrice_loi_binomiale_subset.csv", sep=",", header=FALSE)

Simu_all_neutral_subt <- read.table(file = "out_simu_stop_1000000_generations_10000_repeats_1e-08_mu_500.0_ind_0_migra_0_migrants_0_mig_0_h_0_s.csv", sep=",", header=FALSE)


#extract column of 27 pseudogenes (28 as first column is 0)
Simu_all_neutral_subt <- Simu_all_neutral_subt %>% 
  dplyr::select(V28) %>%
  mutate(V28 = V28/10000) %>%
  mutate(generation = 1:n())%>%
  mutate(Model = "Simulation - 62 neutral genes")

Binomial_all_neutral_subt <-  Binomial_all_neutral_subt %>%
  dplyr::select(V28) %>%
  mutate(generation = 1:n()) %>%
  mutate(Model = "Binomial")

Merged_Models <- bind_rows(Simu_all_neutral_subt, Binomial_all_neutral_subt)



Simulation_graph_subt <- Merged_Models %>% ggplot(., aes(x=generation, y=(V28), color=Model)) +
  geom_line() +
  theme_minimal() +
  xlab(label = "Number of generations (x10)") +
  ylab(label = "Probability") +
  ggtitle(label = "Probability to find 27 pseudogenes among 62 genes") +
  geom_hline(yintercept = 0.05, color="red") +
  annotate("text", x=52000, y=0.103, label="569,566", col="red") +
  annotate("text", x=63000, y=0.115, label="614,470", col="blue") +
  theme(legend.position="bottom", legend.box = "horizontal") 



Simu_all_neutral_subt %>% arrange(desc(V28))
Binomial_all_neutral_subt %>% arrange(desc(V28))

Binomial_all_neutral_subt %>% filter(V28 > 0.05) %>%
  arrange(desc(generation))



#### Simulation for datation Lamprologus lethops  ---------------------------------

setwd("~/Desktop/These_Supplemental_Analysis/Datations/Simulations_lethops/")

#Simulations taking account gene lengths
#system("/Users/maxime/miniconda3/bin/python Script_simu_lethops.py  -g 1000000 -r 10000 -m 1e-8 -s 0.1278851 -N 500.0")


#Analytical models assuming small populations sizes and equal gene lengths
#system("/Users/maxime/miniconda3/bin/python Analatycal_gen_lethops.py")
#system("python2 Simulation_Analytique_Lethops.py")


#Analytical models assuming small populations sizes and only genes found as pseudogene in 1 species


nb_gene_neutre_lethops <- nrow(Lamprologus_lethops_table %>% filter(Gene %in% Pseudogene_list))
moy_length_neutre_lethops <- Lamprologus_lethops_table %>% filter(Gene %in% Pseudogene_list) %>% pull(CDS.length) %>% mean(.)

#system("/Users/maxime/miniconda3/bin/python Analatycal_gen_lethops_subset.py")
#system("python2 Simulation_Analytique_Lethops_subset.py")


#Graphs

Binomial_all_neutral_lethops <- read.table(file = "Matrice_loi_binomiale.csv", sep=",", header=FALSE)
Binomial_subset_neutral_lethops <- read.table(file = "Matrice_loi_binomiale_subset.csv", sep=",", header=FALSE)
Simu_all_neutral_lethops <- read.table(file = "out_simu_stop_1000000_generations_10000_repeats_1e-08_mu_500.0_ind_0_migra_0_migrants_0_mig_0_h_0_s.csv", sep=",", header=FALSE)


#extract column of 2 pseudogenes (3 as first column is 0)
Simu_all_neutral_lethops <- Simu_all_neutral_lethops %>% 
  dplyr::select(V3) %>%
  mutate(V3 = V3/10000) %>%
  mutate(generation = 1:n())%>%
  mutate(Model = "Simulation - 79 neutral genes")


##
Ngeneration <- rownames(Binomial_all_neutral_lethops)

Binomial_all_neutral_lethops <-  Binomial_all_neutral_lethops %>%
  dplyr::select(V3) 

Binomial_all_neutral_lethops <- cbind(generation=Ngeneration, Binomial_all_neutral_lethops) 
Binomial_all_neutral_lethops <- Binomial_all_neutral_lethops %>% mutate(Model = "Binomial")


#Binomial_all_neutral_lethops <-  Binomial_all_neutral_lethops %>%
#  dplyr::select(V3) %>%
#  mutate(generation = 1:n()) %>%
  mutate(Model = "Binomial")

Merged_Models <- bind_rows(Simu_all_neutral_lethops, Binomial_all_neutral_lethops)




Simulation_graph_lethops <- Merged_Models %>% ggplot(., aes(x=generation, y=(V3), color=Model)) +
  geom_line() +
  theme_minimal() +
  xlab(label = "Number of generations (x10)") +
  xlim(0, 20000) +
  ylab(label = "Probability") +
  ggtitle(label = "Probability to find 2 pseudogenes among 79 genes") +
  geom_hline(yintercept = 0.05, color="red") +
  annotate("text", x=1630, y=0.28, label="19,187", col="red") +
  annotate("text", x=2650, y=0.285, label="21,020", col="blue") +
  theme(legend.position="bottom", legend.box = "horizontal") 



Simu_all_neutral_lethops %>% arrange(desc(V3))
Binomial_all_neutral_lethops %>% arrange(desc(V3))


head(Binomial_all_neutral_lethops %>% filter(V3 > 0.05) %>%
       arrange(desc(generation)))

head(Simu_all_neutral_lethops %>% filter(V3 > 0.05) %>%
       arrange(generation))




(Simulation_graph_subt | Simulation_graph_lethops)





## Other graph ##

Binomial_all_neutral_lethops <-  Binomial_all_neutral_lethops %>%
  dplyr::select(V3) %>%
  mutate(generation = 1:n()) %>%
  mutate(Model = "Binomial_all")

Binomial_all_neutral_subt <-  Binomial_all_neutral_subt %>%
  dplyr::select(V28) %>%
  mutate(generation = 1:n()) %>%
  mutate(Model = "Binomial_all")


Binomial_subset_neutral_lethops <-  Binomial_subset_neutral_lethops %>%
  dplyr::select(V3) %>%
  mutate(generation = 1:n()) %>%
  mutate(Model = "Binomial_subset")

Binomial_subset_neutral_subt <-  Binomial_subset_neutral_subt %>%
  dplyr::select(V14) %>%
  mutate(generation = 1:n()) %>%
  mutate(Model = "Binomial_subset")



Merged_Models_lethops <- bind_rows(Binomial_subset_neutral_lethops, Binomial_all_neutral_lethops)
Merged_Models_subt <- bind_rows(Binomial_subset_neutral_subt, Binomial_all_neutral_subt)

 



Merged_Models_lethops %>% ggplot(., aes(x=generation, y=(V3), color=Model)) +
  geom_line() +
  theme_minimal() +
  xlab(label = "Number of generations (x10)") +
  xlim(0, 50000) +
  ylab(label = "Probability") +
  ggtitle(label = "Probability to find 2 pseudogenes among 79 genes") +
  geom_hline(yintercept = 0.05, color="red") +
  #annotate("text", x=1630, y=0.28, label="19,187", col="red") +
  #annotate("text", x=2650, y=0.285, label="21,020", col="blue") +
  theme(legend.position="bottom", legend.box = "horizontal") 

Merged_Models_subt %>% ggplot(., aes(x=generation, y=(V28), color=Model)) +
  geom_line() +
  theme_minimal() +
  xlab(label = "Number of generations (x10)") +
  #xlim(0, 100000) +
  ylab(label = "Probability") +
  ggtitle(label = "Probability to find 27 pseudogenes among 62 genes") +
  geom_hline(yintercept = 0.05, color="red") +
  #annotate("text", x=1630, y=0.28, label="19,187", col="red") +
  #annotate("text", x=2650, y=0.285, label="21,020", col="blue") +
  theme(legend.position="bottom", legend.box = "horizontal") 



colnames(Merged_Models_lethops) <- c("Probability", "generation", "Model")
colnames(Merged_Models_subt) <- c("Probability", "generation", "Model")


Merged_Models_lethops <- Merged_Models_lethops %>% mutate(Species = "Lamprologus_lethops")
Merged_Models_subt <- Merged_Models_subt %>% mutate(Species = "Typhlichthys_subterraneus")


Merged_Models_lethops <- Merged_Models_lethops %>% mutate(Model_species = paste(Model, Species, sep=""))
Merged_Models_subt <- Merged_Models_subt %>% mutate(Model_species = paste(Model, Species, sep=""))


Merged_Models_subt %>% pull(Model) %>% unique(.)



Merged_Models_all <- bind_rows(Merged_Models_lethops, Merged_Models_subt)



Merged_Models_all %>% ggplot(., aes(x=generation, y=(Probability), color=Model_species)) +
  geom_line() +
  theme_minimal() +
  xlab(label = "Number of generations (x10)") +
  xlim(0, 280000) +
  ylim(0.000001, 0.3) +
  ylab(label = "Probability") +
  #ggtitle(label = "Probability to find 27 pseudogenes among 62 genes") +
  geom_hline(yintercept = 0.05, color="red") +
  #annotate("text", x=1630, y=0.28, label="19,187", col="red") +
  #annotate("text", x=2650, y=0.285, label="21,020", col="blue") +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_color_manual(values = c("cornflowerblue", "red", "blue", "brown4"))




###

Binomial_subset_neutral_lethops <- Binomial_subset_neutral_lethops %>% mutate(Species = "Lamprologus lethops")
Binomial_all_neutral_lethops <- Binomial_all_neutral_lethops %>% mutate(Species = "Lamprologus lethops")
Binomial_subset_neutral_subt <- Binomial_subset_neutral_subt %>% mutate(Species = "Typhlichthys subterraneus")
Binomial_all_neutral_subt <- Binomial_all_neutral_subt %>% mutate(Species = "Typhlichthys subterraneus")

colnames(Binomial_subset_neutral_lethops) <- c("Probability", "generation", "Model", "Species")
colnames(Binomial_all_neutral_lethops) <- c("Probability", "generation", "Model", "Species")
colnames(Binomial_subset_neutral_subt) <- c("Probability", "generation", "Model", "Species")
colnames(Binomial_all_neutral_subt) <- c("Probability", "generation", "Model", "Species")


Merged_Models_subset <- bind_rows(Binomial_subset_neutral_lethops, Binomial_subset_neutral_subt)
Merged_Models_all <- bind_rows(Binomial_all_neutral_lethops, Binomial_all_neutral_subt)



Binomial_all_neutral_subt %>% ggplot(., aes(x=generation, y=(Probability))) +
  geom_line(linetype = "solid", colour="red")  +
  geom_line(data=Binomial_subset_neutral_subt, aes(x=generation, y=(Probability)), colour="red",  linetype = "longdash") +
  geom_line(data=Binomial_all_neutral_lethops, aes(x=generation, y=(Probability)), linetype = "solid", colour="blue") +
  geom_line(data=Binomial_subset_neutral_lethops, aes(x=generation, y=(Probability)), colour="blue",  linetype = "longdash") +
  xlim(0, 2e+5) +
  ylim(0.001, 0.3) +
  theme_tufte()
  



Binomial_all_neutral_subt %>% arrange(desc(Probability))
Binomial_subset_neutral_subt %>% arrange(desc(Probability))
Binomial_all_neutral_lethops %>% arrange(desc(Probability))
Binomial_subset_neutral_lethops %>% arrange(desc(Probability))



head(Binomial_subset_neutral_subt %>% filter(Probability >= 0.05) %>%
       arrange(desc(generation)))

head(Binomial_subset_neutral_subt %>% filter(Probability >= 0.05) %>%
       arrange(generation))






##

#### Graph mars 2021 for datation curves ---------------------------------

## dentata


setwd("~/Desktop/These_Supplemental_Analysis/Datations/Old_data_Policarpo_2020")

Dentata_table <- read.table("Dentata_table.tsv", header=FALSE, sep="\t")
colnames(Dentata_table) <- c("Species", "Gene", "genomic_pos", "strand", "cds_length", "Gene_type", "Lof", "Exon_nb", "splice", "Sequence")

Astyanax_table <- read.table("Astyanax_Table.tsv", header=FALSE, sep="\t")
colnames(Astyanax_table) <- c("Species", "Gene", "genomic_pos", "strand", "cds_length", "Gene_type", "Lof", "Exon_nb", "splice", "Sequence", "blank")


Pseudogene_list <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn4m3", "opn6b", "opn7a", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "tmt3b", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1", "gnb3b", "gc2", "gc3", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")
Pseudogene_list_wo_dent <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn6b", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1", "gc2", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")

nb_gene_neutre_dent <- nrow(Dentata_table %>% filter(Gene %in% Pseudogene_list_wo_dent))
moy_length_neutre_dent <- Dentata_table %>% filter(Gene %in% Pseudogene_list_wo_dent) %>% pull(cds_length) %>% mean(.)
Nb_pseudo_subset_dent <- nrow(Dentata_table %>% filter(Gene %in% Pseudogene_list_wo_dent) %>% filter(Gene_type == "Pseudogene"))


setwd("~/Desktop/These_Supplemental_Analysis/Datations/Simulations_old_Lucifuga_Astyanax")
system("python2 Analatycal_gen_dentata.py")
system("/Users/maxime/miniconda3/bin/python Analatycal_gen_dentata_subset.py")

Binomial_all_neutral_dent <- read.table(file = "Matrice_loi_binomiale_dent.csv", sep=",", header=FALSE)
Binomial_subset_neutral_dent <- read.table(file = "Matrice_loi_binomiale_subset_dent.csv", sep=",", header=FALSE)




## astyanax

Pseudogene_list <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn4m3", "opn6b", "opn7a", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "tmt3b", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1.1","rh1.2", "gnb3b", "gc2", "gc3", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")
Pseudogene_list_wo_astyanax <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn4m3", "opn6b", "opn7a", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "tmt3b", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1.1","rh1.2", "gnb3b", "gc2", "gc3", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")
nb_gene_neutre_astyanax <- nrow(Astyanax_table %>% filter(Gene %in% Pseudogene_list_wo_astyanax))
moy_length_neutre_astyanax <- Astyanax_table %>% filter(Gene %in% Pseudogene_list_wo_astyanax) %>% pull(cds_length) %>% mean(.)
Nb_pseudo_subset_astyanax <- nrow(Astyanax_table %>% filter(Gene %in% Pseudogene_list_wo_dent) %>% filter(Gene_type == "Pseudogene"))

setwd("~/Desktop/These_Supplemental_Analysis/Datations/Simulations_old_Lucifuga_Astyanax")
system("python2 Analatycal_gen_astyanax.py")
system("/Users/maxime/miniconda3/bin/python Analatycal_gen_astyanax_subset.py")

Binomial_all_neutral_asty <- read.table(file = "Matrice_loi_binomiale_asty.csv", sep=",", header=FALSE)
Binomial_subset_neutral_asty <- read.table(file = "Matrice_loi_binomiale_subset_asty.csv", sep=",", header=FALSE)






setwd("~/Desktop/These_Supplemental_Analysis/Datations/Simulations_Subterraneus/")

Binomial_all_neutral_subt <- read.table(file = "Matrice_loi_binomiale.csv", sep=",", header=FALSE)
Binomial_subset_neutral_subt <- read.table(file = "Matrice_loi_binomiale_subset.csv", sep=",", header=FALSE)


setwd("~/Desktop/These_Supplemental_Analysis/Datations/Simulations_lethops/")

Binomial_all_neutral_lethops <- read.table(file = "Matrice_loi_binomiale.csv", sep=",", header=FALSE)
Binomial_subset_neutral_lethops <- read.table(file = "Matrice_loi_binomiale_subset.csv", sep=",", header=FALSE)



## Extract interesting columns


Binomial_all_neutral_lethops <-  Binomial_all_neutral_lethops %>%
  dplyr::select(V3) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(Binomial_all_neutral_lethops))
Binomial_all_neutral_lethops <- cbind(generation=Ngeneration, Binomial_all_neutral_lethops) 

Binomial_all_neutral_subt <-  Binomial_all_neutral_subt %>%
  dplyr::select(V28) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(Binomial_all_neutral_subt))
Binomial_all_neutral_subt <- cbind(generation=Ngeneration, Binomial_all_neutral_subt) 

Binomial_all_neutral_dent <-  Binomial_all_neutral_dent %>%
  dplyr::select(V20) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(Binomial_all_neutral_dent))
Binomial_all_neutral_dent <- cbind(generation=Ngeneration, Binomial_all_neutral_dent) 

Binomial_all_neutral_asty <-  Binomial_all_neutral_asty %>%
  dplyr::select(V2) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(Binomial_all_neutral_asty))
Binomial_all_neutral_asty <- cbind(generation=Ngeneration, Binomial_all_neutral_asty) 


Binomial_subset_neutral_lethops <-  Binomial_subset_neutral_lethops %>%
  dplyr::select(V3) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Binomial_subset_neutral_lethops))
Binomial_subset_neutral_lethops <- cbind(generation=Ngeneration, Binomial_subset_neutral_lethops) 

Binomial_subset_neutral_subt <-  Binomial_subset_neutral_subt %>%
  dplyr::select(V14) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Binomial_subset_neutral_subt))
Binomial_subset_neutral_subt <- cbind(generation=Ngeneration, Binomial_subset_neutral_subt) 


Binomial_subset_neutral_dent <-  Binomial_subset_neutral_dent %>%
  dplyr::select(V14) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Binomial_subset_neutral_dent))
Binomial_subset_neutral_dent <- cbind(generation=Ngeneration, Binomial_subset_neutral_dent) 


Binomial_subset_neutral_asty <-  Binomial_subset_neutral_asty %>%
  dplyr::select(V2) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Binomial_subset_neutral_asty))
Binomial_subset_neutral_asty <- cbind(generation=Ngeneration, Binomial_subset_neutral_asty) 



## Graph creation

Binomial_subset_neutral_lethops <- Binomial_subset_neutral_lethops %>% mutate(Species = "Lamprologus lethops")
Binomial_all_neutral_lethops <- Binomial_all_neutral_lethops %>% mutate(Species = "Lamprologus lethops")
Binomial_subset_neutral_subt <- Binomial_subset_neutral_subt %>% mutate(Species = "Typhlichthys subterraneus")
Binomial_all_neutral_subt <- Binomial_all_neutral_subt %>% mutate(Species = "Typhlichthys subterraneus")
Binomial_subset_neutral_dent <- Binomial_subset_neutral_dent %>% mutate(Species = "Lucifuga dentata")
Binomial_all_neutral_dent <- Binomial_all_neutral_dent %>% mutate(Species = "Lucifuga dentata")
Binomial_subset_neutral_asty <- Binomial_subset_neutral_asty %>% mutate(Species = "Astyanax mexicanus")
Binomial_all_neutral_asty <- Binomial_all_neutral_asty %>% mutate(Species = "Astyanax mexicanus")




colnames(Binomial_subset_neutral_lethops) <- c("generation", "Probability", "Model", "Species")
colnames(Binomial_all_neutral_lethops) <- c("generation", "Probability", "Model", "Species")
colnames(Binomial_subset_neutral_subt) <- c("generation", "Probability", "Model", "Species")
colnames(Binomial_all_neutral_subt) <- c("generation", "Probability", "Model", "Species")
colnames(Binomial_subset_neutral_dent) <- c("generation", "Probability", "Model", "Species")
colnames(Binomial_all_neutral_dent) <- c("generation", "Probability", "Model", "Species")
colnames(Binomial_subset_neutral_asty) <- c("generation", "Probability", "Model", "Species")
colnames(Binomial_all_neutral_asty) <- c("generation", "Probability", "Model", "Species")



#library(ggthemes)


Binomial_all_neutral_subt %>% ggplot(., aes(x=generation, y=(Probability))) +
  geom_line(linetype = "solid", colour="red")  +
  geom_line(data=Binomial_subset_neutral_subt, aes(x=generation, y=(Probability)), colour="red",  linetype = "longdash") +
  geom_line(data=Binomial_all_neutral_lethops, aes(x=generation, y=(Probability)), linetype = "solid", colour="blue") +
  geom_line(data=Binomial_subset_neutral_lethops, aes(x=generation, y=(Probability)), colour="blue",  linetype = "longdash") +
  xlim(0, 2e+5) +
  ylim(0.001, 0.3) +
  theme_tufte()


## astyanax and dent added

Binomial_all_neutral_subt %>% ggplot(., aes(x=generation, y=(Probability))) +
  geom_line(linetype = "solid", colour="red")  +
  #geom_line(data=Binomial_subset_neutral_subt, aes(x=generation, y=(Probability)), colour="red",  linetype = "longdash") +
  geom_line(data=Binomial_all_neutral_lethops, aes(x=generation, y=(Probability)), linetype = "solid", colour="blue") +
  #geom_line(data=Binomial_subset_neutral_lethops, aes(x=generation, y=(Probability)), colour="blue",  linetype = "longdash") +
  geom_line(data=Binomial_all_neutral_dent, aes(x=generation, y=(Probability)), linetype = "solid", colour="black") +
  #geom_line(data=Binomial_subset_neutral_dent, aes(x=generation, y=(Probability)), colour="black",  linetype = "longdash") +
  geom_line(data=Binomial_all_neutral_asty, aes(x=generation, y=(Probability)), linetype = "solid", colour="gray") +
  #geom_line(data=Binomial_subset_neutral_asty, aes(x=generation, y=(Probability)), colour="gray",  linetype = "longdash") +
  #xlim(0, 2e+5) +
  xlim(0, 1e+5) +
  ylim(0.001, 0.4) +
  theme_tufte()









Binomial_all_neutral_subt %>% arrange(desc(Probability))
Binomial_subset_neutral_subt %>% arrange(desc(Probability))
Binomial_all_neutral_lethops %>% arrange(desc(Probability))
Binomial_subset_neutral_lethops %>% arrange(desc(Probability))
Binomial_all_neutral_dent %>% arrange(desc(Probability))
Binomial_subset_neutral_dent %>% arrange(desc(Probability))



head(Binomial_subset_neutral_subt %>% filter(Probability >= 0.05) %>%
       arrange(desc(generation)))

head(Binomial_subset_neutral_subt %>% filter(Probability >= 0.05) %>%
       arrange(generation))





#### dN/dS on vision genes concatenations that are conserved between species for datation  --------------------------------


setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Vision")

#Align genes that are found in L. lethops, L. tigripictilis and N. brichardi


as.vector(Vision_Sequences %>%
            filter(Species == "Lamprologus_tigripic" | 
                     Species == "Lamprologus_lethops" | 
                     Species == "Neolamprologus_brichardi") %>%
            group_by(Gene) %>%
            summarise(value = n()) %>%
            filter(value == 3) %>%
            pull(Gene) %>%
            unique(.))

system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i cryba1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cryba1l1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cryba2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cryba4_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybb1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybb1l1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybb1l3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybgx_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crygm5_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crygn2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gc2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gc3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gcap1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gcap2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gcap3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gcap4_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnat1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnat2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnb1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnb1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnb3a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gngt1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gngt2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk7a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk7b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gucy2f_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6c_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6ga_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6gb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6ha_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6hb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rcv1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rcv1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rcv2a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rcv2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rgr1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rgr2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rh1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rh2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rpe65a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rpe65b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta saga_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sagb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Cryaa_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Exorh_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Lws1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn4m1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn4m2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn4m3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn4x1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn4x2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn5_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn6a_1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn6a_2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn7a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn7b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn7d_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn8b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn8c_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Parapinopsin_2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Parietopsin_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Rrh_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Sws1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Sws2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Tmt1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Tmt1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Tmt2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Tmt3a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Tmt3b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Va2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta arr3a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Concatenation_lamprologus_danio_vision.fa")
system("mv partitions.txt Concatenation_lamprologus_danio_vision.partition")

#Align genes that are found in Gadus mrohua, Percopsis transmontana and T. subterraneus

as.vector(Vision_Sequences %>%
            filter(Species == "Gadus_morhua" | 
                     Species == "Percopsis_transmontana" |
                     Species == "Typhlichthys_subterraneus") %>%
            group_by(Gene) %>%
            summarise(value = n()) %>%
            filter(value == 3) %>%
            pull(Gene) %>%
            unique(.))


system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i Va2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Tmt3a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Tmt1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rh1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rgr2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rgr1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rcv2a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rcv1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6hb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6gb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pde6a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta parapinopsin_1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Parietopsin_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn8b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn7d_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn7b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn6a_2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn6a_1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn5_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn4x1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn4m3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Opn3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta arr3a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta arr3b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cryba1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cryba1l1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cryba2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cryba4_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybb1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybb1l1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybb1l3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crybgx_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crygn2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Exorh_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gc2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gc3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gcap1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gcap3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnat1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnat2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnb1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnb1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnb3a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gngt1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gngt2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk7a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta grk7b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gucy2f_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rpe65a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rpe65b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta saga_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta Rrh_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sagb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Concatenation_subterraneus_percopsis_gadus_danio_vision.fa")
system("mv partitions.txt Concatenation_subterraneus_percopsis_gadus_danio_vision.partition")


#Extract species of interest
system("mv Concatenation_lamprologus_danio_vision.partition ../Concatenated_Alignments_Subset/Cichliformes/")
system("mv Concatenation_lamprologus_danio_vision.fa ../Concatenated_Alignments_Subset/Cichliformes/")


system("mv Concatenation_subterraneus_percopsis_gadus_danio_vision.fa ../Concatenated_Alignments_Subset/Percopsiformes/")
system("mv Concatenation_subterraneus_percopsis_gadus_danio_vision.partition ../Concatenated_Alignments_Subset/Percopsiformes/")


#Now, we can run codeml with CodonFreq = 2 and cleandata = 1 with different branch models




#### Datations L. lethops with concatenated vision genes  --------------------------------



#L. lethops <- SPA
#L. tigripictilis <- SPB
#N. brichardi <- SPC

#p-distances at each site computed on MEGA (Complete deletion model)
#file : Alignements/Concatenated_Alignment/All_species/Vision_Genes_Concatenated_Alignment.fa

p_AB_1 <- 0.0028363
p_AB_2 <- 0.0021616
p_AB_3 <- 0.0049973
p_AC_1 <- 0.0031064
p_AC_2 <- 0.0027020
p_AC_3 <- 0.0091842
p_BC_1 <- 0.0040519
p_BC_2 <- 0.0022967
p_BC_3 <- 0.0095894



d_AB_1 <- (-3/4)*log(1-(4/3)*p_AB_1)
d_AB_2 <- (-3/4)*log(1-(4/3)*p_AB_2)
d_AB_3 <- (-3/4)*log(1-(4/3)*p_AB_3)
d_AC_1 <- (-3/4)*log(1-(4/3)*p_AC_1)
d_AC_2 <- (-3/4)*log(1-(4/3)*p_AC_2)
d_AC_3 <- (-3/4)*log(1-(4/3)*p_AC_3)
d_BC_1 <- (-3/4)*log(1-(4/3)*p_BC_1)
d_BC_2 <- (-3/4)*log(1-(4/3)*p_BC_2)
d_BC_3 <- (-3/4)*log(1-(4/3)*p_BC_3)

#Divergence time between tigripic and brichardi (ref Aardema et al. 2020)
T_divergence <- 2355000


L1 <- (d_AB_1+d_AC_1-d_BC_1)/2
L2 <- (d_AB_2+d_AC_2-d_BC_2)/2
L3 <- (d_AB_3+d_AC_3-d_BC_3)/2

M1 <- (d_AB_1-d_AC_1+d_BC_1)/2
M2 <- (d_AB_2-d_AC_2+d_BC_2)/2
M3 <- (d_AB_3-d_AC_3+d_BC_3)/2

N1 <- (d_AB_1+d_AC_1+d_BC_1)/2
N2 <- (d_AB_2+d_AC_2+d_BC_2)/2
N3 <- (d_AB_3+d_AC_3+d_BC_3)/2

a1 <- d_BC_1/(2 * T_divergence)
a2 <- d_BC_2/(2 * T_divergence) 
a3 <- d_BC_3/(2 * T_divergence)

y1 <- d_AC_1 - d_BC_1
y2 <- d_AC_2 - d_BC_2
y3 <- d_AC_3 - d_BC_3


y12 <- (y1 + y2)/2
a12 <- (a1 + a2)/2

Tn_lethops <- (y12 - y3) / (a3 - a12)
Td_lamprologus <- ((d_AB_1 + d_AB_2 + d_AB_3) - (y1 + y2 + y3)) / (2 * (a1 + a2 + a3))

Td_over_Tn_lethops <- Td_lamprologus / Tn_lethops
Tn_over_Td_lethops <- Tn_lethops / Td_lamprologus




## Meredith

#dn/ds extracted from codeml (Concatenated_Alignments_Subset/Cichliformes/)


dn_ds_vision_tigripic <- 0.212706
dn_ds_vision_lethops <- 0.368392

TF_lethops <- (Td_lamprologus * (dn_ds_vision_lethops - 1)) / (dn_ds_vision_tigripic -1)
TP_meredith_lethops <- Td_lamprologus - TF_lethops



#### Datations T. subterrnaeus with concatenated vision genes  --------------------------------

#p-distances at each site computed on MEGA (Complete deletion model, folder Concatenated_Alignment)
#file : Alignements/Concatenated_Alignment/All_species/Vision_Genes_Concatenated_Alignment.fa

#species A       Typhlichthys subterraneus
#species B       Percopsis transmontana
#species C       Gadus morhua


p_AB_1 <- 0.0729335
p_AB_2 <- 0.0451229
p_AB_3 <- 0.2101567
p_AC_1 <- 0.1391140
p_AC_2 <- 0.0884896
p_AC_3 <- 0.3899244
p_BC_1 <- 0.1257428
p_BC_2 <- 0.0790327
p_BC_3 <- 0.3847920



d_AB_1 <- (-3/4)*log(1-(4/3)*p_AB_1)
d_AB_2 <- (-3/4)*log(1-(4/3)*p_AB_2)
d_AB_3 <- (-3/4)*log(1-(4/3)*p_AB_3)
d_AC_1 <- (-3/4)*log(1-(4/3)*p_AC_1)
d_AC_2 <- (-3/4)*log(1-(4/3)*p_AC_2)
d_AC_3 <- (-3/4)*log(1-(4/3)*p_AC_3)
d_BC_1 <- (-3/4)*log(1-(4/3)*p_BC_1)
d_BC_2 <- (-3/4)*log(1-(4/3)*p_BC_2)
d_BC_3 <- (-3/4)*log(1-(4/3)*p_BC_3)

#Divergence time between gadus morhua and percopsis estimated on TimeTree
T_divergence <- 132000000

L1 <- (d_AB_1+d_AC_1-d_BC_1)/2
L2 <- (d_AB_2+d_AC_2-d_BC_2)/2
L3 <- (d_AB_3+d_AC_3-d_BC_3)/2

M1 <- (d_AB_1-d_AC_1+d_BC_1)/2
M2 <- (d_AB_2-d_AC_2+d_BC_2)/2
M3 <- (d_AB_3-d_AC_3+d_BC_3)/2

N1 <- (d_AB_1+d_AC_1+d_BC_1)/2
N2 <- (d_AB_2+d_AC_2+d_BC_2)/2
N3 <- (d_AB_3+d_AC_3+d_BC_3)/2

a1 <- d_BC_1/(2 * T_divergence)
a2 <- d_BC_2/(2 * T_divergence) 
a3 <- d_BC_3/(2 * T_divergence)

y1 <- d_AC_1 - d_BC_1
y2 <- d_AC_2 - d_BC_2
y3 <- d_AC_3 - d_BC_3


y12 <- (y1 + y2)/2
a12 <- (a1 + a2)/2

Tn_subt <- (y12 - y3) / (a3 - a12)
Td_subt <- ((d_AB_1 + d_AB_2 + d_AB_3) - (y1 + y2 + y3)) / (2 * (a1 + a2 + a3))

Tn_over_Td_subt <- Tn_subt / Td_subt 


## Meredith


#dn/ds extracted from codeml (Concatenated_Alignments_Subset/Percopsiformes/)

dn_ds_vision_percopsis <- 0.120395
dn_ds_vision_typhlichthys <- 0.139143


TF_subt <- (Td_subt * (dn_ds_vision_typhlichthys - 1)) / (dn_ds_vision_percopsis -1)
TP_meredith_subt <- Td_subt - TF_subt



#### Datations L. dentata with concatenated vision genes  --------------------------------

#p-distances at each site computed on MEGA (Complete deletion model, folder Concatenated_Alignment)
#file : Alignements/Concatenated_Alignment/All_species/Vision_Genes_Concatenated_Alignment.fa

#species A       Lucifuga dentata
#species B       Lucifuga gibarensis
#species C       Brotula barbata


p_AB_1 <- 0.0085089
p_AB_2 <- 0.0054039
p_AB_3 <- 0.0132361
p_AC_1 <- 0.0614533
p_AC_2 <- 0.0356660
p_AC_3 <- 0.2462183
p_BC_1 <- 0.0586170
p_BC_2 <- 0.0343150
p_BC_3 <- 0.2456780



d_AB_1 <- (-3/4)*log(1-(4/3)*p_AB_1)
d_AB_2 <- (-3/4)*log(1-(4/3)*p_AB_2)
d_AB_3 <- (-3/4)*log(1-(4/3)*p_AB_3)
d_AC_1 <- (-3/4)*log(1-(4/3)*p_AC_1)
d_AC_2 <- (-3/4)*log(1-(4/3)*p_AC_2)
d_AC_3 <- (-3/4)*log(1-(4/3)*p_AC_3)
d_BC_1 <- (-3/4)*log(1-(4/3)*p_BC_1)
d_BC_2 <- (-3/4)*log(1-(4/3)*p_BC_2)
d_BC_3 <- (-3/4)*log(1-(4/3)*p_BC_3)

#Divergence time between Brotula barbata and Lucifuga dentata estimated on TimeTree
T_divergence <- 80000000

L1 <- (d_AB_1+d_AC_1-d_BC_1)/2
L2 <- (d_AB_2+d_AC_2-d_BC_2)/2
L3 <- (d_AB_3+d_AC_3-d_BC_3)/2

M1 <- (d_AB_1-d_AC_1+d_BC_1)/2
M2 <- (d_AB_2-d_AC_2+d_BC_2)/2
M3 <- (d_AB_3-d_AC_3+d_BC_3)/2

N1 <- (d_AB_1+d_AC_1+d_BC_1)/2
N2 <- (d_AB_2+d_AC_2+d_BC_2)/2
N3 <- (d_AB_3+d_AC_3+d_BC_3)/2

a1 <- d_BC_1/(2 * T_divergence)
a2 <- d_BC_2/(2 * T_divergence) 
a3 <- d_BC_3/(2 * T_divergence)

y1 <- d_AC_1 - d_BC_1
y2 <- d_AC_2 - d_BC_2
y3 <- d_AC_3 - d_BC_3


y12 <- (y1 + y2)/2
a12 <- (a1 + a2)/2

Tn_dent <- (y12 - y3) / (a3 - a12)
Td_dent <- ((d_AB_1 + d_AB_2 + d_AB_3) - (y1 + y2 + y3)) / (2 * (a1 + a2 + a3))

Tn_over_Td_dent <- Tn_dent / Td_dent


## Meredith


#dn/ds extracted from codeml (Concatenated_Alignment/All_species/Vision_codeml_CodonFreq2_NoClean_Freeratio.codeml)

dn_ds_vision_gibarensis <- 0.198452
dn_ds_vision_dentata <- 0.376274


TF_dent <- (Td_dent * (dn_ds_vision_dentata - 1)) / (dn_ds_vision_gibarensis -1)
TP_meredith_dent <- Td_dent - TF_dent






#### Mutpred2 on fish mutations   ---------------------------------

setwd("~/Desktop/These_Supplemental_Analysis/Vision_Genes_Mutpred")


#Parse mutpred2 output

system('grep -i -v "[A-Z].*,.*,#.*" Danio_rerio_Vision.rslt | cut -d "," -f1,2,3 > Danio_rerio_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" AGadus_morhua_Vision.rslt | cut -d "," -f1,2,3 > Gadus_morhua_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprologus_lethops_Vision.rslt | cut -d "," -f1,2,3 > Lamprologus_lethops_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Neolamprologus_brichardi_Vision.rslt | cut -d "," -f1,2,3 > Neolamprologus_brichardi_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprologus_tigripic_Vision.rslt | cut -d "," -f1,2,3 > Lamprologus_tigripic_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Percopsis_transmontana_Vision.rslt | cut -d "," -f1,2,3 > Percopsis_transmontana_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Typhlichthys_subterraneus_Vision.rslt | cut -d "," -f1,2,3 > Typhlichthys_subterraneus_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Astyanax_CF_Vision.rslt | cut -d "," -f1,2,3 > Astyanax_SF_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Astyanax_SF_Vision.rslt | cut -d "," -f1,2,3 > Astyanax_CF_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lucifuga_holguinensis_Vision.rslt | cut -d "," -f1,2,3 > Lucifuga_holguinensis_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprogrammus_exutus_Vision.rslt | cut -d "," -f1,2,3 > Lamprogrammus_exutus_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Brotula_barbata_Vision.rslt | cut -d "," -f1,2,3 > Brotula_barbata_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Carapus_acus_Vision.rslt | cut -d "," -f1,2,3 > Carapus_acus_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Pygocentrus_nattereri_Vision.rslt | cut -d "," -f1,2,3 > Pygocentrus_nattereri_Vision_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lucifuga_dentata_Vision.rslt | cut -d "," -f1,2,3 > Lucifuga_dentata_Vision_Scores.csv')




#Import tables

Danio_rerio_Vision_Scores <- read.table("Danio_rerio_Vision_Scores.csv", header=TRUE, sep=",")
Gadus_morhua_Vision_Scores <- read.table("Gadus_morhua_Vision_Scores.csv", header=TRUE, sep=",")
Lamprologus_lethops_Vision_Scores <- read.table("Lamprologus_lethops_Vision_Scores.csv", header=TRUE, sep=",")
Neolamprologus_brichardi_Vision_Scores <- read.table("Neolamprologus_brichardi_Vision_Scores.csv", header=TRUE, sep=",")
Lamprologus_tigripic_Vision_Scores <- read.table("Lamprologus_tigripic_Vision_Scores.csv", header=TRUE, sep=",")
Percopsis_transmontana_Vision_Scores <- read.table("Percopsis_transmontana_Vision_Scores.csv", header=TRUE, sep=",")
Typhlichthys_subterraneus_Vision_Scores <- read.table("Typhlichthys_subterraneus_Vision_Scores.csv", header=TRUE, sep=",")
Astyanax_CF_Vision_Scores <- read.table("Astyanax_CF_Vision_Scores.csv", header=TRUE, sep=",")
Astyanax_SF_Vision_Scores <- read.table("Astyanax_SF_Vision_Scores.csv", header=TRUE, sep=",")
Lucifuga_holguinensis_Vision_Scores <- read.table("Lucifuga_holguinensis_Vision_Scores.csv", header=TRUE, sep=",")
Lamprogrammus_exutus_Vision_Scores <- read.table("Lamprogrammus_exutus_Vision_Scores.csv", header=TRUE, sep=",")
Brotula_barbata_Vision_Scores <- read.table("Brotula_barbata_Vision_Scores.csv", header=TRUE, sep=",")
Carapus_acus_Vision_Scores <- read.table("Carapus_acus_Vision_Scores.csv", header=TRUE, sep=",")
Pygocentrus_nattereri_Vision_Scores <- read.table("Pygocentrus_nattereri_Vision_Scores.csv", header=TRUE, sep=",")
Lucifuga_dentata_Vision_Scores<- read.table("Lucifuga_dentata_Vision_Scores.csv", header=TRUE, sep=",")




Vision_mutpred2 <- bind_rows(
  Danio_rerio_Vision_Scores %>% mutate(Species = "Danio_rerio"),
  Gadus_morhua_Vision_Scores %>% mutate(Species = "Gadus_morhua"),
  Lamprologus_lethops_Vision_Scores %>% mutate(Species = "Lamprologus_lethops"),
  Neolamprologus_brichardi_Vision_Scores %>% mutate(Species = "Neolamprologus_brichardi"),
  Lamprologus_tigripic_Vision_Scores %>% mutate(Species = "Lamprologus_tigripic"),
  Percopsis_transmontana_Vision_Scores %>% mutate(Species = "Percopsis_transmontana"),
  Typhlichthys_subterraneus_Vision_Scores %>% mutate(Species = "Typhlichthys_subterraneus"),
  Astyanax_CF_Vision_Scores %>% mutate(Species = "Astyanax_CF"),
  Astyanax_SF_Vision_Scores %>% mutate(Species = "Astyanax_SF"),
  Lucifuga_holguinensis_Vision_Scores %>% mutate(Species = "Lucifuga_holguinensis"),
  Lamprogrammus_exutus_Vision_Scores %>% mutate(Species = "Lamprogrammus_exutus"),
  Brotula_barbata_Vision_Scores %>% mutate(Species = "Brotula_barbata"),
  Carapus_acus_Vision_Scores %>% mutate(Species = "Carapus_acus"),
  Pygocentrus_nattereri_Vision_Scores %>% mutate(Species = "Pygocentrus_nattereri"),
  Lucifuga_dentata_Vision_Scores %>% mutate(Species = "Lucifuga_dentata")
)



# Kernels densities

gene_to_keep <- Vision_mutpred2 %>% 
  filter(Species == "Percopsis_transmontana") %>%
  pull(ID) %>%
  unique()



Vision_mutpred2 %>%
  filter(Species == "Typhlichthys_subterraneus") %>%
  filter(ID %in% gene_to_keep) %>%
  pull(MutPred2.score) %>%
  mean()

Vision_mutpred2 %>%
  filter(Species == "Typhlichthys_subterraneus") %>%
  pull(MutPred2.score) %>%
  mean()






Vision_mutpred2 %>%
  ggplot(., aes(x=MutPred2.score, colour=Species)) +
  geom_density()

Vision_mutpred2 %>%
  ggplot(., aes(x=MutPred2.score, colour=Species)) +
  stat_ecdf(geom = "step", pad = FALSE)










#### Mutpred2 on neutral simulations   ---------------------------------


setwd("~/Desktop/These_Supplemental_Analysis/Vision_Genes_Mutpred/Simulations")

R_ratio_from_paml <- 1.85056  #extracted from Concatenated alignment
Number_mutations <- as.numeric((Vision_mutpred2 %>% group_by(Species) %>%
                                  summarise(value = n()) %>%
                                  arrange(value))[1,2])  #=65

File_Sequences_Danio <- Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
  mutate(seq_wo_stop = gsub('.{3}$', '', Sequence)) %>% dplyr::select(Gene, seq_wo_stop)



#Run100 simulations 



write.table(File_Sequences_Danio, file="File_Sequences_Danio.tsv", sep="\t", col.names = FALSE, row.names = FALSE)
system('sed -i "" "s/\"//g" File_Sequences_Danio.tsv')
system("for i in {1..100} ; do /Users/maxime/miniconda3/bin/python Neutral_evolution_for_mutpred.py -f File_Sequences_Danio.tsv -r 1.85056 -n 73 ; done")


#Run mutpred2 on every files


system("rm File_Sequences_Danio.tsv")
system("for i in * ; do run_mutpred2.sh -i $i -o $i.rslt ; done")
system("cat *.rslt > Simulations_Visions_MP2_Results")


#Perform graphs with observed mutations and simulated mutations

system('grep -v "[1-9].*,.*,#.*" Simulations_Visions_MP2_Results | cut -f1,2,3 -d "," | sed "s/---/,/g" | sed "s/^ID,/Simulation,ID,/g" | grep -v "Substitution" > Simulations_Visions_MP2_Results_Table.tsv')
Simulations_Vision_Scores <- read.table("Simulations_Visions_MP2_Results_Table.tsv", header=FALSE, sep=",")
colnames(Simulations_Vision_Scores) <- c("Simulation", "ID", "Substitution", "MutPred2.score")
Simulations_Vision_Scores <- Simulations_Vision_Scores %>% mutate(Species = "Simulation")

Simulations_Vision_Scores$MutPred2.score <- as.numeric(as.character(Simulations_Vision_Scores$MutPred2.score))
Simulations_Vision_Scores$ID <- as.character(Simulations_Vision_Scores$ID)
Simulations_Vision_Scores$Substitution <- as.character(Simulations_Vision_Scores$Substitution)


Simulations_Vision_Scores_v1 <- Simulations_Vision_Scores %>% dplyr::select(ID, Substitution, MutPred2.score, Species)
Vision_mutpred2 <- bind_rows(Vision_mutpred2, Simulations_Vision_Scores_v1)



Vision_mutpred2_new_sp <- Vision_mutpred2 %>% filter(
  Species == "Danio_rerio" | Species == "Gadus_morhua" | Species == "Lamprologus_lethops" | Species == "Neolamprologus_brichardi" | Species == "Lamprologus_tigripic" | Species == "Percopsis_transmontana" | Species == "Typhlichthys_subterraneus" | Species == "Simulation"
)

Vision_mutpred2_old_sp <- Vision_mutpred2 %>% filter(
  Species == "Astyanax_CF" | Species == "Astyanax_SF" | Species == "Lucifuga_holguinensis" | Species == "Lamprogrammus_exutus" | Species == "Brotula_barbata" | Species == "Carapus_acus" | Species == "Pygocentrus_nattereri" | Species == "Lucifuga_dentata"
)


pmutpred_vision <- Simulations_Vision_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  stat_ecdf(show.legend = FALSE, geom = "step", pad = FALSE) +
  stat_ecdf(data = Vision_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "solid") +
  stat_ecdf(data = Vision_mutpred2_old_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "dashed") +
  theme( panel.background = element_rect(fill = "white",
                                         colour = "black",
                                         size = 0.5, linetype = "solid"),
         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                         colour = "white"), 
         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                         colour = "white")) +
  scale_color_manual(values = c(rep("gray93", 100), "plum1", "deepskyblue", "darkolivegreen1", "darkolivegreen4", "darkviolet", "lawngreen", "khaki", "red", "darkorchid1", "sienna4", "sandybrown", "royalblue1", "paleturquoise1", "green", "black", "yellow")) +
  xlab("Mutpred2 score") +
  ylab("Cumulative frequency") 


pmutpred_vision_2 <- Simulations_Vision_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  stat_ecdf(show.legend = FALSE, geom = "step", pad = FALSE) +
  stat_ecdf(data = Vision_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "solid") +
  theme( panel.background = element_rect(fill = "white",
                                         colour = "black",
                                         size = 0.5, linetype = "solid"),
         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                         colour = "white"), 
         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                         colour = "white")) +
  scale_color_manual(values = c(rep("gray93", 100), "darkviolet", "lawngreen", "red", "darkorchid1", "royalblue1", "paleturquoise1", "black", "yellow")) +
  xlab("Mutpred2 score") +
  ylab("Cumulative frequency") 





#pmutpred_vision <- Simulations_Vision_Scores %>%
#  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
#  stat_ecdf(show.legend = FALSE, geom = "step", pad = FALSE) +
#  stat_ecdf(data = Vision_mutpred2, aes(x=MutPred2.score, colour=Species, linetype=linetype_sp), geom = "step", pad = FALSE, show.legend = FALSE) +
#  theme_minimal() +
#  scale_color_manual(values = c(rep("gray93", 100), "plum1", "deepskyblue", "darkolivegreen1", "darkolivegreen4", "darkviolet", "lawngreen", "khaki", "red", "darkorchid1", "sienna4", "sandybrown", "royalblue1", "paleturquoise1", "green", "black", "yellow")) +
#  xlab("Mutpred2 score") +
#  ylab("Cumulative frequency") #+
#  #scale_linetype_manual(values = c(rep("solid", 100), "dashed", "dashed", "dashed", "dashed", "solid", "solid", "dashed", "solid", "solid", "dashed", "dashed", "solid", "solid", "dashed", "solid", "solid"))





p2 <- Simulations_Vision_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  geom_density(show.legend = FALSE, trim=TRUE) +
  geom_density(data = Vision_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), show.legend = FALSE, trim=TRUE) +
  theme_minimal() +
  scale_color_manual(values = c(rep("gray", 100), "darkviolet", "lawngreen", "red", "darkorchid1", "royalblue1", "paleturquoise1", "black", "yellow")) 






Species_list <- Vision_mutpred2 %>%
  pull(Species) %>%
  unique(.)




ks.test(Vision_mutpred2 %>% filter(Species == "Lamprologus_lethops") %>% pull(MutPred2.score),
        Vision_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score))

ks.test(Vision_mutpred2 %>% filter(Species == "Lucifuga_dentata") %>% pull(MutPred2.score),
        Vision_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score))

ks.test(Vision_mutpred2 %>% filter(Species == "Astyanax_CF") %>% pull(MutPred2.score),
        Vision_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score))





Vision_mutpred2 %>% ggplot(., aes(x=MutPred2.score, colour=Species)) +
  stat_ecdf(geom = "step", pad = FALSE) +
  theme_minimal()




Vision_mutpred2 %>%
  group_by(Species) %>%
  summarise(value = n())

Gene_present_in_percopsis <- Vision_mutpred2 %>% filter(Species == "Percopsis_transmontana") %>% pull(ID) %>% unique()



#### Datations using mutpred2 results ---------------------------------
Neutral_mutations <- Vision_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score)
Danio_mutations <- Vision_mutpred2 %>% filter(Species == "Danio_rerio") %>% pull(MutPred2.score)
Lethops_mutations <- Vision_mutpred2 %>% filter(Species == "Lamprologus_lethops") %>% pull(MutPred2.score)
Subterraneus_mutations <- Vision_mutpred2 %>% filter(Species == "Typhlichthys_subterraneus") %>% filter(ID %in% Gene_present_in_percopsis) %>% pull(MutPred2.score)  #To be more precise, take only genes also presnent in P.transmontana, even if the diff in results is marginal
Percopsis_mutations <- Vision_mutpred2 %>% filter(Species == "Percopsis_transmontana") %>% pull(MutPred2.score)
Brichardi_mutations <- Vision_mutpred2 %>% filter(Species == "Neolamprologus_brichardi") %>% pull(MutPred2.score)
Morhua_mutations <- Vision_mutpred2 %>% filter(Species == "Gadus_morhua") %>% pull(MutPred2.score)
Tigripic_mutations <- Vision_mutpred2 %>% filter(Species == "Lamprologus_tigripic") %>% pull(MutPred2.score)
astyanaxCF_mutations <- Vision_mutpred2 %>% filter(Species == "Astyanax_CF") %>% pull(MutPred2.score)
astyanaxSF_mutations <- Vision_mutpred2 %>% filter(Species == "Astyanax_SF") %>% pull(MutPred2.score)
holguinensis_mutations <- Vision_mutpred2 %>% filter(Species == "Lucifuga_holguinensis") %>% pull(MutPred2.score)
lamprogrammus_mutations <- Vision_mutpred2 %>% filter(Species == "Lamprogrammus_exutus") %>% pull(MutPred2.score)
brotula_mutations <- Vision_mutpred2 %>% filter(Species == "Brotula_barbata") %>% pull(MutPred2.score)
carapus_mutations <- Vision_mutpred2 %>% filter(Species == "Carapus_acus") %>% pull(MutPred2.score)
pygocentrus_mutations <- Vision_mutpred2 %>% filter(Species == "Pygocentrus_nattereri") %>% pull(MutPred2.score)
dentata_mutations <- Vision_mutpred2 %>% filter(Species == "Lucifuga_dentata") %>% pull(MutPred2.score)



Vision_admixture_simulations <- data.frame(NULL)
Vision_admixture_lethops <- data.frame(NULL)
Vision_admixture_subterraneus <- data.frame(NULL)
Vision_admixture_percopsis <- data.frame(NULL)
Vision_admixture_brichardi <- data.frame(NULL)
Vision_admixture_morhua <- data.frame(NULL)
Vision_admixture_tigripic <- data.frame(NULL)
Vision_admixture_Astyanax_CF <- data.frame(NULL)
Vision_admixture_Astyanax_SF <- data.frame(NULL)
Vision_admixture_holguinensis <- data.frame(NULL)
Vision_admixture_lamprogrammus <- data.frame(NULL)
Vision_admixture_brotula <- data.frame(NULL)
Vision_admixture_carapus <- data.frame(NULL)
Vision_admixture_pygocentrus <- data.frame(NULL)
Vision_admixture_dentata <- data.frame(NULL)
Vision_admixture_danio <- data.frame(NULL)
iter=0
for(j in 1:1000){
  my_vector_test_lethops  <- c()
  my_vector_test_subterraneus  <- c()
  my_vector_test_percopsis  <- c() 
  my_vector_test_brichardi  <- c() 
  my_vector_test_morhua  <- c() 
  my_vector_test_tigripic  <- c()
  my_vector_test_Astyanax_CF  <- c()
  my_vector_test_Astyanax_SF  <- c()
  my_vector_test_holguinensis  <- c()
  my_vector_test_lamprogrammus  <- c()
  my_vector_test_brotula  <- c()
  my_vector_test_carapus  <- c()
  my_vector_test_pygocentrus  <- c()
  my_vector_test_dentata <- c()
  my_vector_test_danio <- c()
  my_vector_name <- c()
  for(i in 1:1000){
    my_mutations_final <- c(sample(Danio_mutations, 1000-i, replace = FALSE, prob = NULL), sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
    my_vector_test_lethops <- c(my_vector_test_lethops, ks.test(Lethops_mutations, my_mutations_final)$p)
    my_vector_test_subterraneus <- c(my_vector_test_subterraneus, ks.test(Subterraneus_mutations, my_mutations_final)$p)
    my_vector_test_percopsis <- c(my_vector_test_percopsis, ks.test(Percopsis_mutations, my_mutations_final)$p)
    my_vector_test_brichardi <- c(my_vector_test_brichardi, ks.test(Brichardi_mutations, my_mutations_final)$p)
    my_vector_test_morhua <- c(my_vector_test_morhua, ks.test(Morhua_mutations, my_mutations_final)$p)
    my_vector_test_tigripic <- c(my_vector_test_tigripic, ks.test(Tigripic_mutations, my_mutations_final)$p)
    my_vector_test_Astyanax_CF <- c(my_vector_test_Astyanax_CF, ks.test(astyanaxCF_mutations, my_mutations_final)$p)
    my_vector_test_Astyanax_SF <- c(my_vector_test_Astyanax_SF, ks.test(astyanaxSF_mutations, my_mutations_final)$p)
    my_vector_test_holguinensis <- c(my_vector_test_holguinensis, ks.test(holguinensis_mutations, my_mutations_final)$p)
    my_vector_test_lamprogrammus <- c(my_vector_test_lamprogrammus, ks.test(lamprogrammus_mutations, my_mutations_final)$p)
    my_vector_test_brotula <- c(my_vector_test_brotula, ks.test(brotula_mutations, my_mutations_final)$p)
    my_vector_test_carapus <- c(my_vector_test_carapus, ks.test(carapus_mutations, my_mutations_final)$p)
    my_vector_test_pygocentrus <- c(my_vector_test_pygocentrus, ks.test(pygocentrus_mutations, my_mutations_final)$p)
    my_vector_test_dentata <- c(my_vector_test_dentata, ks.test(dentata_mutations, my_mutations_final)$p)
    my_vector_test_danio <- c(my_vector_test_danio, ks.test(Danio_mutations, my_mutations_final)$p)
    my_vector_name <- c(my_vector_name, as.character(j))
  }
  
  iter=iter+1
  print(iter)
  
  Vision_admixture_lethops <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_lethops), as.data.frame(rep("lethops", 1000)))
  colnames(Vision_admixture_lethops) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Vision_admixture_subterraneus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_subterraneus), as.data.frame(rep("subterraneus", 1000)))
  colnames(Vision_admixture_subterraneus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Vision_admixture_percopsis <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_percopsis), as.data.frame(rep("percopsis", 1000)))
  colnames(Vision_admixture_percopsis) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Vision_admixture_brichardi <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_brichardi), as.data.frame(rep("brichardi", 1000)))
  colnames(Vision_admixture_brichardi) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Vision_admixture_morhua <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_morhua), as.data.frame(rep("morhua", 1000)))
  colnames(Vision_admixture_morhua) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Vision_admixture_tigripic <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_tigripic), as.data.frame(rep("tigripic", 1000)))
  colnames(Vision_admixture_tigripic) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Vision_admixture_Astyanax_CF <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_Astyanax_CF), as.data.frame(rep("Astyanax_CF", 1000)))
  colnames(Vision_admixture_Astyanax_CF) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Vision_admixture_Astyanax_SF <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_Astyanax_SF), as.data.frame(rep("Astyanax_SF", 1000)))
  colnames(Vision_admixture_Astyanax_SF) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Vision_admixture_holguinensis <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_holguinensis), as.data.frame(rep("holguinensis", 1000)))
  colnames(Vision_admixture_holguinensis) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Vision_admixture_lamprogrammus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_lamprogrammus), as.data.frame(rep("lamprogrammus", 1000)))
  colnames(Vision_admixture_lamprogrammus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Vision_admixture_brotula <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_brotula), as.data.frame(rep("brotula", 1000)))
  colnames(Vision_admixture_brotula) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Vision_admixture_carapus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_carapus), as.data.frame(rep("carapus", 1000)))
  colnames(Vision_admixture_carapus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Vision_admixture_pygocentrus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_pygocentrus), as.data.frame(rep("pygocentrus", 1000)))
  colnames(Vision_admixture_pygocentrus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Vision_admixture_dentata <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_dentata), as.data.frame(rep("dentata", 1000)))
  colnames(Vision_admixture_dentata) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Vision_admixture_danio <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_danio), as.data.frame(rep("danio", 1000)))
  colnames(Vision_admixture_danio) <- c("Simulation_nb", "KS_test_p", "Species")  
  

  
  
  
  Vision_admixture_simulations <- bind_rows(Vision_admixture_simulations,
                                            Vision_admixture_lethops,
                                            Vision_admixture_subterraneus, 
                                            Vision_admixture_percopsis,
                                            Vision_admixture_brichardi,
                                            Vision_admixture_morhua,
                                            Vision_admixture_tigripic,
                                            Vision_admixture_Astyanax_CF,
                                            Vision_admixture_Astyanax_SF,
                                            Vision_admixture_holguinensis,
                                            Vision_admixture_lamprogrammus,
                                            Vision_admixture_brotula,
                                            Vision_admixture_carapus,
                                            Vision_admixture_pygocentrus,
                                            Vision_admixture_dentata,
                                            Vision_admixture_danio
  )
  
}


Vision_admixture_simulations <- Vision_admixture_simulations %>% mutate(Neutral_number = rep((1:1000), 15000))


Vision_admixture_simulations_old_species <- Vision_admixture_simulations %>% filter(Species == "dentata"|
                                                                                    Species == "pygocentrus"|
                                                                                    Species == "carapus"|
                                                                                    Species == "brotula"|
                                                                                    Species == "lamprogrammus"|
                                                                                    Species == "holguinensis"|
                                                                                    Species == "Astyanax_SF"|
                                                                                    Species == "Astyanax_CF")

Vision_admixture_simulations_new_species <- Vision_admixture_simulations %>% filter(Species == "brichardi"| 
                                                                                    Species == "danio"|
                                                                                    Species == "lethops"|
                                                                                    Species == "morhua"|
                                                                                    Species == "percopsis"|
                                                                                    Species == "subterraneus"|
                                                                                    Species == "tigripic")




admix_vision <- Vision_admixture_simulations_new_species %>% 
  group_by(Species, Neutral_number) %>%
  summarise(moyenne = mean(KS_test_p)) %>%
  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
  geom_line(show.legend = FALSE) +
  geom_line(data = (Vision_admixture_simulations_old_species %>% group_by(Species, Neutral_number) %>% summarise(moyenne = mean(KS_test_p))), aes(x=Neutral_number, y=moyenne, color=Species), show.legend = FALSE, linetype="dashed") +
  xlab("Neutral mutations proportion") +
  ylab("KS-test pvalue") +
  theme_minimal() +
  scale_x_continuous(labels = c("0", "25", "50", "75", "100"))


admix_vision <- Vision_admixture_simulations_new_species %>% 
  group_by(Species, Neutral_number) %>%
  summarise(moyenne = mean(KS_test_p)) %>%
  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
  geom_line() +
  geom_line(data = (Vision_admixture_simulations_old_species %>% group_by(Species, Neutral_number) %>% summarise(moyenne = mean(KS_test_p))), aes(x=Neutral_number, y=moyenne, color=Species), linetype="dashed") +
  xlab("Neutral mutations proportion") +
  ylab("KS-test pvalue") +
  theme_minimal() +
  scale_x_continuous(labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(values = c("plum1", "deepskyblue","royalblue1", "darkolivegreen1", "darkolivegreen4", "darkviolet", "sienna4", "sandybrown", "khaki", "red", "lawngreen", "paleturquoise1", "green", "yellow", "darkorchid1" ))



admix_vision_2 <- Vision_admixture_simulations_new_species %>% 
  group_by(Species, Neutral_number) %>%
  summarise(moyenne = mean(KS_test_p)) %>%
  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
  geom_line(show.legend = FALSE) +
  xlab("Neutral mutations proportion") +
  ylab("KS-test pvalue") +
  theme_minimal() +
  scale_x_continuous(labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(values = c("royalblue1", "darkviolet", "red", "lawngreen", "paleturquoise1", "yellow", "darkorchid1" ))



#admix_vision <- Vision_admixture_simulations %>% 
#  group_by(Species, Neutral_number) %>%
#  summarise(moyenne = mean(KS_test_p)) %>%
#  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
#  geom_line(show.legend = FALSE) +
#  xlab("Neutral mutations proportion") +
#  ylab("KS-test pvalue") +
#  theme_minimal() +
#  scale_x_continuous(labels = c("0", "25", "50", "75", "100"))




Proportion_lethops <- head(Vision_admixture_simulations %>%
                             filter(Species == "lethops") %>%
                             group_by(Neutral_number) %>%
                             summarise(moyenne = mean(KS_test_p)) %>%
                             arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000

Proportion_subterraneus <- head(Vision_admixture_simulations %>%
                                  filter(Species == "subterraneus") %>%
                                  group_by(Neutral_number) %>%
                                  summarise(moyenne = mean(KS_test_p)) %>%
                                  arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000


Proportion_dentata <- head(Vision_admixture_simulations %>%
                                  filter(Species == "dentata") %>%
                                  group_by(Neutral_number) %>%
                                  summarise(moyenne = mean(KS_test_p)) %>%
                                  arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000


rslt_propor <- c()
species_curr_list <- Vision_admixture_simulations %>% pull(Species) %>% unique(.)
for (species_curr in species_curr_list){
  rslt_propor <- c(rslt_propor, head(Vision_admixture_simulations %>%
         filter(Species == species_curr) %>%
         group_by(Neutral_number) %>%
         summarise(moyenne = mean(KS_test_p)) %>%
         arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000)
}



cbind(species_curr_list, rslt_propor)






Time_div_subt <- Td_subt
Time_div_lethops <- Td_lamprologus
Time_div_dent <- Td_dent

dn_ds_vision_tigripic <- 0.212706
dn_ds_vision_percopsis <- 0.120395
dn_ds_vision_gibarensis <- 0.198452

dn_ds_vision_tigripic <- 0.08


datation_lethops_mp2 = (dn_ds_vision_tigripic/1) * (0.697 / (1-0.697)) * Time_div_lethops
datation_subterraneus_mp2 = (dn_ds_vision_percopsis/1) * (0.399 / (1-0.399)) * Time_div_subt
datation_dentata_mp2 = (dn_ds_vision_gibarensis/1) * (0.653 / (1-0.653)) * Time_div_dent


datation_lethops_mp2_min = (dn_ds_vision_tigripic/1) * (0.361 / (1-0.361)) * Time_div_lethops
datation_subterraneus_mp2_min = (dn_ds_vision_percopsis/1) * (0.283 / (1-0.283)) * Time_div_subt
datation_dentata_mp2_min = (dn_ds_vision_gibarensis/1) * (0.446 / (1-0.446)) * Time_div_dent



datation_lethops_mp2_max = (dn_ds_vision_tigripic/1) * (1 / (1-1)) * Time_div_lethops
datation_subterraneus_mp2_max = (dn_ds_vision_percopsis/1) * (0.482 / (1-0.482)) * Time_div_subt
datation_dentata_mp2_max = (dn_ds_vision_gibarensis/1) * (0.888 / (1-0.888)) * Time_div_dent




#### Summary datations results ---------------------------------



datation_lethops_mp2
datation_subterraneus_mp2

datation_meredith_subt <- TP_meredith_subt
datation_meredith_lethops <- TP_meredith_lethops

datation_simu_subt <- 569566
datation_simu_lethops <- 19187


summary_lethops <- c(datation_meredith_lethops, datation_lethops_mp2, datation_simu_lethops)
summary_subt <- c(datation_meredith_subt, datation_subterraneus_mp2, datation_simu_subt)














############### SECTION CIRCADIAN ############### 


#### Import datas  ---------------------------------


setwd("~/Desktop/New_Cavefishes_Analysis/Analysis_Scripts/Circadian")
Circadian_Sequences <- read.table("Circadian_Sequences_Table.tsv", header=TRUE, sep="\t")
Circadian_LoFs <- read.table("LoF_Table_Circadian.tsv", header=TRUE, sep="\t")


Circadian_Sequences %>% filter(Species == "Lamprologus_tigripic") %>%
  filter(Gene != "Not_Found") %>%
  group_by(Gene_Type) %>%
  summarise(n())


Circadian_Sequences %>% 
  filter(Gene != "Not_Found") %>%
  group_by(Species, Gene_Type) %>%
  summarise(n())



#### Perform alignments and run aaaml  ---------------------------------

#Perform the same commands performed for vision genes


#### Gene concatenations  ---------------------------------

setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Circadian_Clock")


system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i aanat1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta aanat2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bmal1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bmal2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ck1da_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ck1db_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clock1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clock2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clock3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry2a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dbpb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hlfa_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hlfb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nfil3_2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nfil3_6_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nfil3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d2a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d4a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d4b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta per1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta per2_1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rorab_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rorb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rorca_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rorcb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tefa_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tefb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Concatenation_subterraneus_percopsis_gadus_danio_circadian.fa")
system("mv partitions.txt Concatenation_subterraneus_percopsis_gadus_danio_circadian.partition")


system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i aanat1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta aanat2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bmal1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bmal2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ck1da_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ck1db_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clock1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clock2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clock3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry-dash_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry1a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry2a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cry5_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dbpa_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dbpb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hlfa_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hlfb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nfil3_2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nfil3_4_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nfil3_6_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nfil3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d2a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d2b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d4a_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nr1d4b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta per1b_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta per2_1_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta per2_2_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta per3_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta roraa_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rorb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rorca_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rorcb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tefa_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tefb_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Concatenation_lamprologus_danio_circadian.fa")
system("mv partitions.txt Concatenation_lamprologus_danio_circadian.partition")


system("mv Concatenation_lamprologus_danio_circadian.partition ../Concatenated_Alignments_Subset/Cichliformes/")
system("mv Concatenation_lamprologus_danio_circadian.fa ../Concatenated_Alignments_Subset/Cichliformes/")
system("mv Concatenation_subterraneus_percopsis_gadus_danio_circadian.fa ../Concatenated_Alignments_Subset/Percopsiformes/")
system("mv Concatenation_subterraneus_percopsis_gadus_danio_circadian.partition ../Concatenated_Alignments_Subset/Percopsiformes/")

system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i *.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out All_Species_circadian_concat.fa")
system("mv All_Species_circadian_concat.fa ../Concatenated_Alignment/")
system("mv partitions.txt ../Concatenated_Alignment/Circadian_concat.partition")


#### Mutpred2 on fish substitutions   ---------------------------------


setwd("~/Desktop/These_Supplemental_Analysis/Circadian_clock_Mutpred")


#Parse mutpred2 output

system('grep -i -v "[A-Z].*,.*,#.*" Danio_rerio_Vision.rslt | cut -d "," -f1,2,3 > Danio_rerio_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" AGadus_morhua_Vision.rslt | cut -d "," -f1,2,3 > Gadus_morhua_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprologus_lethops_Vision.rslt | cut -d "," -f1,2,3 > Lamprologus_lethops_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Neolamprologus_brichardi_Vision.rslt | cut -d "," -f1,2,3 > Neolamprologus_brichardi_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprologus_tigripic_Vision.rslt | cut -d "," -f1,2,3 > Lamprologus_tigripic_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Percopsis_transmontana_Vision.rslt | cut -d "," -f1,2,3 > Percopsis_transmontana_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Typhlichthys_subterraneus_Vision.rslt | cut -d "," -f1,2,3 > Typhlichthys_subterraneus_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Astyanax_CF_Vision.rslt | cut -d "," -f1,2,3 > Astyanax_SF_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Astyanax_SF_Vision.rslt | cut -d "," -f1,2,3 > Astyanax_CF_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lucifuga_holguinensis_Vision.rslt | cut -d "," -f1,2,3 > Lucifuga_holguinensis_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprogrammus_exutus_Vision.rslt | cut -d "," -f1,2,3 > Lamprogrammus_exutus_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Brotula_barbata_Vision.rslt | cut -d "," -f1,2,3 > Brotula_barbata_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Carapus_acus_Vision.rslt | cut -d "," -f1,2,3 > Carapus_acus_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Pygocentrus_nattereri_Vision.rslt | cut -d "," -f1,2,3 > Pygocentrus_nattereri_Circadian_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lucifuga_dentata_Vision.rslt | cut -d "," -f1,2,3 > Lucifuga_dentata_Circadian_Scores.csv')




#Import tables

Danio_rerio_Circadian_Scores <- read.table("Danio_rerio_Circadian_Scores.csv", header=TRUE, sep=",")
Gadus_morhua_Circadian_Scores <- read.table("Gadus_morhua_Circadian_Scores.csv", header=TRUE, sep=",")
Lamprologus_lethops_Circadian_Scores <- read.table("Lamprologus_lethops_Circadian_Scores.csv", header=TRUE, sep=",")
Neolamprologus_brichardi_Circadian_Scores <- read.table("Neolamprologus_brichardi_Circadian_Scores.csv", header=TRUE, sep=",")
Lamprologus_tigripic_Circadian_Scores <- read.table("Lamprologus_tigripic_Circadian_Scores.csv", header=TRUE, sep=",")
Percopsis_transmontana_Circadian_Scores <- read.table("Percopsis_transmontana_Circadian_Scores.csv", header=TRUE, sep=",")
Typhlichthys_subterraneus_Circadian_Scores <- read.table("Typhlichthys_subterraneus_Circadian_Scores.csv", header=TRUE, sep=",")
Astyanax_CF_Circadian_Scores <- read.table("Astyanax_CF_Circadian_Scores.csv", header=TRUE, sep=",")
Astyanax_SF_Circadian_Scores <- read.table("Astyanax_SF_Circadian_Scores.csv", header=TRUE, sep=",")
Lucifuga_holguinensis_Circadian_Scores <- read.table("Lucifuga_holguinensis_Circadian_Scores.csv", header=TRUE, sep=",")
Lamprogrammus_exutus_Circadian_Scores <- read.table("Lamprogrammus_exutus_Circadian_Scores.csv", header=TRUE, sep=",")
Brotula_barbata_Circadian_Scores <- read.table("Brotula_barbata_Circadian_Scores.csv", header=TRUE, sep=",")
Carapus_acus_Circadian_Scores <- read.table("Carapus_acus_Circadian_Scores.csv", header=TRUE, sep=",")
Pygocentrus_nattereri_Circadian_Scores <- read.table("Pygocentrus_nattereri_Circadian_Scores.csv", header=TRUE, sep=",")
Lucifuga_dentata_Circadian_Scores<- read.table("Lucifuga_dentata_Circadian_Scores.csv", header=TRUE, sep=",")




Circadian_mutpred2 <- bind_rows(
  Danio_rerio_Circadian_Scores %>% mutate(Species = "Danio_rerio"),
  Gadus_morhua_Circadian_Scores %>% mutate(Species = "Gadus_morhua"),
  Lamprologus_lethops_Circadian_Scores %>% mutate(Species = "Lamprologus_lethops"),
  Neolamprologus_brichardi_Circadian_Scores %>% mutate(Species = "Neolamprologus_brichardi"),
  Lamprologus_tigripic_Circadian_Scores %>% mutate(Species = "Lamprologus_tigripic"),
  Percopsis_transmontana_Circadian_Scores %>% mutate(Species = "Percopsis_transmontana"),
  Typhlichthys_subterraneus_Circadian_Scores %>% mutate(Species = "Typhlichthys_subterraneus"),
  Astyanax_CF_Circadian_Scores %>% mutate(Species = "Astyanax_CF"),
  Astyanax_SF_Circadian_Scores %>% mutate(Species = "Astyanax_SF"),
  Lucifuga_holguinensis_Circadian_Scores %>% mutate(Species = "Lucifuga_holguinensis"),
  Lamprogrammus_exutus_Circadian_Scores %>% mutate(Species = "Lamprogrammus_exutus"),
  Brotula_barbata_Circadian_Scores %>% mutate(Species = "Brotula_barbata"),
  Carapus_acus_Circadian_Scores %>% mutate(Species = "Carapus_acus"),
  Pygocentrus_nattereri_Circadian_Scores %>% mutate(Species = "Pygocentrus_nattereri"),
  Lucifuga_dentata_Circadian_Scores %>% mutate(Species = "Lucifuga_dentata")
)



# Kernels densities


Circadian_mutpred2 %>%
  ggplot(., aes(x=MutPred2.score, colour=Species)) +
  geom_density()

Circadian_mutpred2 %>%
  ggplot(., aes(x=MutPred2.score, colour=Species)) +
  stat_ecdf(geom = "step", pad = FALSE)




#### Mutpred2 on neutral simulations   ---------------------------------


setwd("~/Desktop/These_Supplemental_Analysis/Circadian_clock_Mutpred/Simulations")

R_ratio_from_paml <- 1.85056  #extracted from Concatenated alignment
Number_mutations <- as.numeric((Circadian_mutpred2 %>% group_by(Species) %>%
                                  summarise(value = n()) %>%
                                  arrange(value))[1,2])  #=65

File_Sequences_Danio <- Circadian_Sequences %>% filter(Species == "Danio_rerio") %>%
  mutate(seq_wo_stop = gsub('.{3}$', '', Sequence)) %>% dplyr::select(Gene, seq_wo_stop)



#Run100 simulations 



write.table(File_Sequences_Danio, file="File_Sequences_Danio.tsv", sep="\t", col.names = FALSE, row.names = FALSE)
system('sed -i "" "s/\"//g" File_Sequences_Danio.tsv')
system("for i in {1..100} ; do /Users/maxime/miniconda3/bin/python Neutral_evolution_for_mutpred.py -f File_Sequences_Danio.tsv -r 1.85056 -n 39 ; done")


#Run mutpred2 on every files


system("rm File_Sequences_Danio.tsv")
system("for i in * ; do run_mutpred2.sh -i $i -o $i.rslt ; done")
system("cat *.rslt > Simulations_Circadian_MP2_Results")


#Perform graphs with observed mutations and simulated mutations



#system('grep -v "[1-9].*,.*,#.*" Simulations_Circadian_MP2_Results | cut -f1,2,3 -d "," | sed "s/---/,/g" | sed "s/^ID,/Simulation,ID,/g" | grep -v "Substitution" > Simulations_Circadian_MP2_Results_Table.tsv')
Simulations_Circadian_Scores <- read.table("Simulations_Circadian_MP2_Results_Table.tsv", header=FALSE, sep=",")
colnames(Simulations_Circadian_Scores) <- c("Simulation", "ID", "Substitution", "MutPred2.score")
Simulations_Circadian_Scores <- Simulations_Circadian_Scores %>% mutate(Species = "Simulation")

Simulations_Circadian_Scores$MutPred2.score <- as.numeric(as.character(Simulations_Circadian_Scores$MutPred2.score))
Simulations_Circadian_Scores$ID <- as.character(Simulations_Circadian_Scores$ID)
Simulations_Circadian_Scores$Substitution <- as.character(Simulations_Circadian_Scores$Substitution)


Simulations_Circadian_Scores_v1 <- Simulations_Circadian_Scores %>% dplyr::select(ID, Substitution, MutPred2.score, Species)
Circadian_mutpred2 <- bind_rows(Circadian_mutpred2, Simulations_Circadian_Scores_v1)



Circadian_mutpred2_new_sp <- Circadian_mutpred2 %>% filter(
  Species == "Danio_rerio" | Species == "Gadus_morhua" | Species == "Lamprologus_lethops" | Species == "Neolamprologus_brichardi" | Species == "Lamprologus_tigripic" | Species == "Percopsis_transmontana" | Species == "Typhlichthys_subterraneus" | Species == "Simulation"
)

Circadian_mutpred2_old_sp <- Circadian_mutpred2 %>% filter(
  Species == "Astyanax_CF" | Species == "Astyanax_SF" | Species == "Lucifuga_holguinensis" | Species == "Lamprogrammus_exutus" | Species == "Brotula_barbata" | Species == "Carapus_acus" | Species == "Pygocentrus_nattereri" | Species == "Lucifuga_dentata"
)


pmutpred_circadian <- Simulations_Circadian_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  stat_ecdf(show.legend = FALSE, geom = "step", pad = FALSE) +
  stat_ecdf(data = Circadian_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "solid") +
  stat_ecdf(data = Circadian_mutpred2_old_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "dashed") +
  theme( panel.background = element_rect(fill = "white",
                                         colour = "black",
                                         size = 0.5, linetype = "solid"),
         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                         colour = "white"), 
         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                         colour = "white")) +
  scale_color_manual(values = c(rep("gray93", 100), "plum1", "deepskyblue", "darkolivegreen1", "darkolivegreen4", "darkviolet", "lawngreen", "khaki", "red", "darkorchid1", "sienna4", "sandybrown", "royalblue1", "paleturquoise1", "green", "black", "yellow")) +
  xlab("Mutpred2 score") +
  ylab("Cumulative frequency") 


pmutpred_circadian_2 <- Simulations_Circadian_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  stat_ecdf(show.legend = FALSE, geom = "step", pad = FALSE) +
  stat_ecdf(data = Circadian_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "solid") +
  theme( panel.background = element_rect(fill = "white",
                                         colour = "black",
                                         size = 0.5, linetype = "solid"),
         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                         colour = "white"), 
         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                         colour = "white")) +
  scale_color_manual(values = c(rep("gray93", 100), "darkviolet", "lawngreen", "red", "darkorchid1", "royalblue1", "paleturquoise1", "black", "yellow")) +
  xlab("Mutpred2 score") +
  ylab("Cumulative frequency") 



ks.test(Circadian_mutpred2 %>% filter(Species == "Lamprologus_lethops") %>% pull(MutPred2.score),
        Circadian_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score))

ks.test(Circadian_mutpred2 %>% filter(Species == "Lucifuga_dentata") %>% pull(MutPred2.score),
        Circadian_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score))

ks.test(Circadian_mutpred2 %>% filter(Species == "Astyanax_CF") %>% pull(MutPred2.score),
        Circadian_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score))






Circadian_mutpred2 %>%
  group_by(Species) %>%
  summarise(value = n())


Simulations_Circadian_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  geom_density(show.legend = FALSE) +
  geom_density(data = Circadian_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), show.legend = FALSE) +
  theme_minimal() +
  scale_color_manual(values = c(rep("gray", 100), "darkviolet", "lawngreen", "red", "darkorchid1", "royalblue1", "paleturquoise1", "black", "yellow")) 








#### Mutpred2 mixture results ---------------------------------
Neutral_mutations <- Circadian_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score)
Danio_mutations <- Circadian_mutpred2 %>% filter(Species == "Danio_rerio") %>% pull(MutPred2.score)
Lethops_mutations <- Circadian_mutpred2 %>% filter(Species == "Lamprologus_lethops") %>% pull(MutPred2.score)
Subterraneus_mutations <- Circadian_mutpred2 %>% filter(Species == "Typhlichthys_subterraneus") %>% pull(MutPred2.score)
Percopsis_mutations <- Circadian_mutpred2 %>% filter(Species == "Percopsis_transmontana") %>% pull(MutPred2.score)
Brichardi_mutations <- Circadian_mutpred2 %>% filter(Species == "Neolamprologus_brichardi") %>% pull(MutPred2.score)
Morhua_mutations <- Circadian_mutpred2 %>% filter(Species == "Gadus_morhua") %>% pull(MutPred2.score)
Tigripic_mutations <- Circadian_mutpred2 %>% filter(Species == "Lamprologus_tigripic") %>% pull(MutPred2.score)
astyanaxCF_mutations <- Circadian_mutpred2 %>% filter(Species == "Astyanax_CF") %>% pull(MutPred2.score)
astyanaxSF_mutations <- Circadian_mutpred2 %>% filter(Species == "Astyanax_SF") %>% pull(MutPred2.score)
holguinensis_mutations <- Circadian_mutpred2 %>% filter(Species == "Lucifuga_holguinensis") %>% pull(MutPred2.score)
lamprogrammus_mutations <- Circadian_mutpred2 %>% filter(Species == "Lamprogrammus_exutus") %>% pull(MutPred2.score)
brotula_mutations <- Circadian_mutpred2 %>% filter(Species == "Brotula_barbata") %>% pull(MutPred2.score)
carapus_mutations <- Circadian_mutpred2 %>% filter(Species == "Carapus_acus") %>% pull(MutPred2.score)
pygocentrus_mutations <- Circadian_mutpred2 %>% filter(Species == "Pygocentrus_nattereri") %>% pull(MutPred2.score)
dentata_mutations <- Circadian_mutpred2 %>% filter(Species == "Lucifuga_dentata") %>% pull(MutPred2.score)



Circadian_admixture_simulations <- data.frame(NULL)
Circadian_admixture_lethops <- data.frame(NULL)
Circadian_admixture_subterraneus <- data.frame(NULL)
Circadian_admixture_percopsis <- data.frame(NULL)
Circadian_admixture_brichardi <- data.frame(NULL)
Circadian_admixture_morhua <- data.frame(NULL)
Circadian_admixture_tigripic <- data.frame(NULL)
Circadian_admixture_Astyanax_CF <- data.frame(NULL)
Circadian_admixture_Astyanax_SF <- data.frame(NULL)
Circadian_admixture_holguinensis <- data.frame(NULL)
Circadian_admixture_lamprogrammus <- data.frame(NULL)
Circadian_admixture_brotula <- data.frame(NULL)
Circadian_admixture_carapus <- data.frame(NULL)
Circadian_admixture_pygocentrus <- data.frame(NULL)
Circadian_admixture_dentata <- data.frame(NULL)
Circadian_admixture_danio <- data.frame(NULL)
iter=0
for(j in 1:1000){
  my_vector_test_lethops  <- c()
  my_vector_test_subterraneus  <- c()
  my_vector_test_percopsis  <- c() 
  my_vector_test_brichardi  <- c() 
  my_vector_test_morhua  <- c() 
  my_vector_test_tigripic  <- c()
  my_vector_test_Astyanax_CF  <- c()
  my_vector_test_Astyanax_SF  <- c()
  my_vector_test_holguinensis  <- c()
  my_vector_test_lamprogrammus  <- c()
  my_vector_test_brotula  <- c()
  my_vector_test_carapus  <- c()
  my_vector_test_pygocentrus  <- c()
  my_vector_test_dentata <- c()
  my_vector_test_danio <- c()
  my_vector_name <- c()
  for(i in 1:1000){
    my_mutations_final <- c(sample(Danio_mutations, 1000-i, replace = FALSE, prob = NULL), sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
    my_vector_test_lethops <- c(my_vector_test_lethops, ks.test(Lethops_mutations, my_mutations_final)$p)
    my_vector_test_subterraneus <- c(my_vector_test_subterraneus, ks.test(Subterraneus_mutations, my_mutations_final)$p)
    my_vector_test_percopsis <- c(my_vector_test_percopsis, ks.test(Percopsis_mutations, my_mutations_final)$p)
    my_vector_test_brichardi <- c(my_vector_test_brichardi, ks.test(Brichardi_mutations, my_mutations_final)$p)
    my_vector_test_morhua <- c(my_vector_test_morhua, ks.test(Morhua_mutations, my_mutations_final)$p)
    my_vector_test_tigripic <- c(my_vector_test_tigripic, ks.test(Tigripic_mutations, my_mutations_final)$p)
    my_vector_test_Astyanax_CF <- c(my_vector_test_Astyanax_CF, ks.test(astyanaxCF_mutations, my_mutations_final)$p)
    my_vector_test_Astyanax_SF <- c(my_vector_test_Astyanax_SF, ks.test(astyanaxSF_mutations, my_mutations_final)$p)
    my_vector_test_holguinensis <- c(my_vector_test_holguinensis, ks.test(holguinensis_mutations, my_mutations_final)$p)
    my_vector_test_lamprogrammus <- c(my_vector_test_lamprogrammus, ks.test(lamprogrammus_mutations, my_mutations_final)$p)
    my_vector_test_brotula <- c(my_vector_test_brotula, ks.test(brotula_mutations, my_mutations_final)$p)
    my_vector_test_carapus <- c(my_vector_test_carapus, ks.test(carapus_mutations, my_mutations_final)$p)
    my_vector_test_pygocentrus <- c(my_vector_test_pygocentrus, ks.test(pygocentrus_mutations, my_mutations_final)$p)
    my_vector_test_dentata <- c(my_vector_test_dentata, ks.test(dentata_mutations, my_mutations_final)$p)
    my_vector_test_danio <- c(my_vector_test_danio, ks.test(Danio_mutations, my_mutations_final)$p)
    my_vector_name <- c(my_vector_name, as.character(j))
    

  }
  Circadian_admixture_lethops <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_lethops), as.data.frame(rep("lethops", 1000)))
  colnames(Circadian_admixture_lethops) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Circadian_admixture_subterraneus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_subterraneus), as.data.frame(rep("subterraneus", 1000)))
  colnames(Circadian_admixture_subterraneus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Circadian_admixture_percopsis <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_percopsis), as.data.frame(rep("percopsis", 1000)))
  colnames(Circadian_admixture_percopsis) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Circadian_admixture_brichardi <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_brichardi), as.data.frame(rep("brichardi", 1000)))
  colnames(Circadian_admixture_brichardi) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Circadian_admixture_morhua <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_morhua), as.data.frame(rep("morhua", 1000)))
  colnames(Circadian_admixture_morhua) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Circadian_admixture_tigripic <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_tigripic), as.data.frame(rep("tigripic", 1000)))
  colnames(Circadian_admixture_tigripic) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_Astyanax_CF <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_Astyanax_CF), as.data.frame(rep("Astyanax_CF", 1000)))
  colnames(Circadian_admixture_Astyanax_CF) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Circadian_admixture_Astyanax_SF <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_Astyanax_SF), as.data.frame(rep("Astyanax_SF", 1000)))
  colnames(Circadian_admixture_Astyanax_SF) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_holguinensis <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_holguinensis), as.data.frame(rep("holguinensis", 1000)))
  colnames(Circadian_admixture_holguinensis) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_lamprogrammus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_lamprogrammus), as.data.frame(rep("lamprogrammus", 1000)))
  colnames(Circadian_admixture_lamprogrammus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_brotula <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_brotula), as.data.frame(rep("brotula", 1000)))
  colnames(Circadian_admixture_brotula) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_carapus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_carapus), as.data.frame(rep("carapus", 1000)))
  colnames(Circadian_admixture_carapus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_pygocentrus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_pygocentrus), as.data.frame(rep("pygocentrus", 1000)))
  colnames(Circadian_admixture_pygocentrus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_dentata <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_dentata), as.data.frame(rep("dentata", 1000)))
  colnames(Circadian_admixture_dentata) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_danio <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_danio), as.data.frame(rep("danio", 1000)))
  colnames(Circadian_admixture_danio) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Circadian_admixture_simulations <- bind_rows(Circadian_admixture_simulations,
                                               Circadian_admixture_lethops,
                                               Circadian_admixture_subterraneus, 
                                               Circadian_admixture_percopsis,
                                               Circadian_admixture_brichardi,
                                               Circadian_admixture_morhua,
                                               Circadian_admixture_tigripic,
                                               Circadian_admixture_Astyanax_CF,
                                               Circadian_admixture_Astyanax_SF,
                                               Circadian_admixture_holguinensis,
                                               Circadian_admixture_lamprogrammus,
                                               Circadian_admixture_brotula,
                                               Circadian_admixture_carapus,
                                               Circadian_admixture_pygocentrus,
                                               Circadian_admixture_dentata,
                                               Circadian_admixture_danio
  )
  
  iter=iter+1
  print(iter)
  
}


Circadian_admixture_simulations <- Circadian_admixture_simulations %>% mutate(Neutral_number = rep((1:1000), 15000))


Circadian_admixture_simulations_old_species <- Circadian_admixture_simulations %>% filter(Species == "dentata"|
                                                                                      Species == "pygocentrus"|
                                                                                      Species == "carapus"|
                                                                                      Species == "brotula"|
                                                                                      Species == "lamprogrammus"|
                                                                                      Species == "holguinensis"|
                                                                                      Species == "Astyanax_SF"|
                                                                                      Species == "Astyanax_CF")

Circadian_admixture_simulations_new_species <- Circadian_admixture_simulations %>% filter(Species == "brichardi"| 
                                                                                      Species == "danio"|
                                                                                      Species == "lethops"|
                                                                                      Species == "morhua"|
                                                                                      Species == "percopsis"|
                                                                                      Species == "subterraneus"|
                                                                                      Species == "tigripic")




admix_Circadian <- Circadian_admixture_simulations_new_species %>% 
  group_by(Species, Neutral_number) %>%
  summarise(moyenne = mean(KS_test_p)) %>%
  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
  geom_line(show.legend = FALSE) +
  geom_line(data = (Circadian_admixture_simulations_old_species %>% group_by(Species, Neutral_number) %>% summarise(moyenne = mean(KS_test_p))), aes(x=Neutral_number, y=moyenne, color=Species), show.legend = FALSE, linetype="dashed") +
  xlab("Neutral mutations proportion") +
  ylab("KS-test pvalue") +
  theme_minimal() +
  scale_x_continuous(labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(values = c("plum1", "deepskyblue","royalblue1", "darkolivegreen1", "darkolivegreen4", "darkviolet", "sienna4", "sandybrown", "khaki", "red", "lawngreen", "paleturquoise1", "green", "yellow", "darkorchid1" ))


admix_circadian_2 <- Circadian_admixture_simulations_new_species %>% 
  group_by(Species, Neutral_number) %>%
  summarise(moyenne = mean(KS_test_p)) %>%
  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
  geom_line(show.legend = FALSE) +
  xlab("Neutral mutations proportion") +
  ylab("KS-test pvalue") +
  theme_minimal() +
  scale_x_continuous(labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(values = c("royalblue1", "darkviolet", "red", "lawngreen", "paleturquoise1", "yellow", "darkorchid1" ))



#admix_Circadian <- Circadian_admixture_simulations %>% 
#  group_by(Species, Neutral_number) %>%
#  summarise(moyenne = mean(KS_test_p)) %>%
#  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
#  geom_line(show.legend = FALSE) +
#  xlab("Neutral mutations proportion") +
#  ylab("KS-test pvalue") +
#  theme_minimal() +
#  scale_x_continuous(labels = c("0", "25", "50", "75", "100"))




Proportion_lethops <- head(Circadian_admixture_simulations %>%
                             filter(Species == "lethops") %>%
                             group_by(Neutral_number) %>%
                             summarise(moyenne = mean(KS_test_p)) %>%
                             arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000

Proportion_subterraneus <- head(Circadian_admixture_simulations %>%
                                  filter(Species == "subterraneus") %>%
                                  group_by(Neutral_number) %>%
                                  summarise(moyenne = mean(KS_test_p)) %>%
                                  arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000


rslt_propor <- c()
species_curr_list <- Circadian_admixture_simulations %>% pull(Species) %>% unique(.)
for (species_curr in species_curr_list){
  rslt_propor <- c(rslt_propor, head(Circadian_admixture_simulations %>%
                                       filter(Species == species_curr) %>%
                                       group_by(Neutral_number) %>%
                                       summarise(moyenne = mean(KS_test_p)) %>%
                                       arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000)
}



cbind(species_curr_list, rslt_propor)






############### SECTION PIGMENTATION ############### 

#### Import datas  ---------------------------------


setwd("~/Desktop/New_Cavefishes_Analysis/Analysis_Scripts/Pigmentation/")
Pigmentation_Sequences <- read.table("Pigmentation_Table.tsv", header=TRUE, sep="\t")
Pigmentation_LoFs <- read.table("LoF_Table_Pigmentation.tsv", header=TRUE, sep="\t")


Pigmentation_Sequences %>% filter(Species == "Lamprologus_tigripic") %>%
  filter(Gene != "Not_Found") %>%
  group_by(Gene_Type) %>%
  summarise(n())


Pigmentation_Sequences %>% 
  filter(Gene != "Not_Found") %>%
  group_by(Species, Gene_Type) %>%
  summarise(n())


#### Perform alignments and run aaml  ---------------------------------

#Same commands than those performed for vision and circadian clock genes





#### Gene concatenations  ---------------------------------

setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Pigmentation")

system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i adam17a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta adam17b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta adamts20.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta adrb2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ankrd27.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap1g1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap1m1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap3b1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap3d1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta apc.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta arcn1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta asip1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta atp6v0b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta atrn.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bcl2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s6.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bnc2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta brsk2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta brsk2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta c10orf11.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cd63.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cdh2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta chm.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cited1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clcn7.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta creb1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta creb1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crhb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta csf1a.seqlist.fa.NT_ALN.fasta.align_noFS_NT.fasta csf1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta csf1ra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ctns.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dct.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dctn1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dctn1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dctn2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dock7.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta drd2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dtnbp1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ece1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ece2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ece2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta eda.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta edar.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta edn3b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ednrb2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ednrba.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta egfra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta en1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta erbb3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta erbb3b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fbxw4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fgfr2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fhl2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fig4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta foxd3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta frem2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta frem2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fzd4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gart.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gata3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gch1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gch2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gchfr.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gfpt1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ghra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ghrb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gja5a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gli3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gna11a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnaq.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpc3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpnmb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpr143.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpr161.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hdac1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps6.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hsd3b1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta igsf11.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ikbkg.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta impdh1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta impdh1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ippk.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta irf4a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta irf4b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta itgb1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta itgb1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kcnj13.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kif13a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kita.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kitb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kitlga.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta lef1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta leo1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ltk.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta lyst.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta map2k1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mbtps1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mc1r.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mcoln3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta med12.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mef2ca.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mef2cb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mgrn1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mib1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mib2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mitfa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mitfb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mlpha.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mpv17.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mreg.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mtrex.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mycb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mycbp2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myg1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myo5aa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myo5ab.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myo7aa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myo7ab.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nf1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nrarpa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nsfb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta oca2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta oprm1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ostm1.seqlist.fa.NT_ALN.fasta.align_noFS_NT.fasta ovol1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pah.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta paics.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pax3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pax3b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pax7a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pcbd1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pmela.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pmelb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pnp4a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pnp4b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pomca.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta psen1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pts.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta qdpra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab11a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab11ba.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab11bb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab1aa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab1ab.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab27a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab32a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab38a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab38b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab3ip.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab8a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rabggta.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta recql4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ric8b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rnf41.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta scarb2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta scarb2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta scg2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sf3b1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sfxn1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta shroom2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta shroom2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc24a4b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc24a5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a11a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a11l.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a15a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a15b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc45a2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc7a11.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta smtla.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta snai2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta snapin.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox10.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox9a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox9b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta spra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tfap2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tfap2e.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta th.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tmem33.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tpcn2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trappc6bl.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trim33.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trpm1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trpm1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trpm7.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tspan10.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta txndc5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tyr.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tyrp1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta usp13.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vat1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps11.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps18.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps33a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps39.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta wnt1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta wnt3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta xdh.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta zeb2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta zeb2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta zic2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Concatenation_subterraneus_percopsis_gadus_danio_pigmentation.fa")
system("mv partitions.txt Concatenation_subterraneus_percopsis_gadus_danio_pigmentation.partition")


system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i adam17a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta adam17b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta adamts20.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta adrb2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ankrd27.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap1g1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap1m1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap3b1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ap3d1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta apc.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta arcn1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta asip1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta atp6v0b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta atrn.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bcl2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bloc1s6.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta bnc2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta brsk2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta brsk2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta c10orf11.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cd63.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cdh2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta chm.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta cited1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta clcn7.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta corin.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta creb1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta creb1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta crhb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta csf1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta csf1ra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ctns.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dct.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dctn1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dctn1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dctn2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dock7.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dtnbp1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta dtnbp1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ece1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ece2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ece2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta eda.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta edar.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta edn3b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ednrb2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ednrba.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta egfra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta en1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta en1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta erbb3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta erbb3b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fbxw4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fgfr2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fhl2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fhl2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fig4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta foxd3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta frem2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta frem2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta fzd4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gart.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gata3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gch1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gch2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gchfr.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gfpt1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ghra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ghrb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gja5a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gli3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gna11a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gna11b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gnaq.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpc3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpnmb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpr143.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta gpr161.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hdac1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps3.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hps6.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta hsd3b1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta igsf11.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ikbkg.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta impdh1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta impdh1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ippk.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta irf4a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta irf4b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta itgb1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta itgb1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kcnj13.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kif13a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kita.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kitb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta kitlga.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta lef1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta leo1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta lmx1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ltk.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta lyst.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta map2k1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mbtps1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mc1r.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mcoln3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mcoln3b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta med12.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mef2cb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mgrn1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mib1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mib2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mitfa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mitfb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mlpha.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mlphb.seqlist.fa.NT_ALN.fasta.align_noFS_NT.fasta mpv17.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mreg.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mtrex.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mycb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta mycbp2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myg1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myo5ab.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myo7aa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta myo7ab.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nf1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nrarpa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nsfa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta nsfb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta oca2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta oprm1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ostm1.seqlist.fa.NT_ALN.fasta.align_noFS_NT.fasta ovol1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pah.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta paics.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pax3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pax3b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pax7a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pax7b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pcbd1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pmela.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pmelb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pnp4a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pnp4b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pomca.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta psen1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta pts.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta qdpra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab11a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab11ba.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab11bb.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab17.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab1aa.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab1ab.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab27a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab32a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab38a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab38b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab3ip.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rab8a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rabggta.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta recql4.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta ric8b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta rnf41.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta scarb2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta scarb2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta scg2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sf3b1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sfxn1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta shroom2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta shroom2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc24a4b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc24a5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a11a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a11b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a11l.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a15a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc2a15b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta slc45a2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta smtla.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta snai2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta snapin.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox10.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox9a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta sox9b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta spra.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tfap2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tfap2e.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta th.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tmem33.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tpcn2.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trappc6bl.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trim33.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trpm1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trpm1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta trpm7.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tspan10.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta txndc5.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tyr.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tyrp1a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta tyrp1b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta usp13.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vat1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps11.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps18.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps33a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta vps39.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta wnt1.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta wnt3a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta xdh.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta zeb2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta zeb2b.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta zic2a.seqlist_NT.fa.NT_ALN.fasta.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Concatenation_lamprologus_danio_pigmentation.fa")
system("mv partitions.txt Concatenation_lamprologus_danio_pigmentation.partition")

system("mv Concatenation_lamprologus_danio_pigmentation.partition ../Concatenated_Alignments_Subset/Cichliformes/")
system("mv Concatenation_lamprologus_danio_pigmentation.fa ../Concatenated_Alignments_Subset/Cichliformes/")
system("mv Concatenation_subterraneus_percopsis_gadus_danio_pigmentation.fa ../Concatenated_Alignments_Subset/Percopsiformes/")
system("mv Concatenation_subterraneus_percopsis_gadus_danio_pigmentation.partition ../Concatenated_Alignments_Subset/Percopsiformes/")

system("python /Users/maxime/Downloads/AMAS-master/amas/AMAS.py concat -i *.align_noFS_NT.fasta -f fasta -d dna -u fasta --concat-out Pigmentation_concatenation_all_species.fa")
system("mv Pigmentation_concatenation_all_species.fa ../Concatenated_Alignment/")
system("mv partitions.txt ../Concatenated_Alignment/Pigmentation_concatenation.partition")



#### Mutpred2 on fish substitutions   ---------------------------------


setwd("~/Desktop/These_Supplemental_Analysis/Pigmentation_Genes_Mutpred")


#Parse mutpred2 output

system('grep -i -v "[A-Z].*,.*,#.*" Danio_rerio_Vision.rslt | cut -d "," -f1,2,3 > Danio_rerio_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" AGadus_morhua_Vision.rslt | cut -d "," -f1,2,3 > Gadus_morhua_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprologus_lethops_Vision.rslt | cut -d "," -f1,2,3 > Lamprologus_lethops_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Neolamprologus_brichardi_Vision.rslt | cut -d "," -f1,2,3 > Neolamprologus_brichardi_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprologus_tigripic_Vision.rslt | cut -d "," -f1,2,3 > Lamprologus_tigripic_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Percopsis_transmontana_Vision.rslt | cut -d "," -f1,2,3 > Percopsis_transmontana_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Typhlichthys_subterraneus_Vision.rslt | cut -d "," -f1,2,3 > Typhlichthys_subterraneus_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Astyanax_CF_Vision.rslt | cut -d "," -f1,2,3 > Astyanax_SF_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Astyanax_SF_Vision.rslt | cut -d "," -f1,2,3 > Astyanax_CF_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lucifuga_holguinensis_Vision.rslt | cut -d "," -f1,2,3 > Lucifuga_holguinensis_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lamprogrammus_exutus_Vision.rslt | cut -d "," -f1,2,3 > Lamprogrammus_exutus_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Brotula_barbata_Vision.rslt | cut -d "," -f1,2,3 > Brotula_barbata_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Carapus_acus_Vision.rslt | cut -d "," -f1,2,3 > Carapus_acus_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Pygocentrus_nattereri_Vision.rslt | cut -d "," -f1,2,3 > Pygocentrus_nattereri_Pigmentation_Scores.csv')
system('grep -i -v "[A-Z].*,.*,#.*" Lucifuga_dentata_Vision.rslt | cut -d "," -f1,2,3 > Lucifuga_dentata_Pigmentation_Scores.csv')




#Import tables

Danio_rerio_Pigmentation_Scores <- read.table("Danio_rerio_Pigmentation_Scores.csv", header=TRUE, sep=",")
Gadus_morhua_Pigmentation_Scores <- read.table("Gadus_morhua_Pigmentation_Scores.csv", header=TRUE, sep=",")
Lamprologus_lethops_Pigmentation_Scores <- read.table("Lamprologus_lethops_Pigmentation_Scores.csv", header=TRUE, sep=",")
Neolamprologus_brichardi_Pigmentation_Scores <- read.table("Neolamprologus_brichardi_Pigmentation_Scores.csv", header=TRUE, sep=",")
Lamprologus_tigripic_Pigmentation_Scores <- read.table("Lamprologus_tigripic_Pigmentation_Scores.csv", header=TRUE, sep=",")
Percopsis_transmontana_Pigmentation_Scores <- read.table("Percopsis_transmontana_Pigmentation_Scores.csv", header=TRUE, sep=",")
Typhlichthys_subterraneus_Pigmentation_Scores <- read.table("Typhlichthys_subterraneus_Pigmentation_Scores.csv", header=TRUE, sep=",")
Astyanax_CF_Pigmentation_Scores <- read.table("Astyanax_CF_Pigmentation_Scores.csv", header=TRUE, sep=",")
Astyanax_SF_Pigmentation_Scores <- read.table("Astyanax_SF_Pigmentation_Scores.csv", header=TRUE, sep=",")
Lucifuga_holguinensis_Pigmentation_Scores <- read.table("Lucifuga_holguinensis_Pigmentation_Scores.csv", header=TRUE, sep=",")
Lamprogrammus_exutus_Pigmentation_Scores <- read.table("Lamprogrammus_exutus_Pigmentation_Scores.csv", header=TRUE, sep=",")
Brotula_barbata_Pigmentation_Scores <- read.table("Brotula_barbata_Pigmentation_Scores.csv", header=TRUE, sep=",")
Carapus_acus_Pigmentation_Scores <- read.table("Carapus_acus_Pigmentation_Scores.csv", header=TRUE, sep=",")
Pygocentrus_nattereri_Pigmentation_Scores <- read.table("Pygocentrus_nattereri_Pigmentation_Scores.csv", header=TRUE, sep=",")
Lucifuga_dentata_Pigmentation_Scores<- read.table("Lucifuga_dentata_Pigmentation_Scores.csv", header=TRUE, sep=",")




Pigmentation_mutpred2 <- bind_rows(
  Danio_rerio_Pigmentation_Scores %>% mutate(Species = "Danio_rerio"),
  Gadus_morhua_Pigmentation_Scores %>% mutate(Species = "Gadus_morhua"),
  Lamprologus_lethops_Pigmentation_Scores %>% mutate(Species = "Lamprologus_lethops"),
  Neolamprologus_brichardi_Pigmentation_Scores %>% mutate(Species = "Neolamprologus_brichardi"),
  Lamprologus_tigripic_Pigmentation_Scores %>% mutate(Species = "Lamprologus_tigripic"),
  Percopsis_transmontana_Pigmentation_Scores %>% mutate(Species = "Percopsis_transmontana"),
  Typhlichthys_subterraneus_Pigmentation_Scores %>% mutate(Species = "Typhlichthys_subterraneus"),
  Astyanax_CF_Pigmentation_Scores %>% mutate(Species = "Astyanax_CF"),
  Astyanax_SF_Pigmentation_Scores %>% mutate(Species = "Astyanax_SF"),
  Lucifuga_holguinensis_Pigmentation_Scores %>% mutate(Species = "Lucifuga_holguinensis"),
  Lamprogrammus_exutus_Pigmentation_Scores %>% mutate(Species = "Lamprogrammus_exutus"),
  Brotula_barbata_Pigmentation_Scores %>% mutate(Species = "Brotula_barbata"),
  Carapus_acus_Pigmentation_Scores %>% mutate(Species = "Carapus_acus"),
  Pygocentrus_nattereri_Pigmentation_Scores %>% mutate(Species = "Pygocentrus_nattereri"),
  Lucifuga_dentata_Pigmentation_Scores %>% mutate(Species = "Lucifuga_dentata")
)



# Kernels densities


Pigmentation_mutpred2 %>%
  ggplot(., aes(x=MutPred2.score, colour=Species)) +
  geom_density()

Pigmentation_mutpred2 %>%
  ggplot(., aes(x=MutPred2.score, colour=Species)) +
  stat_ecdf(geom = "step", pad = FALSE)




#### Mutpred2 on neutral simulations   ---------------------------------


setwd("~/Desktop/These_Supplemental_Analysis/Pigmentation_Genes_Mutpred/Simulations")

R_ratio_from_paml <- 1.85056  #extracted from Concatenated alignment
Number_mutations <- as.numeric((Pigmentation_mutpred2 %>% group_by(Species) %>%
                                  summarise(value = n()) %>%
                                  arrange(value))[1,2])  #=65

File_Sequences_Danio <- Pigmentation_Sequences %>% filter(Species == "Danio_rerio") %>%
  mutate(seq_wo_stop = gsub('.{3}$', '', Sequence)) %>% dplyr::select(Gene, seq_wo_stop)



#Run100 simulations 



write.table(File_Sequences_Danio, file="File_Sequences_Danio.tsv", sep="\t", col.names = FALSE, row.names = FALSE)
system('sed -i "" "s/\"//g" File_Sequences_Danio.tsv')
system("for i in {1..100} ; do /Users/maxime/miniconda3/bin/python Neutral_evolution_for_mutpred.py -f File_Sequences_Danio.tsv -r 1.85056 -n 39 ; done")


#Run mutpred2 on every files


system("rm File_Sequences_Danio.tsv")
system("for i in * ; do run_mutpred2.sh -i $i -o $i.rslt ; done")
system("cat *.rslt > Simulations_Pigmentation_MP2_Results")


#Perform graphs with observed mutations and simulated mutations



system('grep -v "[1-9].*,.*,#.*" Simulations_Pigmentation_MP2_Results | cut -f1,2,3 -d "," | sed "s/---/,/g" | sed "s/^ID,/Simulation,ID,/g" | grep -v "Substitution" > Simulations_Pigmentation_MP2_Results_Table.tsv')
Simulations_Pigmentation_Scores <- read.table("Simulations_Pigmentation_MP2_Results_Table.tsv", header=FALSE, sep=",")
colnames(Simulations_Pigmentation_Scores) <- c("Simulation", "ID", "Substitution", "MutPred2.score")
Simulations_Pigmentation_Scores <- Simulations_Pigmentation_Scores %>% mutate(Species = "Simulation")

Simulations_Pigmentation_Scores$MutPred2.score <- as.numeric(as.character(Simulations_Pigmentation_Scores$MutPred2.score))
Simulations_Pigmentation_Scores$ID <- as.character(Simulations_Pigmentation_Scores$ID)
Simulations_Pigmentation_Scores$Substitution <- as.character(Simulations_Pigmentation_Scores$Substitution)


Simulations_Pigmentation_Scores_v1 <- Simulations_Pigmentation_Scores %>% dplyr::select(ID, Substitution, MutPred2.score, Species)
Pigmentation_mutpred2 <- bind_rows(Pigmentation_mutpred2, Simulations_Pigmentation_Scores_v1)



Pigmentation_mutpred2_new_sp <- Pigmentation_mutpred2 %>% filter(
  Species == "Danio_rerio" | Species == "Gadus_morhua" | Species == "Lamprologus_lethops" | Species == "Neolamprologus_brichardi" | Species == "Lamprologus_tigripic" | Species == "Percopsis_transmontana" | Species == "Typhlichthys_subterraneus" | Species == "Simulation"
)

Pigmentation_mutpred2_old_sp <- Pigmentation_mutpred2 %>% filter(
  Species == "Astyanax_CF" | Species == "Astyanax_SF" | Species == "Lucifuga_holguinensis" | Species == "Lamprogrammus_exutus" | Species == "Brotula_barbata" | Species == "Carapus_acus" | Species == "Pygocentrus_nattereri" | Species == "Lucifuga_dentata"
)


pmutpred_Pigmentation <- Simulations_Pigmentation_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  stat_ecdf(show.legend = FALSE, geom = "step", pad = FALSE) +
  stat_ecdf(data = Pigmentation_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "solid") +
  stat_ecdf(data = Pigmentation_mutpred2_old_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "dashed") +
  theme( panel.background = element_rect(fill = "white",
                                         colour = "black",
                                         size = 0.5, linetype = "solid"),
         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                         colour = "white"), 
         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                         colour = "white")) +
  scale_color_manual(values = c(rep("gray93", 100), "plum1", "deepskyblue", "darkolivegreen1", "darkolivegreen4", "darkviolet", "lawngreen", "khaki", "red", "darkorchid1", "sienna4", "sandybrown", "royalblue1", "paleturquoise1", "green", "black", "yellow")) +
  xlab("Mutpred2 score") +
  ylab("Cumulative frequency") 

pmutpred_pigmentation_2 <- Simulations_Pigmentation_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  stat_ecdf(show.legend = FALSE, geom = "step", pad = FALSE) +
  stat_ecdf(data = Pigmentation_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), geom = "step", pad = FALSE, show.legend = FALSE, linetype = "solid") +
  theme( panel.background = element_rect(fill = "white",
                                         colour = "black",
                                         size = 0.5, linetype = "solid"),
         panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                         colour = "white"), 
         panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                         colour = "white")) +
  scale_color_manual(values = c(rep("gray93", 100), "darkviolet", "lawngreen", "red", "darkorchid1", "royalblue1", "paleturquoise1", "black", "yellow")) +
  xlab("Mutpred2 score") +
  ylab("Cumulative frequency") 



ks.test(Pigmentation_mutpred2 %>% filter(Species == "Lamprologus_lethops") %>% pull(MutPred2.score),
        Pigmentation_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score))





Pigmentation_mutpred2 %>%
  group_by(Species) %>%
  summarise(value = n())


Simulations_Pigmentation_Scores %>%
  ggplot(., aes(x = MutPred2.score, colour=Simulation)) +
  geom_density(show.legend = FALSE) +
  geom_density(data = Pigmentation_mutpred2_new_sp, aes(x=MutPred2.score, colour=Species), show.legend = FALSE) +
  theme_minimal() +
  scale_color_manual(values = c(rep("gray", 100), "darkviolet", "lawngreen", "red", "darkorchid1", "royalblue1", "paleturquoise1", "black", "yellow")) 









#### Mutpred2 mixture results ---------------------------------
Neutral_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Simulation") %>% pull(MutPred2.score)
Danio_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Danio_rerio") %>% pull(MutPred2.score)
Lethops_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Lamprologus_lethops") %>% pull(MutPred2.score)
Subterraneus_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Typhlichthys_subterraneus") %>% pull(MutPred2.score)
Percopsis_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Percopsis_transmontana") %>% pull(MutPred2.score)
Brichardi_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Neolamprologus_brichardi") %>% pull(MutPred2.score)
Morhua_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Gadus_morhua") %>% pull(MutPred2.score)
Tigripic_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Lamprologus_tigripic") %>% pull(MutPred2.score)
astyanaxCF_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Astyanax_CF") %>% pull(MutPred2.score)
astyanaxSF_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Astyanax_SF") %>% pull(MutPred2.score)
holguinensis_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Lucifuga_holguinensis") %>% pull(MutPred2.score)
lamprogrammus_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Lamprogrammus_exutus") %>% pull(MutPred2.score)
brotula_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Brotula_barbata") %>% pull(MutPred2.score)
carapus_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Carapus_acus") %>% pull(MutPred2.score)
pygocentrus_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Pygocentrus_nattereri") %>% pull(MutPred2.score)
dentata_mutations <- Pigmentation_mutpred2 %>% filter(Species == "Lucifuga_dentata") %>% pull(MutPred2.score)



Pigmentation_admixture_simulations <- data.frame(NULL)
Pigmentation_admixture_lethops <- data.frame(NULL)
Pigmentation_admixture_subterraneus <- data.frame(NULL)
Pigmentation_admixture_percopsis <- data.frame(NULL)
Pigmentation_admixture_brichardi <- data.frame(NULL)
Pigmentation_admixture_morhua <- data.frame(NULL)
Pigmentation_admixture_tigripic <- data.frame(NULL)
Pigmentation_admixture_Astyanax_CF <- data.frame(NULL)
Pigmentation_admixture_Astyanax_SF <- data.frame(NULL)
Pigmentation_admixture_holguinensis <- data.frame(NULL)
Pigmentation_admixture_lamprogrammus <- data.frame(NULL)
Pigmentation_admixture_brotula <- data.frame(NULL)
Pigmentation_admixture_carapus <- data.frame(NULL)
Pigmentation_admixture_pygocentrus <- data.frame(NULL)
Pigmentation_admixture_dentata <- data.frame(NULL)
Pigmentation_admixture_danio <- data.frame(NULL)
iter=0
for(j in 1:1000){
  my_vector_test_lethops  <- c()
  my_vector_test_subterraneus  <- c()
  my_vector_test_percopsis  <- c() 
  my_vector_test_brichardi  <- c() 
  my_vector_test_morhua  <- c() 
  my_vector_test_tigripic  <- c()
  my_vector_test_Astyanax_CF  <- c()
  my_vector_test_Astyanax_SF  <- c()
  my_vector_test_holguinensis  <- c()
  my_vector_test_lamprogrammus  <- c()
  my_vector_test_brotula  <- c()
  my_vector_test_carapus  <- c()
  my_vector_test_pygocentrus  <- c()
  my_vector_test_dentata <- c()
  my_vector_test_danio <- c()
  my_vector_name <- c()
  for(i in 1:1000){
    my_mutations_final <- c(sample(Danio_mutations, 1000-i, replace = FALSE, prob = NULL), sample(Neutral_mutations, i, replace = FALSE, prob = NULL))
    my_vector_test_lethops <- c(my_vector_test_lethops, ks.test(Lethops_mutations, my_mutations_final)$p)
    my_vector_test_subterraneus <- c(my_vector_test_subterraneus, ks.test(Subterraneus_mutations, my_mutations_final)$p)
    my_vector_test_percopsis <- c(my_vector_test_percopsis, ks.test(Percopsis_mutations, my_mutations_final)$p)
    my_vector_test_brichardi <- c(my_vector_test_brichardi, ks.test(Brichardi_mutations, my_mutations_final)$p)
    my_vector_test_morhua <- c(my_vector_test_morhua, ks.test(Morhua_mutations, my_mutations_final)$p)
    my_vector_test_tigripic <- c(my_vector_test_tigripic, ks.test(Tigripic_mutations, my_mutations_final)$p)
    my_vector_test_Astyanax_CF <- c(my_vector_test_Astyanax_CF, ks.test(astyanaxCF_mutations, my_mutations_final)$p)
    my_vector_test_Astyanax_SF <- c(my_vector_test_Astyanax_SF, ks.test(astyanaxSF_mutations, my_mutations_final)$p)
    my_vector_test_holguinensis <- c(my_vector_test_holguinensis, ks.test(holguinensis_mutations, my_mutations_final)$p)
    my_vector_test_lamprogrammus <- c(my_vector_test_lamprogrammus, ks.test(lamprogrammus_mutations, my_mutations_final)$p)
    my_vector_test_brotula <- c(my_vector_test_brotula, ks.test(brotula_mutations, my_mutations_final)$p)
    my_vector_test_carapus <- c(my_vector_test_carapus, ks.test(carapus_mutations, my_mutations_final)$p)
    my_vector_test_pygocentrus <- c(my_vector_test_pygocentrus, ks.test(pygocentrus_mutations, my_mutations_final)$p)
    my_vector_test_dentata <- c(my_vector_test_dentata, ks.test(dentata_mutations, my_mutations_final)$p)
    my_vector_test_danio <- c(my_vector_test_danio, ks.test(Danio_mutations, my_mutations_final)$p)
    my_vector_name <- c(my_vector_name, as.character(j))
    
    
  }
  Pigmentation_admixture_lethops <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_lethops), as.data.frame(rep("lethops", 1000)))
  colnames(Pigmentation_admixture_lethops) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Pigmentation_admixture_subterraneus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_subterraneus), as.data.frame(rep("subterraneus", 1000)))
  colnames(Pigmentation_admixture_subterraneus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Pigmentation_admixture_percopsis <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_percopsis), as.data.frame(rep("percopsis", 1000)))
  colnames(Pigmentation_admixture_percopsis) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Pigmentation_admixture_brichardi <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_brichardi), as.data.frame(rep("brichardi", 1000)))
  colnames(Pigmentation_admixture_brichardi) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Pigmentation_admixture_morhua <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_morhua), as.data.frame(rep("morhua", 1000)))
  colnames(Pigmentation_admixture_morhua) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Pigmentation_admixture_tigripic <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_tigripic), as.data.frame(rep("tigripic", 1000)))
  colnames(Pigmentation_admixture_tigripic) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_Astyanax_CF <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_Astyanax_CF), as.data.frame(rep("Astyanax_CF", 1000)))
  colnames(Pigmentation_admixture_Astyanax_CF) <- c("Simulation_nb", "KS_test_p", "Species")
  
  Pigmentation_admixture_Astyanax_SF <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_Astyanax_SF), as.data.frame(rep("Astyanax_SF", 1000)))
  colnames(Pigmentation_admixture_Astyanax_SF) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_holguinensis <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_holguinensis), as.data.frame(rep("holguinensis", 1000)))
  colnames(Pigmentation_admixture_holguinensis) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_lamprogrammus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_lamprogrammus), as.data.frame(rep("lamprogrammus", 1000)))
  colnames(Pigmentation_admixture_lamprogrammus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_brotula <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_brotula), as.data.frame(rep("brotula", 1000)))
  colnames(Pigmentation_admixture_brotula) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_carapus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_carapus), as.data.frame(rep("carapus", 1000)))
  colnames(Pigmentation_admixture_carapus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_pygocentrus <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_pygocentrus), as.data.frame(rep("pygocentrus", 1000)))
  colnames(Pigmentation_admixture_pygocentrus) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_dentata <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_dentata), as.data.frame(rep("dentata", 1000)))
  colnames(Pigmentation_admixture_dentata) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_danio <- bind_cols(as.data.frame(my_vector_name), as.data.frame(my_vector_test_danio), as.data.frame(rep("danio", 1000)))
  colnames(Pigmentation_admixture_danio) <- c("Simulation_nb", "KS_test_p", "Species")
  
  
  Pigmentation_admixture_simulations <- bind_rows(Pigmentation_admixture_simulations,
                                                  Pigmentation_admixture_lethops,
                                                  Pigmentation_admixture_subterraneus, 
                                                  Pigmentation_admixture_percopsis,
                                                  Pigmentation_admixture_brichardi,
                                                  Pigmentation_admixture_morhua,
                                                  Pigmentation_admixture_tigripic,
                                                  Pigmentation_admixture_Astyanax_CF,
                                                  Pigmentation_admixture_Astyanax_SF,
                                                  Pigmentation_admixture_holguinensis,
                                                  Pigmentation_admixture_lamprogrammus,
                                                  Pigmentation_admixture_brotula,
                                                  Pigmentation_admixture_carapus,
                                                  Pigmentation_admixture_pygocentrus,
                                                  Pigmentation_admixture_dentata,
                                                  Pigmentation_admixture_danio
  )
  
  iter=iter+1
  print(iter)
  
}


Pigmentation_admixture_simulations <- Pigmentation_admixture_simulations %>% mutate(Neutral_number = rep((1:1000), 15000))


Pigmentation_admixture_simulations_old_species <- Pigmentation_admixture_simulations %>% filter(Species == "dentata"|
                                                                                                  Species == "pygocentrus"|
                                                                                                  Species == "carapus"|
                                                                                                  Species == "brotula"|
                                                                                                  Species == "lamprogrammus"|
                                                                                                  Species == "holguinensis"|
                                                                                                  Species == "Astyanax_SF"|
                                                                                                  Species == "Astyanax_CF")

Pigmentation_admixture_simulations_new_species <- Pigmentation_admixture_simulations %>% filter(Species == "brichardi"| 
                                                                                                  Species == "danio"|
                                                                                                  Species == "lethops"|
                                                                                                  Species == "morhua"|
                                                                                                  Species == "percopsis"|
                                                                                                  Species == "subterraneus"|
                                                                                                  Species == "tigripic")




admix_Pigmentation <- Pigmentation_admixture_simulations_new_species %>% 
  group_by(Species, Neutral_number) %>%
  summarise(moyenne = mean(KS_test_p)) %>%
  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
  geom_line(show.legend = FALSE) +
  geom_line(data = (Pigmentation_admixture_simulations_old_species %>% group_by(Species, Neutral_number) %>% summarise(moyenne = mean(KS_test_p))), aes(x=Neutral_number, y=moyenne, color=Species), show.legend = FALSE, linetype="dashed") +
  xlab("Neutral mutations proportion") +
  ylab("KS-test pvalue") +
  theme_minimal() +
  scale_x_continuous(labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(values = c("plum1", "deepskyblue","royalblue1", "darkolivegreen1", "darkolivegreen4", "darkviolet", "sienna4", "sandybrown", "khaki", "red", "lawngreen", "paleturquoise1", "green", "yellow", "darkorchid1" ))



admix_pigmentation_2 <- Pigmentation_admixture_simulations_new_species %>% 
  group_by(Species, Neutral_number) %>%
  summarise(moyenne = mean(KS_test_p)) %>%
  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
  geom_line(show.legend = FALSE) +
  xlab("Neutral mutations proportion") +
  ylab("KS-test pvalue") +
  theme_minimal() +
  scale_x_continuous(labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(values = c("royalblue1", "darkviolet", "red", "lawngreen", "paleturquoise1", "yellow", "darkorchid1" ))



#admix_Pigmentation <- Pigmentation_admixture_simulations %>% 
#  group_by(Species, Neutral_number) %>%
#  summarise(moyenne = mean(KS_test_p)) %>%
#  ggplot(., aes(x=Neutral_number, y=moyenne, color=Species)) +
#  geom_line(show.legend = FALSE) +
#  xlab("Neutral mutations proportion") +
#  ylab("KS-test pvalue") +
#  theme_minimal() +
#  scale_x_continuous(labels = c("0", "25", "50", "75", "100"))




Proportion_lethops <- head(Pigmentation_admixture_simulations %>%
                             filter(Species == "lethops") %>%
                             group_by(Neutral_number) %>%
                             summarise(moyenne = mean(KS_test_p)) %>%
                             arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000

Proportion_subterraneus <- head(Pigmentation_admixture_simulations %>%
                                  filter(Species == "subterraneus") %>%
                                  group_by(Neutral_number) %>%
                                  summarise(moyenne = mean(KS_test_p)) %>%
                                  arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000


rslt_propor <- c()
species_curr_list <- Pigmentation_admixture_simulations %>% pull(Species) %>% unique(.)
for (species_curr in species_curr_list){
  rslt_propor <- c(rslt_propor, head(Pigmentation_admixture_simulations %>%
                                       filter(Species == species_curr) %>%
                                       group_by(Neutral_number) %>%
                                       summarise(moyenne = mean(KS_test_p)) %>%
                                       arrange(desc(moyenne)), 1) %>% pull(Neutral_number)/1000)
}



cbind(species_curr_list, rslt_propor)

############### SECTION GLOBAL RESULTS ############### 


#### Random position of LoF along CDS sequences  ---------------------------------


LoFs_Summary <- bind_rows(Vision_LoFs, Circadian_LoFs, Pigmentation_LoFs)

LoFs_Summary_stop_fs <- LoFs_Summary %>% 
  filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops") %>%
  filter(LoF == "Deletion" | LoF == "Insertion" | LoF == "Stop_gain") %>%
  mutate(Perc_cds = Pos_CDS/Danio_rerio_Cds_length) %>%
  mutate(y_pos_graph = 1)

p0 <- LoFs_Summary_stop_fs %>% ggplot(., aes(y=y_pos_graph, x=Perc_cds, color=LoF)) +
  geom_hline(yintercept = 1) +
  geom_point(aes(shape = LoF, color=Gene_Set), size = 2.5) +
  scale_shape_manual(values = c(18, 17, 16)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 



LoFs_Summary_stop_fs %>% filter(Gene_Set == "Vision") %>% ggplot(., aes(y=y_pos_graph, x=Perc_cds, color=LoF)) +
  geom_hline(yintercept = 1) +
  geom_point(aes(color = LoF), size = 2.5) +
  scale_shape_manual(values = c(18, 17, 16)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


LoFs_Summary_stop_fs %>% filter(Gene_Set == "Circadian") %>% ggplot(., aes(y=y_pos_graph, x=Perc_cds, color=LoF)) +
  geom_hline(yintercept = 1) +
  geom_point(aes(color = LoF), size = 2.5) +
  scale_shape_manual(values = c(18, 17, 16)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


LoFs_Summary_stop_fs %>% filter(Gene_Set == "Pigmentation") %>% ggplot(., aes(y=y_pos_graph, x=Perc_cds, color=LoF)) +
  geom_hline(yintercept = 1) +
  geom_point(aes(color = LoF), size = 2.5) +
  scale_shape_manual(values = c(18, 17, 16)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

#Make simulations to test for random position of stop codons

stop_position <- LoFs_Summary_stop_fs %>% 
  filter(LoF == "Stop_gain") %>%
  arrange(Perc_cds) %>%
  mutate(Perc_cds_100 = Perc_cds * 100) %>%
  pull(Perc_cds_100)


Nes_stop <- c()
for(j in 1:100000){
  random_segmentation <- runif(length(stop_position),0,100)           #random segmentation of our genes X times (X = number of observed stop codons)
  random_segmentation_sorted <- sort(random_segmentation)            
  random_segmentation_sorted <- append(0,random_segmentation_sorted)
  random_segmentation_sorted <- append(random_segmentation_sorted,100)
  length_random_segments <- c()
  i <- 1
  while(i < length(random_segmentation_sorted)){
    length_random_segments <- append(length_random_segments, random_segmentation_sorted[i+1]-random_segmentation_sorted[i])
    i=i+1
  }
  Relative_length <- length_random_segments/sum(length_random_segments)    #Calculation of the simulated effective segments sizes
  Ne <- 1/sum(Relative_length^2)
  Nes_stop <- append(Nes_stop, Ne)                                                  #Regroupment of simulations results
}


stop_position <- append(0, stop_position)
stop_position <- append(stop_position,100)

length_real_segments=c()
i <- 1
while(i < length(stop_position)){
  length_real_segments = append(length_real_segments, stop_position[i+1]-stop_position[i])  
  i <- i+1
}

relative_length <- length_real_segments/sum(length_real_segments)
ne_stop <- 1/sum(relative_length^2)

p1 <- qplot(as.data.frame(Nes_stop)$Nes_stop, geom="histogram",  
            main = "Histogram for random segmentation (stop codons)", 
            xlab = "Effective segmentation size",  
            ylab= "Number of values",
            fill=I("white"), 
            col=I("black")) +
  theme_minimal() +
  geom_vline(xintercept = ne_stop, color="red", size=1.5)



#Make simulations to test for random position of frameshifts


fs_position <- LoFs_Summary_stop_fs %>% 
  filter(LoF == "Deletion" | LoF == "Insertion") %>%
  arrange(Perc_cds) %>%
  mutate(Perc_cds_100 = Perc_cds * 100) %>%
  pull(Perc_cds_100)


Nes_fs <- c()
for(j in 1:100000){
  random_segmentation <- runif(length(fs_position),0,100)           #random segmentation of our genes X times (X = number of observed fs codons)
  random_segmentation_sorted <- sort(random_segmentation)            
  random_segmentation_sorted <- append(0,random_segmentation_sorted)
  random_segmentation_sorted <- append(random_segmentation_sorted,100)
  length_random_segments <- c()
  i <- 1
  while(i < length(random_segmentation_sorted)){
    length_random_segments <- append(length_random_segments, random_segmentation_sorted[i+1]-random_segmentation_sorted[i])
    i=i+1
  }
  Relative_length <- length_random_segments/sum(length_random_segments)    #Calculation of the simulated effective segments sizes
  Ne <- 1/sum(Relative_length^2)
  Nes_fs <- append(Nes_fs, Ne)                                                  #Regroupment of simulations results
}


fs_position <- append(0, fs_position)
fs_position <- append(fs_position,100)

length_real_segments=c()
i <- 1
while(i < length(fs_position)){
  length_real_segments = append(length_real_segments, fs_position[i+1]-fs_position[i])  
  i <- i+1
}

relative_length <- length_real_segments/sum(length_real_segments)
ne_fs <- 1/sum(relative_length^2)

p2 <- qplot(as.data.frame(Nes_fs)$Nes_fs, geom="histogram",  
            main = "Histogram for random segmentation (frameshifts)", 
            xlab = "Effective segmentation size",  
            ylab= "Number of values",
            fill=I("white"), 
            col=I("black")) +
  theme_minimal() +
  geom_vline(xintercept = ne_fs, color="blue", size=1.5)






(p0 / (p1 | p2))




#### LoF distributions  ---------------------------------


labels_lof <- LoFs_Summary %>% 
  filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops") %>%
  group_by(LoF) %>%
  summarise(value = n()) %>%
  unite(n_string, c(value, LoF), sep=" ") %>%
  pull(n_string)


LoF_distrib <- LoFs_Summary %>% 
  filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops") %>%
  group_by(LoF) %>%
  summarise(value = n()) %>%
  ggplot(., aes(x="", y=value, fill=LoF)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)  +
  theme_minimal() +
  theme_void() +
  scale_fill_manual(name=paste("LoF Type", "(",nrow(LoFs_Summary %>% filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops")) ,"mutations )"),
                    values=c("gray", "black", "purple", "orange", "red", "yellow"),
                    labels=labels_lof) + 
  ggtitle(label = "LoF mutations distribution in cavefishes")


#### Frameshift distributions  ---------------------------------

Fs_distrib <- LoFs_Summary %>% filter(LoF == "Deletion" | LoF == "Insertion") %>%
  filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops") %>%
  group_by(LoF, Nb_Bases) %>%
  summarise(value=n()) %>%
  ggplot(., aes(x=Nb_Bases, y=value, fill=LoF)) +
  geom_bar(stat='identity', position='dodge') + theme_minimal() +
  xlab("Indel size") +
  ylab("Number") +
  ggtitle("Distribution of indel size in cavefishes") +
  #scale_fill_discrete(name="Indel type", labels=c("Deletions", "Insertions")) +
  scale_fill_manual(name="Indel type", labels=c("Deletions", "Insertions"), values=c("gray", "black")) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=18,face="bold")) +
  scale_x_continuous("Indel size", c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22))



LoF_distrib / Fs_distrib

#### dNdS delta between cf and sf  ---------------------------------



df2 <- data.frame(species=rep(c("T. subterraneus", "P. transmontana", "L. lethops", "L. tigripictilis"), each=3),
                  gene_set=rep(c("Vision", "Circadian", "Pigmentation"),4),
                  omega=c(0.139552, 0.0956292, 0.0943482, 0.131664, 0.0850416, 0.0855979, 0.277382, 0.137322, 0.146352, 0.223823, 0.0515404, 0.0891504))


ggplot(data=df2, aes(x=gene_set, y=omega, fill=species)) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  theme_minimal() +
  scale_fill_manual("Species", values=c("L. lethops" = "mediumpurple4", "L. tigripictilis" = "mediumpurple1", "P. transmontana" = "lightcoral", "T. subterraneus" = "indianred4"))


## per gene

ds_dN_merged_vision_v1  <- ds_dN_merged_vision %>% dplyr::select(Gene, omega_percopsis, omega_subterraneus) %>%
  filter(omega_percopsis != "NA") %>%
  filter(omega_subterraneus != "NA") %>%
  filter(Gene != "Pseudogenes_subterraneus_Conatenated_Alignlent.codeml") %>%
  filter(Gene != "Vision_Genes_Concatenated_Alignment.codeml") %>%
  filter(Gene != "tmt1b" & Gene != "pde6a") %>%
  mutate(delta_cf_sf = omega_subterraneus - omega_percopsis)


ds_dN_merged_vision_v1 %>%
  filter(delta_cf_sf > 0)

ds_dN_merged_vision_v1 %>%
  ggplot(., aes(y=delta_cf_sf)) +
  geom_boxplot() +
  ylim(-1, 1)





#### Mutpred2 curves  ---------------------------------

(pmutpred_vision | pmutpred_circadian) / (pmutpred_pigmentation + guide_area())


#### Mutpred2 admixtures  ---------------------------------

(admix_vision | admix_circ) / (admix_pigm + guide_area())




#### Plot dn/ds tree and distributions  ---------------------------------


#Complete tree

library(treeio)
library(ggtree)

setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Concatenated_Alignment/All_species")

codeml_tree <- read.codeml_mlc("Vision_codeml_CodonFreq2_Clean_Freeratio.codeml")
ggtree(codeml_tree, aes(color=dN_vs_dS), size=2, branch.length="branch.length")+
  scale_color_continuous(low="black", high="red") +
  theme(legend.position="bottom") + 
  geom_tiplab() + 
  xlim(0, 2)


codeml_tree <- read.codeml_mlc("Circadian_codeml_CodonFreq2_Clean_Freeratio.codeml")
ggtree(codeml_tree, aes(color=dN_vs_dS), size=2, branch.length="branch.length")+
  scale_color_continuous(low="black", high="red") +
  theme(legend.position="bottom") + 
  geom_tiplab() + 
  xlim(0, 2)


codeml_tree <- read.codeml_mlc("Pigmentation_codeml_CodonFreq2_Clean_Freeratio.codeml")
ggtree(codeml_tree, aes(color=dN_vs_dS), size=2, branch.length="branch.length")+
  scale_color_continuous(low="black", high="red") +
  theme(legend.position="bottom") + 
  geom_tiplab() + 
  xlim(0, 2)




# Distributions with studied species



species_list <- c("Carapus_acus", "Lamprogrammus", "Brotula", "Dentata", "Gibarensis","Lamprologus_tigripic","Lamprologus_lethops", "Neolamprologus_brichardi", "Percopsis_transmontana", "Typhlichthys_subterraneus","Gadus","Astyanax_CF", "Astyanax_SF","Pygocentrus", "Danio_rerio")
morph_list <- c("Surface", "Surface", "Surface", "Lucifuga dentata", "Lucifuga gibarensis", "Surface", "Lamprologus lethops", "Surface", "Surface", "Typhlichthys subterraneus","Surface", "Astyanax mexicanus Cave", "Surface", "Surface", "Surface")


omega_vision <- c(0.062, 0.135, 0.080, 0.376 ,0.198, 0.213, 0.368, 0.146, 0.120, 0.139, 0.076, 0.197, 0.272, 0.068, 0.053)
omega_circadian <- c(0.043, 0.063, 0.052, 0.272, 0.092, 0.151, 0.107, 0.041, 0.093, 0.082, 0.057, 0.178, 0.472, 0.031, 0.041)
omega_pigmentation <- c(0.048, 0.1, 0.071, 0.168, 0.090, 0.094, 0.145, 0.095, 0.082, 0.1, 0.054, 0.376, 0.120, 0.056, 0.050)


df_1 <- data.frame(species_list, omega_vision, morph_list) %>%
  mutate(Class = "Vision")
colnames(df_1) <- c("Species", "Omega", "Morph","Dataset")

df_2 <- data.frame(species_list, omega_circadian, morph_list) %>%
  mutate(Class = "Circadian")
colnames(df_2) <- c("Species", "Omega", "Morph","Dataset")


df_3 <- data.frame(species_list, omega_pigmentation, morph_list) %>%
  mutate(Class = "Pigmentation")
colnames(df_3) <- c("Species", "Omega", "Morph","Dataset")

ma_paml_df <- bind_rows(df_1, df_2, df_3)

ma_paml_df <- ma_paml_df %>% mutate(Omega_class = case_when(
  Omega > 0 & Omega <= 0.05 ~ "C1",
  Omega > 0.05 & Omega <= 0.10 ~ "C2",
  Omega > 0.10 & Omega <= 0.15 ~ "C3",
  Omega > 0.15 & Omega <= 0.20 ~ "C4",
  Omega > 0.20 & Omega <= 0.25 ~ "C5",
  Omega > 0.25 & Omega <= 0.30 ~ "C6",
  Omega > 0.30 & Omega <= 0.35 ~ "C7",
  Omega > 0.35 & Omega <= 0.40 ~ "C8",
  Omega > 0.40 & Omega <= 0.45 ~ "C9",
  Omega > 0.45 & Omega <= 0.50 ~ "C10"
  ))


ma_paml_df$Omega_class = factor(ma_paml_df$Omega_class, levels=c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11"))


Vision_paml <- ma_paml_df %>% filter(Dataset == "Vision") %>% 
  group_by(Omega_class, Morph, .drop=FALSE) %>% 
  summarise(value=length(Omega_class)) %>%
  ggplot(.,  aes(y=value, x=Omega_class, fill=Morph)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("maroon1", "red", "rosybrown1", "orange3", "forestgreen", "khaki1")) + 
  scale_x_discrete(breaks=c("C1","C2","C3","C4","C5","C6","C7","C8","C9", "C10"),
                   labels=c("[0,0.05]", "]0.05,0.10]", "]0.10,0.15]","]0.15,0.20]","]0.20,0.25]","]0.25,0.30]","]0.30,0.35]", 
                            "]0.35,0.40]", "]0.40,0.45]", "]0.45,0.50]")) +
  xlab("Omega") +
  ylab("Number of species") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8))



Circadian_paml <- ma_paml_df %>% filter(Dataset == "Circadian") %>% 
  group_by(Omega_class, Morph, .drop=FALSE) %>% 
  summarise(value=length(Omega_class)) %>%
  ggplot(.,  aes(y=value, x=Omega_class, fill=Morph)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("maroon1", "red", "rosybrown1", "orange3", "forestgreen", "khaki1")) + 
  scale_x_discrete(breaks=c("C1","C2","C3","C4","C5","C6","C7","C8","C9", "C10"),
                   labels=c("[0,0.05]", "]0.05,0.10]", "]0.10,0.15]","]0.15,0.20]","]0.20,0.25]","]0.25,0.30]","]0.30,0.35]", 
                            "]0.35,0.40]", "]0.40,0.45]", "]0.45,0.50]")) +
  xlab("Omega") +
  ylab("Number of species") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8))



Pigmentation_paml <- ma_paml_df %>% filter(Dataset == "Pigmentation") %>% 
  group_by(Omega_class, Morph, .drop=FALSE) %>% 
  summarise(value=length(Omega_class)) %>%
  ggplot(.,  aes(y=value, x=Omega_class, fill=Morph)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("maroon1", "red", "rosybrown1", "orange3", "forestgreen", "khaki1")) + 
  scale_x_discrete(breaks=c("C1","C2","C3","C4","C5","C6","C7","C8","C9", "C10"),
                   labels=c("[0,0.05]", "]0.05,0.10]", "]0.10,0.15]","]0.15,0.20]","]0.20,0.25]","]0.25,0.30]","]0.30,0.35]", 
                            "]0.35,0.40]", "]0.40,0.45]", "]0.45,0.50]")) +
  xlab("Omega") +
  ylab("Number of species") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8))



Vision_paml / Circadian_paml / Pigmentation_paml
# WIthout species with a ds < 0.01

Vision_paml <- ma_paml_df %>% filter(Dataset == "Vision") %>% 
  filter(!Species %in% c("Lamprologus_tigripic","Lamprologus_lethops", "Neolamprologus_brichardi", "Astyanax_CF", "Astyanax_SF")) %>%
  group_by(Omega_class, Morph, .drop=FALSE) %>% 
  summarise(value=length(Omega_class)) %>%
  ggplot(.,  aes(y=value, x=Omega_class, fill=Morph)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("maroon1", "red", "rosybrown1", "orange3", "forestgreen", "khaki1")) + 
  scale_x_discrete(breaks=c("C1","C2","C3","C4","C5","C6","C7","C8","C9", "C10"),
                   labels=c("[0,0.05]", "]0.05,0.10]", "]0.10,0.15]","]0.15,0.20]","]0.20,0.25]","]0.25,0.30]","]0.30,0.35]", 
                            "]0.35,0.40]", "]0.40,0.45]", "]0.45,0.50]")) +
  xlab("Omega") +
  ylab("Number of species") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8))



Circadian_paml <- ma_paml_df %>% filter(Dataset == "Circadian") %>% 
  filter(!Species %in% c("Lamprologus_lethops", "Neolamprologus_brichardi", "Astyanax_CF", "Astyanax_SF", "Gibarensis")) %>%
  group_by(Omega_class, Morph, .drop=FALSE) %>% 
  summarise(value=length(Omega_class)) %>%
  ggplot(.,  aes(y=value, x=Omega_class, fill=Morph)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("maroon1", "red", "rosybrown1", "orange3", "forestgreen", "khaki1")) + 
  scale_x_discrete(breaks=c("C1","C2","C3","C4","C5","C6","C7","C8","C9", "C10"),
                   labels=c("[0,0.05]", "]0.05,0.10]", "]0.10,0.15]","]0.15,0.20]","]0.20,0.25]","]0.25,0.30]","]0.30,0.35]", 
                            "]0.35,0.40]", "]0.40,0.45]", "]0.45,0.50]")) +
  xlab("Omega") +
  ylab("Number of species") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8))




Pigmentation_paml <- ma_paml_df %>% filter(Dataset == "Pigmentation") %>% 
  filter(!Species %in% c("Lamprologus_tigripic","Lamprologus_lethops", "Neolamprologus_brichardi", "Astyanax_CF", "Astyanax_SF", "Gibarensis")) %>%
  group_by(Omega_class, Morph, .drop=FALSE) %>% 
  summarise(value=length(Omega_class)) %>%
  ggplot(.,  aes(y=value, x=Omega_class, fill=Morph)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("maroon1", "red", "rosybrown1", "orange3", "forestgreen", "khaki1")) + 
  scale_x_discrete(breaks=c("C1","C2","C3","C4","C5","C6","C7","C8","C9", "C10"),
                   labels=c("[0,0.05]", "]0.05,0.10]", "]0.10,0.15]","]0.15,0.20]","]0.20,0.25]","]0.25,0.30]","]0.30,0.35]", 
                            "]0.35,0.40]", "]0.40,0.45]", "]0.45,0.50]")) +
  xlab("Omega") +
  ylab("Number of species") +
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8))

Vision_paml / Circadian_paml / Pigmentation_paml

############### SECTION TEST RER METHOD ############### 

## RER METHOD ##

setwd("~/Desktop/These_Supplemental_Analysis/RERconverge_data")


library(devtools)
library("RERconverge")


#Import trees

system('for i in ../../Alignements/Pigmentation_genes/*.aaml ; do x=`echo "$i"` ; y=`tail $i | grep "Danio"` ; echo "$x      $y" ; done > Pigmentation_Trees')
system("sed 's/.*Pigmentation_genes.//g' Pigmentation_Trees | sed 's/_NT.fasta.NT_ALN.fasta.align_noFS_AA.fasta.aaml//g' > Pigmentation_Trees_new")

system('for i in ../../Alignements/Vision_genes/*.aaml ; do x=`echo "$i"` ; y=`tail $i | grep "Danio"` ; echo "$x      $y" ; done > Vision_Trees')
system("sed 's/.*Vision_genes.//g' Vision_Trees | sed 's/_NT.fa.NT_ALN.fasta.align_noFS_AA.fasta.aaml//g' > Vision_Trees_new")


system('for i in ../../Alignements/Circadian_genes/*.aaml ; do x=`echo "$i"` ; y=`tail $i | grep "Danio"` ; echo "$x      $y" ; done > Circadian_Trees')
system("sed 's/.*Circadian_genes.//g' Circadian_Trees | sed 's/_NT.fa.NT_ALN.fasta.align_noFS_AA.fasta.aaml//g' > Circadian_Trees_new")


system("cat Pigmentation_Trees_new Circadian_Trees_new Vision_Trees_new BUSCO_Trees > Concatenated_Tab_Trees.txt")

mytreefiles = "Concatenated_Tab_Trees.txt"
myTrees=readTrees(mytreefiles)


#Estimating relative evolutionary rates


species_analyzed = c("Danio_rerio", "Gadus_morhua", "Lamprologus_lethops", "Lamprologus_tigripic", 
                     "Neolamprologus_brichardi", "Percopsis_transmontana", "Typhlichthys_subterraneus", 
                     "Astyanax_mexicanus_Cave", "Astyanax_mexicanus_Surface", "Pygocentrus_nattereri", 
                     "Carapus_acus", "Lamprogrammus_exutus", "Lucifuga_dentata", "Lucifuga_holguinensis", 
                     "Brotula_barbata")


mamRERw=RERconverge::getAllResiduals(myTrees,useSpecies=species_analyzed,
                                     transform = "sqrt", weighted = T, scale = T, min.sp = 10)


#Plot average and gene trees


outgroup_name <- c("Danio_rerio", "Astyanax_mexicanus_Surface", "Pygocentrus_nattereri", "Astyanax_mexicanus_Cave")
par(mfrow=c(1,2))
avgtree=plotTreeHighlightBranches(myTrees$masterTree, outgroup=outgroup_name,
                                  hlspecies=c("Lucifuga_dentata","Typhlichthys_subterraneus"), hlcols=c("red","red"),
                                  main="Average tree") #plot average tree
pde6b_tree=plotTreeHighlightBranches(myTrees$trees$pde6b, outgroup=outgroup_name,
                                     hlspecies=c("Lucifuga_dentata","Typhlichthys_subterraneus"), hlcols=c("red","red"),
                                     main="pde6b tree") #plot individual gene tree


## test with every cf

avgtree=plotTreeHighlightBranches(myTrees$masterTree, outgroup=outgroup_name,
                                  hlspecies=c("Lucifuga_dentata","Typhlichthys_subterraneus", "Lucifuga_holguinensis", "Lamprologus_lethops", "Astyanax_mexicanus_Cave"), hlcols=c("red","red"),
                                  main="Average tree") #plot average tree
pde6b_tree=plotTreeHighlightBranches(myTrees$trees$pde6b, outgroup=outgroup_name,
                                     hlspecies=c("Lucifuga_dentata","Typhlichthys_subterraneus", "Lucifuga_holguinensis", "Lamprologus_lethops", "Astyanax_mexicanus_Cave"), hlcols=c("red","red"),
                                     main="pde6b tree") #plot individual gene tree


#Plot RERs


foreground_branches = c("Lucifuga_dentata", "Typhlichthys_subterraneus")
foreground_branches = c("Lucifuga_dentata","Typhlichthys_subterraneus", "Lamprologus_lethops", "Astyanax_mexicanus_Cave")
foreground_branches = c("Lucifuga_dentata","Typhlichthys_subterraneus", "Lucifuga_holguinensis", "Lamprologus_lethops", "Astyanax_mexicanus_Cave")
foreground_branches = c("Lucifuga_dentata","Typhlichthys_subterraneus", "Lucifuga_holguinensis", "Lamprologus_lethops")


marineb2b = foreground2Tree(foreground_branches, myTrees, clade="terminal")


phenvMarine2=tree2Paths(marineb2b, myTrees)

corMarine=correlateWithBinaryPhenotype(mamRERw, phenvMarine2, weighted = "auto")


hist(corMarine$P, breaks=15, xlab="Kendall P-value",
     main="P-values for correlation between 200 genes and marine environment")


corMarine_name <- tibble::rownames_to_column(corMarine)


#Vision genes accelerated in cavefishes

corMarine_name %>% 
  filter(Rho > 0) %>%
  filter(P < 0.05)
  
  
#Circadian clock genes accelerated in cavefishes
  



#Pigmentation genes accelerated in cavefishes
  





###



corMarine_name %>% 
  filter(Rho > 0) %>%
  filter(P < 0.05) %>%
  filter(N == 27) %>%
  arrange(P) 




corMarine_name %>% 
  filter(Rho < 0) %>%
  filter(P < 0.05)






plotRers(mamRERw,"hdac1",phenv=phenvMarine2)
plotRers(mamRERw,"ednrba",phenv=phenvMarine2)
plotRers(mamRERw,"grk1a",phenv=phenvMarine2)
plotRers(mamRERw,"cryba1l1",phenv=phenvMarine2)
plotRers(mamRERw,"pde6b",phenv=phenvMarine2)
plotRers(mamRERw,"Opn7d",phenv=phenvMarine2)
plotRers(mamRERw,"per2_1",phenv=phenvMarine2)


Vision_mutpred2 %>% filter(ID == "cryba4")





#small graph test


lethops_vision <- Vision_Sequences %>% filter(Gene != "Not_Found") %>%
  filter(Species == "Lamprologus_lethops") %>% pull(Gene)

typhlichthys_vision <- Vision_Sequences %>% filter(Gene != "Not_Found") %>%
  filter(Species == "Typhlichthys_subterraneus") %>% pull(Gene)


lethops_circadian <- Circadian_Sequences %>% filter(Gene != "Not_Found") %>%
  filter(Species == "Lamprologus_lethops") %>% pull(Gene)

typhlichthys_circadian <- Circadian_Sequences %>% filter(Gene != "Not_Found") %>%
  filter(Species == "Typhlichthys_subterraneus") %>% pull(Gene)

lethops_pigmentation <- Pigmentation_Sequences %>% filter(Gene != "Not_Found") %>%
  filter(Species == "Lamprologus_lethops") %>% pull(Gene)

typhlichthys_pigmentation <- Pigmentation_Sequences %>% filter(Gene != "Not_Found") %>%
  filter(Species == "Typhlichthys_subterraneus") %>% pull(Gene)


nb_genes_common_vision <- length(intersect(lethops_vision, typhlichthys_vision))
nb_genes_common_circadian <- length(intersect(lethops_circadian, typhlichthys_circadian))
nb_genes_common_pigmentation <- length(intersect(lethops_pigmentation, typhlichthys_pigmentation))

accelerated_genes <- corMarine_name %>% 
  filter(Rho > 0) %>%
  filter(P < 0.05) %>%
  pull(rowname)

accel_pigm <- nrow(Pigmentation_Sequences %>% filter(Species == "Danio_rerio") %>%
                     filter(Gene %in% accelerated_genes))

accel_vision <- nrow(Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
                       filter(Gene %in% accelerated_genes))

accel_circadian <- nrow(Circadian_Sequences %>% filter(Species == "Danio_rerio") %>%
                          filter(Gene %in% accelerated_genes))


df2 <- data.frame(supp=rep(c("Normal", "Accelerated"), each=3),
                  gene_set=rep(c("Vision", "Circadian", "Pigmentation"),2),
                  len=c(nb_genes_common_vision, nb_genes_common_circadian, nb_genes_common_pigmentation, accel_vision, accel_circadian, accel_pigm))


ggplot(data=df2, aes(x=gene_set, y=len, fill=supp)) +
  geom_bar(stat="identity")






#### contrast-fel ###

setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Hyphy_relax_test")

vision_fel_rslt <- read.csv("datamonkey-table-vision-FEL.csv")

#also test with qvalue
test_fel <- vision_fel_rslt %>%
  filter(P.value_overall < 0.05) %>%
  mutate(Sens = if_else((beta_TEST > beta_background),
         "Increased_omega",
         "Decreased_omega"))
  
  
test_fel <- vision_fel_rslt %>%
  filter(Q.value_overall < 0.05) %>%
  mutate(Sens = if_else((beta_TEST > beta_background),
                        "Increased_omega",
                        "Decreased_omega"))

test_fel %>% group_by(Sens) %>% summarise(n())





##

############### TRASH ############### 

#Value of p-distances for Ldentata


## test busco ##
library(VennDiagram)

setwd("~/Desktop/These_Supplemental_Analysis/BUSCO_Subt_Perc")

percopsis_busco <- scan("list_percopsis", character(), quote = "")
subterr_busco <- scan("list_subt", character(), quote = "")


missing_common <- length(intersect(percopsis_busco,subterr_busco))
perco_missing <- 714
subt_missing <- 1004


grid.newpage()
d1 <- draw.pairwise.venn(area1 = perco_missing, 
                   area2 = subt_missing, 
                   cross.area = missing_common, 
                   category = c("    Percopsis", "Typhlichthys       "))


#vision genes missing in T.subterraneus and P. transmontana

list_missing_percopsis <- c("sws1", "sws2", "opn4m2", "opn4x2", "opn6b", "opn7c", "opn8a", "opn8c", "opn9", "para-2", "tmt2a", "tmt2b", "va1", "cryba2a", "crybb1l2", "crygm5", "gcap4", "gcap5", "gcap7", "pde6c", "rcv1b", "gnb3b", "gngt2a")
list_missing_subterr <- c("sws1", "sws2", "opn4m1", "opn4m2", "opn4x2", "opn6b", "opn7a", "opn7c", "opn8a", "opn9", "para-2", "tmt1a", "tmt2a", "tmt2b", "tmt3b", "va1", "cryba2a", "crybb1l2", "crygm5", "gcap4", "gcap5", "gcap7", "pde6ga", "pde6ha_1", "pde6ha_2", "pde6ha_3", "rcv1b", "rcv2b", "gnb3b", "gngt2a")


perco_vision_missing <- 23 #length(list_missing_percopsis)
subt_vision_missing <- 30 #length(list_missing_subterr)
missing_common_vision <- 21 #length(intersect(list_missing_percopsis, list_missing_subterr))

grid.newpage()
d2 <- draw.pairwise.venn(area1 = perco_vision_missing, 
                   area2 = subt_vision_missing, 
                   cross.area = missing_common_vision, 
                   category = c("    Percopsis", "Typhlichthys       "))



fisher.test(cbind(c(447, 557, 267),c(21,9,2)))

  


## add mp2 ##

setwd("~/Desktop/New_Cavefishes_Analysis/Pigmentation_Genes_MutPred/Simulations")
Simulations_Vision_Scores <- read.table("Simulations_Pigmentation_MP2_Results_Table.tsv", header=FALSE, sep=",")

setwd("~/Desktop/These_Supplemental_Analysis/Vision_Genes_Mutpred/Simulations/TEMP_MORE")
Simulations_Vision_Scores_to_add <- read.table("Simulations_Visions_MP2_Results_Table.tsv", header=FALSE, sep=",")




to_add_df <- sample_n(Simulations_Vision_Scores_to_add, size = 600, replace=FALSE)




vector_write <- rep(as.vector(Simulations_Vision_Scores %>% pull(V1) %>% unique(.)), 6)
write(vector_write, "TEST.txt")

write.table(to_add_df, file="NEW_DF", col.names=FALSE, row.names=FALSE, sep=",")



### dentata ##

setwd("~/Downloads/Info_SDV_Vulgarisation/Simulations_Subterraneus")

system("python2 Simulation_Analytique_dentata.py")



#Graphs

Binomial_all_neutral_dent <- read.table(file = "Matrice_loi_binomiale.csv", sep=",", header=FALSE)


#extract column of 19 pseudogenes (20 as first column is 0)


Binomial_all_neutral_dent_t <-  Binomial_all_neutral_dent %>%
  dplyr::select(V20) %>%
  mutate(generation = 1:n()) %>%
  mutate(Model = "Binomial")




Binomial_all_neutral_dent_t %>% ggplot(., aes(x=generation, y=(V20), color=Model)) +
  geom_line(show.legend = FALSE) +
  theme_minimal() +
  xlab(label = "Nombre de generations") +
  ylab(label = "Probabilite") +
  ggtitle(label = "Probabilite de trouver 19 pseudogenes") 



