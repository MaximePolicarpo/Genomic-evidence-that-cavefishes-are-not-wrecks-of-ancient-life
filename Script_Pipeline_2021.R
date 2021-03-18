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


#Data from Policarpo et al. 2020

Dentata_Table  <- read.table("Dentata_Table.tsv", header=FALSE, sep="\t")
Astyanax_Table  <- read.table("Astyanax_Table.tsv", header=FALSE, sep="\t")




#### Perform alignments and create phylogenies  ---------------------------------

#setwd("~/Desktop/These_Supplemental_Analysis/Alignements/Vision")

#For each gene, create a fasta file (.fa) with sequences of each species (The easiest way to do it from the sequences tables is to use tsv2fasta online or with python)
#.fa files are available in supplementary data of the paper !


#trimNonHomologousFragments with MACSE
system("for i in *.fa ; do java -jar macse_v2.03.jar -prog trimNonHomologousFragments -seq $i -min_homology_to_keep_seq 0 ; done")
system("rm *_AA.fa ; rm *_detail_NT.fa ; rm *.csv")

#alignSequences with MACSE
system("for i in *_NT.fa ; do java -jar macse_v2.03.jar -prog alignSequences -seq $i -out_NT $i.NT_ALN.fasta -out_AA $i.AA_ALN.fasta ; done")
system("rm *.AA_ALN.fasta")

#exportAlignment with MACSE and rename sequences to replace '!' integrated by MACSE  
system("for i in *.NT_ALN.fasta ; do java -jar macse_v2.03.jar -prog exportAlignment -align $i -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS - -out_NT $i.align_noFS_NT.fasta -out_AA $i.align_noFS_AA.fasta ; done")
system("for i in *_noFS_NT.fasta ; do sed -i '' 's/!/-/g' $i ; done")
system("for i in *_noFS_AA.fasta ; do sed -i '' 's/!/-/g' $i ; done")


#Make phylogenies with IQ-TREE with 1000 ultra fast bootstrap

system("for i in *.align_noFS_AA.fasta ; do iqtree -s $i -st AA -bb 1000 ; done")


#Now, you can visualize and root trees on iTOL

#### Gene number in each species  ---------------------------------

Vision_Sequences %>% filter(Gene != "Not_Found") %>%
  group_by(Species, Gene.Type) %>%
  summarise(n())


#### LoF stats  ---------------------------------


#LoF distribution

#Prepare the legend for the figure
labels_lof <- Vision_LoFs %>% 
  filter(Species == "Typhlichthys_subterraneus" | Species == "Lamprologus_lethops") %>%
  group_by(LoF) %>%
  summarise(value = n()) %>%
  unite(n_string, c(value, LoF), sep=" ") %>%
  pull(n_string)

#Draw figure
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
  ggtitle(label = "LoF mutations distribution in vision genes")


#Frameshift distribution

Vision_LoFs %>% filter(LoF == "Deletion" | LoF == "Insertion") %>%
  group_by(LoF, Nb_Bases) %>%
  summarise(value=n()) %>%
  ggplot(., aes(x=Nb_Bases, y=value, fill=LoF)) +
  geom_bar(stat='identity', position='dodge') + theme_minimal() +
  xlab("Indel size") +
  ylab("Number") +
  ggtitle("Distribution of indel size in vision genes") +
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

#Total number of genes in T.subt
Total_Number_Genes <- nrow(Vision_Sequences %>% 
                             filter(Species == "Typhlichthys_subterraneus") %>%
                             filter(Gene != "Not_Found")) 
#Number of pseudogenes
Total_Number_Pseudo <- nrow(Vision_Sequences %>% 
                              filter(Species == "Typhlichthys_subterraneus") %>%
                              filter(Gene != "Not_Found") %>% 
                              filter(Gene.Type == "Pseudogene"))



#Number of Loss-of-function mutations
Number_Lof <- nrow(Vision_LoFs %>% filter(Species == "Typhlichthys_subterraneus"))

#Number of splice site mutations
Nb_lof_intron <- nrow(Vision_LoFs %>% filter(Species == "Typhlichthys_subterraneus") %>% filter(LoF == "Splice"))

#Number of LoF other than splice site mutations
Nb_lof_cds <- Number_Lof - Nb_lof_intron



#Lets start the simulation

#62, 56, 50, 43, 37 and 31 genes correspond to 1, 0.9, 0.8, 0.7, 0.6 and 0.5 in gene proportion


Table_rslt_simu_final <- as.data.frame(NULL)
nb_simu <- 10000
for(NB_neutral_gene in c(62 ,56, 50, 43, 37, 31)){
  
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




#Draw the graphic with T. subt only

graph_LoF_sitrib <- Table_rslt_simu_final %>% 
  ggplot(., aes(x=LoF_nb_vector, y=LoF_count_vector, group=as.factor(Nb_neutral_gene_vector))) +
  geom_bar(stat="identity", position="dodge", aes(fill= as.factor(Nb_neutral_gene_vector))) +
  theme_minimal() +
  scale_fill_discrete(name = "Neutral gene number", labels = c("62" ,"56", "50", "43", "37", "31")) +
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

#Exact same thing that was done with T. subt

Total_Number_Genes_dent <- nrow(Dentata_Table)

Total_Number_Pseudo_dent <- nrow(Dentata_Table %>% 
                                   filter(V6 == "Pseudogene"))




Number_Lof_dent <- 22
Nb_lof_intron_dent <- 6


Nb_lof_cds_dent <- Number_Lof_dent - Nb_lof_intron_dent



## Simulated numbers of LoF ##

c(0.5, 0.6, 0.7, 0.8, 0.9, 1) * 76 # Get the number of gene of each proportion

Table_rslt_simu_final_dent <- as.data.frame(NULL)
nb_simu <- 10000
for(NB_neutral_gene in c(76 ,68, 61, 53, 46, 38)){
  
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


#Create the graphic for L.dent only

Table_rslt_simu_final_dent %>% 
  ggplot(., aes(x=LoF_nb_vector_dent, y=LoF_count_minus, group=as.factor(Nb_neutral_gene_vector_dent))) +
  geom_bar(stat="identity", position="dodge", aes(fill= as.factor(Nb_neutral_gene_vector_dent)), show.legend = FALSE) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line=element_blank(), panel.background=element_blank()) +
  geom_bar(data=Table_rslt_simu_final, aes(x=LoF_nb_vector, y=LoF_count_vector, group=as.factor(Nb_neutral_gene_vector), fill= as.factor(Nb_neutral_gene_vector)), stat="identity", position="dodge",  show.legend = FALSE) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#009E73","#0072B2", "#F0E442", "#D55E00", "#0072B2", "#D55E00")) +
  scale_x_discrete(name ="Number of LoF mutations per gene", limits=c("0","1","2", "3", "4", "5")) +
  geom_segment(aes(x = -0.45, y = obs_0_lof, xend = 0.5, yend = obs_0_lof, colour = "red"), size=0.1) +
  geom_segment(aes(x = 0.55, y = obs_1_lof, xend = 1.4, yend = obs_1_lof, colour = "red"), size=0.1) +
  geom_segment(aes(x = 1.55, y = obs_2_lof, xend = 2.4, yend = obs_2_lof, colour = "red"), size=0.1) +
  geom_segment(aes(x = 2.55, y = obs_3_lof, xend = 3.4, yend = obs_3_lof, colour = "red"), size=0.1) +
  geom_segment(aes(x = 3.55, y = obs_4_lof, xend = 4.4, yend = obs_4_lof, colour = "red"), size=0.1) +
  geom_segment(aes(x = 4.55, y = obs_5_lof, xend = 5.4, yend = obs_5_lof, colour = "red"), size=0.1) +
  geom_segment(aes(x = -0.45, y = -obs_0_lof_dent, xend = 0.5, yend = -obs_0_lof_dent, colour = "red"), size=0.1) +
  geom_segment(aes(x = 0.55, y = -obs_1_lof_dent, xend = 1.4, yend = -obs_1_lof_dent, colour = "red"), size=0.1) +
  geom_segment(aes(x = 1.55, y = -obs_2_lof_dent, xend = 2.4, yend = -obs_2_lof_dent, colour = "red"), size=0.1) +
  geom_segment(aes(x = 2.55, y = -obs_3_lof_dent, xend = 3.4, yend = -obs_3_lof_dent, colour = "red"), size=0.1) +
  geom_segment(aes(x = 3.55, y = -obs_4_lof_dent, xend = 4.4, yend = -obs_4_lof_dent, colour = "red"), size=0.1) +
  geom_segment(aes(x = 4.55, y = -obs_5_lof_dent, xend = 5.4, yend = -obs_5_lof_dent, colour = "red"), size=0.1) +
  coord_flip()





#### Intersection Dentata - Subterraneus  ---------------------------------

# Prepare data #

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

#List of common pseudogenes
pseudo_in_common <- intersect(Dentata_Table %>% filter(V6 == "Pseudogene") %>% pull(V2), Typhlichthys_subterraneus_table %>% filter(Gene.Type == "Pseudogene") %>% pull(Gene))

#List of common genes
genes_in_common <- intersect(dentata_genes, typhlichthys_genes)

#Only retain commmon genes from the table
subset_table <- Typhlichthys_subterraneus_table %>% filter(Gene %in% genes_in_common)

#Number of pseudogenes among common genes in T.subt
nb_subset_pseudo_subt <- nrow(Typhlichthys_subterraneus_table %>% filter(Gene %in% genes_in_common) %>% filter(Gene.Type == "Pseudogene"))

#Number of pseudogenes among common genes in L.dent
nb_subset_pseudo_dent <- nrow(Dentata_Table %>% filter(V2 %in% genes_in_common) %>% filter(V6 == "Pseudogene"))


#Lets start the simulation without taking into account gene length
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


pseudo_in_common <- intersect(Dentata_Table %>% filter(V6 == "Pseudogene") %>% pull(V2), Typhlichthys_subterraneus_table %>% filter(Gene.Type == "Pseudogene") %>% pull(Gene))
genes_in_common <- intersect(dentata_genes, typhlichthys_genes)


#some genes dont have the same annotation, so change them individually
genes_in_common[6] <- "opn6a"
genes_in_common[15] <- "rh1.1"
genes_in_common[16] <- "rh2.1"


#Compute mean length of pseudos and non-pseudo taking Danio rerio as reference

Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
  filter(Gene %in% pseudo_in_common) %>% pull(CDS.length) %>% mean()


Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
  filter(Gene %in% genes_in_common) %>% pull(CDS.length) %>% mean()


subset_table_V2 <- Vision_Sequences %>% filter(Species == "Danio_rerio") %>%
  filter(Gene %in% genes_in_common)




#Lets start the simulations that take genes lengths into account

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




# Now, lets take the number of intron into account


subset_table_V2 <- subset_table_V2 %>%
  mutate(Intron_nb = Exon.Count - 1)

#Normalize cds length and intron number by the LoF proba rate in exon or in intron 


subset_table_V2 <- subset_table_V2 %>%
  mutate(exon_weight = (proba_stop_gain + proba_frameshift + proba_stop_loss + proba_start_loss) * CDS.length) %>%
  mutate(intron_weight = proba_splice * Intron_nb) %>%
  mutate(Proba_mut = intron_weight+exon_weight)


#Lets start the simulations
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


#Lets make a graphic with the 3 simulations results 

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





#### Datation of cavefishes ---------------------------------



#Compute the number of pseudogene per species
Number_pseudogenes_Tsubt <- nrow(Vision_Sequences %>% filter(Species == "Typhlichthys_subterraneus") %>% filter(Gene != "Not_Found") %>% filter(Gene.Type == "Pseudogene"))
Number_pseudogenes_Lethops <- nrow(Vision_Sequences %>% filter(Species == "Lamprologus_lethops") %>% filter(Gene != "Not_Found") %>% filter(Gene.Type == "Pseudogene"))
Number_pseudogenes_Dentata <- nrow(Dentata_Table %>% filter(V6 == "Pseudogene"))
Number_pseudogenes_Astyanax <- 1 #pde6b 

#Compute the number of genes per species
Number_genes_Tsubt <- nrow(Vision_Sequences %>% filter(Species == "Typhlichthys_subterraneus") %>% filter(Gene != "Not_Found"))
Number_genes_Lethops <- nrow(Vision_Sequences %>% filter(Species == "Lamprologus_lethops") %>% filter(Gene != "Not_Found"))
Number_genes_Dentata <- nrow(Dentata_Table)
Number_genes_Astyanax <- nrow(Astyanax_Table)


#Compute the mean length of genes per species
Mean_gene_length_Tsubt <- Vision_Sequences %>% filter(Species == "Typhlichthys_subterraneus") %>% filter(Gene != "Not_Found") %>% pull(CDS.length) %>% mean(.)
Mean_gene_length_Lethops <- Vision_Sequences %>% filter(Species == "Lamprologus_lethops") %>% filter(Gene != "Not_Found") %>% pull(CDS.length) %>% mean(.)
Mean_gene_length_Dentata <- Dentata_Table %>% pull(V5) %>% mean(.)
Mean_gene_length_Astyanax <- Astyanax_Table %>% pull(V5) %>% mean(.)




#Lets make a list with every genes that are pseudogenes in atleast 1 cave species
Pseudogene_list <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn4m3", "opn4x1","opn6b", "opn7a", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "tmt3b", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1", "gnb3b", "gc2", "gc3", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")

#Lets make a list with every genes that are pseudogenes in atleast 1 cave species but remove those unique to T.subterraneus
Pseudogene_list_wo_Tsubt_uniq <- c("opn4m1", "opn4m2", "opn4m3", "opn6b", "opn4x1","opn7a", "opn7b", "opn7c", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt2a", "tmt3a", "tmt3b","cryaa", "crybb1", "crybgx", "crygm5", "crygn2", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh1", "gnb3b", "gc2", "gc3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")

#Lets make a list with every genes that are pseudogenes in atleast 1 cave species but remove those unique to L.lethops
Pseudogene_list_wo_Llethops_uniq <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn4m3","opn6b", "opn7a", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "tmt3b", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1", "gnb3b", "gc2", "gc3", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")

#Lets make a list with every genes that are pseudogenes in atleast 1 cave species but remove those unique to L.dentata
Pseudogene_list_wo_Ldentata_uniq <- c("gnat2", "opn6a_1", "opn6a_2", "opn4x1","opn4m1", "opn4m2", "opn6b", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1", "gc2", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b")

#Lets make a list with every genes that are pseudogenes in atleast 1 cave species but remove those unique to A.mexicanus
Pseudogene_list_wo_Amexicanus_uniq <- c("gnat2", "opn6a_1", "opn6a_2", "opn4m1", "opn4m2", "opn4m3", "opn6b", "opn7a", "opn7b", "opn7c", "opn7d", "parapinopsin-1", "parietopsin", "rgr1", "tmt1a", "tmt1b", "tmt2a", "tmt3a", "tmt3b", "rpe65b", "cryaa", "cryba1l1", "cryba4", "crybb1", "crybgx", "crygm5", "crygn2", "rcv2a", "saga", "arr3a", "pde6ga", "pde6hb", "pde6b", "pde6c", "rh2", "lws1", "rh1.1","rh1.2", "gnb3b", "gc2", "gc3", "gcap3", "gcap1", "gcap2", "grk1b", "grk7a", "grk7b", "opn4x1")



#Values for datations with only the pseudogene subset in T. subterraneus

Nb_gene_subset_Tsubt <- nrow(Typhlichthys_subterraneus_table %>% filter(Gene != "Not_Found") %>% filter(Gene %in% Pseudogene_list_wo_Tsubt_uniq))
Nb_pseudo_subset_Tsubt <- nrow(Typhlichthys_subterraneus_table %>% filter(Gene != "Not_Found") %>% filter(Gene %in% Pseudogene_list_wo_Tsubt_uniq) %>% filter(Gene.Type == "Pseudogene"))
Moy_length_subset_Tsubt <- Typhlichthys_subterraneus_table %>% filter(Gene != "Not_Found") %>% filter(Gene %in% Pseudogene_list_wo_Tsubt_uniq) %>% pull(CDS.length) %>% mean(.)


#Values for datations with only the pseudogene subset in L. lethops

Nb_gene_subset_Lethops <- nrow(Lamprologus_lethops_table %>% filter(Gene != "Not_Found") %>% filter(Gene %in% Pseudogene_list_wo_Llethops_uniq))
Nb_pseudo_subset_Lethops <- nrow(Lamprologus_lethops_table %>% filter(Gene != "Not_Found") %>% filter(Gene %in% Pseudogene_list_wo_Llethops_uniq) %>% filter(Gene.Type == "Pseudogene"))
Moy_length_subset_Lethops <- Lamprologus_lethops_table %>% filter(Gene != "Not_Found") %>% filter(Gene %in% Pseudogene_list_wo_Llethops_uniq) %>% pull(CDS.length) %>% mean(.)


#Values for datations with only the pseudogene subset in L. dentata


Nb_gene_subset_Ldentata <- nrow(Dentata_Table %>% filter(V2 %in% Pseudogene_list_wo_Ldentata_uniq))
Nb_pseudo_subset_Ldentata <- nrow(Dentata_Table %>% filter(V2 %in% Pseudogene_list_wo_Ldentata_uniq) %>% filter(V6 == "Pseudogene"))
Moy_length_subset_Ldentat <- Dentata_Table %>% filter(V2 %in% Pseudogene_list_wo_Ldentata_uniq) %>% pull(V5) %>% mean(.)


#Values for datations with only the pseudogene subset in A. mexicanus


Nb_gene_subset_Amex <- nrow(Astyanax_Table %>% filter(V2 %in% Pseudogene_list_wo_Amexicanus_uniq))
Nb_pseudo_subset_Amex <- nrow(Astyanax_Table %>% filter(V2 %in% Pseudogene_list_wo_Amexicanus_uniq) %>% filter(V6 == "Pseudogene"))
Moy_length_subset_Amex <- Astyanax_Table %>% filter(V2 %in% Pseudogene_list_wo_Amexicanus_uniq) %>% pull(V5) %>% mean(.)




#Now, lets put these values in simulations script and run them with python : 

setwd("~/Desktop/These_Supplemental_Analysis/Datations/Mars_2021_datation/")
system("for i in *.py ; do python2 $i ; done")




#Now, lets import result tables

Binomial_all_neutral_subt <- read.table(file = "Matrice_loi_binomiale_subt.csv", sep=",", header=FALSE)
Binomial_subset_neutral_subt <- read.table(file = "Matrice_loi_binomiale_subset_subt.csv", sep=",", header=FALSE)

Binomial_all_neutral_lethops <- read.table(file = "Matrice_loi_binomiale_lethops.csv", sep=",", header=FALSE)
Binomial_subset_neutral_lethops <- read.table(file = "Matrice_loi_binomiale_lethops_subset.csv", sep=",", header=FALSE)

Binomial_all_neutral_dent <- read.table(file = "Matrice_loi_binomiale_dent.csv", sep=",", header=FALSE)
Binomial_subset_neutral_dent <- read.table(file = "Matrice_loi_binomiale_subset_dent.csv", sep=",", header=FALSE)

Binomial_all_neutral_asty <- read.table(file = "Matrice_loi_binomiale_asty.csv", sep=",", header=FALSE)
Binomial_subset_neutral_asty <- read.table(file = "Matrice_loi_binomiale_subset_asty.csv", sep=",", header=FALSE)


#Now lets make the graph 


#Start by extracting columns corresponding to the number of pseudogenes in each species (N+1 as the first column correspond to 0 pseudogenes)

All_neutral_lethops <-  Binomial_all_neutral_lethops %>%
  dplyr::select(V4) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(All_neutral_lethops))
All_neutral_lethops <- cbind(generation=Ngeneration, All_neutral_lethops) 

All_neutral_subt <-  Binomial_all_neutral_subt %>%
  dplyr::select(V28) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(All_neutral_subt))
All_neutral_subt <- cbind(generation=Ngeneration, All_neutral_subt) 

All_neutral_dent <-  Binomial_all_neutral_dent %>%
  dplyr::select(V20) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(All_neutral_dent))
All_neutral_dent <- cbind(generation=Ngeneration, All_neutral_dent) 

All_neutral_asty <-  Binomial_all_neutral_asty %>%
  dplyr::select(V2) %>%
  mutate(Model = "Binomial_all")
Ngeneration <- as.numeric(rownames(All_neutral_asty))
All_neutral_asty <- cbind(generation=Ngeneration, All_neutral_asty) 


Subset_neutral_lethops <-  Binomial_subset_neutral_lethops %>%
  dplyr::select(V3) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Subset_neutral_lethops))
Subset_neutral_lethops <- cbind(generation=Ngeneration, Subset_neutral_lethops) 

Subset_neutral_subt <-  Binomial_subset_neutral_subt %>%
  dplyr::select(V14) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Subset_neutral_subt))
Subset_neutral_subt <- cbind(generation=Ngeneration, Subset_neutral_subt) 


Subset_neutral_dent <-  Binomial_subset_neutral_dent %>%
  dplyr::select(V14) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Subset_neutral_dent))
Subset_neutral_dent <- cbind(generation=Ngeneration, Subset_neutral_dent) 


Subset_neutral_asty <-  Binomial_subset_neutral_asty %>%
  dplyr::select(V2) %>%
  mutate(Model = "Binomial_subset")
Ngeneration <- as.numeric(rownames(Subset_neutral_asty))
Subset_neutral_asty <- cbind(generation=Ngeneration, Subset_neutral_asty) 



Subset_neutral_lethops <- Subset_neutral_lethops %>% mutate(Species = "Lamprologus lethops")
All_neutral_lethops <- All_neutral_lethops %>% mutate(Species = "Lamprologus lethops")
Subset_neutral_subt <- Subset_neutral_subt %>% mutate(Species = "Typhlichthys subterraneus")
All_neutral_subt <- All_neutral_subt %>% mutate(Species = "Typhlichthys subterraneus")
Subset_neutral_dent <- Subset_neutral_dent %>% mutate(Species = "Lucifuga dentata")
All_neutral_dent <- All_neutral_dent %>% mutate(Species = "Lucifuga dentata")
Subset_neutral_asty <- Subset_neutral_asty %>% mutate(Species = "Astyanax mexicanus")
All_neutral_asty <- All_neutral_asty %>% mutate(Species = "Astyanax mexicanus")



colnames(Subset_neutral_lethops) <- c("generation", "Probability", "Model", "Species")
colnames(All_neutral_lethops) <- c("generation", "Probability", "Model", "Species")
colnames(Subset_neutral_subt) <- c("generation", "Probability", "Model", "Species")
colnames(All_neutral_subt) <- c("generation", "Probability", "Model", "Species")
colnames(Subset_neutral_dent) <- c("generation", "Probability", "Model", "Species")
colnames(All_neutral_dent) <- c("generation", "Probability", "Model", "Species")
colnames(Subset_neutral_asty) <- c("generation", "Probability", "Model", "Species")
colnames(All_neutral_asty) <- c("generation", "Probability", "Model", "Species")

#Graphic creation

All_neutral_subt %>% ggplot(., aes(x=generation, y=(Probability))) +
  geom_line(linetype = "solid", colour="red")  +
  geom_line(data=Subset_neutral_subt, aes(x=generation, y=(Probability)), colour="red",  linetype = "longdash") +
  geom_line(data=All_neutral_lethops, aes(x=generation, y=(Probability)), linetype = "solid", colour="blue") +
  geom_line(data=Subset_neutral_lethops, aes(x=generation, y=(Probability)), colour="blue",  linetype = "longdash") +
  geom_line(data=All_neutral_dent, aes(x=generation, y=(Probability)), linetype = "solid", colour="black") +
  geom_line(data=Subset_neutral_dent, aes(x=generation, y=(Probability)), colour="black",  linetype = "longdash") +
  geom_line(data=All_neutral_asty, aes(x=generation, y=(Probability)), linetype = "solid", colour="gray") +
  geom_line(data=Subset_neutral_asty, aes(x=generation, y=(Probability)), colour="gray",  linetype = "longdash") +
  xlim(0, 2e+5) +
  #xlim(0, 1e+5) +
  ylim(0.001, 0.4) +
  theme_tufte()





#Get the 0.05 confidence intervals for all datations



#Get the maximum generation and it's probability 

head(All_neutral_subt %>% arrange(desc(Probability)), 1)
head(All_neutral_lethops %>% arrange(desc(Probability)), 1)
head(All_neutral_dent %>% arrange(desc(Probability)), 1)
head(All_neutral_asty %>% arrange(desc(Probability)), 1)

head(Subset_neutral_subt %>% arrange(desc(Probability)), 1)
head(Subset_neutral_lethops %>% arrange(desc(Probability)), 1)
head(Subset_neutral_dent %>% arrange(desc(Probability)), 1)
head(Subset_neutral_asty %>% arrange(desc(Probability)), 1)


#Get the 0.05 confidence intervals for all datations

head(All_neutral_subt %>% filter(Probability >= 0.05), 1)
tail(tail(All_neutral_subt %>% filter(Probability >= 0.05)), 1)
head(Subset_neutral_subt %>% filter(Probability >= 0.05), 1)
tail(Subset_neutral_subt %>% filter(Probability >= 0.05), 1)


head(All_neutral_lethops %>% filter(Probability >= 0.05), 1)
tail(tail(All_neutral_lethops %>% filter(Probability >= 0.05)), 1)
head(Subset_neutral_lethops %>% filter(Probability >= 0.05), 1)
tail(Subset_neutral_lethops %>% filter(Probability >= 0.05), 1)


head(All_neutral_dent %>% filter(Probability >= 0.05), 1)
tail(tail(All_neutral_dent %>% filter(Probability >= 0.05)), 1)
head(Subset_neutral_dent %>% filter(Probability >= 0.05), 1)
tail(Subset_neutral_dent %>% filter(Probability >= 0.05), 1)


head(All_neutral_asty %>% filter(Probability >= 0.05), 1)
tail(tail(All_neutral_asty %>% filter(Probability >= 0.05)), 1)
head(Subset_neutral_asty %>% filter(Probability >= 0.05), 1)
tail(Subset_neutral_asty %>% filter(Probability >= 0.05), 1)






############### SECTION CIRCADIAN ############### 


#### Import datas  ---------------------------------


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



#### Perform alignments   ---------------------------------

#Perform the same commands performed for vision genes




############### SECTION PIGMENTATION ############### 

#### Import datas  ---------------------------------


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


#### Perform alignments  ---------------------------------

#Same commands than those performed for vision and circadian clock genes







############### SECTION GLOBAL RESULTS ############### 

#### Lets see if LoF have random positions  along CDS sequences  ---------------------------------


#Merge LoFs from the three datasets

LoFs_Summary <- bind_rows(Vision_LoFs, Circadian_LoFs, Pigmentation_LoFs)

#Lets plot all these LoF (non splice) on a line in percent of CDS length

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


#Do the same but for vision genes only

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

#Do the same but for circadian clock genes only

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


#Do the same but for pigmentation  genes only

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

p1 | p2


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


