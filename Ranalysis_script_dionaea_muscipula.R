# This program analyzes the relative expression and number of transmembrane helices in 
# genetic transcripts. It requires the Trinity.fasta.transdecoder.pep file and the
# A1-9-matrix-genes-trimm-PwE.counts.matrix

#set the correct working directory and load necessary modules.
setwd("/home/nick/Transcriptome_data/")
library(Biostrings)

# read in peptide sequences and gene IDs from bulk sequences. Note these have selected ORFs with transdecoder
ORF_pep=readAAStringSet("A1-9-group-analysis/Trinity.fasta.transdecoder.pep")
#ORF_dna=readDNAStringSet("A1-9-group-analysis/Trinity.fasta.transdecoder.cds")
#ORF_dna[which(ORF_pep_ID==12161014)]
ORF_pep_ID=as.numeric(gsub("[^0-9]", "", substring(names(ORF_pep), 1, 24)))

# read in expression counts associated with different protein IDs
count_exp=read.table("A1-9-matrix-genes-trimm-PwE.counts.matrix")
count_ID=as.numeric(gsub("[^0-9]", "", row.names(count_exp)))

#declare values for user prompts.
comparison_type='temp'
helices_cutoff=0
rel_exp_cutoff=0
ratio_exp_cutoff=0
min_exp_level=0
#these loops allow the user to specify how many transmembrane helices and what expression levels to use as cutoffs.
print("This script will analyze Venus Fly-trap expression DNA")
while(!(comparison_type == "s" | comparison_type == "d")) {
  comparison_type=as.character(readline(prompt="Press s for subtraction comparison or d for division: "))
}

if (comparison_type=='s') {
  subtracted=TRUE
  divided=FALSE
  
  while(helices_cutoff<=0){
    helices_cutoff=as.integer(readline(prompt="Enter helices cutoff: "))
  }
  while(rel_exp_cutoff<=0){
    rel_exp_cutoff=as.integer(readline(prompt="Enter expression difference cutoff: "))
  }
  
} else if (comparison_type=='d') {
  divided=TRUE
  subtracted=FALSE
  
  while(helices_cutoff<=0){
    helices_cutoff=as.integer(readline(prompt="Enter helices cutoff: "))
  }
  while(ratio_exp_cutoff<=0){
    ratio_exp_cutoff=as.integer(readline(prompt="Enter ratio of expression cutoff: "))
  }
  while(min_exp_level<=0){
    min_exp_level=as.integer(readline(prompt="Enter minimum expression level: "))
  }
}

#combine expression values based on where the samples were harvested.
av_sensory_exp=rowSums(count_exp[,c(1,2,3,4)])/4
av_center_exp=rowSums(count_exp[,c(5,6)])/2
av_neg_exp=rowSums(count_exp[,c(7,8,9)])/3
#center cells are excluded
total_exp=rowSums(count_exp[,c(1,2,3,4,7,8,9)])

#calculate the subtracted or divided values, between conditions, to be used.
sensory_rel_exp=as.numeric(av_sensory_exp-av_neg_exp)
sensory_ratio_exp=as.numeric(av_sensory_exp/av_neg_exp)

#if the user chose to subtract, compare the subtracted values to the cutoff

if (subtracted) {
  exp_peptide=count_exp[which(sensory_rel_exp>rel_exp_cutoff),]
  exp_peptide_sensory=av_sensory_exp[which(sensory_rel_exp>rel_exp_cutoff)]
  exp_peptide_center=av_center_exp[which(sensory_rel_exp>rel_exp_cutoff)]
  exp_peptide_neg=av_neg_exp[which(sensory_rel_exp>rel_exp_cutoff)]
  exp_peptide_names=count_ID[which(sensory_rel_exp>rel_exp_cutoff)]
}

#if the user chose to divide, compare the divided values to the cutoff, and make sure that the total expression levels 
#are above the minimum.

if (divided){
  exp_peptide=count_exp[which(sensory_ratio_exp>ratio_exp_cutoff & total_exp>min_exp_level),]
  exp_peptide_sensory=av_sensory_exp[which(sensory_ratio_exp>ratio_exp_cutoff & total_exp>min_exp_level)]
  exp_peptide_center=av_center_exp[which(sensory_ratio_exp>ratio_exp_cutoff & total_exp>min_exp_level)]
  exp_peptide_neg=av_neg_exp[which(sensory_ratio_exp>ratio_exp_cutoff & total_exp>min_exp_level)]
  exp_peptide_names=count_ID[which(sensory_ratio_exp>ratio_exp_cutoff & total_exp>min_exp_level)]
}

#Find which of the sequences that made the cut are ORFs
exp_ORFs=ORF_pep[which(ORF_pep_ID %in% exp_peptide_names)]
exp_ORFs_names=ORF_pep_ID[which(ORF_pep_ID %in% exp_peptide_names)]

#combine and remove duplicated ORF sequences
exp_ORFs_nodup=exp_ORFs[which(!duplicated(exp_ORFs))]
exp_ORFs_names_nodup=exp_ORFs_names[which(!duplicated(exp_ORFs))]

#add the expression levels to the numeric values we're recorded (there may be a better way to do this)
exp_ORFs_sensory=vector('numeric')
exp_ORFs_center=vector('numeric')
exp_ORFs_neg=vector('numeric')
for (i in 1:length(exp_ORFs_names_nodup)){
  exp_ORFs_sensory=c(exp_ORFs_sensory, exp_peptide_sensory[(exp_peptide_names %in% exp_ORFs_names_nodup[i])])
  exp_ORFs_center=c(exp_ORFs_center, exp_peptide_center[(exp_peptide_names %in% exp_ORFs_names_nodup[i])])
  exp_ORFs_neg=c(exp_ORFs_neg, exp_peptide_neg[(exp_peptide_names %in% exp_ORFs_names_nodup[i])])
}

#create a file in fasta format for tmhmm to process, and then take the output and delete the 
#file. This uses the file tmhmm_single_peptide.txt which calls the tmhmm script
cat(paste("\n>\n", exp_ORFs_nodup), file="exp_ORFs_only_temp.txt")
system("./tmhmm_single_file")
helix_prediction=read.table("tmhmm.pred.temp.out", colClasses="character")
num_helices=as.numeric(gsub("[^0-9]", "", helix_prediction[,4]))
system("rm tmhmm.pred.temp.out")
system("rm exp_ORFs_only_temp.txt")

#compare calculated helix values with user specified cutoff
exp_ORFs_helices=exp_ORFs_nodup[num_helices>=helices_cutoff]
exp_ORFs_helices_names=exp_ORFs_names_nodup[num_helices>=helices_cutoff]
exp_ORFs_sensory_helices=exp_ORFs_sensory[num_helices>=helices_cutoff]
exp_ORFs_center_helices=exp_ORFs_center[num_helices>=helices_cutoff]
exp_ORFs_neg_helices=exp_ORFs_neg[num_helices>=helices_cutoff]
num_helices_cut=num_helices[num_helices>=helices_cutoff]

#create output data vector
output_data=paste("> Gene ID:", paste(substring(exp_ORFs_helices_names, 1, 22)), "   Sensory Expression: ",  
                  as.integer(exp_ORFs_sensory_helices), "   Center Expression: ", as.integer(exp_ORFs_center_helices), 
                  "   Negative Control Expression: ", as.integer(exp_ORFs_neg_helices), "   Number of TM Domains: ",
                  as.character(num_helices_cut), "\n", as.character(exp_ORFs_helices), sep='')
#adjust file name for user preference
if (subtracted){
  file_name=paste(c("predicted_proteins_", helices_cutoff, "_helices_", rel_exp_cutoff, "_subtracted_expression.fasta"), collapse="")
}
if (divided){
  file_name=paste(c("predicted_proteins_", helices_cutoff, "_helices_", min_exp_level, "_min_expression_", ratio_exp_cutoff, 
                    "_ratio_expression.fasta"), collapse="")
}

output_file=paste(length(exp_ORFs_helices_names), " of transcipts identified, with ", 
                  length(exp_ORFs_helices_names[which(!duplicated(exp_ORFs_helices_names))]), " unique genes.")
print(output_file)

#write the data
write(output_data, file=file_name)
