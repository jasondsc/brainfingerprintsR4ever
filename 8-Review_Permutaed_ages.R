# data from SickKids from seq(4,150,0.5) 4 to 150 hz
# 4- 40 Hz
# sickkids 1:73
# CamCan  7-79

# get useful libraries and set colour pallet
library(ggplot2)
library(ggseg)
library(dplyr)
library(colorspace)
library(ggpubr)
library(ggsegDesterieux)

setwd('~/Downloads/SickKidsCode/')

cbbPalette= c("#4f89e0", "#f5ec6c",'#156605',"#76D7C4","#268236", '#4d3d87','#593d80', '#3b49e3')

### Read in CamCan dataset 
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo=read.csv('./approved_data.csv')
full_demo$CCID=paste('sub',full_demo$CCID, sep = '_')
index=full_demo$CCID %in% ids$Ids
full_demo=full_demo[index,]

demo$Group[demo$Age>18 & demo$Age <45]= "young adult"
demo$Group[demo$Age>=45 & demo$Age <65]= "adult"
demo$Group[demo$Age>=65 & demo$Age <90]= "older adult"

# aperiodic corr
ap_corr_psd_rest=read.csv('./Destriuex_corraperioidc_rest1.csv', header = FALSE)
ap_corr_psd_task=read.csv('./Destriuex_corraperiodic_rest2.csv', header = FALSE)

# index to reorder desterieux atlas bc brainstorm has chnaged the order in 2024
index_dest= read.csv('./INDEX_Desteriuex.csv', header = FALSE)

# reorder CamCAN matrix
atlas2=read.csv('./Neuromaps_destrieux_full_atlas.csv')

# save aperiodic corrected data
list=index_dest$V1
roi_index= c(mapply(seq, 79*(list-1)+7, (79*(list-1)+79)))

ap_corr_psd_train= ap_corr_psd_task[,roi_index]
ap_corr_psd_valid= ap_corr_psd_rest[,roi_index]

# read in SickKids data
demoSickkids= read.csv('./ToMove_demo.csv')
SickKid_ap_corr_psd_valid= read.csv('./ToMove_corrspecValidation.csv')
SickKid_ap_corr_psd_train= read.csv('./ToMove_corrspecTraining.csv')

SickKid_ap_corr_psd_valid= SickKid_ap_corr_psd_valid[-1]
SickKid_ap_corr_psd_train= SickKid_ap_corr_psd_train[-1]

list=1:148
roi_index= c(mapply(seq, 93*(list-1)+1, (93*(list-1)+73)))

SickKid_ap_corr_psd_valid= SickKid_ap_corr_psd_valid[,roi_index]
SickKid_ap_corr_psd_train= SickKid_ap_corr_psd_train[,roi_index]

colnames(SickKid_ap_corr_psd_train)=1:10804
colnames(ap_corr_psd_train)=1:10804

colnames(SickKid_ap_corr_psd_valid)=1:10804
colnames(ap_corr_psd_valid)=1:10804

# merge the two datasets into one large one
training= rbind(ap_corr_psd_train, SickKid_ap_corr_psd_train)
validation= rbind(ap_corr_psd_valid, SickKid_ap_corr_psd_valid)

Ages= c(demo$Age, demoSickkids$Age)

rm(SickKid_ap_corr_psd_train, SickKid_ap_corr_psd_valid, ap_corr_psd_task, ap_corr_psd_rest)

temp= array(, c(500,37, 148))
age_mean= array(, c(1,37))

atlas=read.csv('./Neuromaps_destrieux_full_atlas.csv')

for (perm in 1:500){
  
  Ages_temp=sample(Ages)
  age_order= order(Ages_temp)
for (i in 1:37){
  
  # start= 25*(i-1) # sliding window of 50 people with 50% overlap
  # end= (25+(i*25))-1
  # if(i==1){
  #  start=1
  # } #39!!!
  start= 25*(i-1) # sliding window of 100 with 75% overlap
  end= (75+(i*25))-1
  if(i==1){
    start=1
  }
  indsub= age_order[start:end]
  z_target= t(scale(t(training[indsub,])))
  z_database= t(scale(t(validation[indsub,])))
  #age_mean[i]= mean(Ages[indsub])
  icc = c()
  
  n = 100
  k = 2
  df_b = n-1
  df_w = n*(k-1)
  
  
  for (i_edge in 1:length(validation)){
    x= data.frame(unlist(z_target[,i_edge]),unlist(z_database[,i_edge]) )
    x_w_mean =  rowMeans(x)
    x_g_mean = mean(unlist(x))
    ss_t = sum(((x - x_g_mean)^2))
    ss_w = sum((x - (x_w_mean))^2)
    ss_b = ss_t - ss_w
    ms_b = ss_b / df_b
    ms_w = ss_w / df_w
    icc[i_edge] = (ms_b - ms_w) / (ms_b + ((k-1)*ms_w))
    
    
  }
  
  icctemp <- matrix(icc, nrow = 148, byrow = TRUE)
  
  temp[perm,i,] = rowMeans(data.frame(theta= rowMeans(icctemp[,1:8]),
                                 alpha= rowMeans(icctemp[,9:18]),
                                 beta= rowMeans(icctemp[,19:52]), 
                                 gamma= rowMeans(icctemp[,53:73])))
  

}
  
  print(perm)
}



for (i in 1:37){
  
  # start= 25*(i-1) # sliding window of 50 people with 50% overlap
  # end= (25+(i*25))-1
  # if(i==1){
  #  start=1
  # } #39!!!
  start= 25*(i-1) # sliding window of 100 with 75% overlap
  end= (75+(i*25))-1
  if(i==1){
    start=1
  }
  indsub= age_order[start:end]
  z_target= t(scale(t(training[indsub,])))
  z_database= t(scale(t(validation[indsub,])))
  #age_mean[i]= mean(Ages[indsub])
  icc = c()
  
  n = 100
  k = 2
  df_b = n-1
  df_w = n*(k-1)
  
  
  for (i_edge in 1:length(validation)){
    x= data.frame(unlist(z_target[,i_edge]),unlist(z_database[,i_edge]) )
    x_w_mean =  rowMeans(x)
    x_g_mean = mean(unlist(x))
    ss_t = sum(((x - x_g_mean)^2))
    ss_w = sum((x - (x_w_mean))^2)
    ss_b = ss_t - ss_w
    ms_b = ss_b / df_b
    ms_w = ss_w / df_w
    icc[i_edge] = (ms_b - ms_w) / (ms_b + ((k-1)*ms_w))
    
    
  }
  
  icctemp <- matrix(icc, nrow = 148, byrow = TRUE)
  
  temp[perm,i,] = rowMeans(data.frame(theta= rowMeans(icctemp[,1:8]),
                                 alpha= rowMeans(icctemp[,9:18]),
                                 beta= rowMeans(icctemp[,19:52]), 
                                 gamma= rowMeans(icctemp[,53:73])))
  

}
  

save(temp, file = "./Permuted500.rds")

tt=readRDS("./Permuted500.rds")

age_order= order(Ages)
for (i in 1:37){
  
  start= 25*(i-1) # sliding window of 100 with 75% overlap
  end= (75+(i*25))-1
  if(i==1){
    start=1
  }
indsub= age_order[start:end]
age_mean[i]= mean(Ages[indsub])

}


# read in genetic data obtained from abagen
gene_expression= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/Destiruex_abagen/Destrieux_expression_genes.csv')

index= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/Destiruex_abagen/index_destriuex.csv', header = FALSE)

labels= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/Destiruex_abagen/destriuex_atlas_labels.csv')

labels=labels[labels$structure== 'cortex',]

atlas2=read.csv('/Users/jason/Desktop/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

atlas2=atlas2[index$V1,]

gene_names= read.csv('/Users/jason/Documents/HCP_twin_projec/top50_per_gene_loadings.csv')

pos_gene_expression= rowMeans(gene_expression[,which(colnames(gene_expression) %in% gene_names$pos_genes)])

neg_gene_expression= rowMeans(gene_expression[,which(colnames(gene_expression) %in% gene_names$neg_genes)])


effect_pos=c()
effect_neg=c()
for (perm in 1:500){
  
  temp=apply(tt[perm,,index$V1], 1, cor.test, neg_gene_expression)
  neg=unlist(temp)
  temp=apply(tt[perm,,index$V1], 1, cor.test, pos_gene_expression)
  pos=unlist(temp)
  # plot effects
  data4plot= data.frame(value=c(neg[seq(4,370,10)], pos[seq(4,370,10)]),
                        upper=c(neg[seq(10,370,10)], pos[seq(10,370,10)] ),
                        lower= c(neg[seq(9,370,10)], pos[seq(9,370,10)]),
                        gene= rep(c('negative', 'positive'),each=37),
                        group= c(age_mean, age_mean))
  
  data4plotPOS= data4plot[data4plot$gene== 'positive',]
  data4plotPOS$valueZ= DescTools::FisherZ(as.numeric(data4plotPOS$value))
  
  lm3=lm(valueZ ~poly(group, 3), data4plotPOS)
  effect_pos[perm]= coef(lm3)[[4]]
  
}

(sum(0.73< effect_pos)+1)/1001 # 0.000999001


# now let us repeat for the first functional gradient 
atlas=read.csv('/Users/jason/Desktop/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

effect_FC=c()
for (perm in 1:500){
  
  temp=apply(tt[perm,,], 1, cor.test, atlas$fcgradient01)
  neg=unlist(temp)
  # plot effects
  data4plot= data.frame(value=c(neg[seq(4,370,10)]),
                        upper=c(neg[seq(10,370,10)] ),
                        lower= c(neg[seq(9,370,10)]),
                        group= c(age_mean))
  
  data4plot$valueZ= DescTools::FisherZ(as.numeric(data4plot$value))
  
  lm3=lm(valueZ ~poly(group, 3), data4plot)
  effect_FC[perm]= coef(lm3)[[4]]
  
}

(sum(-0.46> effect_FC)+1)/1001 # 0.005994006


