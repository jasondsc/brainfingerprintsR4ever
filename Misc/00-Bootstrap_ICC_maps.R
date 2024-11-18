
library(ggplot2)
library(dplyr)
library(BayesFactor)

# used just to check that the ICC isnt driven by random outlier subjects
# bootstrapped the effects 

############## READ IN DATA  ####################
cbbPalette= c("#4f89e0", "#f54c6c","#76D224", '#e4e000', '#4d1f90')


#ge_colors = c("#000000","#007E7E","#8001FC","#FF8207","#FFFFFF")
#scale_fill_gradientn(colors = ge_colors,oob = scales::squish,limits=c(floor(rng[1]), ceiling(rng[2])))
setwd('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/')

demo= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/demo_clean_no_duplicates.csv')
demo= demo[order(demo$UniqueID),]
demo= demo[1:403,]

corrPSD_training = read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting//outputs/BIC_corrspec_BIC.csv', header = FALSE)
corrPSD_validation = corrPSD_training[seq(2,802,2),]
corrPSD_training= corrPSD_training[seq(1,802,2),]

subj_Id_psd= read.csv('./outputs/subject_brainstorm_ids.csv', header = TRUE)

# note remove this when all subjets have been processed and accounted for
demo=demo[c(1:29, 31:240, 242:404),]
index=order(subj_Id_psd$subjID)
subj_psd_order=subj_Id_psd$subjID[index]

paste('sub-',demo$UniqueID, sep = '') == subj_Id_psd$subjID[index]
sum(paste('sub-',demo$UniqueID, sep = '') == subj_Id_psd$subjID[index])

corrPSD_training= corrPSD_training[index,]
corrPSD_validation= corrPSD_validation[index,]


demo=demo[1:401,]
demo$Group= 'below 12 years old'
demo$Group[demo$Age >=12 & demo$Age <18]= '12-18 years old'
demo$Group[demo$Age >=18]= '18+ years old'



############### KIDS ICC ##################

ind=demo$Group=='below 12 years old'

temp= array(, c(100, 148))

for (i in 1:100){
  
  indsub= sample(which(ind), 50)
  z_target= t(scale(t(corrPSD_validation[indsub,])))
  z_database= t(scale(t(corrPSD_training[indsub,])))
  icc = c()
  
  
  n = 50
  k = 2
  df_b = n-1
  df_w = n*(k-1)
  
  
  for (i_edge in 1:length(corrPSD_training)){
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
  
  temp[i,] = rowMeans(data.frame(theta= rowMeans(icctemp[,1:8]),alpha= rowMeans(icctemp[,9:18]),beta= rowMeans(icctemp[,19:52]), gamma= rowMeans(icctemp[,53:93])))

  }

write.csv(temp, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/BOOTS_ICC_youngkids_aperiodiccorr.csv')

############### TEEN ICC ##################

ind=demo$Group=='12-18 years old'

temp= array(, c(100, 148))

for (i in 1:100){
  
  indsub= sample(which(ind), 50)
  z_target= t(scale(t(corrPSD_validation[indsub,])))
  z_database= t(scale(t(corrPSD_training[indsub,])))
  icc = c()
  
  
  n = 50
  k = 2
  df_b = n-1
  df_w = n*(k-1)
  
  
  for (i_edge in 1:length(corrPSD_training)){
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
  
  temp[i,] = rowMeans(data.frame(theta= rowMeans(icctemp[,1:8]),alpha= rowMeans(icctemp[,9:18]),beta= rowMeans(icctemp[,19:52]), gamma= rowMeans(icctemp[,53:93])))
  
  
}

write.csv(temp, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/BOOTS_ICC_teens_aperiodiccorr.csv')


############### ADULTS ICC ##################

ind=demo$Group=='18+ years old'

temp= array(, c(100, 148))

for (i in 1:100){
  
  indsub= sample(which(ind), 50)
  z_target= t(scale(t(corrPSD_validation[indsub,])))
  z_database= t(scale(t(corrPSD_training[indsub,])))
  icc = c()
  
  
  n = 50
  k = 2
  df_b = n-1
  df_w = n*(k-1)
  
  
  for (i_edge in 1:length(corrPSD_training)){
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
  
  temp[i,] = rowMeans(data.frame(theta= rowMeans(icctemp[,1:8]),alpha= rowMeans(icctemp[,9:18]),beta= rowMeans(icctemp[,19:52]), gamma= rowMeans(icctemp[,53:93])))
  
  
}

write.csv(temp, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/BOOTS_ICC_adults_aperiodiccorr.csv')


############### PLOT EFFECTS  ##################
library(ggplot2)
library(ggseg)
library(ggsegDesterieux)

icc_kids=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/BOOTS_ICC_youngkids_aperiodiccorr.csv')
icc_kids=icc_kids[,-1]


icc_teen=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/BOOTS_ICC_teens_aperiodiccorr.csv')
icc_teen=icc_teen[,-1]

icc_ad=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/BOOTS_ICC_adults_aperiodiccorr.csv')
icc_ad=icc_ad[,-1]



atlas=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

someData6= tidyr::tibble(atlas$region, colMeans(icc_kids), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC < 0.45] = 0.45
someData6$ICC[someData6$ICC > 0.75 ] = 0.75

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  colorspace::scale_fill_continuous_sequential(palette= 'Batlow', rev= FALSE, limits= c(0.45, 0.75)) + 
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/BOOTSTRAP_brain_map_ICC_kids.pdf', device = "pdf")


atlas=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')


someData6= tidyr::tibble(atlas$region, colMeans(icc_teen), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC < 0.45] = 0.45
someData6$ICC[someData6$ICC > 0.75 ] = 0.75

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  colorspace::scale_fill_continuous_sequential(palette= 'Batlow', rev= FALSE, limits= c(0.45, 0.75)) + 
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/BOOTSTRAP_brain_map_ICC_teens.pdf', device = "pdf")



atlas=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

someData6= tidyr::tibble(atlas$region, colMeans(icc_ad), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC < 0.45] = 0.45
someData6$ICC[someData6$ICC > 0.75 ] = 0.75

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  colorspace::scale_fill_continuous_sequential(palette= 'Batlow', rev= FALSE, limits= c(0.45, 0.75)) + 
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/BOOTSTRAP_brain_map_ICC_adults.pdf', device = "pdf")



XCkids= icc_kids-icc_ad

someData6= tidyr::tibble(atlas$region, rowMeans(XCkids), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC < -0.1] = -0.1
someData6$ICC[someData6$ICC > 0.1 ] = 0.1

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  scale_fill_gradient2(low = "#f54c6c", mid = "#FCF7F4", high = "#76D224", midpoint = 0.0, limits=c(-0.1, 0.1 )) +
  theme_void() + scale_color_manual('white')


ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/BOOTSTRAP_brain_map_aperiodic_corr_kids_adults.pdf', device = "pdf")



XCkids= icc_kids-icc_teen
someData6= tidyr::tibble(atlas$region, rowMeans(XCkids), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC < -0.1] = -0.1
someData6$ICC[someData6$ICC > 0.1 ] = 0.1

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  scale_fill_gradient2(low = "#f54c6c", mid = "#FCF7F4", high = "#76D224", midpoint = 0.0, limits=c(-0.1, 0.1 )) +
  theme_void() + scale_color_manual('white')


ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/BOOTSTRAP_brain_map_aperiodic_corr_kids_teen.pdf', device = "pdf")



XCkids= icc_teen- icc_ad
someData6= tidyr::tibble(atlas$region, rowMeans(XCkids), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC < -0.1] = -0.1
someData6$ICC[someData6$ICC > 0.1 ] = 0.1

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  scale_fill_gradient2(low = "#f54c6c", mid = "#FCF7F4", high = "#76D224", midpoint = 0.0, limits=c(-0.1, 0.1 )) +
  theme_void() + scale_color_manual('white')


ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/BOOTSTRAP_brain_map_aperiodic_corr_teen_adults.pdf', device = "pdf")


# check relationship of BOOTstrap ICC and original ICC

icc_sickkidscorr=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_apericorr_difference_youn_adultsDEC2023.csv')

icc_sickkidscorr=icc_sickkidscorr[,-1]

XCkids=data.frame(theta= rowMeans(icc_sickkidscorr[,1:8]),alpha= rowMeans(icc_sickkidscorr[,9:18]),beta= rowMeans(icc_sickkidscorr[,19:52]), gamma= rowMeans(icc_sickkidscorr[,53:93]))

BOOTICC= icc_kids-icc_ad

cor.test(colMeans(BOOTICC), rowMeans(XCkids))

# PERMUTATION TESTS HERE 
permuted_index= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/permuted_indexes_of_destriuex_atlas_SPINs&Twirl.csv')
permuted_index= permuted_index[,-1]
permuted_index=permuted_index+1
#gradient offset
orig=cor.test(colMeans(BOOTICC), rowMeans(XCkids))

permuted_corr=c()
for (i in 1:1000){
  
  cor_temp=cor.test(rowMeans(XCkids), colMeans(BOOTICC)[permuted_index[,i]])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
  
}
(sum(orig$estimate < permuted_corr)+1)/1001 # p 0.0001


cor.test(colMeans(BOOTICC), atlas$fcgradient01) # slightly stronger

orig=cor.test(colMeans(BOOTICC), atlas$fcgradient01)

permuted_corr=c()
for (i in 1:1000){
  
  cor_temp=cor.test(rowMeans(XCkids), atlas$fcgradient01[permuted_index[,i]])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
  
}
(sum(orig$estimate < permuted_corr)+1)/1001 # p 0.0001



gene_expression= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destiruex_abagen/Destrieux_expression_genes.csv')

index= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destiruex_abagen/index_destriuex.csv', header = FALSE)

labels= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destiruex_abagen/destriuex_atlas_labels.csv')

labels=labels[labels$structure== 'cortex',]

atlas2=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

atlas2=atlas2[index$V1,]

gene_names= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destiruex_abagen/top50_per_gene_loadings.csv')

pos_gene_expression= rowMeans(gene_expression[,which(colnames(gene_expression) %in% gene_names$pos_genes)])

neg_gene_expression= rowMeans(gene_expression[,which(colnames(gene_expression) %in% gene_names$neg_genes)])


BOOTICC= BOOTICC[,index$V1]

cor.test(colMeans(BOOTICC), pos_gene_expression) # stronger effect -0,45
cor.test(colMeans(BOOTICC), neg_gene_expression) # 0.44


#gradient offset
orig=cor.test(colMeans(BOOTICC)[order(index$V1)], pos_gene_expression[order(index$V1)])
permuted_corr=c()
for (i in 1:1000){
  
  cor_temp=cor.test(colMeans(BOOTICC)[order(index$V1)], pos_gene_expression[order(index$V1)][permuted_index[,i]])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
  
}
(sum(orig$estimate > permuted_corr)+1)/1001 # p 0.0019

#gradient offset
orig=cor.test(colMeans(BOOTICC)[order(index$V1)], neg_gene_expression[order(index$V1)])
permuted_corr=c()
for (i in 1:1000){
  
  cor_temp=cor.test(colMeans(BOOTICC)[order(index$V1)], neg_gene_expression[order(index$V1)][permuted_index[,i]])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
  
}
(sum(orig$estimate < permuted_corr)+1)/1001 # p 0.0039

