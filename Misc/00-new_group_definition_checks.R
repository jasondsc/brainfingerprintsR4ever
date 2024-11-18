
library(ggplot2)
library(dplyr)
library(BayesFactor)

# repeat prior analysis with new definition of age groups
# see definitions below

## READ IN DATA 
cbbPalette= c("#4f89e0", "#f54c6c","#76D224", '#e4e000', '#4d1f90')
setwd('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/')

demo= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/demo_clean_no_duplicates.csv')
demo= demo[order(demo$UniqueID),]
demo= demo[1:403,]

psd_training = read.csv('./outputs/BIC_PSD_training.csv', header = FALSE)
psd_validation = psd_training[seq(2,802,2),]
psd_training= psd_training[seq(1,802,2),]


aperPSD_training = read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting//outputs/BIC_aperioidc_Rest.csv', header = FALSE)
aperPSD_validation = aperPSD_training[seq(2,802,2),]
aperPSD_training= aperPSD_training[seq(1,802,2),]


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


psd_training= psd_training[index,]
psd_validation= psd_validation[index,]

aperPSD_training= aperPSD_training[index,]
aperPSD_validation= aperPSD_validation[index,]

corrPSD_training= corrPSD_training[index,]
corrPSD_validation= corrPSD_validation[index,]


## CREATED NEW AGE DEFINITIONS FOR DIFFERENTIATION
demo=demo[1:401,]
demo$Group= 'below 11 years old'
demo$Group[demo$Age >=11 & demo$Age <20]= '11-20 years old'
demo$Group[demo$Age >=20]= '20+ years old'


corr_psd=cor(t(psd_training),t(psd_validation)) 
demo$psd_self_corr=diag(corr_psd)
demo$psd_identifiability=diag(apply(corr_psd, 2, scale))

print("psd fingerprinting:")
tt=apply(corr_psd, 1, which.max)
sum(seq(1:401)==tt)/401
tt=apply(corr_psd, 2, which.max)
sum(seq(1:401)==tt)/401


# differentiation accuracies

# look @ age groups now 
# start with below 11
psd_training_age= psd_training[demo$Group== 'below 11 years old',]
psd_validation_age= psd_validation[demo$Group== 'below 11 years old',]
sample_size= dim(psd_validation_age)[1]

fingerprint_boot=matrix(, ncol = 1,nrow=2000)

temp_boot=matrix(, ncol = 2,nrow=1000)

for(i in 1:1000) {
  index_random_par=sample(1:sample_size, 119)
  
  corr_psd=cor(t(psd_training_age[index_random_par,]),t(psd_validation_age[index_random_par,])) 
  tt=apply(corr_psd, 1, which.max)
  temp_boot[i,1]=sum(seq(1:119)==tt)/119
  tt=apply(corr_psd, 2, which.max)
  temp_boot[i,2]=sum(seq(1:119)==tt)/119
  
}

fingerprint_boot[,1]= as.vector(temp_boot)
write.csv(fingerprint_boot, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_fingerprinting_below111.csv')


beepr::beep(2)

# TEENS 11-20
psd_training_age= psd_training[demo$Group== '11-20 years old',]
psd_validation_age= psd_validation[demo$Group== '11-20 years old',]
sample_size= dim(psd_validation_age)[1]

fingerprint_boot=matrix(, ncol = 1,nrow=2000)

temp_boot=matrix(, ncol = 2,nrow=1000)

for(i in 1:1000) {
  index_random_par=sample(1:sample_size, 70)
  
  corr_psd=cor(t(psd_training_age[index_random_par,]),t(psd_validation_age[index_random_par,])) 
  tt=apply(corr_psd, 1, which.max)
  temp_boot[i,1]=sum(seq(1:70)==tt)/70
  tt=apply(corr_psd, 2, which.max)
  temp_boot[i,2]=sum(seq(1:70)==tt)/70
  
}

fingerprint_boot[,1]= as.vector(temp_boot)
write.csv(fingerprint_boot, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_fingerprinting_11-20.csv')



# adults (20 plus) 
psd_training_age= psd_training[demo$Group== '20+ years old',]
psd_validation_age= psd_validation[demo$Group== '20+ years old',]
sample_size= dim(psd_validation_age)[1]

fingerprint_boot=matrix(, ncol = 1,nrow=2000)

temp_boot=matrix(, ncol = 2,nrow=1000)

for(i in 1:1000) {
  index_random_par=sample(1:sample_size, 171)
  
  corr_psd=cor(t(psd_training_age[index_random_par,]),t(psd_validation_age[index_random_par,])) 
  tt=apply(corr_psd, 1, which.max)
  temp_boot[i,1]=sum(seq(1:171)==tt)/171
  tt=apply(corr_psd, 2, which.max)
  temp_boot[i,2]=sum(seq(1:171)==tt)/171
  
}

fingerprint_boot[,1]= as.vector(temp_boot)
write.csv(fingerprint_boot, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_fingerprinting_20.csv')


beepr::beep(2)

# plot effects 

data1=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_fingerprinting_fullcohort2.csv')

data1=data1[,-1]
data1=data1[,1]
data4plot= data.frame(mean_differen=mean(data1), lower=quantile(data1, probs=c(0.025,0.975))[1], upper= quantile(data1, probs=c(0.025,0.975))[2])

data4plot$band= c("full band")
data4plot$group= 'full cohort'

data1=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_fingerprinting_below111.csv')

data1=data1[,-1]

data4plot2= data.frame(mean_differen=mean(data1), lower=quantile(data1, probs=c(0.025,0.975))[1], upper= quantile(data1, probs=c(0.025,0.975))[2])

data4plot2$band= c("full band")
data4plot2$group= 'below 11'

data4plot= rbind(data4plot, data4plot2)


data1=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_fingerprinting_11-20.csv')

data1=data1[,-1]

data4plot2= data.frame(mean_differen=mean(data1), lower=quantile(data1, probs=c(0.025,0.975))[1], upper= quantile(data1, probs=c(0.025,0.975))[2])

data4plot2$band= c("full band")
data4plot2$group= '11-20 years old'

data4plot= rbind(data4plot, data4plot2)


data1=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_fingerprinting_20.csv')

data1=data1[,-1]

data4plot2= data.frame(mean_differen=mean(data1), lower=quantile(data1, probs=c(0.025,0.975))[1], upper= quantile(data1, probs=c(0.025,0.975))[2])

data4plot2$band= c("full band")
data4plot2$group= '20+'

data4plot= rbind(data4plot, data4plot2)


data4plot[,c(1,2,3)]= data4plot[,c(1,2,3)]*100

data4plot$band= factor(data4plot$band, levels = c("full band"))

data4plot$group= factor(data4plot$group, levels = c("full cohort", "below 11", "11-20 years old", "20+"))


ggplot(data=data4plot, aes(x=band, y=mean_differen, fill=group, width=.5)) +  coord_cartesian(ylim=c(0,100))+geom_bar(stat = "identity", position=position_dodge(width = 0.5)) + geom_errorbar(aes(ymin=lower, ymax=upper), width=.2 ,position=position_dodge(.5))  +theme_minimal() + scale_fill_manual(values=cbbPalette) + labs(y="Accuracy (%)", x="",title="") + theme(text = element_text(size=40), legend.title=element_blank()) + ggpubr::theme_classic2()


ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/differentiation_accuracy_bootstrapped_across_age_groups_newdef_NOV2023.pdf', device = "pdf", width=25, height=4.51, dpi=800)

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/differentiation_accuracy_bootstrapped_across_age_groupsnewdef_NOV2023.jpeg', device = "jpeg", width=25, height=4.51, dpi=800)



######################## aperiodic data ######################################


# look @ age groups now 
# start with below 11
psd_training_age= aperPSD_training[demo$Group== 'below 11 years old',]
psd_validation_age= aperPSD_validation[demo$Group== 'below 11 years old',]
sample_size= dim(psd_validation_age)[1]

fingerprint_boot=matrix(, ncol = 1,nrow=2000)

temp_boot=matrix(, ncol = 2,nrow=1000)

for(i in 1:1000) {
  index_random_par=sample(1:sample_size, 133)
  
  corr_psd=cor(t(psd_training_age[index_random_par,]),t(psd_validation_age[index_random_par,])) 
  tt=apply(corr_psd, 1, which.max)
  temp_boot[i,1]=sum(seq(1:133)==tt)/133
  tt=apply(corr_psd, 2, which.max)
  temp_boot[i,2]=sum(seq(1:133)==tt)/133
  
}

fingerprint_boot[,1]= as.vector(temp_boot)
write.csv(fingerprint_boot, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_aperiodic_fingerprinting_below11.csv')

beepr::beep(2)

# TEENS 12-18
psd_training_age= aperPSD_training[demo$Group== '11-20 years old',]
psd_validation_age= aperPSD_validation[demo$Group== '11-20 years old',]
sample_size= dim(psd_validation_age)[1]

fingerprint_boot=matrix(, ncol = 1,nrow=2000)

temp_boot=matrix(, ncol = 2,nrow=1000)

for(i in 1:1000) {
  index_random_par=sample(1:sample_size, 51)
  
  corr_psd=cor(t(psd_training_age[index_random_par,]),t(psd_validation_age[index_random_par,])) 
  tt=apply(corr_psd, 1, which.max)
  temp_boot[i,1]=sum(seq(1:51)==tt)/51
  tt=apply(corr_psd, 2, which.max)
  temp_boot[i,2]=sum(seq(1:51)==tt)/51
  
}

fingerprint_boot[,1]= as.vector(temp_boot)
write.csv(fingerprint_boot, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_aperiodic_fingerprinting_11-20.csv')

beepr::beep(2)



# adults
psd_training_age= aperPSD_training[demo$Group== '20+ years old',]
psd_validation_age= aperPSD_validation[demo$Group== '20+ years old',]
sample_size= dim(psd_validation_age)[1]

fingerprint_boot=matrix(, ncol = 1,nrow=2000)

temp_boot=matrix(, ncol = 2,nrow=1000)

for(i in 1:1000) {
  index_random_par=sample(1:sample_size, 171)
  
  corr_psd=cor(t(psd_training_age[index_random_par,]),t(psd_validation_age[index_random_par,])) 
  tt=apply(corr_psd, 1, which.max)
  temp_boot[i,1]=sum(seq(1:171)==tt)/171
  tt=apply(corr_psd, 2, which.max)
  temp_boot[i,2]=sum(seq(1:171)==tt)/171
  
}

fingerprint_boot[,1]= as.vector(temp_boot)
write.csv(fingerprint_boot, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/bootstrapped_aperiodic_fingerprinting_20.csv')

beepr::beep(2)



###################### Look at ICC of new age groups just in case ###############################
atlas= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_atlas.csv')

# look at age groups 
ind=demo$Group=='below 11 years old'
z_target= t(scale(t(psd_validation[ind,])))
z_database= t(scale(t(psd_training[ind,])))
icc = c()


n = length(ind)
k = 2
df_b = n-1
df_w = n*(k-1)


for (i_edge in 1:length(psd_training)){
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

icc_mat_young <- matrix(icc, nrow = 148, byrow = TRUE)

write.csv(icc_mat_young, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_NEWGROUPDEFF_BELOW11.csv')

someData1= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,1:8]), atlas$hemi, band= 'theta')
colnames(someData1)= c('region', 'ICC', 'hemi', 'band')

someData2= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,9:18]), atlas$hemi, band= 'alpha')
colnames(someData2)=  c('region', 'ICC', 'hemi', 'band')

someData3= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,19:52]), atlas$hemi, band= 'beta')
colnames(someData3)= c('region', 'ICC', 'hemi', 'band')

someData4= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,53:92]), atlas$hemi, band= 'gamma')
colnames(someData4)= c('region', 'ICC', 'hemi', 'band')

someData5= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,93:293]), atlas$hemi, band= 'high gamma')
colnames(someData5)= c('region', 'ICC', 'hemi', 'band')

someData= rbind(someData1, someData2, someData3, someData4, someData5)

someData$ICC[someData$ICC<0.7]=0.7
someData$ICC[someData$ICC>1]=1


someData$band= ordered(someData$band, levels= c("theta", "alpha", "beta", "gamma", "high gamma"))


someData %>%
  group_by(band) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =ICC)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ band, nrow = 1) +  colorspace::scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  colorspace::scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_NOV2023_below11.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_NOV2023_below11.jpeg', device = "jpeg")


# look at age groups 
ind=demo$Group=='11-20 years old'
z_target= t(scale(t(psd_validation[ind,])))
z_database= t(scale(t(psd_training[ind,])))
icc = c()


n = length(ind)
k = 2
df_b = n-1
df_w = n*(k-1)


for (i_edge in 1:length(psd_training)){
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

icc_mat_teens <- matrix(icc, nrow = 148, byrow = TRUE)

write.csv(icc_mat_teens, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_NEWGROUPDEFF_11-20.csv')

someData1= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,1:8]), atlas$hemi, band= 'theta')
colnames(someData1)= c('region', 'ICC', 'hemi', 'band')

someData2= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,9:18]), atlas$hemi, band= 'alpha')
colnames(someData2)=  c('region', 'ICC', 'hemi', 'band')

someData3= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,19:52]), atlas$hemi, band= 'beta')
colnames(someData3)= c('region', 'ICC', 'hemi', 'band')

someData4= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,53:92]), atlas$hemi, band= 'gamma')
colnames(someData4)= c('region', 'ICC', 'hemi', 'band')

someData5= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,93:293]), atlas$hemi, band= 'high gamma')
colnames(someData5)= c('region', 'ICC', 'hemi', 'band')

someData= rbind(someData1, someData2, someData3, someData4, someData5)

someData$ICC[someData$ICC<0.7]=0.7
someData$ICC[someData$ICC>1]=1


someData$band= ordered(someData$band, levels= c("theta", "alpha", "beta", "gamma", "high gamma"))


someData %>%
  group_by(band) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =ICC)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ band, nrow = 1) +  colorspace::scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  colorspace::scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_NOV2023_11-20.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_NOV2023_11-20.jpeg', device = "jpeg")

# look at age groups 
ind=demo$Group=='20+ years old'
z_target= t(scale(t(psd_validation[ind,])))
z_database= t(scale(t(psd_training[ind,])))
icc = c()


n = length(ind)
k = 2
df_b = n-1
df_w = n*(k-1)


for (i_edge in 1:length(psd_training)){
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

icc_mat_adults <- matrix(icc, nrow = 148, byrow = TRUE)

write.csv(icc_mat_adults, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_NEWGROUPDEFF_20.csv')

someData1= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,1:8]), atlas$hemi, band= 'theta')
colnames(someData1)= c('region', 'ICC', 'hemi', 'band')

someData2= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,9:18]), atlas$hemi, band= 'alpha')
colnames(someData2)=  c('region', 'ICC', 'hemi', 'band')

someData3= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,19:52]), atlas$hemi, band= 'beta')
colnames(someData3)= c('region', 'ICC', 'hemi', 'band')

someData4= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,53:92]), atlas$hemi, band= 'gamma')
colnames(someData4)= c('region', 'ICC', 'hemi', 'band')

someData5= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,93:293]), atlas$hemi, band= 'high gamma')
colnames(someData5)= c('region', 'ICC', 'hemi', 'band')

someData= rbind(someData1, someData2, someData3, someData4, someData5)

someData$ICC[someData$ICC<0.7]=0.7
someData$ICC[someData$ICC>1]=1


someData$band= ordered(someData$band, levels= c("theta", "alpha", "beta", "gamma", "high gamma"))


someData %>%
  group_by(band) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =ICC)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ band, nrow = 1) +  colorspace::scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  colorspace::scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_NOV2023_20.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_NOV2023_20.jpeg', device = "jpeg")



# difference between kids and adults
icc_diff= icc_mat_young- icc_mat_adults

write.csv(icc_diff, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_difference_youn_adults_NOV2023.csv')



# check how similar the two age definitions are
icc_sickkids=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_difference_youn_adults.csv')

icc_sickkids=icc_sickkids[,-1]

Xkids=data.frame(theta= rowMeans(icc_diff[,1:8]),alpha= rowMeans(icc_diff[,9:18]),beta= rowMeans(icc_diff[,19:52]), gamma= rowMeans(icc_diff[,53:92]), hgamma= rowMeans( icc_diff[,93:293]))

OrigXkids=data.frame(theta= rowMeans(icc_sickkids[,1:8]),alpha= rowMeans(icc_sickkids[,9:18]),beta= rowMeans(icc_sickkids[,19:52]), gamma= rowMeans(icc_sickkids[,53:92]), hgamma= rowMeans( icc_sickkids[,93:293]))


cor.test(rowMeans(Xkids), rowMeans(OrigXkids))

table(demo$Group)




# plot effect 
atlas=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

someData6= tidyr::tibble(atlas$region, rowMeans(Xkids), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC >0.12]=0.12
someData6$ICC[someData6$ICC < -0.12]= -0.12

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  scale_fill_gradient2(low = "#f54c6c", mid = "#FCF7F4", high = "#76D224", midpoint = 0.0, limits=c(-0.12, 0.12 )) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_difference_10years.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_difference_10years.jpeg', device = "jpeg")


ggplot(someData6, aes(x=rowMeans(OrigXkids), y = rowMeans(Xkids), fill='yes', colour='yes')) + 
  geom_jitter() + 
  stat_smooth(method = "lm",fullrange = T) + scale_fill_manual(values=cbbPalette[5])  + scale_colour_manual(values=cbbPalette[5]) +
  theme_minimal() +   xlab("original ICC") + ylab("ICC new cohort")+ ggpubr::theme_classic2() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/scatter_ICC_new_10years.pdf', device = "pdf", width = 5, height = 5)



# PERMUTATION TESTS HERE 
permuted_index= read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/permuted_indexes_of_destriuex_atlas_SPINs&Twirl.csv')
permuted_index= permuted_index[,-1]
permuted_index=permuted_index+1
#gradient offset
orig=cor.test(rowMeans(Xkids), rowMeans(OrigXkids))
permuted_corr=c()
for (i in 1:1000){
  
  cor_temp=cor.test(rowMeans(Xkids), rowMeans(OrigXkids)[permuted_index[,i]])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
  
}
sum(orig$estimate < permuted_corr)/1000 # p 0.00




################################################## aperiodic corr ICC compare ######################


# look at age groups 
ind=demo$Group=='below 11 years old'
z_target= t(scale(t(corrPSD_validation[ind,])))
z_database= t(scale(t(corrPSD_training[ind,])))
icc = c()


n = sum(ind)
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

icc_mat_young <- matrix(icc, nrow = 148, byrow = TRUE)

write.csv(icc_mat_young, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_apericorr_youngkids_analysis_NOV2023.csv')

someData1= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,1:8]), atlas$hemi, band= 'theta')
colnames(someData1)= c('region', 'ICC', 'hemi', 'band')

someData2= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,9:18]), atlas$hemi, band= 'alpha')
colnames(someData2)=  c('region', 'ICC', 'hemi', 'band')

someData3= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,19:52]), atlas$hemi, band= 'beta')
colnames(someData3)= c('region', 'ICC', 'hemi', 'band')

someData4= tidyr::tibble(atlas$region, rowMeans(icc_mat_young[,53:92]), atlas$hemi, band= 'gamma')
colnames(someData4)= c('region', 'ICC', 'hemi', 'band')

someData= rbind(someData1, someData2, someData3, someData4)

someData$ICC[someData$ICC<0.7]=0.7
someData$ICC[someData$ICC>1]=1


someData$band= ordered(someData$band, levels= c("theta", "alpha", "beta", "gamma"))


someData %>%
  group_by(band) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =ICC)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ band, nrow = 1) +  colorspace::scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  colorspace::scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_apericorr_ICC_youngkids_NOV2023.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_apericorr_ICC_youngkids_NOV2023.jpeg', device = "jpeg")


# look at age groups 
ind=demo$Group=='11-20 years old'
z_target= t(scale(t(corrPSD_validation[ind,])))
z_database= t(scale(t(corrPSD_training[ind,])))
icc = c()


n = sum(ind)
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

icc_mat_teens <- matrix(icc, nrow = 148, byrow = TRUE)

write.csv(icc_mat_teens, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_apericorr_11-20_analysis_NOV2023.csv')

someData1= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,1:8]), atlas$hemi, band= 'theta')
colnames(someData1)= c('region', 'ICC', 'hemi', 'band')

someData2= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,9:18]), atlas$hemi, band= 'alpha')
colnames(someData2)=  c('region', 'ICC', 'hemi', 'band')

someData3= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,19:52]), atlas$hemi, band= 'beta')
colnames(someData3)= c('region', 'ICC', 'hemi', 'band')

someData4= tidyr::tibble(atlas$region, rowMeans(icc_mat_teens[,53:92]), atlas$hemi, band= 'gamma')
colnames(someData4)= c('region', 'ICC', 'hemi', 'band')

someData= rbind(someData1, someData2, someData3, someData4)

someData$ICC[someData$ICC<0.7]=0.7
someData$ICC[someData$ICC>1]=1


someData$band= ordered(someData$band, levels= c("theta", "alpha", "beta", "gamma"))


someData %>%
  group_by(band) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =ICC)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ band, nrow = 1) +  colorspace::scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  colorspace::scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_apericorr_ICC_11-20NOV2023.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_apericorr_ICC_11-20NOV2023.jpeg', device = "jpeg")

# look at age groups 
ind=demo$Group=='20+ years old'
z_target= t(scale(t(corrPSD_validation[ind,])))
z_database= t(scale(t(corrPSD_training[ind,])))
icc = c()


n = sum(ind)
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

icc_mat_adults <- matrix(icc, nrow = 148, byrow = TRUE)

write.csv(icc_mat_adults, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_apericorr_20+_analysis_NOV2023.csv')

someData1= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,1:8]), atlas$hemi, band= 'theta')
colnames(someData1)= c('region', 'ICC', 'hemi', 'band')

someData2= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,9:18]), atlas$hemi, band= 'alpha')
colnames(someData2)=  c('region', 'ICC', 'hemi', 'band')

someData3= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,19:52]), atlas$hemi, band= 'beta')
colnames(someData3)= c('region', 'ICC', 'hemi', 'band')

someData4= tidyr::tibble(atlas$region, rowMeans(icc_mat_adults[,53:92]), atlas$hemi, band= 'gamma')
colnames(someData4)= c('region', 'ICC', 'hemi', 'band')

someData= rbind(someData1, someData2, someData3, someData4)

someData$ICC[someData$ICC<0.7]=0.7
someData$ICC[someData$ICC>1]=1


someData$band= ordered(someData$band, levels= c("theta", "alpha", "beta", "gamma"))


someData %>%
  group_by(band) %>% 
  brain_join(desterieux) %>% 
  reposition_brain(hemi ~ side) %>% 
  ggplot(aes(fill =ICC)) + 
  geom_sf(show.legend = TRUE) + 
  facet_wrap( ~ band, nrow = 1) +  colorspace::scale_fill_continuous_sequential(palette= 'Sunset', rev= FALSE) +  colorspace::scale_color_continuous_sequential(palette= 'Sunset', rev= FALSE) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_apericorr_ICC_20+NOV2023.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_apericorr_ICC_20+NOV2023.jpeg', device = "jpeg")



# difference between kids and adults


icc_diff= icc_mat_young- icc_mat_adults

write.csv(icc_diff, 'C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_apericorr_difference_youn_adults_NOV2023.csv')

icc_sickkidscorr=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_output/ICC_apericorr_difference_youn_adultsDEC2023.csv')

icc_sickkidscorr=icc_sickkidscorr[,-1]

OrigCkids=data.frame(theta= rowMeans(icc_sickkidscorr[,1:8]),alpha= rowMeans(icc_sickkidscorr[,9:18]),beta= rowMeans(icc_sickkidscorr[,19:52]), gamma= rowMeans(icc_sickkidscorr[,53:93]))

NEWCkids=data.frame(theta= rowMeans(icc_diff[,1:8]),alpha= rowMeans(icc_diff[,9:18]),beta= rowMeans(icc_diff[,19:52]), gamma= rowMeans(icc_diff[,53:93]))


cor.test(rowMeans(OrigCkids), rowMeans(NEWCkids))



atlas=read.csv('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

someData6= tidyr::tibble(atlas$region, rowMeans(NEWCkids), atlas$hemi, band= 'broadband')
colnames(someData6)= c('region', 'ICC', 'hemi', 'band')

someData6$ICC[someData6$ICC >0.12]=0.12
someData6$ICC[someData6$ICC < -0.12]= -0.12

ggplot(someData6) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = ICC)) + 
  scale_fill_gradient2(low = "#f54c6c", mid = "#FCF7F4", high = "#76D224", midpoint = 0.0, limits=c(-0.12, 0.12 )) +
  theme_void() + scale_color_manual('white')

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_CORRSPEC_difference_10years.pdf', device = "pdf")

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/brain_map_ICC_CORRSPEC_difference_10years.jpeg', device = "jpeg")


ggplot(someData6, aes(x=rowMeans(OrigCkids), y = rowMeans(NEWCkids), fill='yes', colour='yes')) + 
  geom_jitter() + 
  stat_smooth(method = "lm",fullrange = T) + scale_fill_manual(values=cbbPalette[5])  + scale_colour_manual(values=cbbPalette[5]) +
  theme_minimal() +   xlab("original ICC") + ylab("ICC new cohort")+ ggpubr::theme_classic2() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 

ggsave('C:/Users/J Da Silva Castanhei/Documents/SickKidz_fingerprinting/R_figures/scatter_aperiodiccorr_ICC_new_10years.pdf', device = "pdf", width = 5, height = 5)


#gradient offset
orig=cor.test(rowMeans(OrigCkids), rowMeans(NEWCkids))
permuted_corr=c()
for (i in 1:1000){
  
  cor_temp=cor.test(rowMeans(OrigCkids), rowMeans(NEWCkids)[permuted_index[,i]])
  permuted_corr= c(permuted_corr, cor_temp$estimate)
  
}
sum(orig$estimate < permuted_corr)/1000 # p 0.00

