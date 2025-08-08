# data from SickKids from seq(4,150,0.5) 4 to 150 hz
# 4- 40 Hz
# sickkids 1:73
# CamCan  7-79

library(ggplot2)
library(ggseg)
library(dplyr)
library(colorspace)
library(ggpubr)
library(ggsegDesterieux)

cbbPalette= c("#4f89e0", "#f5ec6c",'#156605',"#76D7C4","#268236", '#4d3d87','#593d80', '#3b49e3')

### Read in CamCan dataset 
demo=read.csv('~/Documents/CAMCAN_outputs/CAMCAN_Analysis/standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('~/Documents/CAMCAN_outputs/CAMCAN_Analysis/subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo=read.csv('~/Documents/CAMCAN_outputs/approved_data.csv')
full_demo$CCID=paste('sub',full_demo$CCID, sep = '_')
index=full_demo$CCID %in% ids$Ids
full_demo=full_demo[index,]

demo$Group[demo$Age>18 & demo$Age <45]= "young adult"
demo$Group[demo$Age>=45 & demo$Age <65]= "adult"
demo$Group[demo$Age>=65 & demo$Age <90]= "older adult"

# aperiodic corr
ap_corr_psd_rest=read.csv('~/Documents/CAMCAN_outputs/CAMCAN_destriuex_fingerprintPCA/Destriuex_corraperioidc_rest1.csv', header = FALSE)
ap_corr_psd_task=read.csv('~/Documents/CAMCAN_outputs/CAMCAN_destriuex_fingerprintPCA/Destriuex_corraperiodic_rest2.csv', header = FALSE)

# index to reorder desterieux atlas bc brainstorm has chnaged the order in 2024
index_dest= read.csv('~/Documents/CAMCAN_outputs/INDEX_Desteriuex.csv', header = FALSE)

# reorder CamCAN matrix
atlas2=read.csv('/Users/jason/Desktop/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

# save aperiodic corrected data
list=index_dest$V1
roi_index= c(mapply(seq, 79*(list-1)+7, (79*(list-1)+79)))

ap_corr_psd_train= ap_corr_psd_task[,roi_index]
ap_corr_psd_valid= ap_corr_psd_rest[,roi_index]

demoSickkids= read.csv('/Users/jason/Documents/SickKids/Data2Move/ToMove_demo.csv')
SickKid_ap_corr_psd_valid= read.csv('/Users/jason/Documents/SickKids/Data2Move/ToMove_corrspecValidation.csv')
SickKid_ap_corr_psd_train= read.csv('/Users/jason/Documents/SickKids/Data2Move/ToMove_corrspecTraining.csv')

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

# combine datasets
training= rbind(ap_corr_psd_train, SickKid_ap_corr_psd_train)
validation= rbind(ap_corr_psd_valid, SickKid_ap_corr_psd_valid)

Ages= c(demo$Age, demoSickkids$Age)

rm(SickKid_ap_corr_psd_train, SickKid_ap_corr_psd_valid, ap_corr_psd_task, ap_corr_psd_rest)

# Compute ICC from aperiodic corrected data 

age_order= order(Ages)

temp= as.matrix(array(, c(37, 148)))
age_mean= array(, c(1,37))

atlas=read.csv('/Users/jason/Desktop/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')


for (i in 1:37){
  
  #start= 25*(i-1) # sliding window of 50 people with 50% overlap
  #end= (25+(i*25))-1
  #if(i==1){
  #  start=1
  #} 39!!!
  start= 25*(i-1) # sliding window of 100 with 75% overlap
  end= (75+(i*25))-1
  if(i==1){
    start=1
  }
  indsub= age_order[start:end]
  z_target= t(scale(t(training[indsub,])))
  z_database= t(scale(t(validation[indsub,])))
  age_mean[i]= mean(Ages[indsub])
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
  
  temp[i,] = rowMeans(data.frame(theta= rowMeans(icctemp[,1:8]),
                                 alpha= rowMeans(icctemp[,9:18]),
                                 beta= rowMeans(icctemp[,19:52]), 
                                 gamma= rowMeans(icctemp[,53:73])))

}

# look at alignmnent to FC acorss all age bins
tt=apply(temp, 1, cor.test, atlas$fcgradient01)
neg=unlist(tt)


# plot effects
data4plot= data.frame(value=c(neg[seq(4,370,10)]),
                      upper=c(neg[seq(10,370,10)]),
                      lower= c(neg[seq(9,370,10)]),
                      group= c(age_mean))

data4plot$value= as.numeric(data4plot$value)
data4plot$upper= as.numeric(data4plot$upper)
data4plot$lower= as.numeric(data4plot$lower)

ggplot(data4plot, aes(x=group, y=value, color= 'fill', fill='fill')) + geom_point(size=5) +
  stat_smooth(method = "glm",formula = y ~ poly(x, 3),fullrange = T)+ 
  ggpubr::theme_classic2() + scale_colour_manual(values=c('#268236', '#e4e000')) +
  scale_fill_manual(values=c('#268236', '#e4e000')) +
  labs(y="Pearson correlation", x="age",title="") + coord_cartesian(ylim=c(-0.5, 0.1)) +
  theme(text = element_text(size=40), legend.title=element_blank())

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_FCGRAD01_continuous_mergedCohorts_POLY3rd.pdf', device = "pdf", width= 10, height =5)


# run non linear 3rd order model
lm0=lm(DescTools::FisherZ(value) ~1, data4plot)
lm1=lm(DescTools::FisherZ(value) ~poly(group, 1), data4plot)
lm2=lm(DescTools::FisherZ(value) ~poly(group, 2), data4plot)
lm3=lm(DescTools::FisherZ(value) ~poly(group, 3), data4plot)

anova(lm0,lm1,lm2, lm3)
summary(lm3)
sjPlot::tab_model(lm3)

# spatial autocorrelation permutation tests
permuted_index= read.csv('~/Documents/SickKids/abagen_analysis/csv_data4MATLAB/permuted_indexes_of_destriuex_atlas_SPINs&Twirl.csv')
permuted_index= permuted_index[,-1]
permuted_index=permuted_index+1

index= read.csv('/Users/jason/Documents/HCP_twin_projec/abagen_analysis/Destiruex_abagen/index_destriuex.csv', header = FALSE)

beta1=c()
beta2=c()
beta3=c()
rsq= c()
FCGrad01= atlas$fcgradient01[order(index$V1)]

for (i in 1:1000){
  
  ttt=apply(temp, 1, cor.test, FCGrad01[permuted_index[,i]])
  neg=unlist(ttt)
  
  # plot effects
  tempdata4plot= data.frame(value=c(neg[seq(4,370,10)]),
                            upper=c(neg[seq(10,370,10)]),
                            lower= c(neg[seq(9,370,10)]),
                            group= c(age_mean))
  tempdata4plot$value= as.numeric(tempdata4plot$value)
  
  lm3=lm(DescTools::FisherZ(value) ~poly(group, 3), tempdata4plot)
  
  
  beta1[i]=lm3$coefficients[2]
  beta2[i]=lm3$coefficients[3]
  beta3[i]=lm3$coefficients[4]
  rsq[i]=summary(lm3)$r.squared
  
}

(sum(-0.46> beta3)+1)/1001 # 0.00399
