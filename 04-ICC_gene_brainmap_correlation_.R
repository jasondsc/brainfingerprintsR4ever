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

# read in SickKids data
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

# merge the two datasets into one large one
training= rbind(ap_corr_psd_train, SickKid_ap_corr_psd_train)
validation= rbind(ap_corr_psd_valid, SickKid_ap_corr_psd_valid)

Ages= c(demo$Age, demoSickkids$Age)

rm(SickKid_ap_corr_psd_train, SickKid_ap_corr_psd_valid, ap_corr_psd_task, ap_corr_psd_rest)

# Compute ICC from aperiodic corrected data 
age_order= order(Ages)

temp= as.matrix(array(, c(37, 148)))
age_mean= array(, c(1,37))

atlas=read.csv('/Users/jason/Desktop/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')

# loop through age windows
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
  
  atlas$ICC= temp[i,]
  #atlas$ICC[atlas$ICC >0.8] =0.8
  #atlas$ICC[atlas$ICC <0.3] =0.3

  ggplot(atlas) +
    geom_brain(atlas = desterieux,
               position = position_brain(hemi ~ side),
               aes(fill = ICC)) +  colorspace::scale_fill_continuous_sequential(palette= 'Viridis', rev= FALSE, limits= c(min(atlas$ICC),max(atlas$ICC)))  + theme_void()
  filename=paste('~/Documents/SickKids/figures/ICC/MinMaxNEW_brain_topo_ICC_merged_', round(mean(Ages[indsub])), '_group.pdf', sep='')
  ggsave(filename, device = "pdf")
  
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

# look at alignment to the graident of neg and pos genes across all age windows
tt=apply(temp[,index$V1], 1, cor.test, neg_gene_expression)
neg=unlist(tt)

tt=apply(temp[,index$V1], 1, cor.test, pos_gene_expression)
pos=unlist(tt)

# plot topo of negative and positive genes
atlas=read.csv('/Users/jason/Desktop/Destrieux_neuromaps/Neuromaps_destrieux_full_atlas.csv')
someData= data.frame(region=atlas$region[index$V1], FC= neg_gene_expression, atlas$hemi[index$V1])

someData$FC= scale(someData$FC)[,1]
someData$FC[someData$FC <= 0] = NA
someData$FC[someData$FC > 0] = -3

ggplot(someData) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = FC)) +  colorspace::scale_fill_continuous_sequential('YlOrBr', limits=c(-3, 3 ), rev= FALSE)  +
  theme_void() 

ggsave('~/Documents/CAMCAN_outputs/figures/NEG_genes_threshold_brain.pdf', device = "pdf", width= 5, height =5)

someData= data.frame(region=atlas$region[index$V1], FC= pos_gene_expression, atlas$hemi[index$V1])

someData$FC= scale(someData$FC)[,1]
someData$FC[someData$FC > 0] = 3
someData$FC[someData$FC <= 0] = NA

ggplot(someData) +
  geom_brain(atlas = desterieux, 
             position = position_brain(hemi ~ side),
             aes(fill = FC)) +  colorspace::scale_fill_continuous_sequential('YlOrBr', limits=c(-3, 3 ), rev= FALSE)  +
  theme_void() 

ggsave('~/Documents/CAMCAN_outputs/figures/POS_genes_threshold_brain.pdf', device = "pdf", width= 5, height =5)



# plot effects
data4plot= data.frame(value=c(neg[seq(4,370,10)], pos[seq(4,370,10)]),
                      upper=c(neg[seq(10,370,10)], pos[seq(10,370,10)] ),
                      lower= c(neg[seq(9,370,10)], pos[seq(9,370,10)]),
                      gene= rep(c('negative', 'positive'),each=37),
                      group= c(age_mean, age_mean))

data4plot$value= as.numeric(data4plot$value)
data4plot$upper= as.numeric(data4plot$upper)
data4plot$lower= as.numeric(data4plot$lower)

#write.csv(data4plot, '~/Documents/SickKids/abagen_analysis/Gene_corrs.csv')

data4plot$value[data4plot$gene=='negative']= -1*data4plot$value[data4plot$gene=='negative']
data4plot$valueZ=DescTools::FisherZ(data4plot$value)

ggplot(data4plot, aes(x=group, y=valueZ, group=gene, color = gene, fill=gene)) + geom_point(size=5) +
  stat_smooth(method = "glm",formula = y ~ poly(x, 3),fullrange = F)+
  ggpubr::theme_classic2() + scale_colour_manual(values=c('#9D11C6', '#e4e000')) +
  scale_fill_manual(values=c('#9D11C6', '#e4e000')) +
  labs(y="Pearson correlation", x="age",title="") + coord_cartesian(ylim=c(0, 0.9)) +
  theme(text = element_text(size=40), legend.title=element_blank()) + facet_wrap(~gene, nrow = 2)

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_genecorr_continuous_mergedCohorts_POLY3rd.pdf', device = "pdf", width= 10, height =10)


# plot effects
data4plotNEG= data4plot[data4plot$gene== 'negative',]
data4plotNEG$valueZ= DescTools::FisherZ(data4plotNEG$value)

# fit polynomial model
lm3=lm(DescTools::FisherZ(value) ~poly(group, 3), data4plotNEG)
sjPlot::tab_model(lm3)


ggplot(data4plotNEG, aes(x=group, y=valueZ, group=gene, color = gene, fill=gene)) + geom_point(size=5) + 
  geom_line(data=data4plotNEG, 
            aes(y = pred.quad, x=group), size = 1, col="black") +
  ggpubr::theme_classic2() + scale_colour_manual(values=c('#9D11C6', '#e4e000')) +
  scale_fill_manual(values=c('#9D11C6', '#e4e000')) +
  labs(y="Pearson correlation", x="age",title="") + coord_cartesian(ylim=c(-0.9, 0.0)) +
  theme(text = element_text(size=40), legend.title=element_blank()) + facet_wrap(~gene, nrow = 2)

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_genecorr_continuous_merged_piecewise_NEG.pdf', device = "pdf", width= 10, height =5)




data4plotPOS= data4plot[data4plot$gene== 'positive',]
data4plotPOS$valueZ= DescTools::FisherZ(data4plotPOS$value)

# fit polynomial model
lm3=lm(DescTools::FisherZ(value) ~poly(group, 3), data4plotPOS)
sjPlot::tab_model(lm3)

ggplot(data4plotPOS, aes(x=group, y=valueZ, group=gene, color = gene, fill=gene)) + geom_point(size=5) + 
  geom_line(data=data4plotPOS, 
            aes(y = pred.quad, x=group), size = 1, col="black") +
  ggpubr::theme_classic2() + scale_colour_manual(values=c( '#e4e000')) +
  scale_fill_manual(values=c( '#e4e000')) +
  labs(y="Pearson correlation", x="age",title="") + coord_cartesian(ylim=c(0, 0.9)) +
  theme(text = element_text(size=40), legend.title=element_blank()) 

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_genecorr_continuous_merged_piecewise_POS.pdf', device = "pdf", width= 10, height =5)



# fit polynomial function for negative
lm0=lm(DescTools::FisherZ(value) ~1, data4plotNEG)
lm1=lm(DescTools::FisherZ(value) ~group, data4plotNEG)
lm2=lm(DescTools::FisherZ(value) ~group, data4plotNEG)
lm3=lm(DescTools::FisherZ(value) ~poly(group, 3), data4plotNEG)

anova(lm0,lm1,lm2, lm3)
summary(lm3)
sjPlot::tab_model(lm3)

sjPlot::tab_model(lm1)
sjPlot::tab_model(lm2)

## Run spatial permutations over the entire set of relationships 
# permute map, then run all corrs with ICC, then fit 3rd order poly 

permuted_index= read.csv('~/Documents/SickKids/abagen_analysis/csv_data4MATLAB/permuted_indexes_of_destriuex_atlas_SPINs&Twirl.csv')
permuted_index= permuted_index[,-1]
permuted_index=permuted_index+1

beta1=c()
beta2=c()
beta3=c()
rsq= c()
neg_genes= neg_gene_expression[order(index$V1)]

for (i in 1:1000){
  
  ttt=apply(temp, 1, cor.test, neg_genes[permuted_index[,i]])
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

(sum(-0.82> beta3)+1)/1001 # 0.001998002



data4plotPOS= data4plot[data4plot$gene== 'positive',]
lm0=lm(DescTools::FisherZ(value) ~1, data4plotPOS)
lm1=lm(DescTools::FisherZ(value) ~poly(group, 1), data4plotPOS)
lm2=lm(DescTools::FisherZ(value) ~poly(group, 2), data4plotPOS)
lm3=lm(DescTools::FisherZ(value) ~poly(group, 3), data4plotPOS)

anova(lm0,lm1,lm2, lm3)
summary(lm3)
sjPlot::tab_model(lm3)



permuted_index= read.csv('~/Documents/SickKids/abagen_analysis/csv_data4MATLAB/permuted_indexes_of_destriuex_atlas_SPINs&Twirl.csv')
permuted_index= permuted_index[,-1]
permuted_index=permuted_index+1

beta1=c()
beta2=c()
beta3=c()
rsq= c()
pos_genes= pos_gene_expression[order(index$V1)]

for (i in 1:1000){
  
  ttt=apply(temp, 1, cor.test, pos_genes[permuted_index[,i]])
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

(sum(0.73< beta3)+1)/1001 # 0.00399


############# Plot results of PLS #################
dataPLS= read.csv('~/Documents/SickKids/abagen_analysis/PLS_genes_across_lifespan_2.csv')
dataPLS$group=t(age_mean)

dataPLS$SIG= dataPLS$pvalFDR <0.05

# covarience explained basically stays the same 
ggplot(dataPLS, aes(x=group, y=varexplained, fill='var', colour= 'var')) + geom_point(size=5, aes(alpha= SIG)) + 
  geom_line(aes(x=group, y=varexplained), size=1.3) +
  geom_errorbar(aes(ymin=CIlower,ymax=CIupper),width=0.0, size=1.3) +
  ggpubr::theme_classic2() + scale_colour_manual(values=c('#76D7C4', '#e4e000')) +
  scale_fill_manual(values=c('#76D7C4', '#e4e000')) +
  labs(y="covarience explained", x="age",title="") + coord_cartesian(ylim=c(0, 1)) +
  theme(text = element_text(size=40), legend.title=element_blank(), legend.position = 'none')

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_PLS_var_explained.pdf', device = "pdf", width= 5, height =5)


# covarience explained basically stays the same 
ggplot(dataPLS, aes(x=group, y=pvalFDR, fill='var', colour= 'var')) + geom_point(size=5, aes(alpha= SIG)) + 
  geom_hline(yintercept = 0.05, color= 'Grey', linetype='dashed') +
  ggpubr::theme_classic2() + scale_colour_manual(values=c('#76D7C4', '#e4e000')) +
  scale_fill_manual(values=c('#76D7C4', '#e4e000')) +
  labs(y="p values (FDR)", x="age",title="") + coord_cartesian(ylim=c(0, 0.2)) +
  theme(text = element_text(size=40), legend.title=element_blank(), legend.position = 'none')

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_spun_pvaluesFDR.pdf', device = "pdf", width= 5, height =5)




# covarience explained basically stays the same 
ggplot(dataPLS, aes(x=group, y=overlapGeneLoad, fill='var', colour= 'var')) + geom_point(size=5, aes(alpha= SIG)) + 
  stat_smooth(method = "glm",formula = y ~ poly(x, 2),fullrange = T)+ 
  ggpubr::theme_classic2() + scale_colour_manual(values=c('#76D7C4', '#e4e000')) +
  scale_fill_manual(values=c('#76D7C4', '#e4e000')) +
  labs(y="gene loading similarity", x="age",title="") + coord_cartesian(ylim=c(0, 1)) +
  theme(text = element_text(size=40), legend.title=element_blank(), legend.position = 'none')

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_PLS_GeneLoadings.pdf', device = "pdf", width= 5, height =5)



lm0=lm(DescTools::FisherZ(overlapGeneLoad) ~1, dataPLS)
lm1=lm(DescTools::FisherZ(overlapGeneLoad) ~poly(group, 1), dataPLS)
lm2=lm(DescTools::FisherZ(overlapGeneLoad) ~poly(group, 2), dataPLS)
lm3=lm(DescTools::FisherZ(overlapGeneLoad) ~poly(group, 3), dataPLS)

anova(lm0,lm1,lm2, lm3)
anova(lm0,lm1)
summary(lm2)
sjPlot::tab_model(lm2)

stacked_data=cbind(stack(dataPLS[,7:10]), dataPLS[,11], rep(factor(c('theta', 'alpha', 'beta', 'gamma'), levels = c('theta', 'alpha', 'beta', 'gamma'),), each=37))
colnames(stacked_data)[3:4]= c('group','band')

ggplot(stacked_data, aes(x=group, y=-1*values, fill=band, colour= band)) + geom_point()+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 4), se = TRUE) + 
  #stat_smooth(method = "glm",formula = y ~ poly(x, 2),fullrange = T)+
  ggpubr::theme_classic2() + scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  labs(y="Neurophysiology loadings", x="age",title="") + coord_cartesian(ylim=c(-0.75, 0.75)) + facet_wrap(~band, nrow=1)

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_PLS_NeuroPhysloadings_over_time.pdf', device = "pdf", width= 12, height =4)



cbbPalette= c("#4f89e0",'#4d3d87', '#156605',"#76D7C4",'#4d3d87', '#AC11C6')

stacked_data=cbind(stack(dataPLS[dataPLS$pvalFDR<0.05,7:10]), dataPLS[dataPLS$pvalFDR<0.05,11], rep(factor(c('theta', 'alpha', 'beta', 'gamma'), levels = c('theta', 'alpha', 'beta', 'gamma'),), each=28))
colnames(stacked_data)[3:4]= c('group','band')

ggplot(stacked_data, aes(x=group, y=-1*values, fill=band, colour= band)) + geom_point()+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 4), se = TRUE) + 
  #stat_smooth(method = "glm",formula = y ~ poly(x, 2),fullrange = T)+
  ggpubr::theme_classic2() + scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  labs(y="Neurophysiology loadings", x="age",title="") + coord_cartesian(ylim=c(-0.75, 0.75)) + facet_wrap(~band, nrow=1)

ggsave('~/Documents/CAMCAN_outputs/figures/HundredPpl_Aperiodiccorr_PLS_NeuroPhysloadings_over_timeONLYSIGEFFECTS.pdf', device = "pdf", width= 12, height =4)



