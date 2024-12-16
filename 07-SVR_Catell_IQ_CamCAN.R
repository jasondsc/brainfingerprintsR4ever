
# load desterieux CamCan data -- corrected spectra
## Use PCA to reduce features and then run SVR
setwd('~/Desktop/Predict_Catel_desterieux/')
demo=read.csv('./demographics_4_Guilia.csv')


psd_rest1=read.csv('./Destriuex_corraperioidc_rest1.csv', header = FALSE)
psd_rest2=read.csv('./Destriuex_corraperiodic_rest2.csv', header = FALSE)

index_d= read.csv('./INDEX_Desteriuex.csv')

psd=0.5*(psd_rest1+psd_rest2)

list=index_d$X17
roi_index= c(mapply(seq, 79*(list-1)+7, (79*(list-1)+79)))

# note to self
# theta 1-8
# alpha 9-18
# beta 19-52
# gamma 53-73

target_LSP= psd[!is.na(demo$Catell_score),roi_index]
demo= demo[!is.na(demo$Catell_score),]
predictions= c()

rm(psd_rest1, psd_rest2)

# run PCA
library(factoextra)
X_scale= scale(target_LSP, center = TRUE, scale = TRUE)
pca_result <- prcomp(X_scale, scale = TRUE, rank=68)

target_LSP= as.data.frame(pca_result$x)

predictions=matrix(, nrow =1000 , ncol =603 )

# run 1000 iterations of SVR and then see if it predicts held out sample 
for (i in 1:1000){
  
  ids_pos=sample(603,603)
  
  datatrain= target_LSP[ids_pos[1:482], ] # use 80% of the data
  datatrain$Catell= demo$Catell_score[ids_pos[1:482]]
  datatest= target_LSP[ ids_pos[483:603], ]
  datatest$Catell= demo$Catell_score[ids_pos[483:603]]

  svm_model=e1071::svm(Catell ~., data= datatrain, kernel = "linear", cost = 1)
  
  pred <- predict(svm_model,datatest)
  predictions[i,ids_pos[483:603]] = pred
  write.csv(predictions, '~/Desktop/Predict_Catel_desterieux/decode_SVR_Catell_80_20_from_all_feat_PCA.csv')
  
}

# correlate predictions for each iteration
apply(predictions, 2, cor.test)

# take col means and see how it predicts subjects acorss all iterations
cor.test(colMeans(predictions, na.rm=TRUE), demo$Catell_score)

corr_pred=c()
for (i in 1:1000) {
  
  temp=predictions[i,]
  indd= which(!is.na(temp))
  corr_pred[i]= cor(temp[indd], demo$Catell_score[indd])
  
}

hist(corr_pred)
