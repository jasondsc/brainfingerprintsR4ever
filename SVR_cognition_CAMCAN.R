setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

for (j in 1:68){

roi_index= (300*(j-1)+1):(300*(j-1)+300)

target_LSP= target[,roi_index]
database_LSP= database[,roi_index]
target_LSP= (target_LSP+database_LSP)/2
target_LSP=log(target_LSP)
#target_LSP=log(target_LSP+(2*abs(min(target_LSP))))

# USE BETA BAND
target_LSP=target_LSP[, 27:60]


prediction=c()

  for (i in 1:1000){

      ids_pos=sample(nrow(target),605)

      datatrain= target_LSP[ids_pos[1:484], ]
      datatrain$ACER= full_demo$additional_acer[ids_pos[1:484]]
      datatest= target_LSP[ ids_pos[485:605], ]
      datatest$ACER= full_demo$additional_acer[ids_pos[485:605]]


      svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)

      pred <- predict(svm_model,datatest)
      cor(pred, datatest$ACER)
      prediction[i]=cor(pred, datatest$ACER)
  }

predictions= cbind(predictions,prediction)
print(j)
predictions= as.data.frame(predictions)

write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_FROMBETA_80_20_new_log.csv')

}




### random 

setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  target_LSP=log(target_LSP)
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$additional_acer[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$additional_acer[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20.csv')
  
}




## 70 30 

setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  target_LSP=log(target_LSP)
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:424], ]
    datatrain$ACER= full_demo$additional_acer[ids_pos[1:424]]
    datatest= target_LSP[ ids_pos[425:605], ]
    datatest$ACER= full_demo$additional_acer[ids_pos[425:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_70_30.csv')
  
}



### use artifact corrected spectra

setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./psd_rest_residuals.csv', header = FALSE)
psd_task=read.csv('./psd_task_residuals.csv', header = FALSE)

psd_rest= psd_rest[-1,]
psd_task= psd_task[-1,]

psd_rest <- as.data.frame(sapply(psd_rest, as.numeric))
psd_task <- as.data.frame(sapply(psd_task, as.numeric))

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  target_LSP=log(target_LSP+(2*abs(min(target_LSP))))
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$additional_acer[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$additional_acer[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20_residuals_artifacts.csv')
  
}




##### remove 

setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

full_demo$ACER_new= rowSums(full_demo[,39:42])

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  #target_LSP=log(target_LSP)
  target_LSP=log(target_LSP+(2*abs(min(target_LSP))))
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$ACER_new[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$ACER_new[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20_noattention.csv')
  
}



setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

full_demo$ACER_new= rowSums(full_demo[,c(38, 40:42)])

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  #target_LSP=log(target_LSP)
  target_LSP=log(target_LSP+(2*abs(min(target_LSP))))
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$ACER_new[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$ACER_new[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20_nomemory.csv')
  
}


setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

full_demo$ACER_new= rowSums(full_demo[,c(38:39, 41:42)])

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  #target_LSP=log(target_LSP)
  target_LSP=log(target_LSP+(2*abs(min(target_LSP))))
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$ACER_new[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$ACER_new[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20_nofluency.csv')
  
}




setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

full_demo$ACER_new= rowSums(full_demo[,c(38:40, 42)])

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  #target_LSP=log(target_LSP)
  target_LSP=log(target_LSP+(2*abs(min(target_LSP))))
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$ACER_new[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$ACER_new[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20_nolanguage.csv')
  
}




setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./PSD_Rest.csv', header = FALSE)
psd_task=read.csv('./PSD_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

full_demo$ACER_new= rowSums(full_demo[,c(38:41)])

for (j in 1:68){
  
  roi_index= (300*(j-1)+1):(300*(j-1)+300)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  #target_LSP=log(target_LSP)
  target_LSP=log(target_LSP+(2*abs(min(target_LSP))))
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$ACER_new[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$ACER_new[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20_novisuo.csv')
  
}



### aperiodic corrected

setwd('~/Documents/CAMCAN_outputs/NPM/')
demo=read.csv('./standard_data.csv')
demo$CCID=paste('sub',demo$CCID, sep = '_')
ids=read.csv('./subject_codes.csv')
index=demo$CCID %in% ids$Ids
demo=demo[index,]

full_demo= read.csv('./approved_data.tsv', sep = '\t')

index=full_demo$CCID %in% substr(ids$Ids, 5,20)
full_demo=full_demo[index,]

psd_rest=read.csv('./corrspec_Rest.csv', header = FALSE)
psd_task=read.csv('./corrspec_Task.csv', header = FALSE)

database= psd_rest[!is.na(full_demo$additional_acer),]
target= psd_task[!is.na(full_demo$additional_acer),]

full_demo= full_demo[!is.na(full_demo$additional_acer),]

predictions= c()

for (j in 1:68){
  
  roi_index= (78*(j-1)+1):(78*(j-1)+78)
  
  target_LSP= target[,roi_index]
  database_LSP= database[,roi_index]
  target_LSP= (target_LSP+database_LSP)/2
  #target_LSP=log(target_LSP)
  target_LSP=log(target_LSP+(2*abs(min(target_LSP))))
  
  
  prediction=c()
  
  for (i in 1:1000){
    
    ids_pos=sample(nrow(target),605)
    
    datatrain= target_LSP[ids_pos[1:484], ]
    datatrain$ACER= full_demo$additional_acer[ids_pos[1:484]]
    datatest= target_LSP[ ids_pos[485:605], ]
    datatest$ACER= full_demo$additional_acer[ids_pos[485:605]]
    
    
    svm_model=e1071::svm(ACER ~., data= datatrain, kernel = "linear", cost = 1)
    
    pred <- predict(svm_model,datatest)
    cor(pred, datatest$ACER)
    prediction[i]=cor(pred, datatest$ACER)
  }
  
  predictions= cbind(predictions,prediction)
  print(j)
  predictions= as.data.frame(predictions)
  
  write.csv(predictions, '~/Documents/CAMCAN_outputs/decode_SVR_ACER_80_20_APERIODIC_CORR.csv')
  
}

