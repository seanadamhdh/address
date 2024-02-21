
source("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/packages.R")
source("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/evaluate_model_adjusted.R")

#simple function for plotting
plot_spc<-function(spc){
  matplot(x=spc%>%colnames(),
          y=t(spc),
          type="l",
          col=rgb(0,0,0,.1),
          lty="solid",
          #ylim=c(-1,10),
          xlab="wavelength [nm]",
          ylab="absorbance [freedom /sq inch]") #look unit up
}


# load datasets ####

# load spc
spc_data_raw <- read_excel("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/data/processed/allspecoriginal_Oct2023.xlsx")
spc_data<-tibble(spc_data_raw[1],spc=spc_data_raw[-1])

# plot raw spc
plot_spc(spc_data$spc)

# rm very off spc (outliers, bad scans in aggregates, etc.)
filter(spc_data,rowSums(spc<=0)==0&rowSums(spc>4.5)==0)->spc_data_clean
plot_spc(spc_data_clean$spc)

# No. of removed spc
cat(
  "\ntotal spectra:\t\t",nrow(spc_data),
  "\nspectra removed:\t",nrow(spc_data)-nrow(spc_data_clean),
  "\nremaining spectra:\t",nrow(spc_data_clean)
)

################################
# IMPLEMENTED SOME PREPROCESSING #
################################

# smoothing
spc_data_clean$spc_sg11<-
  savitzkyGolay(spc_data_clean$spc,
                p=3,
                m=0,
                w=11)
plot_spc(spc_data_clean$spc_sg11)

# SNV
spc_data_clean$spc_sg11_snv<-standardNormalVariate(spc_data_clean$spc_sg11)
plot_spc(spc_data_clean$spc_sg11_snv)


# load manual samples
Regular <- read_excel("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/data/Regularsampling_A12_clean.xlsx")
Regular_TOP <- filter(Regular,!str_detect(Sample_ID,"BH1"))
Regular_TOP_spc<-inner_join(Regular_TOP,spc_data_clean,by="date")

Regular_BOT <- filter(Regular,str_detect(Sample_ID,"BH1"))
Regular_BOT_spc<-inner_join(Regular_BOT,spc_data_clean,by="date")



# load autosampler
Auto <- read_csv("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/data/Autosampler_A12_clean.csv")[-1]
Auto_spc<-inner_join(Auto,spc_data_clean,by="date")




######################################
# Calibration with Auto Sets ####
{
  #select spc set (when loop inactive)
  set<-"spc_sg11_snv"
  
  # Variable selection for model fitting
  variables<-names(Auto)[-c(1:4,7,10)]
  
  # Cubist prarams
  trainGrid<-expand.grid(committees=c(1,2,5,10,20,50),
                         neighbors=c(0:9))
  
  fitControl <- trainControl(
    ## 10-fold CV
    method = "cv",
    number = 10,
    ## print progress
    verboseIter = TRUE)
  
  # initialize loop
  for (set in c("spc","spc_sg11","spc_sg11_snv")){
    for (i in variables){
      for (trans in c("none","log1p")){
        cat("\n\n\n\n started ",set,"-",trans,"-",i,"\n\n")
        tmp<-na.omit(Auto_spc[,c(i,set)])
        
        if(trans=="log1p"){ 
          tmp[[i]]<-log(tmp[[i]]+1)
        }
        
        # 75 % train, 25 % test, based on variable
        inTrain<-createDataPartition(tmp[[i]],p = .75)[[1]]
        train<-tmp[inTrain,]
        
        test<-tmp[-inTrain,]
        
        
        out<-train(
          x=train[[set]],
          y=train[[i]],
          method="cubist",
          tuneGrid=trainGrid,
          trControl=fitControl
        )
        # appending test-set and metadata
        out$testingData<-test
        out$partition<-inTrain
        out$documentation<-list(
          variable=i,
          trans=trans,
          set=set,
          n_train=nrow(train),
          n_test=nrow(test),
          test_eval_finalModel=
            evaluate_model_adjusted(data.frame(obs=test[[i]],pred=predict(out$finalModel,test[[set]])),obs="obs",pred="pred")
          #case_when(trans=="none"~evaluate_model_adjusted(data.frame(obs=test[[i]],pred=predict(out$finalModel,test[[set]])),obs="obs",pred="pred"),
          #          trans=="log1p"~evaluate_model_adjusted(data.frame(obs=exp(test[[i]])-1,pred=exp(predict(out$finalModel,test[[set]]))-1),obs="obs",pred="pred")
          #),
          #note="test_eval was carried out with de-logged predictions"
        )
        
        # saving output to /R_main/temp/ as rds. Named cubist-auto_#spc_set#-#trans#-#variable#
        saveRDS(out,paste0("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/temp/cubist-auto_",set,"-",trans,"-",i))
        cat("\n\n finished ",set,"-",trans,"-",i,"\n\n")
        print(out$documentation$test_eval_finalModel)
      }
      
    }
    
  }
}

# aggregating evaluation for test sets ####
plotlist<-list()
cubist_model_eval_table<-c()
for (i in list.files("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/temp/",full.names = T)){
  model<-read_rds(i)
  model_name<-basename(i)
  # calculating obs / pred
  if (model_name%>%str_detect("spc-")){
    if(model$documentation$trans=="none"){
      test_ObsPred<-data.frame(pred=predict(model,model$testingData$spc),obs=model$testingData[[1]])
    }else if (model$documentation$trans=="log1p"){
      test_ObsPred<-data.frame(pred=exp(predict(model,model$testingData$spc))-1,obs=exp(model$testingData[[1]])-1)
    }
    set<-"spc"
  }else if (model_name%>%str_detect("spc_sg11-")){
    if(model$documentation$trans=="none"){
      test_ObsPred<-data.frame(pred=predict(model,model$testingData$spc_sg11),obs=model$testingData[[1]])
    }else if (model$documentation$tr=="log1p"){
      test_ObsPred<-data.frame(pred=exp(predict(model,model$testingData$spc_sg11))-1,obs=exp(model$testingData[[1]])-1)
    }
    set<-"spc_sg11"
  }else if (model_name%>%str_detect("spc_sg11_snv-")){
    if(model$documentation$trans=="none"){
      test_ObsPred<-data.frame(pred=predict(model,model$testingData$spc_sg11_snv),obs=model$testingData[[1]])
    }else if (model$documentation$trans=="log1p"){
      test_ObsPred<-data.frame(pred=exp(predict(model,model$testingData$spc_sg11_snv))-1,obs=exp(model$testingData[[1]])-1)
    }
    set<-"spc_sg11_snv"
  }
  #plotting
  ggplot(test_ObsPred,aes(x=obs,y=pred))+
    geom_point()+
    geom_abline(intercept=0,slope=1,linetype="dotted")+
    ggtitle(basename(i))+
    theme_pubr()->plotlist[[basename(i)]]
  
  
  #calculating stats
  cubist_model_eval_table<-bind_rows(cubist_model_eval_table,
                                     data.frame(set=set,
                                                trans=model$documentation$trans,
                                                variable=basename(i)%>%strsplit("-")%>%last%>%last,
                                                data.frame(evaluate_model_adjusted(test_ObsPred,obs="obs",pred="pred"))))
}


list(eval=cubist_model_eval_table,plots=plotlist)->Cubist_test_evaluation
#saveRDS(Cubist_test_evaluation,"~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/temp/Cubist_evaluation")
# choose ther dir. code breaks with other objects in /temp

var_sel<-"doc_mgL"
ggplotly(
{
var_mods<-filter(Cubist_test_evaluation$eval,variable==var_sel)
best_mod<-var_mods[which.min(var_mods$rmse),]
mod<-read_rds(paste0("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/temp/cubist-auto_",
                     best_mod$set,
                     "-",
                     best_mod$trans,
                     "-",
                     var_sel))


if(best_mod$trans=="log1p"){
predictions_train<-data.frame(date=Auto_spc$date,
                        pred=exp(predict(mod,Auto_spc[[best_mod$set]]))-1

  )
} else {
    predictions<-data.frame(date=Auto_spc$date,
                            pred=predict(mod,Auto_spc[[best_mod$set]])
    )
}


  ggplot(data=Auto,aes(x=as.POSIXct(date),y=.data[[var_sel]]))+
    geom_point(size=.1)+
  geom_point(data=Auto[-c(mod$partition),],shape=3,col="red")+
#  geom_line(linewidth=.2,alpha=.5)+
  geom_line(data=predictions,aes(x=date,y=pred))+
  ggtitle(paste0("cubist-auto_",
                     best_mod$set,
                     "-",
                     best_mod$trans,
                     "-",
                     var_sel))+
  ylab(var_sel)+
  xlab("date")+
  theme_pubr()
}
)




### spectral importance plots
ggplotly({
mod$trainingData%>%
  select(-c(`.outcome`))%>%
  summarise_all(mean)%>%
  t%>%
  data.frame%>%
  rownames_to_column->avg_spc
names(avg_spc)<-c("nm","avg_absorbance")

ggplot(avg_spc,aes(x=as.numeric(nm),
                   y=avg_absorbance
                   ))+
  geom_line()+
  geom_col(data=mod$finalModel$usage,aes(x=as.numeric(Variable),
                                         y=Conditions/25,
                                         fill="Conditions"))+
  geom_col(data=mod$finalModel$usage,aes(x=as.numeric(Variable),
                                         y=Model/25,
                                         fill="Model"))+
  ggtitle(paste0("cubist-auto_",
                 best_mod$set,
                 "-",
                 best_mod$trans,
                 "-",
                 var_sel))+
  theme_pubr()+
  scale_y_continuous(sec.axis = sec_axis("Usage",trans = ~.*25))
})

















