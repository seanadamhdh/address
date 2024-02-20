
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
# IMPLEMENT SOME PREPROCESSING #
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

#select spc set (when loop inactive)
set<-"spc_sg11_snv"

# Variable selection for model fitting
variables<-c(
  "doc_mgL",
  "Fe_mgL"
)
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
      out$documentation<-list(
        variable=i,
        trans=trans,
        n_train=nrow(train),
        n_test=nrow(test),
        test_eval_finalModel=
          evaluate_model_adjusted(data.frame(obs=test[[i]],pred=predict(out$finalModel,test[[set]])),obs="obs",pred="pred")
          #case_when(trans=="none"~evaluate_model_adjusted(data.frame(obs=test[[i]],pred=predict(out$finalModel,test[[set]])),obs="obs",pred="pred"),
          #          trans=="log1p"~evaluate_model_adjusted(data.frame(obs=exp(test[[i]])-1,pred=exp(predict(out$finalModel,test[[set]]))-1),obs="obs",pred="pred")
          #),
        #note="test_eval was carried out with de-logged predictions"
      )
      
      
      
      
      saveRDS(out,paste0("~/Documents/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/temp/cubist-auto_",set,"-",trans,"-",i))
      break()
    }
    
  }
  
}



























