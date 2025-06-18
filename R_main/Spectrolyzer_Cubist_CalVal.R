



####################################################################################################################################################### #
# Loading Packages, sourcing code, loading and pre-processing spectra and reference data ####

{
  
  
  if(stringr::str_detect(osVersion,"Windows")){
    #workpc/win-sean
    data_dir="//zfs1.hrz.tu-freiberg.de/fak3biogeochemie/03 Projects - Projekte/ADDRESS/ADDRESS/"
    code_dir="C:/Users/adam/Documents/GitHub/" #<- set user
  }else if(stringr::str_detect(osVersion,"Ubuntu")){
    #ubuntu                       
    data_dir="/run/user/1000/gvfs/smb-share:server=zfs1.hrz.tu-freiberg.de,share=fak3biogeochemie/03 Projects - Projekte/ADDRESS/ADDRESS/"
    code_dir="/home/hydropedo/Documents/GitHub/" #<- set user
  }else{# e.g. macos
    data_dir=""#...set
    code_dir=""#...set
  }
  
  # sourcing some scripts from R_main
  source(paste0(code_dir,"/address/R_main/packages.R"))
  source(paste0(code_dir,"/address/R_main/evaluate_model_adjusted.R"))
  
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
}  
  
  
  ###
  # load datasets ####
  
  ## load spc ####
  spc_data_raw <- read_excel(paste0(data_dir,"/data/processed/allspecoriginal_Oct2023.xlsx"))
  spc_data<-tibble(spc_data_raw[1],spc=spc_data_raw[-1])
  
  ### plot raw spc ####
  plot_spc(spc_data$spc)
  
  ### remove bad scans and plot again ####
  # rm very off spc (outliers, bad scans in aggregates, etc.) 
  filter(spc_data,rowSums(spc<=0)==0&rowSums(spc>4.5)==0)->spc_data_clean
  plot_spc(spc_data_clean$spc)
  
  # No. of removed spc
  cat(
    "\ntotal spectra:\t\t",nrow(spc_data),
    "\nspectra removed:\t",nrow(spc_data)-nrow(spc_data_clean),
    "\nremaining spectra:\t",nrow(spc_data_clean)
  )
  
  ################################## #
  # IMPLEMENTED SOME PREPROCESSING #
  ################################## #
  

  ### smoothing ####
  spc_data_clean$spc_sg11<-
    savitzkyGolay(spc_data_clean$spc,
                  p=3,
                  m=0,
                  w=11)
  plot_spc(spc_data_clean$spc_sg11)
  
  ### SNV ####
  spc_data_clean$spc_sg11_snv<-standardNormalVariate(spc_data_clean$spc_sg11)
  plot_spc(spc_data_clean$spc_sg11_snv)
  
  
  ## load manual samples ####
  Regular <- read_excel(paste0(data_dir,"/data/Regularsampling_A12_clean.xlsx"))
  Regular_TOP <- filter(Regular,!str_detect(Sample_ID,"BH1"))
  Regular_TOP_spc<-inner_join(Regular_TOP,spc_data_clean,by="date")
  
  Regular_BOT <- filter(Regular,str_detect(Sample_ID,"BH1"))
  Regular_BOT_spc<-inner_join(Regular_BOT,spc_data_clean,by="date")
  
  
  
  ## load autosampler ####
  Auto <- read_csv(paste0(data_dir,"/data/Autosampler_A12_clean.csv"))[-1]
  Auto_spc<-inner_join(Auto,spc_data_clean,by="date")
}
####################################################################################################################################################### #


################################################################################################################# #
# Calibration with Auto Sets ####
if(F){ # safety... set T to run !!! Legacy paths do not match
  #select spc set manually (when loop inactive)
  # set<-"spc_sg11_snv"
  
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
  # for all pre-processing techniques run above
  for (set in c("spc","spc_sg11","spc_sg11_snv")){
    # for all (numerical) variables
    for (i in variables){
      # both without and with log(1+p) transformation
      for (trans in c("none","log1p")){
        # print progress
        cat("\n\n\n\n started ",set,"-",trans,"-",i,"\n\n")
        
        # subset current variable + spc set
        tmp<-na.omit(Auto_spc[,c(i,set)])
        
        # if TRUE, log transformation
        if(trans=="log1p"){ 
          tmp[[i]]<-log(tmp[[i]]+1)
        }
        
        # 75 % train, 25 % test, based on variable
        inTrain<-createDataPartition(tmp[[i]],p = .75)[[1]]
        
        train<-tmp[inTrain,]
        test<-tmp[-inTrain,]
        
        # train Cubist models on tuning grid
        out<-train(
          x=train[[set]],
          y=train[[i]],
          method="cubist",
          tuneGrid=trainGrid,
          trControl=fitControl
        )
        
        # appending test-set and metadata to model object for easier evaluation and documentation
        out$testingData<-test
        out$partition<-inTrain
        out$documentation<-list(
          variable=i,
          trans=trans,
          set=set,
          n_train=nrow(train),
          n_test=nrow(test),
          test_eval_finalModel=
            # !!!currently no back-transformation here, maybe fix at some point !!!
            evaluate_model_adjusted(data.frame(obs=test[[i]],pred=predict(out$finalModel,test[[set]])),obs="obs",pred="pred")
          #case_when(trans=="none"~evaluate_model_adjusted(data.frame(obs=test[[i]],pred=predict(out$finalModel,test[[set]])),obs="obs",pred="pred"),
          #          trans=="log1p"~evaluate_model_adjusted(data.frame(obs=exp(test[[i]])-1,pred=exp(predict(out$finalModel,test[[set]]))-1),obs="obs",pred="pred")
          #),
          #note="test_eval was carried out with de-logged predictions"
        )
        
        # saving output to /R_main/temp/ as rds. Named cubist-auto_#spc_set#-#trans#-#variable#
        ############ LEGACY DO NOT RUN   ############ #
        saveRDS(out,paste0(root_dir,"/GitHub/ADDRESS-adit_drainage_solute_source_control/R_main/temp/cubist-auto_",set,"-",trans,"-",i))
        cat("\n\n finished ",set,"-",trans,"-",i,"\n\n")
        print(out$documentation$test_eval_finalModel)
      }
      
    }
    
  }
}
####################################################################################################################################################### #

########################################################################################################################################################## #
#  Model evaluation using test sets ####
# init
if(F){ #safety, set to T if you want to rerun this chunk
  plotlist<-list()
  cubist_model_eval_table<-c()
  
  # set folder name
  model_folder<-"/models/Cubist_2024-02-21"
  
  # loop for all models in dir
  for (i in list.files(paste0(root_dir,"/GitHub/ADDRESS-adit_drainage_solute_source_control",model_folder),full.names = T,
                       pattern = "cubist") # small "c" cubist is used as pattern for model class objects in /temp -> results etc should be named Cubist with capital "C"
  ){
    model<-read_rds(i)
    model_name<-basename(i)
    
    # calculating obs / pred based on transformation and spc-set. Un-transforms log(1+p) data
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
      ggtitle(model_name)+
      theme_pubr()->plotlist[[model_name]]
    
    
    #calculating stats
    cubist_model_eval_table<-bind_rows(cubist_model_eval_table,
                                       data.frame(set=set,
                                                  trans=model$documentation$trans,
                                                  variable=basename(i)%>%strsplit("-")%>%last%>%last,
                                                  data.frame(evaluate_model_adjusted(test_ObsPred,obs="obs",pred="pred"))))
  }
  
  # aggregating
  list(eval=cubist_model_eval_table,plots=plotlist)->Cubist_test_evaluation
  
  # save results 
  saveRDS(Cubist_test_evaluation,paste0(root_dir,"/GitHub/ADDRESS-adit_drainage_solute_source_control",model_folder,"/Cubist_evaluation"))
}

  
########################################################################################################################################################## #
# displaying results ####  
# Assuming root_dir and model_folder variables are defined and hold the desired path

## re-loading results ####
# ...originally created in chunk above
Cubist_test_evaluation <- readRDS(paste0(root_dir,"/GitHub/ADDRESS-adit_drainage_solute_source_control",model_folder,"/Cubist_evaluation"))

# check  change in case it is different from above
# set folder name
model_folder<-"/models/Cubist_2024-02-21"

# available variables
Cubist_test_evaluation$eval$variable%>%unique
  ## plot timeseries for variable:___ ####
  var_sel<-"Cd_mgL"
  # best combination of spc-set and transformation is choosen based on lowest RMSEP (test-data)
  
  var_mods<-filter(Cubist_test_evaluation$eval,variable==var_sel)
  best_mod<-var_mods[which.min(var_mods$rmse),]
  mod<-read_rds(paste0(root_dir,"/GitHub/ADDRESS-adit_drainage_solute_source_control",model_folder,"/cubist-auto_",
                       best_mod$set,
                       "-",
                       best_mod$trans,
                       "-",
                       var_sel))
  
  
  
  
  ## main plot 1 - timeseries ####
  ggplotly(
    {
      if(best_mod$trans=="log1p"){
        predictions<-data.frame(Auto_spc,
                                pred=exp(predict(mod,Auto_spc[[best_mod$set]]))-1
        )
      } else {
        predictions<-data.frame(Auto_spc,
                                pred=predict(mod,Auto_spc[[best_mod$set]])
        )
      }
      
      # !!! there should not be NAs except in var_sel !!!
      ref_data<-na.omit(select(Auto_spc,all_of(c("Sample_ID","campaign","site_id","date",var_sel))))
      
      ggplot(data=ref_data,aes(x=as.POSIXct(date),
                               y=.data[[var_sel]]))+
        geom_point(aes(size="training",
                       shape="training",
                       color="training"
        ))+
        geom_point(data=ref_data[-c(mod$partition),],
                   aes(size="testing",
                       shape="testing",
                       color="testing"
                   ),stroke=1)+
        #  geom_line(linewidth=.2,alpha=.5)+
        geom_line(data=predictions,aes(
          x=date,
          y=pred))+
        ggtitle(paste0("cubist-auto_",
                       best_mod$set,
                       "-",
                       best_mod$trans,
                       "-",
                       var_sel))+
        ylab(var_sel)+
        xlab("date")+
        ylab("Cd (mg/L)")+
        scale_size_manual("Set",breaks=c("training","testing"),values=c(.5,2))+
        scale_shape_manual("Set",breaks=c("training","testing"),values=c(16,3))+
        scale_color_manual("Set",breaks=c("training","testing"),values=c("black","red3"))+
        theme_pubr()
    }
  )
  
## main plot 2 - timeseries cropped ####
  #modifing time frame   
  ggplotly(
    {
      if(best_mod$trans=="log1p"){
        predictions<-data.frame(Auto_spc,
                                pred=exp(predict(mod,Auto_spc[[best_mod$set]]))-1
        )
      } else {
        predictions<-data.frame(Auto_spc,
                                pred=predict(mod,Auto_spc[[best_mod$set]])
        )
      }
      
      # !!! there should not be NAs except in var_sel !!!
      ref_data<-na.omit(select(Auto_spc,all_of(c("Sample_ID","campaign","site_id","date",var_sel))))
      
      ggplot(data=ref_data,aes(x=as.POSIXct(date),
                               y=.data[[var_sel]]))+
        geom_point(aes(size="training",
                       shape="training",
                       color="training"
        ))+
        geom_point(data=ref_data[-c(mod$partition),],
                   aes(size="testing",
                       shape="testing",
                       color="testing"
                   ),stroke=1)+
        #  geom_line(linewidth=.2,alpha=.5)+
        geom_line(data=predictions,aes(
          x=date,
          y=pred))+
        ggtitle(paste0("cubist-auto_",
                       best_mod$set,
                       "-",
                       best_mod$trans,
                       "-",
                       var_sel))+
        ylab(var_sel)+
        xlab("date")+
        ylab(var_sel)+
        scale_size_manual("Set",breaks=c("training","testing"),values=c(.5,2))+
        scale_shape_manual("Set",breaks=c("training","testing"),values=c(16,3))+
        scale_color_manual("Set",breaks=c("training","testing"),values=c("black","red3"))+
        scale_x_datetime(limits = as.POSIXct(c("2022-05-16", "2022-09-20"))) +
        theme_pubr()
    }
  )
  
  ## main plot 3 - spectral importance ####
  # spectral importance plots
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
      scale_y_continuous(sec.axis = sec_axis("Usage",transform = ~.*25))
  })

####################################################################################################################################################### #

## plot 4 - raw model testing obs vs predictions ####
plot(mod$testingData[[var_sel]],predict(mod,mod$testingData[[best_mod$set]]))
 abline(0,1)












