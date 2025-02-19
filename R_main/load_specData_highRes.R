

####################################################################################################################################################### #
# Loading Packages, sourcing code, loading and pre-processing spectra and reference data ####
{
  # depending on OS root is different. Please adjust here if necessary
  
  if(stringr::str_detect(osVersion,"Windows")){
    #workpc/win-sean
    data_dir="//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/"
    code_dir="C:/Users/adam/Documents/GitLab/" #<- set user
  }else if(stringr::str_detect(osVersion,"Ubuntu")){
    #ubuntu                       
    data_dir="/run/user/1000/gvfs/smb-share:server=zfs1.hrz.tu-freiberg.de,share=fak3ibf/Hydropedo/"
    code_dir="/home/hydropedo/Documents/GitLab/" #<- set user
  }else{# e.g. macos
    data_dir=""#...set
    code_dir=""#...set
  }
  ################################## #
  # X-SCALING: conversion Lambda365-abs to spctrolyzer-abs... log10 -> ln conversion and cm-1 -> m-1 
  ################################## #
  X_scaling=230.3
  ################################### #
  # sourcing some scripts from R_main
  source(paste0(code_dir,"/address/R_main/packages.R"))
  source(paste0(code_dir,"/address/R_main/evaluate_model_adjusted.R"))
  source("./R_main/Spectrolyzer_load_good.R") #updated
  
  # simple function for plotting ####
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

# Spectrolyzers in dataset
# SerialNo 21470203
# SerialNo 22130200 <- used to be in Sosa 


############################# #

# If wavelength range is < max, no of cols remain the same but are filled with nan --> keep seq(200,750,2.5), even though data only to 720

############################# #

############################################################## #
#' function for chemometric model application
#' 
#' currenrly not really generaliszed, tailored to existing models
#' 
#' @param X spectra set, tibble with colnames == spc wavelengths
#' @param model_dir path to folders with models
#' @param variable target variable
#' @param trans target variable transformation (e.g. log1p / none)
#' @param set spc preprocessing (e.g., spc / spc_sg11  /spc_sg11_snv)
#' @param prefix currently only cubist-auto, but for generalisation
#' @returns Vector of predictions with length(Yu)=nrow(X). Vector indices correspond to rownumbers of X.
#' 
#' 
predict_spectrolyzer=function(
  X, 
  model_dir="//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/projects/ADDRESS/models/Cubist_2024-02-21/",
  variable="Suva254",
  trans="log1p",
  set="spc",
  prefix="cubist-auto" 
){
  print(paste("Using", paste0(model_dir,"/",prefix,"_",set,"-",trans,"-",variable)))
  #read model
  mod=readRDS(paste0(model_dir,"/",prefix,"_",set,"-",trans,"-",variable))
  
  # check spc match requriements (larger than needed is ok - i think)
  if(min(names(X)%>%as.numeric()) > min(mod$trainingData%>%select(-.outcome)%>%names%>%as.numeric)|
     max(names(X)%>%as.numeric()) < max(mod$trainingData%>%select(-.outcome)%>%names%>%as.numeric)
  ){
    errorCondition("Xu does not cover wavenumber range of Xr")
    return()
  }
  ##########  IMPLEMENT MORE ERROR CHECKING / SAFETY FEATURES HERE ################### #
  # ... if i'll ever have the time...
  
  
  # Pre-Processing, if required
  ## currently only implemented for the preprocessing types used in original calibration runs (default model_dir)
  if(str_detect(set,"sg11")){
    X=savitzkyGolay(X,m=0,p=3,w=11)
  }
  if(str_detect(set,"snv")){
    X=standardNormalVariate(X)
  }
  
  if(trans=="none"){
    Y=predict(mod,X)
    return(Y)
  }else if(trans=="log1p"){
    Y=expm1(predict(mod,X))
    return(Y)
  }else{
    errorCondition("trans error")
    return()
  }
}




########################################################################################################### #
# load and process ####


## load spc from raw data ####
# takes a while, if no changes, use serialized rds data (loaded below)
if(F){ 
  all_mine_spc=Spectro_batch_load(parent_dir = paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/Spectro_data/"),
                                  
                                  wavelengths = seq(200,750,2.5), #keep as is
                                  parameters = c("DOCeq",
                                                 "TOCeq",
                                                 "Turbidity",
                                                 "Temperature"
                                  ),
                                  zip = F,
                                  read_param = T)
  
  
  
  saveRDS(all_mine_spc,paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_inclParam"))
}

## reload raw data ####
all_mine_spc=readRDS(paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_inclParam"))

## nesting spc for easier handling ####
all_mine_prep=tibble(
  # non spc cols (not nested)
  all_mine_spc%>%select(Date_Time,Serial_No,DOCeq,Flags_DOCeq,TOCeq,Flags_TOCeq,
                     Turbidity,Flags_Turbidity,Temperature,Flags_Temperature,identifier),
  # spc cols
  spc=all_mine_spc%>%select(-c(Date_Time,Serial_No,DOCeq,Flags_DOCeq,TOCeq,Flags_TOCeq,
                            Turbidity,Flags_Turbidity,Temperature,Flags_Temperature,identifier)))


## checking and cleaning data ####
# spc rows containing NAs
which(is.na(all_mine_prep$spc)%>%rowSums()>0)
# ...cols
which(is.na(all_mine_prep$spc)%>%colSums()>0)
#'
#' currently 25 measurements (rows) with 
#' only lower end of wavelengths are affected (200-375 nm), total of 71 cols

# rm rows with missing spc values
all_mine=filter(all_mine_prep,rowSums(is.na(spc))==0)

# using same cutoff as for Lambda. Converting to SPectrolyzer abs unit with scaling factor
# aside from bad scans, there are multiple duplicated date_times... zip-folders duplicated?
filter(all_mine,rowSums(spc<=0)==0&rowSums(spc>(4.5*X_scaling))==0)%>%filter(!duplicated(Date_Time))->all_mine_clean

## save clean data ####
saveRDS(all_mine_clean,paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_clean"))

## reload clean data ####
# if available, steps above may be skipped (but are not too slow, so not omitted with if(F){...})
all_mine_clean=readRDS(paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_clean"))

## load model evaluation statistics ####
# set model dir. If no access to fak3/Hydropedo, this must link to local dir
model_dir=paste0(data_dir,"/projects/ADDRESS/models/Cubist_2024-02-21/")

# load model evaluation data (testset validation)
eval=readRDS(paste0(model_dir,"Cubist_evaluation"))


# save eval for Conrad as csv
if(F){
  write_excel_csv(eval$eval,paste0(data_dir,"projects/ADDRESS/models/Cubist_2024_02_21_eval.csv"))
}


# note eval$plots are deprecated and will not plot properly
View(eval$eval)

# Predictions ####
## example predictions 1 (manual) ####
if(F){
variable="Zn_mgL"
best_mod_eval=filter(eval$eval,variable==variable)%>%filter(rmse==min(rmse))

Zn_pred=predict_spectrolyzer(all_mine$spc,variable=variable,
                     trans=best_mod_eval$trans,
                     set = best_mod_eval$set,
                     model_dir = model_dir,
                     prefix = "cubist-auto"
                     )

}



## example 2 (loop for best models and save) ####

#' loop fetches best model for each variable Y, predicts Yu and saves predictions with timestamp as csv,
#' as well as a textfile with relevant metadata about the model params 
for (var_name in unique(eval$eval$variable)){
  
  # find best model (test-eval)
  best_mod_eval=filter(eval$eval,variable==var_name)%>%filter(rmse==min(rmse))
  trans=best_mod_eval$trans
  set=best_mod_eval$set
  
  # apply model / predict Yu
  Y_pred=predict_spectrolyzer(
    X=all_mine_clean$spc/X_scaling,
    model_dir = model_dir,
    variable = var_name,
    #trans="none",
    trans=trans,
    #set="spc"
    set=set
  )
  
  # compile and save metadata (simple df, saved as txt)
  meta=data.frame(Variable=var_name,
              X_scaling=X_scaling,
              Y_tranformation=trans,
              X_transformation=set,
              model_origin=model_dir,
              notes="X_scaling used to convert from Lambda365 to Spectrolyzer absorbance cm-m and ln-log10. X_transformation=spc_preprocessing.")
  write_delim(meta,paste0(data_dir,"projects/ADDRESS/model_out/",var_name,"_meta.txt"))
  
  # add timestamp to Yu pred. Save as csv
  tibble(Date_Time=all_mine_clean$Date_Time,Y_pred)%>%write_excel_csv(paste0(data_dir,"projects/ADDRESS/model_out/",var_name,".csv"))
}


## reload all predictions from example 2 ####
# init with all_mine_clean data
all_pred=all_mine_clean
# loop loads and appends predictions
pb=progress::progress_bar$new(total=length(list.files(paste0(data_dir,"projects/ADDRESS/model_out/"),pattern=".csv",full.names = T)))
for (i in list.files(paste0(data_dir,"projects/ADDRESS/model_out/"),pattern=".csv",full.names = T)){
  pred_i=read_csv(i,col_types = cols(Date_Time = col_datetime(format = "%Y/%m/%d %H:%M:%S")),progress = F)
  names(pred_i)=c("Date_Time",basename(i)%>%str_replace(".csv","_pred"))
  pred_i$Date_Time=(pred_i$Date_Time%>%as.POSIXct())
  #print(head(pred_i)) #debug
  all_pred=left_join(all_pred,pred_i,by="Date_Time")
  pb$tick()
}



# Plotting predicted concentrations ####
# Daily averages to remove noise.
# list of predicted variables
var_list=list.files(paste0(data_dir,"projects/ADDRESS/model_out/"),pattern=".csv")%>%str_replace(".csv","_pred")
var_list_obs=var_list%>%str_remove("_pred")

autosampler_data=read_csv(paste0(data_dir,"/projects/ADDRESS/data/Autosampler_A12_clean.csv"))
# note camp_date==date, redundant col?
manual_data=read_excel(paste0(data_dir,"/projects/ADDRESS/data/Regularsampling_A12_clean.xlsx"))
manual_data$date=as.Date(manual_data$date)

#interactive ggplotly-plot
ggplotly(
all_pred%>%
  group_by(date=date(Date_Time))%>%
  summarise_all(.funs = ~mean(.,na.rm=T))%>%
  #mutate_at(.vars = all_of(var_list),.funs = scale)%>%
  pivot_longer(cols=all_of(var_list))%>%
  ggplot(aes(x=date,y=value,col=str_remove(name,"_pred")))+
  geom_line()+
  geom_point(
    data=autosampler_data%>%
      select(all_of(c("date","site_id",var_list_obs)))%>%
      pivot_longer(cols=var_list_obs),
    aes(x=date,y=value,group=name),shape=3)+
  geom_point(
    data=manual_data%>%
      select(all_of(c("date","site_id",var_list_obs[-7])))%>% # no Durchfluss
      pivot_longer(cols=var_list_obs[-7]),
    aes(x=date,y=value,group=name),shape=13)+
  theme_minimal()+
  ylab("Predicted value [mgL, muS/cm, l/s, abs-ratio,...]")+
  ggtitle("Timeseries of predicted variables")+
  scale_color_discrete("")+
  theme(axis.title.x = element_blank())
  )
# Fe, Al are baaaad... possibly snv problem .. yup bug fixed

# legacy code
# ggplot(all_mine_clean,
#        aes(x=Date_Time,
#            y=spc$`425.0`))+
#   geom_line(linewidth=.1)+
#   geom_line(aes(y=spc$`250.0`),linewidth=.1,col="blue4")+
#   geom_line(data=allmine  spc_data,aes(x=date,y=spc$`250.0`*X_scaling),linewidth=.1,col="red4")+
#   geom_line(data=spc_data,aes(x=date,y=spc$`425.0`*X_scaling),linewidth=.1,col="red2")+
#   geom_point(data=Regular,aes(x=date,y=.data[[var_name]]*Y_scaling,shape=site_id),col="green4")+
#   geom_point(data=Auto,aes(x=as.POSIXct(camp_date),y=.data[[var_name]]*Y_scaling,shape="auto"),col="green4")+
#   scale_shape_manual(values=c(11,4,1))+
#   geom_line(data=tibble(dt=as.POSIXct(all_mine_clean$Date_Time),Y_pred),aes(x=dt,Y_pred*Y_scaling),col="green")+ggtitle(var_name)
# 
# 
# bind_rows(
# tibble(set="spectro",id=all_mine_clean$Date_Time,all_mine_clean$spc/230)%>%
#   pivot_longer(cols = names(all_mine_clean$spc)),
# tibble(set="lambda",id=spc_data_clean$date,spc_data_clean$spc)%>%
#   pivot_longer(cols = names(spc_data_clean$spc))
# )%>%ggplot(aes(x=as.numeric(name),y=value,col=set))+geom_line(linewidth=.1,alpha=.25)












