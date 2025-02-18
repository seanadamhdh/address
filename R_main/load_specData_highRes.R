

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
  ##################################
  # X-SCALING: conversion Lambda365-abs to spctrolyzer-abs... empiric
  ##################################
  X_scaling=230.3
  ###################################
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

# Spectrolyzers in dataset
# SerialNo 21470203
# SerialNo 22130200 <- used to be in Sosa 


##############################

# If wavelength range is < max, no of cols remain the same but are filled with nan --> keep seq(200,750,2.5), even though data only to 720

##############################
source("./R_main/packages.R")

source("./R_main/Spectrolyzer_load_good.R") #updated
if(F){ #takes a while
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
#reload
all_mine_spc=readRDS(paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_inclParam"))

###############################################################
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
#' 
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
  ##########  IMPLEMENT MORE ERROR CHECKING / SAFETY FEATURES HERE ###################
  
  # Pre-Processing, if required
  if(str_detect(set,"sg11")){
    X=savitzkyGolay(X,m=0,p=3,w=11)
  }
  if(str_detect(set,"snv")){
    standardNormalVariate(X)
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




############################################################
# laod and process ####

#reload
all_mine_spc=readRDS(paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_inclParam"))

# nesting spc for easier handling 
all_mine_prep=tibble(
  # non spc cols (not nested)
  all_mine_spc%>%select(Date_Time,Serial_No,DOCeq,Flags_DOCeq,TOCeq,Flags_TOCeq,
                     Turbidity,Flags_Turbidity,Temperature,Flags_Temperature,identifier),
  # spc cols
  spc=all_mine_spc%>%select(-c(Date_Time,Serial_No,DOCeq,Flags_DOCeq,TOCeq,Flags_TOCeq,
                            Turbidity,Flags_Turbidity,Temperature,Flags_Temperature,identifier)))


# checking data
# spc rows containing NAs
which(is.na(all_mine_prep$spc)%>%rowSums()>0)
# ...cols
which(is.na(all_mine_prep$spc)%>%colSums()>0)
#'
#' currently 25 measurements (rows) with 
#' only lower end of wavelengths are affected (200-375 nm), total of 71 cols

# rm rows with missing spc values
all_mine=filter(all_mine_prep,rowSums(is.na(spc))==0)



#load model evaluation statistics
model_dir=paste0(data_dir,"/projects/ADDRESS/models/Cubist_2024-02-21/")
eval=readRDS(paste0(model_dir,"Cubist_evaluation"))
# note eval$plots are deprecated
View(eval$eval)


# example predictions 1

variable="Zn_mgL"
best_mod_eval=filter(eval$eval,variable==variable)%>%filter(rmse==min(rmse))

Zn_pred=predict_spectrolyzer(all_mine$spc,variable=variable,
                     trans=best_mod_eval$trans,
                     set = best_mod_eval$set,
                     model_dir = model_dir,
                     prefix = "cubist-auto"
                     )


filter(all_mine,rowSums(spc<=0)==0&rowSums(spc>(4.5*X_scaling))==0)->all_mine_clean


saveRDS(all_mine,paste0(data_dir,"/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_clean"))



# example 2


Y_scaling=1 # for better visibility in plot
var_name="Zn_mgL"
best_mod_eval=filter(eval$eval,variable==var_name)%>%filter(rmse==min(rmse))
trans=best_mod_eval$trans
set=best_mod_eval$set


Y_pred=predict_spectrolyzer(
  X=all_mine_clean$spc/X_scaling,
  model_dir = model_dir,
  variable = var_name,
  #trans="none",
  trans=trans,
  #set="spc"
  set=set
)





out=all_mine_clean
out[[paste0(var_name,"_pred")]]=Y_pred

meta=data.frame(Variable=var_name,
            X_scaling=X_scaling,
            Y_tranformation=trans,
            X_transformation=set,
            model_origin=model_dir,
            notes="X_scaling used to convert from Lambda365 to Spectrolyzer absorbance cm-m and ln-log10. X_transformation=spc_preprocessing.")
write_delim(meta,paste0(data_dir,"projects/ADDRESS/model_out/","Zn_pred_meta.txt"))

tibble(Date_Time=all_mine_clean$Date_Time,Y_pred)%>%write_excel_csv(paste0(data_dir,"projects/ADDRESS/model_out/","Zn_pred.csv"))


ggplot(out,aes(x=Date_Time))+
  geom_line(spc)
  geom_point(aes(y=Zn_mgL_pred),col="red3")+
  theme_minimal()


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



















