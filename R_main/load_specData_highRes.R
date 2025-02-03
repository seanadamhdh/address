#SerialNo 21470203


##############################

# If wavelength range is < max, no of cols remain the same but are filled with nan --> keep seq(200,750,2.5), even though data only to 720

##############################
source("./R_main/packages.R")
source("./R_main/old imported scripts/Spectrolyzer_load_good.R")

all_mine_spc=Spectro_batch_load(parent_dir = "//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/field_data/data/ADDRESS/Spectrolyzer/Spectro_data/",
                                wavelengths = seq(200,750,2.5), #keep as is
                                parameters = c("DOCeq",
                                                "TOCeq",
                                                "Turbidity",
                                                "Temperature"
                                               ),
                                zip = F,
                                read_param = T)

saveRDS(all_mine_spc,"//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_inclParam")




###############################################################
# Funktion zur Anwendung der cubist Modelle auf die Spektrolyzer Daten

predict_spectrolyzer=function(
    # spectra set, tibble with colnames == spc wavelengths, 
  X, 
  # path to folders with models
  model_dir="//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/projects/ADDRESS/models/Cubist_2024-02-21/",
  
  # target variable
  variable="Suva254",
  
  # target varaible transformation (log1p / none)
  trans="log1p",
  
  # spc preprocessing (spc / spc_sg11)
  set="spc",
  
  prefix="cubist-auto" # currently only cubist-auto, but for generalisation
  
){
  print(paste("Using", paste0(model_dir,"/",prefix,"_",set,"-",trans,"-",variable)))
  #read model
  mod=readRDS(paste0(model_dir,"/",prefix,"_",set,"-",trans,"-",variable))
  
  # check spc match requriements
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
all_mine_spc=readRDS("//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/field_data/data/ADDRESS/Spectrolyzer/spectrolyzer_all_inclParam")

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
model_dir="//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/projects/ADDRESS/models/Cubist_2024-02-21/"
eval=readRDS(paste0(model_dir,"Cubist_evaluation"))
# note eval$plots are deprecated
View(eval$eval)


# example predictions

predict()






















