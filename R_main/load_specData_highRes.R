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








#################################

## load spc data 
# my_spc_data=readRDS(...path-to-prepared-data) or =read_csv ...


# model folder
model_dir="//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/projects/ADDRESS/models/Cubist_2024-02-21/"

eval=readRDS(paste0(model_dir,"Cubist_evaluation"))

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
  
  mod=readRDS(paste0(model_dir,"/",prefix,"_",set,"-",trans,"-",variable))
  
  
  if(colnames(X)!=colnames(mod$training_data[,-c(".outcome")])){
    errorCondition("Xu does not match Xr")
    
  }
  
  
}














