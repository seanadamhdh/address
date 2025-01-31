#SerialNo 21470203


##############################

# If wavelength range is < max, no of cols remain the same but are filled with nan --> keep seq(200,750,2.5), even though data only to 720

##############################
source("packages.R")
source("./old imported scripts/Spectrolyzer_load_good.R")

all_mine_spc=Spectro_batch_load(parent_dir = "//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/field_data/data/ADRESS/Spectrolyzer/Spectro_data/",
                                wavelengths = seq(200,750,2.5), #keep as is
                                parameters = c("DOCeq",
                                                "Flags",
                                                "TOCeq",
                                                "Flags",
                                                "Turbidity",
                                                "Flags",
                                                "Temperature"
                                               ),
                                zip = F)
saveRDS(all_mine_spc,"//zfs1.hrz.tu-freiberg.de/fak3ibf/Hydropedo/field_data/data/ADRESS/Spectrolyzer/spectrolyzer_all")








#################################

## load spc data 
# my_spc_data=readRDS(...path-to-prepared-data) or =read_csv ...

some_model=readRDS(paste0("C:/Users/adam/Documents/GitHub/",                                      # path to GitHub dir (alt. download from GitHub)
                          "ADDRESS-adit_drainage_solute_source_control/models/Cubist_2024-02-21/",# dir to models
                          "cubist-auto_spc-log1p-Al_mgL"                                          # specific model
                          )
)



predicitions=predict(some_model,my_spc_data)

# vector














