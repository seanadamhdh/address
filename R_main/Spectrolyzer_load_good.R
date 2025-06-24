

# INFO ####
#' Functions for loading raw spectrolyzer data
#'

# lazy loading of packages (make sure package.R is in wd)
#source("./R_main/packages.R")
#------------------------------------------------------------------------------------------------------------------------------------------- -
#' @title generalised fingerprint load-in
#' @description
#' NOTE: Removes flagged spectra.
#' 
#' @param directory filepath of fingerprint-folder
#' @param wavelengths col-names of measured wavelengths. default is 200-750 in 2.5 nm, steps
#' @param sel_wavelengths Range of wavelengths that shall be returned
#' @param id Should serial No be added as column to the tibble (for ID-ing the spectrolyzer later). Default=T
Spectro_load_fingerprint<-function(directory,
                                   wavelengths=seq(200,750,2.5),
                                   sel_wavelengths=seq(200,722.5,2.5),
                                   id=T){
  # load file
  temp<-suppressMessages(#silent read-in
    read_csv(paste0(directory,
                        "/",
                        list.files(directory,pattern = ".csv") #actual filename
                        ),
                 col_types = cols(Timestamp = col_datetime(format = "%Y-%m-%dT%H:%M:%S+0000")))
  )
  # rename cols
  names(temp)<-c("Date_Time",format(wavelengths,nsmall=1),"Flags")
  # POSIXct format
  temp$Date_Time<-as.POSIXct(temp$Date_Time)
  
  # remove flagged cols (e.g. "NO MEDIUM" Error, but also other possible error-rows)
  subset(temp,is.na(Flags))%>%
    #' rm Flags (did it's job)
    select(-c("Flags"))%>%
    #' rm nan cols (no data recorded for these wavelengths) 
    #' using first() because i guess it's faster
    #' if issues come up replace with all()
    select_if(~is.nan(first(.))!=T)%>%
    # select desired range
    select(all_of(c("Date_Time",format(sel_wavelengths,nsmall=1))))->temp
  
  if(id==T){
    Serial_No<-read.delim(paste0(directory,
                                 "/",
                                 list.files(directory,pattern = ".json") #actual filename
                                 ), 
                          header=FALSE)[1,1]%>%substr(str_locate(.,"sensorSerial:")[2]+1, # start of serialNo.
                                                      str_locate(.,"sensorSerial:")[2]+8) # end of serialNo, assuming 8 digits
    temp$Serial_No<-rep(Serial_No,length(temp[[1]]))
  }
  return(temp)
}
#------------------------------------------------------------------------------------------------------------------------------------------- -
#' generalised parameter load-in
#' @param directory filepath of parameter folder
#' @param parameters list of recorded parameters, forming col-names (Default is Temperature only). Must match order in csv. 
#' e.g. ... DOCeq | Flags | TOCeq | Flags would require c("DOCeq","TOCeq")
#' @param id Should serial No be added as column to the tibble (for ID-ing the spectrolyzer later). Default=T
Spectro_load_parameter<-function(directory,
                                 parameters=c("Temp"),
                                 id=T){
  # load file
  temp<-read_csv(paste0(directory,
                 "/",
                 list.files(directory,pattern = ".csv") #actual filename
                 ),
                 col_types = cols(Timestamp = col_datetime(format = "%Y-%m-%dT%H:%M:%S+0000"),
                                  Flags = col_character()))
  # rename cols
  # fix for multiple Flag cols
  namelist=c("Date_Time")
  for (i in parameters){
    namelist=c(namelist,i,paste0("Flags_",i))
  }
  
  print(namelist)
  names(temp)<-namelist
  
  # POSIXct format
  temp$Date_Time<-as.POSIXct(temp$Date_Time)
  return(temp)
}

#### LEGACY CODE
#' used to remove flagged parameter values... Sosa spectrolyzer only recorded Temperature... easy. Mine spectrolyzer has multiple columns with individual flags.
#' Rewrote code above to give unique flag names. No flag columns are handed on to batch load. Easy solution and makes sense to keep flags for future filtering instead of 
#' dumping all rows with any flags.
  #' # remove flagged cols (e.g. "NO MEDIUM" Error, but also other possible error-rows)
  #' temp%>%filter(across(contains("Flags"), ~ is.na(.)) %>% rowSums() == ncol(select(., contains("Flags"))))%>%
  #'   #' rm Flags (did it's job)
  #'   select(-contains("Flags"))%>%
  #'   #' rm nan cols (no data recorded for these wavelengths) 
  #'   #' using first() because i guess it's faster
  #'   #' if issues come up replace with all()
  #'   select_if(~is.nan(first(.))!=T)->temp
  #' if(id==T){
  #' Serial_No<-read.delim(paste0(directory,
  #'                              "/",
  #'                              list.files(directory,pattern = ".csv") #actual filename
  #'                              ), 
  #'                       header=FALSE)[1,1]%>%substr(str_locate(.,"sensorSerial:")[2]+1, # start of serialNo.
  #'                                                   str_locate(.,"sensorSerial:")[2]+8) # end of serialNo, assuming 8 digits
  #' temp$Serial_No<-rep(Serial_No,length(temp[[1]]))
  #' }
  #' return(temp)
  #'}
  #'



#------------------------------------------------------------------------------------------------------------------------------------------- -
#' merging of fingerprint and parameter files
#' @param fingerprint fingerprint tibble, created by Spectro_load_fingerprints
#' @param parameters fingerprint tibble, created by Spectro_load_parameters
#' @note Warning messages: "Unknown or initialised colum: `Serial_No` ... 
#' known warning, function switches automatically to alternate case.
Spectro_merge_col<-function(fingerprint,parameter){ 
  #' case1 both have serial no. fingerprint serial no is extracted and reintroduced after merging... should be fine?
  if(is.null(fingerprint$Serial_No)==F&    
     is.null(parameter$Serial_No)==F){
    Serial_No<-fingerprint$Serial_No
    temp<-left_join(fingerprint,parameter,by=c("Date_Time","Serial_No"))
    temp$Serial_No<-Serial_No
  #' case 2 either only one set contains serial no. -> no .x and .y mess (simply joined according to date)
  #'        or no serial no at all... so no joining of it either 
  }else{
      temp<-left_join(fingerprint,parameter,by=c("Date_Time"))
  }
  return(temp)
}
#------------------------------------------------------------------------------------------------------------------------------------------- -
#' local unzip function (optional new parent folder)
#' @param zip_directory zip-folder that is to be extracted
#' @param new Default is NULL. Input creates new parent directory for the extracte file,
#' e.g. .../new/folder/, where folder/ is named after the unzipped file "folder.zip"
unzip_local<-function(zip_directory,new_dir=NULL){
  if(is.null(new_dir)==F){
    # when new dir is given, this creates the path
    new_dir<-paste0(dirname(zip_directory),"/",new_dir,"/")
    # checks if new_dir already exists, gets rid of annoying warnings
    if(dir.exists(new_dir)==F){
    dir.create(new_dir)
      }
    unzip(zip_directory,exdir = paste0(new_dir,str_remove(basename(zip_directory),".zip")))
  }else{
    # same folder as for zip file 
    unzip(zip_directory,exdir = str_remove(zip_directory,".zip"))
  }
}




######### #
#         #
#  works  #
#         #
######### #
#' dir:   "C:/Users/SeanA/OneDrive/Sosa/Spectrolyzer/"
#'        "C:/Users/SeanA/OneDrive/Sosa/Spectro_new_04_10_22/"
#------------------------------------------------------------------------------------------------------------------------------------------- -
#' enables load-in of multipe scpectrolyzer-zip-files or already un-zipped files to coherent tibble
#' @param parent_dir path of folder containing the zip- or un-zipped files
#' @param wavelengths col-names of measured wavelengths. default is 200-750 in 2.5 m, steps
#' @param sel_wavelengths Range of wavelengths that shall be returned
#' @param parameters list of recorded parameters, forming col-names (Default is Temperature only)
#' @param zip Are zip-files (T) or already unzipped files (F) loaded? Default is zip=T.
#' @param exclude Are there unwanted zip-files (when zip=T) or folders (when zip=F)? 
#' @param read_param If TRUE, parameter fiels are read in, else skipped (...fix issues before setting TRUE)
#' If so, list them in a vector and they will be disregarded.
#' @notes if zip files are loaded, the unzipped folders are put in a new parent folder called "data" inside
#' of the original directory
Spectro_batch_load<-function(parent_dir,
                             wavelengths=seq(200,750,2.5),
                             sel_wavelengths=seq(200,722.5,2.5),
                             parameters=c("Temp"),
                             zip=T,
                             exclude=NA,
                             read_param=FALSE){
  
  # check if directroy string ends with "/"; if not, add it
  if(str_ends(parent_dir,"/")==F){
    parent_dir<-paste0(parent_dir,"/")
  }
  
  # unzippery
  if(zip==T){
    zip_folders<-list.files(parent_dir,pattern=".zip")
    zip_folders<-subset(zip_folders,zip_folders%in%exclude==F)

    for (i in zip_folders){
      unzip_local(paste0(parent_dir,i),new_dir="Spectro_data")
    }
    
    # after unzipping, change parent_dir to new folder with unzipped data
    parent_dir<-paste0(parent_dir,"Spectro_data/")
    
    #' might cause an issue, if the "data" folder is "contaminated" 
    #' with folders that were previously put there.
    #' [-1] because first entry is always parent dir itself
    folders<-list.dirs(strtrim(parent_dir,nchar(parent_dir)-1))[-1]

  
  }else{
    #' when already unzipped data: 
    #' [-1] because first entry is always parent dir itself
    #' also weird: list.dirs wants "path/dir" but NOT "path/dir/" 
    folders<-list.dirs(strtrim(parent_dir,nchar(parent_dir)-1))[-1]
    subset(folders,folders%in%exclude==F)
  }
  
  # initialising
  S_fp<-tibble()
  S_par<-tibble()
  print("List of folders:/n")
  print(folders)
  # load all data into the MOAT (mother of all tibbles (= )
  for (i in folders){
    # fingerprint read-in
    if(str_detect(i,"fingerprint")){
    S_fp<-bind_rows(S_fp,Spectro_load_fingerprint(i, wavelengths = wavelengths, sel_wavelengths = sel_wavelengths, id=T))
    
    # parameter read-in
    }else if(str_detect(i,"parameter")&read_param==TRUE){
    S_par<-bind_rows(S_par,Spectro_load_parameter(i,parameters = parameters,id=F)) # faster without double Serial-No. read-in
    # ERROR case
    }else if(read_param==FALSE){
      print("skipping parameter data")
    }else
      print("ERROR Unknown folder. Only folders containing 'fingerprint' or 'parameter' designation are allowed")
  }
  
  
  # combining fingerprint and parameter tables
  if(read_param){
    Spectro_data<-Spectro_merge_col(S_fp,S_par)
  }else{
    Spectro_data=S_fp
  }
  #' remove duplicates, which may have been caused by overlapping timeframes of teh readouts
  #' which are ultimately the zip files that were read in earlier
  Spectro_data<-distinct(Spectro_data,paste(format(Date_Time,"%d.%m.%Y-%H:%M:%S"),Serial_No),.keep_all = T) #aslo creates unique column
  names(Spectro_data)[length(Spectro_data)]<-"identifier"
  return(Spectro_data)
}



