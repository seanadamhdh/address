source("package_load.R")

##############
#            #
#   works    #
#            #
##############
#' Creates 1st derivative of the absorbtion curves of a Spectrolyzer dataset 
#' @param dataset Input tibble; Created by Sprectro_merge_col or Spectro_batch_load, maybe even Spectro_load_fingerprint 
#' or equally formatted 
#' @param cutoff_long To which wavelength should the derivate be calculated?
#' @param cutoff_short From which wavelength should the derivate be calculated?
#'  
Spectro_derivative<-function(dataset,cutoff_long=NULL,cutoff_short=NULL){
  #' finds the cols with numbers (abs-wavelengths) for subsetting... slightly dodgy, but should work in combination with tibbles created
  #' with functions mentioned above
  #' Then finds min, max, and stepsize of the wavelength spectrum
  spec_range<-summarise(
    tibble(data=as.numeric(
      subset(names(dataset),
             str_detect(names(dataset),"[0-9]+.")
             ))),
    min=min(data),
    max=max(data),
    step=data[2]-data[1])
  
  # initialising output tibble using input tibble (data will be overwritten)
  derivate_data<-dataset 
  
  # setting startpoint
  if(is.null(cutoff_long)){
    cutoff_long<-spec_range$max
  }
  i<-cutoff_long
  
  if(is.null(cutoff_short)){
    cutoff_short<-spec_range$min
  }
  while(i>cutoff_short){
    i_minus_h<-i-spec_range$step
    
    #'calculation of derivative
    #' moving from the higher wavelengths to the lower ones - 
    #' calculates derivative and assigns it to the respective shorter wavelength (i_minus_h)
    #' this decision was made, because then the startpoint remains the shortest wavelength (200 nm in our case)
     
    y<-format(i_minus_h,nsmall=1)
    #print(y)
    di<-transmute(dataset,
                  di=(dataset[format(i,nsmall=1)]-dataset[format(i_minus_h,nsmall=1)])/spec_range$step)[["di"]]
    #print(names(di))
    derivate_data[y]<-di
    #print(names(derivate_data[[y]]))  
    i<-i_minus_h
  }
  derivate_data<-select(derivate_data,-format(seq(cutoff_long,spec_range$max,spec_range$step),nsmall=1))
  
  return(derivate_data)
}



