# automatically installes and loads useful packages.
# some packages might have some issues with dependancies / R versions / OS

# if(!require(Rtools)){
#   install.packages("Rtools")
#  require(Rtools)}
if(!require(devtools)){
  install.packages("devtools")
  require(devtools)}
if(!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)}
if(!require(readxl)){
  install.packages("readxl")
  require(readxl)}


# advanced data handling and model building
if(!require(resemble)){
  install.packages("resemble")
  require(resemble)}
if(!require(caret)){
  install.packages("caret")
  require(caret)}
if(!require(Cubist)){
  install.packages("Cubist")
  require(Cubist)}
if(!require(mdatools)){
  install.packages("mdatools")
  require(mdatools)}
if(!require(simplerspec)){
  install_github("https://github.com/philipp-baumann/simplerspec.git")
  require(simplerspec)}
if(!require(prospectr)){
  #install_packages("prospectr")
  install_github("https://github.com/l-ramirez-lopez/prospectr.git")
  require(prospectr)}

# for fancy plots
if(!require(ggExtra)){
  install.packages("ggExtra")
  require(ggExtra)}
if(!require(ggthemes)){
  install.packages("ggthemes")
  require(ggthemes)}
if(!require(ggpubr)){
  install.packages("ggpubr")
  require(ggpubr)}
if(!require(ggtern)){
  install.packages("ggtern")
  require(ggtern)}
if(!require(plotly)){
  install.packages("plotly")
  require(plotly)}

# geospatial stuff
#if(!require(rgdal)){
#  install.packages("rgdal")
#  require(rgdal)
#}
if(!require(sf)){
  install.packages("sf")
  require(sf)
}
if(!require(crs)){
  install.packages("crs")
  require(crs)
}


# some own functions ####

## basic statistics ####
summarise_metrics<-function(dataset,group=NA,parameters=NA){
  summary_df<-tibble()

  if(is.na(group)){  #run without grouping
    for (i in parameters){
      summary_df<-bind_rows(
        summary_df,
        tibble(
          parameter=i,
          summarise(
            dataset,
            min=min(.data[[i]]),
            q25=quantile(.data[[i]],.25),
            median=median(.data[[i]]),
            mean=mean(.data[[i]]),
            q75=quantile(.data[[i]],.75),
            max=max(.data[[i]]),
            sd=sd(.data[[i]]),
            var=var(.data[[i]])
          )
        )
      )
    }
  }else{    #run with grouping
    for (i in parameters){
      summary_df<-bind_rows(
        summary_df,
        tibble(
          parameter=i,
          summarise(
            dataset%>%group_by(get(group)),  #only works for one grouping argument right now
            min=min(.data[[i]]),
            q25=quantile(.data[[i]],.25),
            median=median(.data[[i]]),
            mean=mean(.data[[i]]),
            q75=quantile(.data[[i]],.75),
            max=max(.data[[i]]),
            sd=sd(.data[[i]]),
            var=var(.data[[i]])
          )
        )
      )
    }
    names(summary_df)[2]<-group  #proprerly name grouping column
  }
  return(summary_df)
}


## wn nm conversions ####
# (eiher one works both ways, technically redundant)
wavenumber_to_wavelength<-function(x){return(10^9/(x*10^2))}

wavelength_to_wavenumber<-function(x){return(10^9/(x*10^2))}

## return logical if duplicates in vector ####
all_duplicates<-function(x){return(duplicated(x)|duplicated(x,fromLast=T))}

## generic spc plotting ####
spectra_plotter<-function(dataset,spc_id="spc_rs",sample_size=nrow(dataset),interactive=F,reverse_x=T){
  dataset<-dataset[sample(nrow(dataset),sample_size),]
  bind_cols(ID=dataset[["sample_id"]],dataset[[spc_id]])%>%
    pivot_longer(
      cols = colnames(.)[-1],
      names_to = "wavenumber",
      names_transform = list("wavenumber"=as.numeric),
      values_to = "absorbance"
    )%>%
    ggplot(aes(x=wavenumber,y=absorbance,col=ID))+
    geom_line(alpha=.1)+
    theme_minimal()+
    theme(legend.position = "none")->plt
  if(reverse_x==T){
    plt<-plt+scale_x_reverse()
  }
  if(interactive==T){
    ggplotly(plt)
    }else{
      plot(plt)
    }
}


## simple validation metrics ####

ME <- function(obs, pred){
  mean(pred - obs, na.rm = TRUE)
}

RMSE <- function(obs, pred){
  sqrt(mean((pred - obs)^2, na.rm = TRUE))
}


R2 <- function(obs, pred){
  # sum of the squared error
  SSE <- sum((pred - obs) ^ 2, na.rm = T)
  # total sum of squares
  SST <- sum((obs - mean(obs, na.rm = T)) ^ 2, na.rm = T)
  R2 <- 1 - SSE/SST
  return(R2)
}

## colorblind safe scale ####
#https://jrnold.github.io/ggthemes/reference/colorblind.html
colorblind_safe_colors<-c("#000000",
                          "#E69F00",
                          "#56B4E9",
                          "#009E73",
                          "#F0E442",
                          "#0072B2",
                          "#D55E00",
                          "#CC79A7")

