##### Spectrolyzer cleanup #####


source("package_load.R")
require(forecastML)
require(zoo)
require(data.table)
# function to remove doublettes
Spectro_duplicate_cleanup<-function(spectro){
  return(
    subset(spectro,!duplicated(spectro,by=c("Date_Time","Serial_No")))
  )
}

###########################################################################################################################################################################
# loding spectro data ####

#' basically 1 big parent folder with all downloaded data would suffice (Spectro_batch_load form there)
#' 
###########################################################################################################
Spectro_batch_load("C:/Users/SeanA/OneDrive/Sosa/Spectrolyzer/alles 06.10/Spectro_data/",zip=F)->spectro_raw1
Spectro_batch_load("C:/Users/SeanA/OneDrive/Sosa/Spectro_new_04_10_22/Spectro_data/",zip=F)->spectro_raw2
spectro_data<-bind_rows(spectro_raw1,spectro_raw2)
spectro_data<-Spectro_duplicate_cleanup(spectro_data)

#new readout from 14.10.22 added
Spectro_batch_load("C:/Users/SeanA/OneDrive/Sosa/Spectrolyzer/14.10/07_JUL/")->spectro_JUL
Spectro_batch_load("C:/Users/SeanA/OneDrive/Sosa/Spectrolyzer/14.10/08_AUG/")->spectro_AUG
Spectro_batch_load("C:/Users/SeanA/OneDrive/Sosa/Spectrolyzer/14.10/09_SEP/")->spectro_SEP
Spectro_batch_load("C:/Users/SeanA/OneDrive/Sosa/Spectrolyzer/14.10/10_OKT/")->spectro_OKT
bind_rows(spectro_data,spectro_JUL,spectro_AUG,spectro_SEP,spectro_OKT)->spectro_new
spectro_new%>%Spectro_duplicate_cleanup()->spectro_new

spectro_data<-spectro_new
###################################################################################################




# add MW col for easier IDing
spectro_data$MW<-transmute(spectro_data,
                           MW=case_when(Serial_No==22130201~"MW1",
                                        Serial_No==22130200~"MW2"))[["MW"]]

# quantile based appraoch
# spectro_data_clean<-subset(spectro_data, spectro_data[["700.0"]]<quantile(spectro_data[["700.0"]],.95))

# seperate cutoffs for MW1 and MW2 - different ranges ####
cutoff_percentile<-.95
qMW1<-quantile(subset(spectro_data,MW=="MW1")[["700.0"]],cutoff_percentile)
qMW2<-quantile(subset(spectro_data,MW=="MW2")[["700.0"]],cutoff_percentile)
#' old, just q95 cutoff
#spectro_data_clean<-subset(spectro_data,
 #                          MW=="MW1"&spectro_data[["700.0"]]<qMW1|
  #                           MW=="MW2"&spectro_data[["700.0"]]<qMW2)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!
# add new timesteps to the interval list when measurement interval is changed
t_steps<-c(2,15) # messy but no other idea how to accomodate for different timesteps
# cutoff amount scaling -> e.g., 4 looks if 4 rows above is NA
n<-5
# !!!!!!!!!!!!!!!!!!!!!!!!!!!
spectro_data_cleaner_MW1<-subset(spectro_data,MW=="MW1"&spectro_data[["700.0"]]<qMW1)
spectro_data_cleaner_MW2<-subset(spectro_data,MW=="MW2"&spectro_data[["700.0"]]<qMW2)

i<-0
while(i<n){
  spectro_data_cleaner_MW1<-subset(spectro_data_cleaner_MW1,
                                   (Date_Time-shift(Date_Time,1))%in%(t_steps))
  
  spectro_data_cleaner_MW2<-subset(spectro_data_cleaner_MW2,
                                   (Date_Time-shift(Date_Time,1))%in%(t_steps))
  i<-i+1

  }
# new, with 
spectro_data_clean<-bind_rows(spectro_data_cleaner_MW1,spectro_data_cleaner_MW2)





# comparison of timeline
# set timeframe and wavelength: 
A<-"300.0"
start<-as.POSIXct("2022-07-01")
end<-as.POSIXct("2022-07-31")
ggplot(subset(spectro_data,Date_Time%within%interval(start = start,end=end)),aes(x=Date_Time))+
  geom_point(aes(y=.data[[A]],col="red"))+
  geom_point(data=subset(spectro_data_clean,Date_Time%within%interval(start = start,end=end)), aes(y=.data[[A]]))+
  geom_point(data=subset(mutate(Wehrproben_R_data,MW=Messwehr),
                         UTC%within%interval(start = start,end=end)),aes(x=UTC,y=600),col="green")+
  facet_wrap(facets = vars(MW))

# cutoff visualisation
# double-log, +5 because of negative absorbtion values
ggplot(spectro_data,aes(x=.data[["700.0"]]+5))+
  geom_histogram()+
  geom_vline(data=filter(spectro_data,MW=="MW1"),aes(xintercept = 5+qMW1),col="red")+
  geom_vline(data=filter(spectro_data,MW=="MW2"),aes(xintercept = 5+qMW1),col="red")+
  scale_y_log10()+
  scale_x_log10(breaks=c(5,15,105,1005),labels=c(0,10,100,1000))+
  facet_wrap(facets = vars(MW))+
  xlab("A700 [Abs/m-1]")

################################################################################################################################################################
#### loading Wehr-Data ####

Wehrproben <- read_excel("C:/Users/SeanA/OneDrive/Sosa/Messfelder/Daten/SOSA_DOC_all.xlsx",
                         sheet = "Wehre", col_types = c("text","text","numeric", "date",
                                                        "date", "date", "numeric","numeric",
                                                        "date", "numeric", "numeric","numeric",
                                                        "numeric", "date", "numeric","numeric",
                                                        "numeric", "numeric", "numeric","text"))
# "mimimi, expecting date got NA"-error is ok, no problems

# setting UTC/MEZ DateTime up correctly
Wehrproben$MEZ<-as.POSIXct(
  paste(
    Wehrproben$Probennahmedatum,
    as_hms(Wehrproben$MEZ)),
  tz="Europe/Berlin"
)
Wehrproben$UTC<-as.POSIXct(
  paste(
    Wehrproben$Probennahmedatum,
    as_hms(Wehrproben$UTC)),
  tz="UTC",
  format="%Y-%m-%d %H:%M:%S"
)

names(Wehrproben)<-c("ID",
                     "Messwehr",
                     "Rep",
                     "Probennahmedatum",
                     "MEZ",
                     "UTC",
                     "EC",
                     "pH",
                     "Analysedatum_DOC_DIC",
                     "VF",
                     "raw_DOC",
                     "DOC",
                     "DIC",
                     "Analysendatum_ICPOES",
                     "Fe",
                     "Al",
                     "Cu",          # some Met-cols are currently characters, because <Best.Gr. entries ... if needed replace with 0 or 1/2 Bestgr and convert to numeric
                     "Pb",
                     "Zn",
                     "Notes")


#######################################################################################################################################################################################################
#### merging with DOC data ####
#' function for merging Spectro data to Wehr data
#' @param Spectro_data - Spectrolyzer data, prepped with MW-column
#' @param Wehr_data - Wehrdata, including UTC column (Date-Time in UTC format) and Messwehr column
#' @param h - search radius in hours
#' @param wavelengths subset of wavelenghts to be selected from spectrolyzer spectra 
Spectro_Wehr_merge<-function(Spectro_data,Wehr_data,h=1,wavelengths=format(seq(200,725,2.5),nsmall=1)){
  Spectro_data<-select(Spectro_data,all_of(c("Date_Time","Serial_No","Temp","MW",wavelengths)))
  names(Spectro_data)[-c(1:4)]<-paste0("A",wavelengths)
  merged<-tibble()
  n<-c()
  for (i in c(1:length(Wehr_data[[1]]))){ 
    # Mean of Spectrolyzer data within search radius
    subset(Spectro_data,
           Date_Time%within%interval(Wehr_data$UTC[i]-hours(h),
                                     Wehr_data$UTC[i]+hours(h))&
             MW==Wehr_data$Messwehr[i])%>%
      summarize_if(is.numeric,
                   mean)->spectro_sub
    
    # No. of obs. within search radius
    subset(Spectro_data,
           Date_Time%within%interval(Wehr_data$UTC[i]-hours(h),
                                     Wehr_data$UTC[i]+hours(h))&
             MW==Wehr_data$Messwehr[i])[[1]]%>%length()->n
    
    merged<-bind_rows(merged,tibble(Wehr_data[i,],spectro_sub,n=n))
  }
  #names(merged)[c(14:225)]<-paste0("A",format(seq(200.0,725.0,5),nsmall=1))
  return(merged)
} 


# using 12 h selection radius -> overall 1 day worth of data
#Spectro_Wehr<-Spectro_Wehr_merge(spectro_data_clean,Wehrproben,h=12)

Spectro_Wehr<-Spectro_Wehr_merge(spectro_data_clean,Wehrproben,h=12)



Spectro_calib<-subset(Spectro_Wehr,is.na(DOC)==F&is.na(A300.0)==F) # A300.0 works just like any other wavelength for checking data availability 

Spectro_calib_icpoes<-subset(Spectro_Wehr,is.na(Fe)==F&is.na(A300.0)==F) # A300.0 works just like any other wavelength for checking data availability 



write_excel_csv(Spectro_calib,"C:/Users/SeanA/OneDrive/Sosa/R/SOSA/tables/Spectro_calib.csv")

Spectro_calib_raw<-Spectro_Wehr_merge(spectro_data,
                   Wehrproben,h=12)
Spectro_calib_raw<-subset(Spectro_calib_raw,
                          is.na(DOC)==F&is.na(A200.0)==F)

######################################################################################################

# filling holes in time series ################

# linear interpolation for MW1
# subsetting for MW1, removing col 213 (727.5 nm), only partially recorded wavelength
Spectro_MW1<-subset(spectro_data_clean,Serial_No==22130201)[-c(213)]

fill_gaps(Spectro_MW1,date_col = 1,frequency = "1 min")->temp # creating timeseries with 1min intervals, equally spaced
# when coarser intervalls are choosen, issues come up, because of 2min and 15 min measurement intervalls in the timeseries - fixable?


# earlier attempt of further cleaning 
if(F){
# have to fill the gaps first, then see where gaps are "NA", then thin out with the subset below, then "refill" the removed cols
# probably there is a much smarter way of doing it - e.g. via Date_Time, but this kind of works as well


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################


#'!!!!!!!!!!!!!!!!!!!!! important: add "is.na(shift(temp[["200.0"]],timestep))==F" with or ("|") argument for each timeinterval of measurements
#' e.g. 2 min and 15 min are currently implemented... and new argument for 30 min (winter)
Spectro_cleaner_MW1<-subset(temp,
                                 is.na(shift(temp[["200.0"]],30))==F|
                                 is.na(shift(temp[["200.0"]],15))==F|
                                   is.na(shift(temp[["200.0"]],2))==F)









# testing 
if(F){
ggplot(subset(spectro_data,MW=="MW1"),aes(x=Date_Time,y=.data[["300.0"]]))+ #raw
  geom_point()+
  geom_point(data=Spectro_MW1,col="red")+        #clean
  geom_point(data=spectro_data_cleaner_MW1,col="blue")#+ # cleaner
  xlim(as.POSIXct(c("2022-08-10","2022-08-25"),format="%Y-%m-%d"))

data.frame(
  "0"=temp[["200.0"]],
  "1"=is.na(shift(temp[["200.0"]],2)),
  "2"=is.na(shift(temp[["200.0"]],15)),
  "3"=is.na(shift(temp[["200.0"]],30)),
  "res"=is.na(shift(temp[["200.0"]],2))==F|
    is.na(shift(temp[["200.0"]],15))==F|
    is.na(shift(temp[["200.0"]],30))==F)%>%
View()

}


fill_gaps(Spectro_cleaner_MW1,date_col = 1,frequency = "1 min")->temp # creating timeseries with 1min intervals, equally spaced


}




# linear interpolation of absorbtion values
mutate_at(temp,.vars=c(format(seq(200,725,2.5),nsmall=1),"Temp"),.funs = ~na.approx(.))->Spectro_MW1_interpolated

# cubic spline interpolation of absorbtion values
# mutate_at(temp,.vars=c(format(seq(200,725,2.5),nsmall=1),"Temp"),.funs = ~na.spline(.))->Spectro_MW1_spline_interpolated
#' that REALLY did not work well


# manually filling in non-numeric cols, also creating data_type col with identifier: original - actual data and interpolated - ... well, interpolated
Spectro_MW1_interpolated$Serial_No<-rep(22130201,nrow(Spectro_MW1_interpolated)) 
Spectro_MW1_interpolated<-mutate(Spectro_MW1_interpolated,data_type=if_else(is.na(MW),"interpolated","original"))
Spectro_MW1_interpolated$identifier<-paste(Spectro_MW1_interpolated$Date_Time,Spectro_MW1_interpolated$Serial_No)
Spectro_MW1_interpolated$MW<-rep("MW1",nrow(Spectro_MW1_interpolated))



# linear interpolation for MW2

# subsetting for MW2, removing col 213 (727.5 nm), only partially recorded wavelength, details see above 
Spectro_MW2<-subset(spectro_data_clean,Serial_No==22130200)[-c(213)]

fill_gaps(Spectro_MW2,date_col = 1,frequency = "1 min")->temp
mutate_at(temp,.vars=c(format(seq(200,725,2.5),nsmall=1),"Temp"),.funs = ~na.approx(.))->Spectro_MW2_interpolated

Spectro_MW2_interpolated$Serial_No<-rep(22130200,nrow(Spectro_MW2_interpolated)) 
Spectro_MW2_interpolated<-mutate(Spectro_MW2_interpolated,data_type=if_else(is.na(MW),"interpolated","original"))
Spectro_MW2_interpolated$identifier<-paste(Spectro_MW2_interpolated$Date_Time,Spectro_MW2_interpolated$Serial_No)
Spectro_MW2_interpolated$MW<-rep("MW2",nrow(Spectro_MW2_interpolated))




#' in theory one could subsample the interpolated data with equidistant time steps to reduce data dnsity - however, this reduces the number of original datapoints in 
#' favour of uniform timesteps
#' e.g. 15 min steps:
Spectro_MW1_interpolated_15min<-subset(Spectro_MW1_interpolated,minute(Spectro_MW1_interpolated$Date_Time)%in%seq(0,45,15))
Spectro_MW2_interpolated_15min<-subset(Spectro_MW2_interpolated,minute(Spectro_MW2_interpolated$Date_Time)%in%seq(0,45,15))

Spectro_MW1_interpolated_10min<-subset(Spectro_MW1_interpolated,minute(Spectro_MW1_interpolated$Date_Time)%in%seq(0,50,10))
Spectro_MW2_interpolated_10min<-subset(Spectro_MW2_interpolated,minute(Spectro_MW2_interpolated$Date_Time)%in%seq(0,50,10))


#### building DOC timeseries data frame ####
require(mdatools)

mod_DOC_1<-read_rds("C:/Users/SeanA/OneDrive/Sosa/R/SOSA/ipls models/200_725_om.rda")
# 200-725 nm, outlier detection based on om

summary(mod_DOC_1$om)
#summary####
#' Info: 
#' Number of selected components: 6
#' Cross-validation: full (leave one out)
#' 
#' Response variable: DOC
#' X cumexpvar Y cumexpvar    R2  RMSE Slope    Bias   RPD
#' Cal    99.99997    99.42434 0.994 1.002 0.994  0.0000 13.34
#' Cv           NA          NA 0.991 1.259 0.995 -0.0329 10.62
#####
predict(mod_DOC_1$om,Spectro_MW1_interpolated[mod_DOC_1$var.selected])->mod_DOC_1_MW1_predicted
mod_DOC_1_MW1_predicted$y.pred%>%data.frame()->predictions_MW1


predict(mod_DOC_1$om,Spectro_MW2_interpolated[mod_DOC_1$var.selected])->mod_DOC_1_MW2_predicted
mod_DOC_1_MW2_predicted$y.pred%>%data.frame()->predictions_MW2

# full spectro dataset with interpolated values and modelled DOC conc. in 1min resolution
DOC_predicted<-bind_rows(
  cbind.data.frame(Spectro_MW1_interpolated,DOC=predictions_MW1$Comp.6.DOC),
  cbind.data.frame(Spectro_MW2_interpolated,DOC=predictions_MW2$Comp.6.DOC))



# thinning out to 15min interval <- this one is used
DOC_predicted_15min<-bind_rows(
  subset(cbind.data.frame(Spectro_MW1_interpolated,DOC=predictions_MW1$Comp.6.DOC),minute(Date_Time)%in%c(0,15,30,45)),
  subset(cbind.data.frame(Spectro_MW2_interpolated,DOC=predictions_MW2$Comp.6.DOC),minute(Date_Time)%in%c(0,15,30,45)))


# mean for intervals ####
#' mean of 1 min interpolated values for 15min bin (easily changable to whatever interval)
#' takes a while
summarise_all(
  group_by(
    select(
      mutate(
        DOC_predicted,
        ts = round_date(
          Date_Time, 
          unit = '15 min')),
      -c("identifier",
         "data_type",
         "Date_Time")),
    ts,MW),
  mean)->DOC_predicted_15min_avg
names(DOC_predicted_15min_avg)[1]<-"Date_Time"
#####

ggplot(DOC_predicted,aes(x=Date_Time,y=DOC))+
  geom_point(aes(col="raw interpolated"))+
  geom_line(aes(col="raw interpolated"))+
  geom_point(data=DOC_predicted_15min,aes(col="single val pick"))+
  geom_line(data=DOC_predicted_15min,aes(col="single val pick"))+
  geom_point(data=DOC_predicted_15min_avg,aes(col="mean over timeinterval"))+
  geom_line(data=DOC_predicted_15min_avg,aes(col="mean over timeinterval"))+
  facet_grid(rows=vars(MW),scales="free_y")+
  ylim(5,55)+
  xlim(as.POSIXct("2022-08-23"),as.POSIXct("2022-08-24"))+
  scale_color_manual(breaks=c("raw interpolated",
                              "single val pick",
                              "mean over timeinterval"),
                     values=c("black",
                              "red",
                              "blue"))
  

# calculation ratios etc ####
# ! runtime allows for a quick coffee break

names(DOC_predicted_15min)[215]<-"ID"
Spectro_predicted_and_indicators<-left_join(mutate(DOC_predicted_15min,
                                             SUVA254=SUVA254(DOC,.data[["255.0"]],abs_unit_length = 1), #using 255, as 254 is not available for Spectro
                                             E2E3=E2_E3(.data[["250.0"]],.data[["365.0"]]),
                                             E2E4=E2_E4(.data[["200.0"]],.data[["400.0"]]),
                                             E4E6=E4_E6(.data[["465.0"]],.data[["665.0"]]),
                                             CDOM_abs=CDOM_absorbance(.data[["440.0"]],.data[["690.0"]],abs_unit_length = 1)
)
,S_R(DOC_predicted_15min,2.5),by="ID")

# save all
write_excel_csv(Spectro_predicted_and_indicators,file="C:/Users/SeanA/OneDrive/Sosa/R/SOSA/tables/Spectro_predicted_and_ratios_ALL.csv")


#exclude spectra
Spectro_indicators<-select(mutate(Spectro_predicted_and_indicators,A300=.data[["300.0"]]),
       -format(seq(200,725,2.5),nsmall=1))
write_excel_csv(Spectro_indicators,
                file="C:/Users/SeanA/OneDrive/Sosa/R/SOSA/tables/Spectro_predicted_and_ratios.csv")


















