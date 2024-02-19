source("package_load.R")
spectro_raw_all$MW<-transmute(spectro_raw_all,
                              MW=case_when(Serial_No==22130201~"MW1",
                                           Serial_No==22130200~"MW2"))[["MW"]]
  



# adding col for Messfeld ID (derived from Serial_No) to Spectro data
  
merged<-tibble()
n<-c()
for (i in c(1:length(Wehrproben_R_data[[1]]))){ 
subset(spectro_raw_all,
       Date_Time%within%interval(Wehrproben_R_data$UTC[i]-hours(1),
                                 Wehrproben_R_data$UTC[i]+hours(1))&
         MW==Wehrproben_R_data$Messwehr[i])%>%
    summarize_if(is.numeric,
                 mean)->spectro_sub
  subset(spectro_raw_all,
         Date_Time%within%interval(Wehrproben_R_data$UTC[i]-hours(1),
                                   Wehrproben_R_data$UTC[i]+hours(1))&
           MW==Wehrproben_R_data$Messwehr[i])[[1]]%>%length()->n
  merged<-bind_rows(merged,tibble(Wehrproben_R_data[i,],spectro_sub,n=n))
  }
names(merged)[c(14:225)]<-paste0("A",format(seq(200.0,725.0,5),nsmall=1))


write_excel_csv(merged,"C:/Users/SeanA/OneDrive/Sosa/R/SOSA/tables/merged.csv")

## 1st der
y<-Spectro_derivative(spectro_raw_all)
## 2nd der
z<-Spectro_derivative(y)





###################

# calibration attempt 1 - simple linear regression using one bandwidth####

raw_lm<-c()
for (i in format(seq(200,727.5,2.5),nsmall=1)){
  mod<-lm(
    data=transmute(merged,DOC=DOC,spectro=.data[[i]]),
    formula=spectro~DOC)
  raw_lm<-rbind.data.frame(
    raw_lm,
    as.numeric(
      c(
        i,
        mod$coefficients[[1]],
        mod$coefficients[[2]],
        summary(mod)$coefficients[2,4],
        summary(mod)$r.squared
        )
      )
  )
}
names(raw_lm)<-c("nm","intercept","slope","p_val","R2")

write_excel_csv(raw_lm,file="C:/Users/SeanA/OneDrive/Sosa/R/SOSA/tables/raw_mod.csv")



#' issue: outliers
#' solution: using high wavelength bandwidth and cutting off data above threshold value
#' all datapoints with A>=100 in the 725 band are removed - lower takes to many precious datapoints away


cap725_lm<-c()
for (i in format(seq(200,727.5,2.5),nsmall=1)){
  mod<-lm(
    data=transmute(subset(merged,merged[["725.0"]]<100),DOC=DOC,spectro=.data[[i]]),
    formula=DOC~spectro)
  cap725_lm<-rbind.data.frame(
    cap725_lm,
    as.numeric(
      c(
        i,
        mod$coefficients[[1]],
        mod$coefficients[[2]],
        summary(mod)$coefficients[2,4],
        summary(mod)$r.squared
      )
    )
  )
}
names(cap725_lm)<-c("nm","intercept","slope","p_val","R2")

write_excel_csv(cap725_lm,file="C:/Users/SeanA/OneDrive/Sosa/R/SOSA/tables/cap725_mod.csv")

# wavelegths with R2 and pval of regression slope better than cutoffs

p_val_max<-0.001
R2_min<-.85
#subset(cap725_lm,cap725_lm$R2>R2_min&p_val<p_val_max)
wavelength_selection<-format(subset(cap725_lm,cap725_lm$R2>R2_min&p_val<p_val_max)$nm,nsmall=1)



# plotting selection with highest fit
temp<-subset(merged,merged[["725.0"]]<100,select=c(wavelength_selection,"DOC","Messwehr"))
temp<-pivot_longer(temp,cols=wavelength_selection,names_to = "nm",values_to = "A")


svg("C:/Users/SeanA/OneDrive/Sosa/R/SOSA/plots/cap725_lm.svg")             
ggplot(temp,aes(x=DOC,y=A))+
  geom_point(aes(col=Messwehr))+
  geom_smooth(method="glm")+
  facet_wrap(facets=vars(nm))+
  theme_test()
dev.off()


# plotting results of cap_lm against wavelengths
svg("C:/Users/SeanA/OneDrive/Sosa/R/SOSA/plots/cap725_results.svg") 
pivot_longer(cap725_lm,cols = names(cap725_lm)[-1])%>%
  ggplot(aes(x=nm,y=value))+
  geom_point()+
  facet_wrap(facets = vars(name),scales="free_y")
dev.off()

# for displaying the timeline, cutoff vals as low as 5 work
ggplot(subset(spectro_raw_all,spectro_raw_all[["725.0"]]<5),aes(x=Date_Time,col=MW))+geom_point(aes(y=.data[["425.0"]]),size=.3)





############################################################################
# multiple linear regression lambda #######


# single wavelength regressions
lambda_lm<-c()
for (i in paste0("A",seq(250,800,1))){
  mod<-lm(
    data=transmute(subset(Messfeld_Lambda_merged,DOC<200&A380<2.5), # manual removal of 2 extreme outliers
                   DOC=DOC,lambda=.data[[i]]),
    formula=DOC~lambda)
  lambda_lm<-rbind.data.frame(
    lambda_lm,
    as.numeric(
      c(
        str_sub(i,2),
        mod$coefficients[[1]],
        mod$coefficients[[2]],
        summary(mod)$coefficients[2,4],
        summary(mod)$r.squared,
        sigma(mod),
        cor.test(Messfeld_Lambda_merged[["DOC"]],Messfeld_Lambda_merged[[i]])$estimate[[1]]
      )
    )
  )
}
names(lambda_lm)<-c("nm","intercept","slope","p_val","R2","RSME","cor_pearson")

pivot_longer(lambda_lm,cols = names(lambda_lm)[-1])%>%
  ggplot(aes(x=nm,y=value))+
  geom_point()+
  facet_wrap(facets = vars(name),scales="free_y")
