
##############################
##### 3. Fourier spectra #####
##############################

### Rainfall ###
Rainfall_spec<-spectrum(Rainfall_monthly_E$Rainfall_standard_mm,demean=T,detrend=T)
Rainfall_spec_df<-data.frame(Freq=Rainfall_spec$freq,Spec=Rainfall_spec$spec,Data="Rainfall",Type="",Site="Savanna")


### Temperature ###
MaxTemp_Saline_spec<-spectrum(Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="Saline"&Temperature_monthly_E$Type=="Max"],demean=T,detrend=T)
MaxTemp_Saline_spec_df<-data.frame(Freq=MaxTemp_Saline_spec$freq,Spec=MaxTemp_Saline_spec$spec,Data="Temp (max)", Type="Max",Site="Forest")

MaxTemp_SEGC_spec<-spectrum(Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="SEGC"&Temperature_monthly_E$Type=="Max"],demean=T,detrend=T)
MaxTemp_SEGC_spec_df<-data.frame(Freq=MaxTemp_SEGC_spec$freq,Spec=MaxTemp_SEGC_spec$spec,Data="Temp (max)", Type="Max",Site="Savanna")

MinTemp_Saline_spec<-spectrum(Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="Saline"&Temperature_monthly_E$Type=="Min"],demean=T,detrend=T)
MinTemp_Saline_spec_df<-data.frame(Freq=MinTemp_Saline_spec$freq,Spec=MinTemp_Saline_spec$spec,Data="Temp (min)", Type="Min",Site="Forest")

MinTemp_SEGC_spec<-spectrum(Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="SEGC"&Temperature_monthly_E$Type=="Min"],demean=T,detrend=T)
MinTemp_SEGC_spec_df<-data.frame(Freq=MinTemp_SEGC_spec$freq,Spec=MinTemp_SEGC_spec$spec,Data="Temp (min)", Type="Min",Site="Savanna")


### Humidity ###
AbsHumidity_Savanna_spec<-spectrum(AbsHumidity_monthly_C$AbsHumidity_standard_gm3[AbsHumidity_monthly_C$Site=="Savanna"],demean=T,detrend=T)
AbsHumidity_Savanna_spec_df<-data.frame(Freq=AbsHumidity_Savanna_spec$freq,Spec=AbsHumidity_Savanna_spec$spec,Data="Humidity", Type="",Site="Savanna")

AbsHumidity_Forest_spec<-spectrum(AbsHumidity_monthly_C$AbsHumidity_standard_gm3[AbsHumidity_monthly_C$Site=="Forest"],demean=T,detrend=T)
AbsHumidity_Forest_spec_df<-data.frame(Freq=AbsHumidity_Forest_spec$freq,Spec=AbsHumidity_Forest_spec$spec,Data="Humidity", Type="",Site="Forest")


### Solar ###
Solar_spec<-spectrum(Solar_monthly_C$Solar_standard_Wm2,demean=T,detrend=T)
Solar_spec_df<-data.frame(Freq=Solar_spec$freq,Spec=Solar_spec$spec,Data="Solar radiation", Type="",Site="Savanna")


### Wind ###
Wind_spec<-spectrum(Wind_monthly_C$Wind_standard_m_s,demean=T,detrend=T)
Wind_spec_df<-data.frame(Freq=Wind_spec$freq,Spec=Wind_spec$spec,Data="Wind spped", Type="",Site="Savanna")


### AOT ###
AOT_spec<-spectrum(Aerosol_monthly_C$AOT_500_standard,demean=T,detrend=T)
AOT_spec_df<-data.frame(Freq=AOT_spec$freq,Spec=AOT_spec$spec,Data="AOD", Type="",Site="Savanna")


## -  Combine all data
Spec_combined<-rbind(Rainfall_spec_df,
                     MaxTemp_Saline_spec_df,
                     MaxTemp_SEGC_spec_df,
                     MinTemp_Saline_spec_df,
                     MinTemp_SEGC_spec_df,
                     AbsHumidity_Savanna_spec_df,
                     AbsHumidity_Forest_spec_df,
                     Solar_spec_df,
                     Wind_spec_df,
                     AOT_spec_df)

## - Make combined plot (Supplementary Figure S1)
ggplot(Spec_combined[Spec_combined$Freq<0.5,],aes(x=Freq,y=Spec))+
  geom_hline(yintercept=-Inf)+
  geom_vline(xintercept=-Inf)+
  geom_line()+
  #scale_linetype_manual(values = c("dotted","solid","dashed"))+
  #scale_colour_manual(values=c("gray","black","darkgray"))+
  #scale_size_manual(values=c(0.8,1,0.5))+
  theme_classic()+
  facet_grid(Data~Site)+
  geom_vline(aes(xintercept=1/12),linetype="dotted")+
  geom_vline(aes(xintercept=1/6),linetype="dotted")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
  ylab("Spectrum")+
  xlab("Frequency (1/months)")


###############################
##### 4. Long-term trends #####
###############################

### Rainfall ###

## - Test for linear trend in total annual rainfall over time
# - Run cpglm

Rain_trend_m1.cpglmm<-cpglmm(Rainfall_mm~Year_r+(1|DOY),
                             data=Rainfall_daily_D,na.action=na.exclude)
Rain_trend_m2.cpglmm<-cpglmm(Rainfall_mm~1+(1|DOY),
                             data=Rainfall_daily_D,na.action=na.exclude)
# - Table 2
Rainfall_trend_AIC<-data.frame(Response="Rainfall",
                               Model=c("Long-term change", "No long-term change"),
                               Predictors=c("Year","Intercept only"),
                               DF= AIC(Rain_trend_m1.cpglmm,Rain_trend_m2.cpglmm)$df,
                               AIC=round(AIC(Rain_trend_m1.cpglmm,Rain_trend_m2.cpglmm)$AIC,digits=1))

Rainfall_trend_AIC$Delta_AIC<-Rainfall_trend_AIC$AIC-min(Rainfall_trend_AIC$AIC)

# - 95%CI
c(round(summary(Rain_trend_m1.cpglmm)$coefs[2,1] - 1.96* summary(Rain_trend_m1.cpglmm)$coefs[2,2],digits=2),
  round(summary(Rain_trend_m1.cpglmm)$coefs[2,1] + 1.96* summary(Rain_trend_m1.cpglmm)$coefs[2,2],digits=2))

# - Predict trend
Rainfall_daily_D$Rainfall_daily_mm_predict<-predict(Rain_trend_m1.cpglmm,Rainfall_daily_D,type="response")

Rainfall_yearly_F$Rainfall_annual_mm_predict<-ddply(Rainfall_daily_D,.(Year),summarise,
                                                    sum(Rainfall_daily_mm_predict))[,2]

# - Change in total annual rainfall (mm) per decade = -75mm per decade
#(Rainfall_yearly_F$Rainfall_annual_mm_predict[Rainfall_yearly_F$Year==2017]-Rainfall_yearly_F$Rainfall_annual_mm_predict[Rainfall_yearly_F$Year==1984])/(2017-1984)*10
#Percentage change compared to mean total annual rainfall = -0.0547
#(Rainfall_yearly_F$Rainfall_annual_mm_predict[Rainfall_yearly_F$Year==2017]-Rainfall_yearly_F$Rainfall_annual_mm_predict[Rainfall_yearly_F$Year==1984])/(2017-1984)*10/mean(Rainfall_yearly_F$Rainfall_annual_mm_predict,na.rm=T)

# - Plot trend and raw data (Figure 3A)
Rainfall_annual_glm_plot<-ggplot()+geom_line(data=Rainfall_yearly_F[c(1,34),],aes(x=Year,y=Rainfall_annual_mm_predict))+
  geom_point(data=Rainfall_yearly_F,aes(x=Year,y=Rainfall_yearly_mm),colour="gray")+
  geom_line(data=Rainfall_yearly_F,aes(x=Year,y=Rainfall_yearly_mm),colour="gray")+
  scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015))+
  theme_classic()+
  ylab("Rainfall (mm/year)")+
  xlab("Year")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

Rainfall_annual_glm_plot

## - Test for linear trend interacting with season

# - CPGLM on daily data including interaction with season
Rain_seasonal_m1.cpglmm<-cpglmm(Rainfall_mm~Season*Year_r+(1|DOY),
                                data=Rainfall_daily_D,na.action=na.exclude)
Rain_seasonal_m2.cpglmm<-cpglmm(Rainfall_mm~Season+Year_r+(1|DOY),
                                data=Rainfall_daily_D,na.action=na.exclude)
# - Table 3

Rainfall_seasonal_trend_AIC<-data.frame(Response="Rainfall",
                                        Model=c("Long-term change by season", "Long-term change not by season"),
                                        Predictors=c("Year * Season","Year + Season"),
                                        DF= AIC(Rain_seasonal_m1.cpglmm,Rain_seasonal_m2.cpglmm)$df,
                                        AIC=round(AIC(Rain_seasonal_m1.cpglmm,Rain_seasonal_m2.cpglmm)$AIC,digits=1))
Rainfall_seasonal_trend_AIC$Delta_AIC<-Rainfall_seasonal_trend_AIC$AIC-min(Rainfall_seasonal_trend_AIC$AIC)

#Remove intercept to directly compare estimates
Rain_seasonal_m1.cpglmm_mod<-cpglmm(Rainfall_mm~-1+Season+Season:Year_r+(1|DOY),
                                    data=Rainfall_daily_D,na.action=na.exclude)

# - Extract estimates
Rain_seasonal_m1.cpglmm_mod_coefs <- data.frame(Response="Rainfall",
                                                Predictor=c("DJF","JJAS","MAM","ON",
                                                            "Year: DJF", "Year: JJAS","Year: MAM","Year: ON"),
                                                Est = round(Rain_seasonal_m1.cpglmm_mod$fixef,digits=2) ,
                                                S.E. = round(summary(Rain_seasonal_m1.cpglmm_mod)$coefs[,2],digits=2),
                                                T.Value= round(summary(Rain_seasonal_m1.cpglmm_mod)$coefs[,3],digits=2))

Rain_seasonal_m1.cpglmm_mod_coefs$LL = round(Rain_seasonal_m1.cpglmm_mod_coefs$Est - 1.96 * Rain_seasonal_m1.cpglmm_mod_coefs$S.E., digits=2)
Rain_seasonal_m1.cpglmm_mod_coefs$UL = round(Rain_seasonal_m1.cpglmm_mod_coefs$Est + 1.96 * Rain_seasonal_m1.cpglmm_mod_coefs$S.E., digits=2)
row.names(Rain_seasonal_m1.cpglmm_mod_coefs)<-NULL

#Table 4
Rainfall_seasonal_trend_glm<-data.frame(Response = "Rainfall",
                                        Predictor = c("DJF","JJAS","MAM","ON","Year: DJF","Year: JJAS","Year: MAM","Year: ON"),
                                        Estimate= round(Rainfall_seasonal_trend_glm[Rainfall_seasonal_trend_glm$group=="fixed","estimate"],digits=2),
                                        SE = round(Rainfall_seasonal_trend_glm[Rainfall_seasonal_trend_glm$group=="fixed","std.error"],digits=2),
                                        "T Value" = round(Rainfall_seasonal_trend_glm[Rainfall_seasonal_trend_glm$group=="fixed","statistic"],digits=2),
                                        "95% CI" =Rainfall_seasonal_trend_glm[Rainfall_seasonal_trend_glm$group=="fixed","95%CI"],
                                        "."=Rainfall_seasonal_trend_glm[Rainfall_seasonal_trend_glm$group=="fixed","sig"])

# - Predict from the GLM
Rainfall_seasonal<-ddply(Rainfall_daily_D,.(Year,Season),summarise,
                         Rainfall_daily_mean=mean(Rainfall_mm,na.rm=T),
                         Rainfall_daily_predict_mm=mean(Rainfall_daily_mm_predict))

# ddply(Rainfall_seasonal,.(Season),summarise,
#       Rainfall_1984=Rainfall_daily_predict_mm[Year==1984],
#       Rainfall_2017=Rainfall_daily_predict_mm[Year==2017],
#      Years=2017-1984,
#      Difference_mm=Rainfall_2017-Rainfall_1984,
#      Difference_mm_yrs=Difference_mm/Years,
#      Difference_mm_dec=Difference_mm_yrs*10,
#      Rainfall_mean_mm=mean(Rainfall_daily_mean,na.rm=T),
#      Difference_per_dec=round((Difference_mm_dec/Rainfall_mean_mm)*100,digits=2))




### Temperature ###

## - Test for linear trend in minimum daily temperature

# - Run LMM to test for linear trend over time
Temp_annual_m1<-lmer(Temperature_c~Year_r+(1|DOY),data=Temperature_daily_F,na.action = na.exclude)
Temp_annual_m2<-lmer(Temperature_c~1+(1|DOY),data=Temperature_daily_F,na.action = na.exclude)

# - Check for temporal autocorrelation
#Temperature_daily_F$Residuals<-NA
#Temperature_daily_F$Residuals<-resid(Temp_annual_m1)
#acf_plot(Temperature_daily_F$Residuals,split_by=list(Temperature_daily_F$DOY),fun=median)
#acf_n_plots(Temperature_daily_F$Residuals,split_by=list(Temperature_daily_F$DOY))

# - Table 2

Temperature_trend_AIC<-data.frame(Response="Temperature",
                                  Model=c("Long-term change", "No long-term change"),
                                  Predictors=c("Year","Intercept only"),
                                  DF= AIC(Temp_annual_m1,Temp_annual_m2)$df,
                                  AIC=round(AIC(Temp_annual_m1,Temp_annual_m2)$AIC,digits=1))
Temperature_trend_AIC$Delta_AIC<-Temperature_trend_AIC$AIC-min(Temperature_trend_AIC$AIC)

# - Predict from LMM
Temperature_daily_F$LMM_predict<-predict(Temp_annual_m1,Temperature_daily_F,re.form=NA)
Temperature_annual_LMM<-ddply(Temperature_daily_F,.(Year),summarise,
                              Temperature_interp_mean_c=mean(Temperature_interp_c,na.rm=T),
                              Temperature_mean_c=mean(Temperature_c,na.rm=T),
                              Temperature_interp_mean_sample=length(Temperature_interp_c[complete.cases(Temperature_interp_c)]),
                              Temperature_mean_sample=length(Temperature_c[complete.cases(Temperature_c)]),
                              Temperature_LMM_c=mean(LMM_predict,na.rm=T))

Temperature_annual_LMM<-Temperature_annual_LMM[complete.cases(Temperature_annual_LMM$Year),]

#(Temperature_annual_LMM$Temperature_LMM_c[Temperature_annual_LMM$Year==2017]-Temperature_annual_LMM$Temperature_LMM_c[Temperature_annual_LMM$Year==1984])/(2017-1984)*10
#0.25 per decade
#(Temperature_annual_LMM$Temperature_LMM_c[Temperature_annual_LMM$Year==2017]-Temperature_annual_LMM$Temperature_LMM_c[Temperature_annual_LMM$Year==1984])/(2017-1984)*10/mean(Temperature_daily_F$Temperature_c,na.rm=T)
#1.1%

# - Figure 3B
Temperature_ym_LMM<-ddply(Temperature_daily_F,.(YearMonth),summarise,
                          Temperature_mean_c=mean(Temperature_c,na.rm=T),
                          Temperature_mean_sample=length(Temperature_c[complete.cases(Temperature_c)]),
                          Temperature_LMM_c=mean(LMM_predict,na.rm=T))

Temperature_ym_LMM$Temperature_mean_c[Temperature_ym_LMM$Temperature_mean_sample<5]<-NA

MinTemp_annual_glm_plot<-ggplot()+
  geom_point(data=Temperature_ym_LMM,aes(x=YearMonth,y=Temperature_mean_c),colour="gray",alpha=0.4)+
  # geom_line(data=Temperature_ym_LMM,aes(x=YearMonth,y=Temperature_mean_c),colour="gray")+
  geom_line(data=Temperature_ym_LMM[month(Temperature_ym_LMM$YearMonth)==6,],aes(x=YearMonth,y=Temperature_LMM_c))+
  scale_x_continuous(breaks=c(1985,1990,1995,2000,2005,2010,2015))+
  theme_classic()+
  ylab("Min temperature (c)")+
  xlab("Year")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

## - Test for linear trend interacting with season
# - GLM on daily data including interaction with season
Temp_seasonal_m1<-lmer(Temperature_c~Year_r*Season+(1|DOY),data=Temperature_daily_F,na.action = na.exclude)
Temp_seasonal_m2<-lmer(Temperature_c~Year_r+Season+(1|DOY),data=Temperature_daily_F,na.action = na.exclude)

# - Table 3
Temperature_seasonal_trend_AIC<-data.frame(Response="Temperature",
                                           Model=c("Long-term change by season", "Long-term change not by season"),
                                           Predictors=c("Year * Season","Year + Season"),
                                           DF= AIC(Temp_seasonal_m1,Temp_seasonal_m2)$df,
                                           AIC=round(AIC(Temp_seasonal_m1,Temp_seasonal_m2)$AIC,digits=1))
Temperature_seasonal_trend_AIC$Delta_AIC<-Temperature_seasonal_trend_AIC$AIC-min(Temperature_seasonal_trend_AIC$AIC)

# - Tabke 4 - Modify best model by removing intercept
Temp_seasonal_m1_mod<-lmer(Temperature_c~-1+Season+Year_r:Season+(1|DOY),data=Temperature_daily_F,na.action = na.exclude)
summary(Temp_seasonal_m1_mod)$coefficients

Temp_seasonal_m1_mod_coefs <- data.frame(Response="Temperature",
                                         Predictor=c("DJF","JJAS","MAM","ON","Year: DJF", "Year: JJAS","Year: MAM","Year: ON"),
                                         Est = round(summary(Temp_seasonal_m1_mod)$coefficients[,1],digits=2),
                                         S.E. = round(summary(Temp_seasonal_m1_mod)$coefficients[,2],digits=2),
                                         T.Value= round(summary(Temp_seasonal_m1_mod)$coefficients[,3],digits=2))


Temp_seasonal_m1_mod_coefs$LL = round(Temp_seasonal_m1_mod_coefs$Est - 1.96 * Temp_seasonal_m1_mod_coefs$S.E., digits=2)
Temp_seasonal_m1_mod_coefs$UL = round(Temp_seasonal_m1_mod_coefs$Est + 1.96 * Temp_seasonal_m1_mod_coefs$S.E., digits=2)
row.names(Temp_seasonal_m1_mod_coefs)<-NULL

# - Predict from the LMM
Temperature_daily_F$Lopee_glm<-predict(Temp_seasonal_m1,Temperature_daily_F,re.form=NA,type="response")

# ddply(Temperature_daily_F,.(Season),summarise,
#      Temp_1984=mean(Lopee_glm[Year==1984],na.rm=T),
#      Temp_2017=mean(Lopee_glm[Year==2017],na.rm=T),
#      Years=2017-1984,
#      Difference_c=Temp_2017-Temp_1984,
#      Difference_c_yrs=Difference_c/Years,
#      Difference_c_dec=Difference_c_yrs*10,
#      Temp_mean_c=mean(Temperature_c,na.rm=T),
#      Difference_per_dec=(Difference_c_dec)/Temp_mean_c)


### - Gridded temperature data for Lopee ###

## - Berkeley
Berkeley_annual_m1<-lmer(Temperature_c~Year_r+(1|DOY),data=Berkeley_daily,na.action = na.exclude)
summary(Berkeley_annual_m1)

# - Predict from GLM
Berkeley_daily$GLM_predict<-predict(Berkeley_annual_m1,Berkeley_daily,re.form=NA,allow.new.levels=TRUE)
Berkeley_yearly<-ddply(Berkeley_daily,.(Year),summarise,
                       Temperature_mean_c=mean(Temperature_c,na.rm=T),Temperature_GLM_c=mean(GLM_predict,na.rm=T))
#(Berkeley_yearly$Temperature_GLM_c[Berkeley_yearly$Year==2017]-Berkeley_yearly$Temperature_GLM_c[Berkeley_yearly$Year==1984])/(2017-1984)*10
# - 0.16 per decade


## - CRU
CRU_annual_m1<-lmer(Temperature_c~Year_r+(1|Month),data=CRU_monthly,na.action = na.exclude)
summary(CRU_annual_m1)

# - Predict from GLM
#CRU_monthly$GLM_predict<-predict(CRU_annual_m1,CRU_monthly,re.form=NA,allow.new.levels=TRUE,type="response")
#CRU_yearly<-ddply(CRU_monthly,.(Year),summarise,Temperature_mean_c=mean(Temperature_c,na.rm=T),Temperature_GLM_c=mean(GLM_predict,na.rm=T))
#(CRU_yearly$Temperature_GLM_c[CRU_yearly$Year==2016]-CRU_yearly$Temperature_GLM_c[CRU_yearly$Year==1984])/(2016-1984)*10
# - 0.19 per decade



####################################
##### 5. Periodicity over time #####
####################################

### Rainfall ###

## - Compute wavelet transform on complete (interpolated) monthly time series for rainfall
Rainfall_wavelet<-wt(Rainfall_monthly_E[,c(2,9)],dj=1/12)
#plot(Rainfall_wavelet,xlab="",ylab="Rainfall period (yrs)") #Figure 3C

# - Extract wavelet trasnform for three components (6-month, 12-month and 2-4 years)
Rainfall_wavelet_extract<-data.frame(YearMonth=Rainfall_wavelet$t,
                                     COI=Rainfall_wavelet$coi,
                                     Power_6months=Rainfall_wavelet$power.corr[which.min(abs(Rainfall_wavelet$period-0.5)),],
                                     Power_12months=Rainfall_wavelet$power.corr[which.min(abs(Rainfall_wavelet$period-1)),],
                                     Power_2_4years=colMeans(Rainfall_wavelet$power.corr[which.min(abs(Rainfall_wavelet$period-2)):which.min(abs(Rainfall_wavelet$period-4)),]))

# - Remove data within the cone of influence (overly influenced by edge effects)
Rainfall_wavelet_extract$Power_6months[Rainfall_wavelet_extract$COI<0.5]<-NA
Rainfall_wavelet_extract$Power_12months[Rainfall_wavelet_extract$COI<1]<-NA
Rainfall_wavelet_extract$Power_2_4years[Rainfall_wavelet_extract$COI<4]<-NA
Rainfall_wavelet_extract$Date<-as.Date(paste(year(Rainfall_wavelet_extract$YearMonth),month(Rainfall_wavelet_extract$YearMonth),"15",sep="-"))

#mean(Rainfall_wavelet_extract$Power_6months,na.rm=T)/mean(Rainfall_wavelet_extract$Power_12months,na.rm=T)
#mean(Rainfall_wavelet_extract$Power_6months,na.rm=T)/mean(Rainfall_wavelet_extract$Power_2_4years,na.rm=T)

# - Plot extracted wavelet components (Figure 3E)
Rainfall_extract_plot<-ggplot(Rainfall_wavelet_extract)+
  geom_line(aes(x=Date,y=Power_6months*0.0001,linetype="Biannual"))+
  geom_line(aes(x=Date,y=Power_12months*0.0001,linetype="Annual"))+
  geom_line(aes(x=Date,y=Power_2_4years*0.0001,linetype="2-4 years"))+
  scale_linetype_manual(values=c("dotted","solid","longdash"))+
  scale_x_date(breaks=c(as.Date("1985-01-01"),as.Date("1990-01-01"),as.Date("1995-01-01"),as.Date("2000-01-01"),as.Date("2005-01-01"),as.Date("2010-01-01"),as.Date("2015-01-01")),date_label="%Y")+
  theme_classic()+
  ylab("Wavelet power (/10000)")+
  xlab("Year")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.title = element_blank(),legend.position=c(0.9,0.9))


### Temperature ###

## - Compute wavelet transform on complete (interpolated) monthly time series for mean minimum daily temperature
Temperature_wavelet<-wt(Temperature_monthly_G[,c(3,8)],dj=1/12)
#plot(Temperature_wavelet,xlab="",ylab="Temperature period (yrs)") # Figure 3D

# - Extract wavelet trasnform for three components (6-month, 12-month and 2-4 years)
Temperature_wavelet_extract<-data.frame(YearMonth=Temperature_wavelet$t,
                                        COI=Temperature_wavelet$coi,
                                        Power_6months=Temperature_wavelet$power.corr[which.min(abs(Temperature_wavelet$period-0.5)),],
                                        Power_12months=Temperature_wavelet$power.corr[which.min(abs(Temperature_wavelet$period-1)),],
                                        Power_2_4years=colMeans(Temperature_wavelet$power.corr[which.min(abs(Temperature_wavelet$period-2)):which.min(abs(Temperature_wavelet$period-4)),]))

# - Remove data within the cone of influence (overly influenced by edge effects)
Temperature_wavelet_extract$Power_6months[Temperature_wavelet_extract$COI<0.5]<-NA
Temperature_wavelet_extract$Power_12months[Temperature_wavelet_extract$COI<1]<-NA
Temperature_wavelet_extract$Power_2_4years[Temperature_wavelet_extract$COI<4]<-NA
Temperature_wavelet_extract$Date<-as.Date(paste(year(Temperature_wavelet_extract$YearMonth),month(Temperature_wavelet_extract$YearMonth),"15",sep="-"))

#mean(Temperature_wavelet_extract$Power_12months,na.rm=T)/mean(Temperature_wavelet_extract$Power_6months,na.rm=T)
#mean(Temperature_wavelet_extract$Power_12months,na.rm=T)/mean(Temperature_wavelet_extract$Power_2_4years,na.rm=T)

# - Plot extracted wavelet components (Figure 3F)
Temperature_extract_plot<-ggplot(Temperature_wavelet_extract)+
  geom_line(aes(x=Date,y=Power_6months,linetype="Biannual"))+
  geom_line(aes(x=Date,y=Power_12months,linetype="Annual"))+
  geom_line(aes(x=Date,y=Power_2_4years,linetype="2-4 years"))+
  scale_linetype_manual(values=c("dotted","solid","longdash"))+
  scale_x_date(breaks=c(as.Date("1985-01-01"),as.Date("1990-01-01"),as.Date("1995-01-01"),as.Date("2000-01-01"),as.Date("2005-01-01"),as.Date("2010-01-01"),as.Date("2015-01-01")),date_label="%Y")+
  theme_classic()+
  ylim(c(0,100))+
  ylab("Wavelet power")+
  xlab("Year")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.title = element_blank(),legend.position=c(0.9,0.9))


###################################
###### 6. Influence of oceans #####
###################################

## - Load teleconnection data
# - Multivariate ENSO Index (MEI; Wolter & Timlin 1993; Wolter & Timlin 1998) sourced from the NOAA website (https://www.esrl.noaa.gov/psd/enso/mei/index.html)
# - Indian Ocean Dipole (IOD) Dipole Mode Index (Saji & Yamagata 2003) sourced from the NOAA website (https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/DMI/)
# - SST anomalies for the tropical north Atlantic (NATL, 5Â°â€“20Â°N, 60Â°â€“30Â°W) and the south equatorial Atlantic (SATL ,0Â°â€“20Â°S, 30Â°Wâ€“10Â°E) sourced from the NOAA National Weather Service Climate Prediction Center (http://www.cpc.ncep.noaa.gov/data/indices/). 

Tele <- read.csv("../Data/Teleconnections_combined")
Tele$YearMonth<-as.yearmon(as.Date(paste(Tele$Year,Tele$Month,"01",sep="-")))

### Rainfall ###

## - Merge rainfall and teleconnection datasets
Rainfall_Tele_monthly<-merge(Rainfall_monthly_E[,c("Month","Year","YearMonth","Season","Rainfall_interp_mm")],Tele[,c("Year","Month","MEI","NATL","NATL_anom","SATL","SATL_anom","IOD","NAO")],c("Year","Month"))
Rainfall_Tele_monthly<-Rainfall_Tele_monthly[order(Rainfall_Tele_monthly$YearMonth),]
Rainfall_Tele_monthly<-Rainfall_Tele_monthly[complete.cases(Rainfall_Tele_monthly$MEI),]

#MEI
#plot(cbind(Rainfall_Tele_monthly[,"YearMonth"],Rainfall_Tele_monthly[,"MEI"]),type="l")
#MEI_wavelet<-wt(cbind(Rainfall_Tele_monthly[,"YearMonth"],Rainfall_Tele_monthly[,"MEI"]),dj=1/12)
#plot(MEI_wavelet,xlab="",ylab="MEI (yrs)")

Rain_MEI_wtc<-wtc(cbind(Rainfall_Tele_monthly[,"YearMonth"],
                        scale(Rainfall_Tele_monthly[,"Rainfall_interp_mm"])),
                  cbind(Rainfall_Tele_monthly[,"YearMonth"],
                        scale(Rainfall_Tele_monthly[,"MEI"])),
                  pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Rain_MEI_wtc$coi)) test<-cbind(test,Rain_MEI_wtc$period>Rain_MEI_wtc$coi[i])
test<-test[,-1]
Rain_MEI_wtc$rsq[test]<-NA

Rain_wtc_global<-data.frame(period=Rain_MEI_wtc$period,
                            global.coherence=rowMeans(Rain_MEI_wtc$rsq,na.rm=T),
                            x="Rainfall",
                            y="MEI")

#NATL_anom_anom
#plot(cbind(Rainfall_Tele_monthly[,"YearMonth"], Rainfall_Tele_monthly[,"NATL_anom"]),type="l")
#NATL_anom_wavelet<-wt(cbind(Rainfall_Tele_monthly[,"YearMonth"],Rainfall_Tele_monthly[,"NATL_anom"]), dj=1/12)
#plot(NATL_anom_wavelet,xlab="",ylab="NATL_anom (yrs)")

Rain_NATL_anom_wtc<-wtc(cbind(Rainfall_Tele_monthly[,"YearMonth"],
                              scale(Rainfall_Tele_monthly[,"Rainfall_interp_mm"])),
                        cbind(Rainfall_Tele_monthly[,"YearMonth"],
                              scale(Rainfall_Tele_monthly[,"NATL_anom"])),
                        pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Rain_NATL_anom_wtc$coi)) test<-cbind(test,Rain_NATL_anom_wtc$period>Rain_NATL_anom_wtc$coi[i])
test<-test[,-1]
Rain_NATL_anom_wtc$rsq[test]<-NA

Rain_wtc_global<-rbind(Rain_wtc_global,
                       data.frame(period=Rain_NATL_anom_wtc$period,
                                  global.coherence=rowMeans(Rain_NATL_anom_wtc$rsq,na.rm=T),
                                  x="Rainfall",
                                  y="NATL_anom"))

#SATL_anom
#plot(cbind(Rainfall_Tele_monthly[,"YearMonth"], Rainfall_Tele_monthly[,"SATL_anom"]),type="l")
#SATL_anom_wavelet<-wt(cbind(Rainfall_Tele_monthly[,"YearMonth"],Rainfall_Tele_monthly[,"SATL_anom"]),dj=1/12)
#plot(SATL_anom_wavelet,xlab="",ylab="SATL_anom (yrs)")

Rain_SATL_anom_wtc<-wtc(cbind(Rainfall_Tele_monthly[,"YearMonth"],
                              scale(Rainfall_Tele_monthly[,"Rainfall_interp_mm"])),
                        cbind(Rainfall_Tele_monthly[,"YearMonth"],
                              scale(Rainfall_Tele_monthly[,"SATL_anom"])),
                        pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Rain_SATL_anom_wtc$coi)) test<-cbind(test,Rain_SATL_anom_wtc$period>Rain_SATL_anom_wtc$coi[i])
test<-test[,-1]
Rain_SATL_anom_wtc$rsq[test]<-NA

Rain_wtc_global<-rbind(Rain_wtc_global,
                       data.frame(period=Rain_SATL_anom_wtc$period,
                                  global.coherence=rowMeans(Rain_SATL_anom_wtc$rsq,na.rm=T),
                                  x="Rainfall",
                                  y="SATL_anom"))

#IOD
#plot(cbind(Rainfall_Tele_monthly[,"YearMonth"],Rainfall_Tele_monthly[,"IOD"]),type="l")
#IOD_wavelet<-wt(cbind(Rainfall_Tele_monthly[,"YearMonth"],Rainfall_Tele_monthly[,"IOD"]),dj=1/12)
#plot(IOD_wavelet,xlab="",ylab="IOD (yrs)")

Rain_IOD_wtc<-wtc(cbind(Rainfall_Tele_monthly[,"YearMonth"],
                        scale(Rainfall_Tele_monthly[,"Rainfall_interp_mm"])),
                  cbind(Rainfall_Tele_monthly[,"YearMonth"],
                        scale(Rainfall_Tele_monthly[,"IOD"])),
                  pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Rain_IOD_wtc$coi)) test<-cbind(test,Rain_IOD_wtc$period>Rain_IOD_wtc$coi[i])
test<-test[,-1]
Rain_IOD_wtc$rsq[test]<-NA

Rain_wtc_global<-rbind(Rain_wtc_global,
                       data.frame(period=Rain_IOD_wtc$period,
                                  global.coherence=rowMeans(Rain_IOD_wtc$rsq,na.rm=T),
                                  x="Rainfall",
                                  y="IOD"))

### Temperature ###

## - Merge temperature and teleconnection datasets
Temperature_Tele_monthly<-merge(Temperature_monthly_G[,c("Month","Year","YearMonth","Season","Temperature_interp_c"),],
                                Tele[,c("Year","Month","MEI","NATL","NATL_anom","SATL","SATL_anom","IOD","NAO")],c("Year","Month"))

Temperature_Tele_monthly<-Temperature_Tele_monthly[order(Temperature_Tele_monthly$YearMonth),]
Temperature_Tele_monthly<-Temperature_Tele_monthly[complete.cases(Temperature_Tele_monthly$MEI),]

#MEI
#plot(cbind(Temperature_Tele_monthly[,"YearMonth"],Temperature_Tele_monthly[,"MEI"]),type="l")
#MEI_wavelet<-wt(cbind(Temperature_Tele_monthly[,"YearMonth"],Temperature_Tele_monthly[,"MEI"]),dj=1/12)
#plot(MEI_wavelet,xlab="",ylab="MEI (yrs)")

Temp_MEI_wtc<-wtc(cbind(Temperature_Tele_monthly[,"YearMonth"],
                        scale(Temperature_Tele_monthly[,"Temperature_interp_c"])),
                  cbind(Temperature_Tele_monthly[,"YearMonth"],
                        scale(Temperature_Tele_monthly[,"MEI"])),
                  pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Temp_MEI_wtc$coi)) test<-cbind(test,Temp_MEI_wtc$period>Temp_MEI_wtc$coi[i])
test<-test[,-1]
Temp_MEI_wtc$rsq[test]<-NA

Temp_wtc_global<-data.frame(period=Temp_MEI_wtc$period,
                            global.coherence=rowMeans(Temp_MEI_wtc$rsq,na.rm=T),
                            x="Temperature",
                            y="MEI")

#NATL_anom
# plot(cbind(Temperature_Tele_monthly[,"YearMonth"], Temperature_Tele_monthly[,"NATL_anom"]),type="l")
# NATL_anom_wavelet<-wt(cbind(Temperature_Tele_monthly[,"YearMonth"], Temperature_Tele_monthly[,"NATL_anom"]),dj=1/12)

Temp_NATL_anom_wtc<-wtc(cbind(Temperature_Tele_monthly[,"YearMonth"],
                              scale(Temperature_Tele_monthly[,"Temperature_interp_c"])),
                        cbind(Temperature_Tele_monthly[,"YearMonth"],
                              scale(Temperature_Tele_monthly[,"NATL_anom"])),
                        pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Temp_NATL_anom_wtc$coi)) test<-cbind(test,Temp_NATL_anom_wtc$period>Temp_NATL_anom_wtc$coi[i])
test<-test[,-1]
Temp_NATL_anom_wtc$rsq[test]<-NA

Temp_wtc_global<-rbind(Temp_wtc_global,
                       data.frame(period=Temp_NATL_anom_wtc$period,
                                  global.coherence=rowMeans(Temp_NATL_anom_wtc$rsq,na.rm=T),
                                  x="Temperature",
                                  y="NATL_anom"))

#SATL_anom
#plot(cbind(Temperature_Tele_monthly[,"YearMonth"],Temperature_Tele_monthly[,"SATL_anom"]),type="l")
#SATL_anom_wavelet<-wt(cbind(Temperature_Tele_monthly[,"YearMonth"], Temperature_Tele_monthly[,"SATL_anom"]),dj=1/12)

Temp_SATL_anom_wtc<-wtc(cbind(Temperature_Tele_monthly[,"YearMonth"],
                              scale(Temperature_Tele_monthly[,"Temperature_interp_c"])),
                        cbind(Temperature_Tele_monthly[,"YearMonth"],
                              scale(Temperature_Tele_monthly[,"SATL_anom"])),
                        pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Temp_SATL_anom_wtc$coi)) test<-cbind(test,Temp_SATL_anom_wtc$period>Temp_SATL_anom_wtc$coi[i])
test<-test[,-1]
Temp_SATL_anom_wtc$rsq[test]<-NA

Temp_wtc_global<-rbind(Temp_wtc_global,
                       data.frame(period=Temp_SATL_anom_wtc$period,
                                  global.coherence=rowMeans(Temp_SATL_anom_wtc$rsq,na.rm=T),
                                  x="Temperature",
                                  y="SATL_anom"))

#IOD
#plot(cbind(Temperature_Tele_monthly[,"YearMonth"],Temperature_Tele_monthly[,"IOD"]),type="l")
#IOD_wavelet<-wt(cbind(Temperature_Tele_monthly[,"YearMonth"],Temperature_Tele_monthly[,"IOD"]),dj=1/12)
#plot(IOD_wavelet,xlab="",ylab="IOD (yrs)")

Temp_IOD_wtc<-wtc(cbind(Temperature_Tele_monthly[,"YearMonth"],
                        scale(Temperature_Tele_monthly[,"Temperature_interp_c"])),
                  cbind(Temperature_Tele_monthly[,"YearMonth"],
                        scale(Temperature_Tele_monthly[,"IOD"])),
                  pad=T,mother="morlet",dj=1/12,nrands=1000)

#Remove wavelet coherence values within the cone of influence
test<-NA
for( i in 1:length(Temp_IOD_wtc$coi)) test<-cbind(test,Temp_IOD_wtc$period>Temp_IOD_wtc$coi[i])
test<-test[,-1]
Temp_IOD_wtc$rsq[test]<-NA

Temp_wtc_global<-rbind(Temp_wtc_global,
                       data.frame(period=Temp_IOD_wtc$period,
                                  global.coherence=rowMeans(Temp_IOD_wtc$rsq,na.rm=T),
                                  x="Temperature",
                                  y="IOD"))


## -  Figure 4 - Plot all wavelet coherency plots together

#Note on phase arrows on wavelet coherence plots
# Arrows pointing to the right mean that x and y are in phase.
# Arrows pointing to the left mean that x and y are in anti-phase.
# Arrows pointing up mean that y leads x by Ï€/2.
# Arrows pointing down mean that x leads y by Ï€/2.


layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = TRUE))
par(oma = c(0, 0, 0, 0), mar = c(2,4,2,2) )

#layout.show()

plot(Rain_MEI_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="A. Rainfall-MEI",adj=0)
plot(Temp_MEI_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="B. Temp-MEI",adj=0)
plot(Rain_NATL_anom_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="C. Rain-NATL",adj=0)
plot(Temp_NATL_anom_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="D. Temp-NATL",adj=0)
plot(Rain_SATL_anom_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="E. Rain-SATL",adj=0)
plot(Temp_SATL_anom_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="F. Temp-SATL",adj=0)
plot(Rain_IOD_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="G. Rain-IOD",adj=0)
plot(Temp_IOD_wtc, plot.cb = FALSE, plot.coi=TRUE,lwd.sig=1.5,plot.phase = TRUE,ylab="Period (yrs)")
title(main="H. Temp-IOD",adj=0)


## - Supplemental Figure 3 - Plot  global wavelet coherency together
Temp_Rain_wtc_global<-rbind(Rain_wtc_global,Temp_wtc_global)

ggplot(Temp_Rain_wtc_global,aes(y=global.coherence,x=period))+
  geom_line(aes(linetype=y,colour=y))+
  scale_color_brewer(type="div")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_log10(breaks=c(0.25,0.5,1,2,4,8))+
  #scale_x_continuous(trans=reverselog_trans(10),breaks=c(0.25,0.5,1,2,4,8))+
  theme_classic()+
  theme(legend.title=element_blank())+
  xlab("Period (years)")+
  ylab("Wavelet Coherence")+
  #coord_flip()+
  facet_wrap(~x,ncol=1)

