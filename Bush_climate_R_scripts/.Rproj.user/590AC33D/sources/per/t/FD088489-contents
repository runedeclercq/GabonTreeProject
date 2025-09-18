##Contents

##### 1. Data upload and preparation - line 26
##### 2. Seasonality - line 384
##### 3. Fourier spectra - line 631
##### 4. Long-term trends - line 706
##### 5. Periodicity over time - line 927
##### 6. Influence of oceans  - line 1003

## - Load libraries
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(timeDate)
library(scales)
library(biwavelet)
library(ggpubr)
library(corrplot)
library(lme4)
library("itsadug")
library("tidyr")
library("broom")
library("gridExtra")
library(cplm)
library(biwavelet)

##########################################
##### 1. Data upload and preparation #####
##########################################

### Rainfall ###

## - Data available for download at the University of Stirling's DataSTORRE (http://hdl.handle.net/11667/133)

## - Upload daily rainfall data (Dataset C - a combination of datasets A and B with adjusted data)
Rainfall_daily_C <- read.csv("Bush_data/Rainfall_daily_v2017-12-31.csv")

Rainfall_daily_C$Date<-as.Date(Rainfall_daily_C$Date)
Rainfall_daily_C<-Rainfall_daily_C[order(Rainfall_daily_C$Date),]

## - Calculate mean value for each day in the time series and interpolate missing values using the 10 day rolling mean (to create Dataset D)
Rainfall_daily_D<-ddply(Rainfall_daily_C,.(Date),summarise,Rainfall_daily_mm=mean(Rainfall_adjusted_mm,na.rm=T))

# - Create full calendar of dates
Rainfall_daily_D<-merge(Rainfall_daily_D,
                        data.frame(Date=seq.Date(from=min(Rainfall_daily_D$Date,na.rm=T),to=max(Rainfall_daily_D$Date,na.rm=T),by="day"))
                        ,"Date",all.y=T)
Rainfall_daily_D$Data<-as.factor(ifelse(is.na(Rainfall_daily_D$Rainfall_daily_mm),"Missing","Data"))

# - Interpolate missing data
Rainfall_daily_D$MeanRain_10day<-rollapply(Rainfall_daily_D$Rainfall_daily_mm,10,mean,na.rm=T,fill=NA)
Rainfall_daily_D$Rainfall_interp_mm<-Rainfall_daily_D$Rainfall_daily_mm
Rainfall_daily_D$Rainfall_interp_mm[is.na(Rainfall_daily_D$Rainfall_interp_mm)]<-Rainfall_daily_D$MeanRain_10day[is.na(Rainfall_daily_D$Rainfall_interp_mm)]
Rainfall_daily_D$Comments<-ifelse(is.na(Rainfall_daily_D$Rainfall_daily_mm),"Interpolated from 10 day running mean", "NA")
Rainfall_daily_D<-data.frame(Date=Rainfall_daily_D$Date, Rainfall_mm=Rainfall_daily_D$Rainfall_daily_mm,Rainfall_interp_mm=Rainfall_daily_D$Rainfall_interp_mm, Comments=Rainfall_daily_D$Comments, Site="SEGC")
Rainfall_daily_D$DOY<-as.integer(format(Rainfall_daily_D$Date,"%j"))
Rainfall_daily_D$Month<-month(Rainfall_daily_D$Date)
Rainfall_daily_D$Season<-mapvalues(Rainfall_daily_D$Month,from=c(1:12),to=c("DJF","DJF","MAM","MAM","MAM","JJAS","JJAS","JJAS","JJAS","ON","ON","DJF"))
Rainfall_daily_D$Year<-year(Rainfall_daily_D$Date)
Rainfall_daily_D$Year_r<-scale(Rainfall_daily_D$Year)

## - Sum rainfall for each month in the time series from interpolated daily data and interpolate missing months (to create Dataset E)
Rainfall_monthly_E<-ddply(Rainfall_daily_D,.(Year,Season,Month,YearMonth=as.yearmon(Date)),summarise, 
                          Rainfall_daily_mm=mean(Rainfall_mm,na.rm=T),
                          Rainfall_mm=sum(Rainfall_interp_mm,na.rm=T),
                          sample=length(Date[complete.cases(Rainfall_interp_mm)]),
                          sample_max=length(Date))

# - Remove monthly values with more than 10% missing days
Rainfall_monthly_E$sample_prop<-Rainfall_monthly_E$sample/Rainfall_monthly_E$sample_max
Rainfall_monthly_E$Rainfall_mm[Rainfall_monthly_E$sample_prop<0.9]<-NA
Rainfall_monthly_E$Rainfall_daily_mm[Rainfall_monthly_E$sample_prop<0.9]<-NA

# - Complete missing values with mean value for corresponding calendar month
Rainfall_monthly_E<-merge(Rainfall_monthly_E,
                          ddply(Rainfall_monthly_E,.(Month=month(as.Date(YearMonth))),summarise,
                                Rainfall_month_mm=mean(Rainfall_mm,na.rm=T)), "Month")
Rainfall_monthly_E$Rainfall_interp_mm<-Rainfall_monthly_E$Rainfall_mm
Rainfall_monthly_E$Rainfall_interp_mm[is.na(Rainfall_monthly_E$Rainfall_interp_mm)]<-Rainfall_monthly_E$Rainfall_month_mm[is.na(Rainfall_monthly_E$Rainfall_interp_mm)]
Rainfall_monthly_E<-Rainfall_monthly_E[order(Rainfall_monthly_E$YearMonth),]
Rainfall_monthly_E$Rainfall_standard_mm<-scale(Rainfall_monthly_E$Rainfall_interp_mm)

## - Sum rainfall for each year in the timeseries from interpolated daily data (to create Dataset F)
Rainfall_yearly_F<-ddply(Rainfall_daily_D,.(Year=year(Rainfall_daily_D$Date)),summarise, 
                         Rainfall_daily_mm=mean(Rainfall_interp_mm,na.rm=T),
                         Rainfall_yearly_mm=sum(Rainfall_interp_mm))
Rainfall_yearly_F$Year_r<-scale(Rainfall_yearly_F$Year)

## - Summary of rainfall dataframes (match flowchart in supporting information)
summary(Rainfall_daily_C) #Lopee daily rainfall dataset including original and adjusted data (Dataset C)
summary(Rainfall_daily_D) #Lopee mean daily rainfall time series with non-interpolated and interpolated data (Dataset D)
summary(Rainfall_monthly_E) #Lopee monthly rainfall time series with non-interpolated and interpolated data (Dataset E)
summary(Rainfall_yearly_F) #Lopee yearly rainfall time series (Dataset F)


### Temperature ####

## - Data available for download at the University of Stirling's DataSTORRE (http://hdl.handle.net/11667/133)

## - Upload daily max/min temperature data (Dataset C - a combination of datasets A and B)
Temperature_daily_C <- read.csv("Bush_data/Temperature_daily_v2018-03-13.csv")

Temperature_daily_C$Date<-as.Date(Temperature_daily_C$Date)

## - Calculate mean value for each day in the time series (to create Dataset D)
Temperature_daily_D<-ddply(Temperature_daily_C,.(Site,Type,Date),summarise,
                           Temperature_c=mean(Temperature_c,na.rm=T))

# - Create full calendar of dates
Temperature_daily_D<-rbind(merge(Temperature_daily_D[Temperature_daily_D$Site=="Saline",],
                                 expand.grid(Date=seq.Date(from=min(Temperature_daily_D$Date[Temperature_daily_D$Site=="Saline"],na.rm=T),to=max(Temperature_daily_D$Date[Temperature_daily_D$Site=="Saline"],na.rm=T),by="day"),
                                             Type=c("Max","Min"),Site=c("Saline")),c("Date","Type","Site"),all=T),
                           merge(Temperature_daily_D[Temperature_daily_D$Site=="SEGC",],
                                 expand.grid(Date=seq.Date(from=min(Temperature_daily_D$Date[Temperature_daily_D$Site=="SEGC"],na.rm=T),to=max(Temperature_daily_D$Date[Temperature_daily_D$Site=="SEGC"],na.rm=T),by="day"),
                                             Type=c("Max","Min"),
                                             Site=c("SEGC")),
                                 c("Date","Type","Site"),all=T))

Temperature_daily_D$DOY<-as.integer(format(Temperature_daily_D$Date,"%j"))
Temperature_daily_D$YearMonth<-as.yearmon(Temperature_daily_D$Date)
Temperature_daily_D$Month<-month(Temperature_daily_D$Date)
Temperature_daily_D$Year<-year(Temperature_daily_D$Date)
Temperature_daily_D$Site2<-ifelse(Temperature_daily_D$Site=="Saline","Forest","Savanna")

## - Calculate monthly mean daily max/min temperature time series for both sites (Dataset E)
Temperature_monthly_E<-ddply(Temperature_daily_D[Temperature_daily_D$Year<2018,],.(Site,Type,Year,Month,YearMonth),summarise,
                             Sample_Lopee=length(Temperature_c[complete.cases(Temperature_c)]),
                             Temperature_c=mean(Temperature_c,na.rm=T))

Temperature_monthly_E$Temperature_c[Temperature_monthly_E$Sample_Lopee<5]<-NA
Temperature_monthly_E$Date<-as.Date(paste(year(as.Date(Temperature_monthly_E$YearMonth)),month(as.Date(Temperature_monthly_E$YearMonth)),"15",sep="-"))

# - Complete missing values with mean value for corresponding calendar month
Temperature_monthly_E<-merge(Temperature_monthly_E,
                             ddply(Temperature_monthly_E,.(Site,Type,Month),summarise,
                                   Temperature_month_mean_c=mean(Temperature_c,na.rm=T)),
                             c("Site","Type","Month"))

Temperature_monthly_E$Temperature_interp_c<-Temperature_monthly_E$Temperature_c
Temperature_monthly_E$Temperature_interp_c[is.na(Temperature_monthly_E$Temperature_interp_c)]<-Temperature_monthly_E$Temperature_month_mean_c[is.na(Temperature_monthly_E$Temperature_interp_c)]
Temperature_monthly_E<-Temperature_monthly_E[order(Temperature_monthly_E$YearMonth),]
Temperature_monthly_E$Temperature_standard_c<-NA
Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="Saline"&Temperature_monthly_E$Type=="Max"]<-scale(Temperature_monthly_E$Temperature_interp_c[Temperature_monthly_E$Site=="Saline"&Temperature_monthly_E$Type=="Max"])
Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="SEGC"&Temperature_monthly_E$Type=="Max"]<-scale(Temperature_monthly_E$Temperature_interp_c[Temperature_monthly_E$Site=="SEGC"&Temperature_monthly_E$Type=="Max"])
Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="Saline"&Temperature_monthly_E$Type=="Min"]<-scale(Temperature_monthly_E$Temperature_interp_c[Temperature_monthly_E$Site=="Saline"&Temperature_monthly_E$Type=="Min"])
Temperature_monthly_E$Temperature_standard_c[Temperature_monthly_E$Site=="SEGC"&Temperature_monthly_E$Type=="Min"]<-scale(Temperature_monthly_E$Temperature_interp_c[Temperature_monthly_E$Site=="SEGC"&Temperature_monthly_E$Type=="Min"])

## - Calculate long-term series of minimum daily temperature combined from both sites and all equipment (Dataset F)
Temperature_daily_F<-ddply(Temperature_daily_D[Temperature_daily_D$Type=="Min",],.(Year, Month=month(Date), YearMonth, DOY,Date),summarise,
                           Temperature_c=mean(Temperature_c,na.rm=T))
Temperature_daily_F<-Temperature_daily_F[Temperature_daily_F$Year<2018,]

# - Interpolate missing data
Temperature_daily_F<-Temperature_daily_F
Temperature_daily_F$MeanTemp_10day<-rollapply(Temperature_daily_F$Temperature_c,10,mean,na.rm=T,fill=NA)
Temperature_daily_F$Temperature_interp_c<-Temperature_daily_F$Temperature_c
Temperature_daily_F$Temperature_interp_c[is.na(Temperature_daily_F$Temperature_c)]<-Temperature_daily_F$MeanTemp_10day[is.na(Temperature_daily_F$Temperature_c)]
Temperature_daily_F$Comments<-ifelse(is.na(Temperature_daily_F$Temperature_c),"Interpolated from 10 day running mean", "NA")

Temperature_daily_F$Year_r<-scale(Temperature_daily_F$Year)
Temperature_daily_F$Season<-mapvalues(Temperature_daily_F$Month,from=c(1:12),to=c("DJF","DJF","MAM","MAM","MAM","JJAS","JJAS","JJAS","JJAS","ON","ON","DJF"))

## - Calculate mean monthly minimum temperaure from combined daily dataset (to create dataset G)
Temperature_monthly_G<-ddply(Temperature_daily_F,.(Year,Season,Month,YearMonth),summarise,
                             Sample_Lopee=length(Temperature_c[complete.cases(Temperature_c)]),
                             Temperature_c=mean(Temperature_c,na.rm=T))
Temperature_monthly_G$Temperature_c[Temperature_monthly_G$Sample_Lopee<5]<-NA
Temperature_monthly_G$Date<-as.Date(paste(year(as.Date(Temperature_monthly_G$YearMonth)),month(as.Date(Temperature_monthly_G$YearMonth)),"15",sep="-"))

# - Complete missing values with mean value for corresponding calendar month
Temperature_monthly_G<-merge(Temperature_monthly_G,
                             ddply(Temperature_monthly_G,.(Month),summarise,
                                   Temperature_month_mean_c=mean(Temperature_c,na.rm=T)),
                             c("Month"))

Temperature_monthly_G$Temperature_interp_c<-Temperature_monthly_G$Temperature_c
Temperature_monthly_G$Temperature_interp_c[is.na(Temperature_monthly_G$Temperature_interp_c)]<-Temperature_monthly_G$Temperature_month_mean_c[is.na(Temperature_monthly_G$Temperature_interp_c)]
Temperature_monthly_G<-Temperature_monthly_G[order(Temperature_monthly_G$YearMonth),]

## - Summary of temperature datasets
summary(Temperature_daily_C) #Lopee daily temperature dataset including original data from each equipment at each site (Dataset C)
summary(Temperature_daily_D)  #Lopee mean daily temperature dataset combining data from different equipment at each site (Dataset D)
summary(Temperature_monthly_E) #Lopee mean monthly temperature dataset for each site including non-interpolated and inteprolated data (Dataset E)
summary(Temperature_daily_F) #Lopee mean daily minimum temperature dataset combining non-interpolated data from each site (Dataset F)
summary(Temperature_monthly_G) #Lopee mean monthly minimum daily temperature dataset for both sites combined including non-interpolated and inteprolated data (Dataset G)


### Berkeley Temp ###

## - download available from http://climexp.knmi.nl/start.cgi for the grid-cell overlapping the SEGC location (0.2N, 11.6E).

Berkeley_daily <- read.csv("Bush_data/iberkeley_tmin_daily_full_11.6E_0.2N_1983-2019_n.dat.csv", 
                           sep = "",        # tab-separated
                           comment.char = "#",# skip metadata lines
                           header = FALSE)    # no column names in the file

names(Berkeley_daily)<-c("Date","Temperature_c")
Berkeley_daily$Date<-as.Date(as.character(Berkeley_daily$Date),"%Y%m%d")
Berkeley_daily$DOY<-format(Berkeley_daily$Date,"%j")
Berkeley_daily$Year<-year(Berkeley_daily$Date)
Berkeley_daily$Month<-month(Berkeley_daily$Date)
Berkeley_daily$Year_r<-scale(Berkeley_daily$Year)
Berkeley_daily<-Berkeley_daily[Berkeley_daily$Year>1983&Berkeley_daily$Year<2018,]


### CRU Temp ###

## - downloaded from http://climexp.knmi.nl/start.cgi for the grid-cell overlapping the SEGC location (0.2N, 11.6E).

CRU_monthly <- read.table("Bush_data/icru4_tmn_11.6E_0.2N_firstyear-lastyear_n.dat.txt", 
                           sep = "",        # tab-separated
                           comment.char = "#",# skip metadata lines
                           header = FALSE)    # no column names in the file

names(CRU_monthly)<-c("Year",1:12)
CRU_monthly<-gather(CRU_monthly,"Month","Temperature_c",2:13)
CRU_monthly$YearMonth<-as.yearmon(paste(CRU_monthly$Year,CRU_monthly$Month),"%Y %m")
CRU_monthly$Year_r<-scale(CRU_monthly$Year)
CRU_monthly<-CRU_monthly[CRU_monthly$Year>1983,]


### Absolute humidity ###

## - Data available for download at the University of Stirling's DataSTORRE (http://hdl.handle.net/11667/133)

## - Upload daily absolute humidity record (Dataset A)
AbsHumidity_daily_A <- read.csv("Bush_data/Humidity_daily_v2018-03-14.csv")

AbsHumidity_daily_A$Date<-as.Date(AbsHumidity_daily_A$Date,format="%d/%m/%y")
AbsHumidity_daily_A<-AbsHumidity_daily_A[complete.cases(AbsHumidity_daily_A$Date),]
AbsHumidity_daily_A$Site2<-ifelse(AbsHumidity_daily_A$Site=="Saline","Forest","Savanna")

AbsHumidity_daily_B<-ddply(AbsHumidity_daily_A[AbsHumidity_daily_A$Night==T,],.(Site=Site2,Date),summarise,AbsoluteHumidity_gm3=mean(AbsoluteHumidity_gm3,na.rm=T))
AbsHumidity_daily_B<-merge(expand.grid(Date=seq(from=min(AbsHumidity_daily_B$Date),to=max(AbsHumidity_daily_B$Date),"day"),
                                       Site=c("Savanna","Forest")),
                           AbsHumidity_daily_B,c("Date","Site"),all.x=T)
AbsHumidity_daily_B$Year<-year(AbsHumidity_daily_B$Date)
AbsHumidity_daily_B$Month<-month(AbsHumidity_daily_B$Date)
AbsHumidity_daily_B$YearMonth<-as.yearmon(AbsHumidity_daily_B$Date)
AbsHumidity_daily_B$DOY<-as.integer(format(AbsHumidity_daily_B$Date,"%j"))

## - Calculate monthly mean time series for both sites (Dataset C)
AbsHumidity_monthly_C<-ddply(AbsHumidity_daily_B[AbsHumidity_daily_B$Year<2018,],.(Site,Year, Month,YearMonth),summarise,
                             Sample_Lopee=length(AbsoluteHumidity_gm3[complete.cases(AbsoluteHumidity_gm3)]),
                             AbsoluteHumidity_gm3=mean(AbsoluteHumidity_gm3,na.rm=T))

AbsHumidity_monthly_C$Date<-as.Date(paste(AbsHumidity_monthly_C$Year,AbsHumidity_monthly_C$Month,"15",sep="-"))
AbsHumidity_monthly_C$AbsoluteHumidity_gm3[AbsHumidity_monthly_C$Sample_Lopee<5]<-NA

#Complete missing values with mean value for corresponding calendar month
AbsHumidity_monthly_C<-merge(AbsHumidity_monthly_C,
                             ddply(AbsHumidity_monthly_C,.(Site,Month),summarise,AbsoluteHumidity_month_mean_gm3=mean(AbsoluteHumidity_gm3,na.rm=T)),
                             c("Site","Month"))

AbsHumidity_monthly_C$AbsHumidity_interp_gm3<-AbsHumidity_monthly_C$AbsoluteHumidity_gm3
AbsHumidity_monthly_C$AbsHumidity_interp_gm3[is.na(AbsHumidity_monthly_C$AbsoluteHumidity_gm3)]<-AbsHumidity_monthly_C$AbsoluteHumidity_month_mean_gm3[is.na(AbsHumidity_monthly_C$AbsoluteHumidity_gm3)]
AbsHumidity_monthly_C$AbsHumidity_standard_gm3<-scale(AbsHumidity_monthly_C$AbsHumidity_interp_gm3)
AbsHumidity_monthly_C<-AbsHumidity_monthly_C[order(AbsHumidity_monthly_C$Site,AbsHumidity_monthly_C$Date),]

summary(AbsHumidity_daily_A)
summary(AbsHumidity_daily_B)
summary(AbsHumidity_monthly_C)


### Wind ###

## - Data available for download at the University of Stirling's DataSTORRE (http://hdl.handle.net/11667/133)

## - Upload daily mean wind speed record (Dataset A)
Wind_daily_A <- read.csv("Bush_data/Wind_daily_v2016-07-01.csv")
names(Wind_daily_A)[3]<-"Wind_m_s"
Wind_daily_A$Date<-as.Date(Wind_daily_A$Date)

Wind_daily_B<-ddply(Wind_daily_A,.(Date),summarise,
                    Wind_m_s=mean(Wind_m_s,na.rm=T))
Wind_daily_B<-merge(data.frame(Date=seq.Date(from=min(Wind_daily_B$Date,na.rm=T),to=max(Wind_daily_B$Date,na.rm=T),by="day")),
                    Wind_daily_B,"Date",all.x=T)
Wind_daily_B$Year<-year(Wind_daily_B$Date)
Wind_daily_B$Month<-month(Wind_daily_B$Date)
Wind_daily_B$YearMonth<-as.yearmon(Wind_daily_B$Date)
Wind_daily_B$DOY<-as.integer(format(Wind_daily_B$Date,"%j"))

## - Calculate monthly mean time series for both sites (Dataset C)
Wind_monthly_C<-ddply(Wind_daily_B[Wind_daily_B$Year<2018,],.(Year, Month,YearMonth),summarise,
                      Sample_Lopee=length(Wind_m_s[complete.cases(Wind_m_s)]),
                      Wind_m_s=mean(Wind_m_s,na.rm=T))
Wind_monthly_C$Date<-as.Date(paste(Wind_monthly_C$Year,Wind_monthly_C$Month,"15",sep="-"))
Wind_monthly_C$Wind_m_s[Wind_monthly_C$Sample_Lopee<5]<-NA

# - Complete missing values with mean value for corresponding calendar month
Wind_monthly_C<-merge(Wind_monthly_C,
                      ddply(Wind_monthly_C,.(Month),summarise,Wind_month_mean_m_s=mean(Wind_m_s,na.rm=T)),
                      c("Month"))

Wind_monthly_C$Wind_interp_m_s<-Wind_monthly_C$Wind_m_s
Wind_monthly_C$Wind_interp_m_s[is.na(Wind_monthly_C$Wind_m_s)]<-Wind_monthly_C$Wind_month_mean_m_s[is.na(Wind_monthly_C$Wind_m_s)]
Wind_monthly_C$Wind_standard_m_s<-scale(Wind_monthly_C$Wind_interp_m_s)
Wind_monthly_C<-Wind_monthly_C[order(Wind_monthly_C$Date),]

summary(Wind_daily_A)
summary(Wind_daily_B)
summary(Wind_monthly_C)


### Solar ###

## - Data available for download at the University of Stirling's DataSTORRE (http://hdl.handle.net/11667/133)

## - Upload daily mean solar radiation record (Dataset A)
Solar_daily_A<- read.csv("Bush_data/Solar_daily_v2015-10-30.csv")

names(Solar_daily_A)[3]<-"Solar_Wm2"
Solar_daily_A$Date<-as.Date(Solar_daily_A$Date)
Solar_daily_A$Month<-month(Solar_daily_A$Date)

Solar_daily_B<-ddply(Solar_daily_A,.(Date),summarise,Solar_Wm2 = mean(Solar_Wm2,na.rm=T))
Solar_daily_B<-merge(data.frame(Date=seq(from=min(Solar_daily_B$Date),to=max(Solar_daily_B$Date),"day")),
                     Solar_daily_B,c("Date"),all.x=T)
Solar_daily_B$Year<-year(Solar_daily_B$Date)
Solar_daily_B$Month<-month(Solar_daily_B$Date)
Solar_daily_B$YearMonth<-as.yearmon(Solar_daily_B$Date)
Solar_daily_B$DOY<-as.integer(format(Solar_daily_B$Date,"%j"))
Solar_daily_B$Solar_Wm2[which(Solar_daily_B$Solar_Wm2<10)]<-NA

## - Calculate monthly mean time series for both sites (Dataset C)
Solar_monthly_C<-ddply(Solar_daily_B[Solar_daily_B$Year<2018,],.(Year, Month,YearMonth),summarise,
                       Sample_Lopee=length(Solar_Wm2[complete.cases(Solar_Wm2)]),
                       Solar_Wm2=mean(Solar_Wm2,na.rm=T))

Solar_monthly_C$Date<-as.Date(paste(Solar_monthly_C$Year,Solar_monthly_C$Month,"15",sep="-"))
Solar_monthly_C$Solar_Wm2[Solar_monthly_C$Sample_Lopee<5]<-NA

# - Complete missing values with mean value for corresponding calendar month
Solar_monthly_C<-merge(Solar_monthly_C,
                       ddply(Solar_monthly_C,.(Month),summarise,Solar_month_mean_Wm2=mean(Solar_Wm2,na.rm=T)),
                       c("Month"))

Solar_monthly_C$Solar_interp_Wm2<-Solar_monthly_C$Solar_Wm2
Solar_monthly_C$Solar_interp_Wm2[is.na(Solar_monthly_C$Solar_Wm2)]<-Solar_monthly_C$Solar_month_mean_Wm2[is.na(Solar_monthly_C$Solar_Wm2)]
Solar_monthly_C$Solar_standard_Wm2<-scale(Solar_monthly_C$Solar_interp_Wm2)
Solar_monthly_C<-Solar_monthly_C[order(Solar_monthly_C$Date),]

summary(Solar_daily_A)
summary(Solar_daily_B)
summary(Solar_monthly_C)


### AOT ###

## - Download available from  the NASA Aerosol Robotic Network (Aeronet; https://aeronet.gsfc.nasa.gov/; Holben et al. 1998). 

library(data.table)
Aerosol_daily_A <- fread("Bush_data/20140101_20181231_SEGC_Lope_Gabon.lev20", 
                         skip = "Date(dd:mm:yyyy)",
                         na.strings = c("-999", "-999."))

Aerosol_daily_A$Date<-Aerosol_daily_A$"Date(dd:mm:yyyy)"
Aerosol_daily_A$Date<-as.Date(as.character(Aerosol_daily_A$Date),"%d:%m:%Y")

Aerosol_daily_B<-merge(data.frame(Date=seq.Date(from=min(Aerosol_daily_A$Date,na.rm=T),to=max(Aerosol_daily_A$Date,na.rm=T),by="day")),
                       Aerosol_daily_A[,c("Date", "Day_of_Year", "AOD_440nm", "AOD_500nm", "AOD_675nm", "Data_Quality_Level")],"Date",all.x=T)

Aerosol_daily_B$Year<-year(Aerosol_daily_B$Date)
Aerosol_daily_B$Month<-month(Aerosol_daily_B$Date)
Aerosol_daily_B$YearMonth<-as.yearmon(Aerosol_daily_B$Date)
Aerosol_daily_B$DOY<-as.integer(format(Aerosol_daily_B$Date,"%j"))

#summary(Aerosol_daily_B$Month[Aerosol_daily_B$Data=="Missing"])
#642 days missing out of 1045 possible days between 2014-04-27 and 2017-03-06 (61% missing).
#More than 70% values missing in months June-November (peaking in August, 85%). Presumably because of cloud cover data removed at quality control.
#Most data available in March (35% missing vales)

## - Calculate monthly mean time series for both sites (Dataset C)
Aerosol_monthly_C<-ddply(Aerosol_daily_B[Aerosol_daily_B$Year<2018,],.(Year, Month,YearMonth),summarise,
                         Sample_Lopee=length(AOD_500nm[complete.cases(AOD_500nm)]),
                         AOT_500=mean(AOD_500nm,na.rm=T))
Aerosol_monthly_C$Date<-as.Date(paste(Aerosol_monthly_C$Year,Aerosol_monthly_C$Month,"15",sep="-"))
Aerosol_monthly_C$AOT_500[Aerosol_monthly_C$Sample_Lopee<5]<-NA

# - Complete missing values with mean value for corresponding calendar month
Aerosol_monthly_C<-merge(Aerosol_monthly_C,
                         ddply(Aerosol_monthly_C,.(Month),summarise,AOT_500_month_mean=mean(AOT_500,na.rm=T)),
                         c("Month"))

Aerosol_monthly_C$AOT_500_interp<-Aerosol_monthly_C$AOT_500
Aerosol_monthly_C$AOT_500_interp[is.na(Aerosol_monthly_C$AOT_500)]<-Aerosol_monthly_C$AOT_500_month_mean[is.na(Aerosol_monthly_C$AOT_500)]
Aerosol_monthly_C$AOT_500_standard<-scale(Aerosol_monthly_C$AOT_500_interp)
Aerosol_monthly_C<-Aerosol_monthly_C[order(Aerosol_monthly_C$Date),]

summary(Aerosol_daily_A)
summary(Aerosol_daily_B)
summary(Aerosol_monthly_C)

##########################
##### 2. Seasonality #####
##########################

### Rainfall ###

## - Calculate rainfall means for each DOY and month
Rainfall_DOY<-ddply(Rainfall_daily_D,.(DOY),summarise,
                    MeanRain_mm=mean(Rainfall_interp_mm,na.rm=T))

Rainfall_DOY$MeanRain_10day_mm<-rollapply(Rainfall_DOY$MeanRain_mm,10,mean,partial=T)
Rainfall_DOY$Dayof2019<-as.Date(as.character(Rainfall_DOY$DOY),"%j")
Rainfall_DOY$Month<-month(Rainfall_DOY$Dayof2019)
Rainfall_DOY<-merge(Rainfall_DOY,ddply(Rainfall_DOY,.(Month),summarise,MeanRain_month_mm=mean(MeanRain_mm,na.rm=T)),"Month")
Rainfall_DOY$month_labels <- mapvalues(Rainfall_DOY$Month,from=1:12,to=c('J','F','M','A','M','J','J','A','S','O','N','D'))

# - Seasonality plot (Figure 2)
Rainfall_seasonal_plot<-ggplot(Rainfall_DOY,aes(x=Dayof2019))+
  geom_line(data=Rainfall_DOY,aes(y=MeanRain_mm),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(data=Rainfall_DOY, aes(y=MeanRain_10day_mm),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(data=Rainfall_DOY, aes(y=MeanRain_month_mm),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  theme_classic()+
  ylab("Rainfall (mm/day)")+
  xlab("Month")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  #ggtitle("A.")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))

Rainfall_seasonal_plot

### Temperature ###

## - Calculate temprature means for each DOY and month
Temperature_DOY<-ddply(Temperature_daily_D,.(Site,Type,Month=month(Date),DOY),summarise,
                       Temperature_c=mean(Temperature_c,na.rm=T))

Temperature_DOY<-data.frame(Temperature_DOY,Temperature_c_10days=ddply(Temperature_DOY,.(Site,Type),summarise,
                                                                       Temperature_c_10days=rollapply(Temperature_c,10,mean,partial=T))[,3])

Temperature_DOY<-merge(Temperature_DOY,ddply(Temperature_DOY,.(Site,Type,Month),summarise,
                                             Temperature_c_month=mean(Temperature_c,na.rm=T)),
                       c("Site","Type","Month"))

Temperature_DOY$Dayof2019<-as.Date(as.character(Temperature_DOY$DOY),"%j")
Temperature_DOY$Site<-mapvalues(Temperature_DOY$Site,from=c("SEGC","Saline"),to=c("Savanna","Forest"))

## - Seasonality plot (Figure 2)
MaxMin_seasonal_plot_forest<-ggplot(Temperature_DOY[Temperature_DOY$Site=="Forest",])+
  geom_line(aes(x=Dayof2019,y=Temperature_c,group=Type),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Temperature_c_10days,group=Type),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Temperature_c_month,group=Type),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  theme_classic()+
  ylab("Temperature (c)")+
  xlab("Month")+
  ylim(c(20,35))+
  #ggtitle("B.")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
  geom_label(label="Forest",x=as.Date("2019-12-01"),y=34.5)
MaxMin_seasonal_plot_forest


MaxMin_seasonal_plot_savanna<-ggplot(Temperature_DOY[Temperature_DOY$Site=="Savanna",])+
  geom_line(aes(x=Dayof2019,y=Temperature_c,group=Type),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Temperature_c_10days,group=Type),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Temperature_c_month,group=Type),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  theme_classic()+
  ylab("Temperature (c)")+
  xlab("Month")+
  ylim(c(20,35))+
  #ggtitle("B.")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
  geom_label(label="Savanna",x=as.Date("2019-12-01"),y=34.5)
MaxMin_seasonal_plot_savanna

### Absolute Humidity ###

Humidity_DOY<-ddply(AbsHumidity_daily_B,.(Site,DOY),summarise,
                    AbsoluteHumidity_gm3=mean(AbsoluteHumidity_gm3,na.rm=T))
Humidity_DOY$Date<-as.Date(as.character(Humidity_DOY$DOY),"%j")
Humidity_DOY$Month<-month(Humidity_DOY$Date)
Humidity_DOY<-data.frame(Humidity_DOY,AbsoluteHumidity_gm3_10days=ddply(Humidity_DOY,.(Site),summarise,
                                                                        AbsoluteHumidity_gm3_10days=rollapply(AbsoluteHumidity_gm3,10,mean,partial=T))[,2])
Humidity_DOY<-merge(Humidity_DOY,ddply(Humidity_DOY,.(Site,Month),summarise,
                                       AbsoluteHumidity_gm3_month=mean(AbsoluteHumidity_gm3,na.rm=T)), c("Site","Month"))

#ddply(Humidity_DOY, .(Site), summarise,
#      Min=min(AbsoluteHumidity_gm3_month),
#      Max=max(AbsoluteHumidity_gm3_month))

Humidity_DOY$Dayof2019<-as.Date(as.character(Humidity_DOY$DOY),"%j")

## - Seasonality plot (Figure 2)
AH_seasonal_plot_forest<-ggplot(Humidity_DOY[Humidity_DOY$Site=="Forest",])+
  geom_line(aes(x=Dayof2019,y=AbsoluteHumidity_gm3),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=AbsoluteHumidity_gm3_10days),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=AbsoluteHumidity_gm3_month),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  theme_classic()+
  ylim(c(17.5,23.5))+
  ylab("Absolute Humidity (g/m3)")+
  xlab("Month")+
  #ggtitle("B.")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
  geom_label(label="Forest",x=as.Date("2019-12-01"),y=23)
AH_seasonal_plot_forest

AH_seasonal_plot_savanna<-ggplot(Humidity_DOY[Humidity_DOY$Site=="Savanna",])+
  geom_line(aes(x=Dayof2019,y=AbsoluteHumidity_gm3),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=AbsoluteHumidity_gm3_10days),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=AbsoluteHumidity_gm3_month),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  theme_classic()+
  ylim(c(17.5,23.5))+
  ylab("Absolute Humidity (g/m3)")+
  xlab("Month")+
  #ggtitle("B.")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
  geom_label(label="Savanna",x=as.Date("2019-12-01"),y=23)
AH_seasonal_plot_savanna

### Solar ###

Solar_DOY<-ddply(Solar_daily_B,.(Month=month(Date),DOY),summarise,
                 Solar_Wm2=mean(Solar_Wm2,na.rm=T))
Solar_DOY$Solar_Wm2_10days<-rollapply(Solar_DOY$Solar_Wm2,10, mean,partial=T)
Solar_DOY<-merge(Solar_DOY,ddply(Solar_DOY,.(Month),summarise,Solar_Wm2_month=mean(Solar_Wm2,na.rm=T)),"Month")
Solar_DOY$Dayof2019<-as.Date(as.character(Solar_DOY$DOY),"%j")

## - Seasonality plot (Figure 2)
Solar_seasonal_plot<-ggplot(Solar_DOY,aes(x=Dayof2019,y=Solar))+
  geom_line(aes(x=Dayof2019,y=Solar_Wm2),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Solar_Wm2_10days),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Solar_Wm2_month),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  theme_classic()+
  ylab("Solar Radiation (W/m2)")+
  xlab("Month")+
  #  ggtitle("D.")+
  theme(legend.title=element_blank(),legend.position="none")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))
Solar_seasonal_plot


### Wind ###

Wind_DOY<-ddply(Wind_daily_B,.(Month=month(Date),DOY),summarise,
                Wind_m_s=mean(Wind_m_s,na.rm=T))
Wind_DOY$Wind_m_s_10days<-rollapply(Wind_DOY$Wind_m_s,10, mean,partial=T)
Wind_DOY<-merge(Wind_DOY,ddply(Wind_DOY,.(Month),summarise,Wind_m_s_month=mean(Wind_m_s,na.rm=T)),"Month")
Wind_DOY$Dayof2019<-as.Date(as.character(Wind_DOY$DOY),"%j")

## - Seasonality plot (Figure 2)
Wind_seasonal_plot<-ggplot(Wind_DOY,aes(x=Dayof2019,y=Solar))+
  geom_line(aes(x=Dayof2019,y=Wind_m_s),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Wind_m_s_10days),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=Wind_m_s_month),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  theme_classic()+
  ylab("Wind Speed (m/s)")+
  xlab("Month")+
  theme(legend.title=element_blank(),legend.position="none")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
  #ggtitle("E.")+
  guides(linetype=guide_legend(nrow=2),alpha=guide_legend(nrow=2),colour=guide_legend(nrow=2))

Wind_seasonal_plot


### Aerosol optical thickness ###

Aerosol_DOY<-ddply(Aerosol_daily_B,.(Month,DOY),summarise,
                   MeanAOT_440=mean(AOD_440nm,na.rm=T),
                   MeanAOT_500=mean(AOD_500nm,na.rm=T),
                   MeanAOT_675=mean(AOD_675nm,na.rm=T))

Aerosol_DOY$MeanAOT_440_10days<-rollapply(Aerosol_DOY$MeanAOT_440,10, mean,na.rm=T,partial=T)
Aerosol_DOY$MeanAOT_500_10days<-rollapply(Aerosol_DOY$MeanAOT_500,10, mean,na.rm=T,partial=T)
Aerosol_DOY$MeanAOT_675_10days<-rollapply(Aerosol_DOY$MeanAOT_675,10, mean,na.rm=T,partial=T)

Aerosol_DOY<-merge(Aerosol_DOY,ddply(Aerosol_DOY,.(Month),summarise,
                                     MeanAOT_440_month=mean(MeanAOT_440,na.rm=T),
                                     MeanAOT_500_month=mean(MeanAOT_500,na.rm=T),
                                     MeanAOT_675_month=mean(MeanAOT_675,na.rm=T)),"Month")
Aerosol_DOY$Dayof2019<-as.Date(as.character(Aerosol_DOY$DOY),"%j")

## - Seasonality plot (Figure 2)
AOT_seasonal_plot<-ggplot(Aerosol_DOY)+
  geom_line(aes(x=Dayof2019,y=MeanAOT_500),alpha=0.7,colour="gray",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=MeanAOT_500_10days),alpha=0.7,colour="black",linewidth=0.4)+
  geom_line(aes(x=Dayof2019,y=MeanAOT_500_month),alpha=1,linewidth=0.7)+
  scale_x_date(breaks=date_breaks("1 month"),labels = function(x) substring(format(as.Date(x),"%b"), 1, 1) )+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
  theme_classic()+
  ylab("Aerosol Optical Depth (500nm)")+
  xlab("Month")+
  theme(legend.title=element_blank(),legend.position="none")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
  #ggtitle("E.")+
  guides(linetype=guide_legend(nrow=2),alpha=guide_legend(nrow=2),colour=guide_legend(nrow=2))

AOT_seasonal_plot

# Supplementary Figure S2
#ggplot(Aerosol_DOY)+
#  geom_line(aes(x=Dayof2019,y=MeanAOT_500_month,linetype="AOT_500"),alpha=1,linewidth=0.7)+
#  geom_line(aes(x=Dayof2019,y=MeanAOT_440_month,linetype="AOT_440"),alpha=1,linewidth=0.7)+
#  geom_line(aes(x=Dayof2019,y=MeanAOT_675_month,linetype="AOT_675"),alpha=1,linewidth=0.7)+
#  scale_x_date(labels=date_format("%b"),breaks=date_breaks("1 month"))+
#  geom_vline(aes(xintercept=as.numeric(as.Date("2019-03-01"))),linetype="dotted")+
#  geom_vline(aes(xintercept=as.numeric(as.Date("2019-06-01"))),linetype="dotted")+
#  geom_vline(aes(xintercept=as.numeric(as.Date("2019-10-01"))),linetype="dotted")+
#  geom_vline(aes(xintercept=as.numeric(as.Date("2019-12-01"))),linetype="dotted")+
#  theme_classic()+
#  ylab("Aerosol Optical Depth")+
#  theme(legend.title=element_blank(),legend.position="right")+
#  xlab("Month")+
#  theme(plot.margin=unit(c(0.5,0.5,0.5,1),"cm"))+
#  guides(linetype=guide_legend(nrow=3),alpha=guide_legend(nrow=3),colour=guide_legend(nrow=2))

## -  Arrange all seasonal plots together  - Figure 2
#ggarrange(Rainfall_seasonal_plot,Wind_seasonal_plot,
#          MaxMin_seasonal_plot_savanna,MaxMin_seasonal_plot_forest,
#          AH_seasonal_plot_savanna,AH_seasonal_plot_forest,
#          Solar_seasonal_plot,AOT_seasonal_plot,
#          ncol=2,nrow=4,align="hv",
#          labels=c("A","B","C","D","E","F","G","H"))

