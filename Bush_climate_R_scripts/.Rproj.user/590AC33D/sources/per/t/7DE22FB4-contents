rm(list = ls())

make_plots <- FALSE



source(file = "Bush2020paper_part1.R")

# List of data frames you want to export
df_list <- list(
  Rainfall_DOY = Rainfall_DOY,
  Temperature_DOY = Temperature_DOY,
  Humidity_DOY = Humidity_DOY,
  Solar_DOY = Solar_DOY,
  Wind_DOY = Wind_DOY,
  Aerosol_DOY = Aerosol_DOY
)

# Make sure the folder exists
if(!dir.exists("processed_data")) dir.create("processed_data")

# Export each data frame to CSV in that folder
for (name in names(df_list)) {
  write.csv(df_list[[name]], file = file.path("processed_data", paste0(name, ".csv")), row.names = FALSE)
}








# ---------------------------------------------------------------------------------

if(make_plots){
  
    # rainfall
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
    
    
    # temperature
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
    
    
    
    # absolute humidity
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
    
    
    
    # solar
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
    
    
    # wind
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
    
    # aerosol
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
}

rm(list = ls())