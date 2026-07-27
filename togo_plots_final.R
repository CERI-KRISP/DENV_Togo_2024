library("readxl")
library(ggplot2)
library(dplyr)
library(geosphere)
library(ggnewscale)
library(rnaturalearth)
library(stringr)
library(sf)
library(geodata)
library(raster)
library(RColorBrewer)
library(tidyr)
library(lubridate)

library(ggtree)
library(tidyverse)
library(tidytree)
library(ape)
library(treeio)
library(cowplot)

#### Figure 1

togo_sequence_metadata<-read_excel('data/samples_metadata.xlsx')
togo_sequence_metadata$date<-as.Date(cut(as.Date(togo_sequence_metadata$`Sample Collection Date (DD-MM-YYYY)`),
                                         breaks = "day",
                                         start.on.monday = FALSE))

#Supplementary Figure S1
Ct_coverage_plot<-ggplot()+
  theme_bw()+
  geom_point(data=togo_sequence_metadata, aes(Ct,coverage),fill='red2',colour='black',shape=21,alpha=0.5,size=2)+
  geom_smooth(data=togo_sequence_metadata, aes(Ct,coverage), colour='red4', alpha=0.2)+
  xlab("Ct Score")+ylab("Genome Coverage (%)")
Ct_coverage_plot


#Fig 1A
Africa_2024_cases_deaths<-read_excel('data/Dengue_2024_Africa_cases.xlsx')

africa_cdc_map<-read.csv('data/Africa_fortified_Name_0.tsv',sep = "\t")

map_data<- africa_cdc_map %>% 
  left_join(Africa_2024_cases_deaths, by = c("id" = "Country"))


centroids <- map_data %>% 
  group_by(id) %>% 
  group_modify(~ data.frame(centroid(cbind(.x$long, .x$lat))))


centroids<- centroids %>% 
  left_join(Africa_2024_cases_deaths, by = c("id" = "Country"))


Africa_map<-ggplot() +
  #geom_polygon()+
  theme_void()+
  geom_map(data=map_data,map=map_data, aes(long, lat,map_id=id), color="grey40", fill='white',size=0.2)+
  
  geom_map(data=map_data,map=map_data, aes(long, lat,map_id=id, fill=log10(as.numeric(TotalCases_September2024))), color="grey40",size=0.2) +
  scale_fill_distiller(palette = "PuBuGn", direction = 1,na.value = "white",
                       breaks = c(0, 1, 2, 3,4), labels = c(1, 10, 100, 1000,10000),
                       
                       name='2024 Dengue Cases') +   
  #scale_fill_distiller(palette = "PuBuGn", direction = 1,na.value = "white",trans = "log",
  #                    breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000), name='Genomes') + 
  coord_fixed()+
  theme(legend.position = c(0.25,0.3),
        legend.key.width = unit(dev.size()[1] / 20, "inches"))+
geom_point(data = centroids,
          aes(x = lon, y = lat,size=TotalDeaths_September2024),
         shape=21, fill='deeppink3',alpha=0.7)+
  scale_size_continuous(name='Deaths')


ggsave(file="figures/Fig1A_Africa_dengue_cases_map.pdf",plot=Africa_map)



### Togo cases
Togo_cases_all_tests<-read_excel('data/Togo_cases_all_tests.xlsx')
Togo_cases_all_tests$date2<-as.Date(cut(Togo_cases_all_tests$`DATE DE PRELEVEMENT`,
                                        breaks = "week",
                                        start.on.monday = FALSE))
Togo_cases_all_tests$date3<-as.Date(cut(Togo_cases_all_tests$`DATE DE PRELEVEMENT`,
                                        breaks = "month",
                                        start.on.monday = FALSE))
#Fig 1B
togo_cases<-ggplot(Togo_cases_all_tests)+
  theme_bw()+
  geom_bar(aes(date2,fill=Classification),colour='black',size=0.3)+
  scale_fill_manual(values=c('dodgerblue3','purple3'),name='',labels=c('Suspected Cases','Confirmed Cases'))+
  geom_rug(data=togo_sequence_metadata, aes(date), colour='red3',alpha=0.5)+
  #geom_line(data=monthly_togo_TP,aes(date_2024-60,IndexP*100))+
  ylab("New Cases")+
  xlab("Sampling Dates")+
  theme(legend.position=c(0.2,0.8),legend.background = element_blank())
togo_cases

#Fig 1B with Rt estimate
rt <- read_csv("../Rt estimates/output/Rt_estimates_Togo_clean.csv")
max_cases <- max(table(Togo_cases_all_tests$date2), na.rm = TRUE) # Scale Rt
sf <- max_cases / 7.5 # Set extent (7.5) to match grid
togo_cases + # Add Rt overlay
  geom_ribbon(data = rt, aes(x = date, ymin = R_low * sf, ymax = R_high * sf),
              fill = "#E69F00", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = rt, aes(x = date, y = R_median * sf),
            colour = "#E69F00", linewidth = 1, inherit.aes = FALSE) +
  geom_hline(yintercept = 1 * sf, linetype = "dashed", linewidth = 0.6, colour = "#E69F00") +
  scale_y_continuous(name = "New Cases", sec.axis = sec_axis(~ . / sf, name = "Rt", breaks = seq(0, 7.5, by = 2))
  )

ggsave("figures/Fig1B_Togo_with_Rt.pdf", width = 8, height = 4,
       units = "in")

#Fig1C
togo_cases_by_district<-ggplot(subset(subset(subset(Togo_cases_all_tests,!is.na(REGION)),REGION!='MARITIME'),REGION!='CENTRALE'))+
  theme_bw()+
  geom_bar(aes(date2,fill=Classification),colour='black',size=0.3)+
  scale_fill_manual(values=c('dodgerblue3','purple3'),name='',labels=c('Suspected Cases','Confirmed Cases'))+
  geom_rug(data=subset(togo_sequence_metadata,REGION!='NOT SPECIFIED'), aes(date), colour='red3',alpha=0.5)+
  ylab("New Cases")+
  xlab("Sampling Dates")+
  theme(legend.position='none',legend.background = element_blank())+
  facet_wrap(~REGION,ncol=2)
togo_cases_by_district


#Fig 1D
##Togo epi / genomes map
dir.create("data/gadm", recursive = TRUE)
Togo_gadm_data_0 <- geodata::gadm(country = "TGO", level = 0,path = "data/gadm")
Togo_gadm_data_1 <- geodata::gadm(country = "TGO", level = 1,path = "data/gadm")
Togo_gadm_data_1_sf <- sf::st_as_sf(Togo_gadm_data_1)

Togo_gadm_data_2 <- geodata::gadm(country = "TGO", level = 2,path = tempdir())
Togo_gadm_data_2_sf <- sf::st_as_sf(Togo_gadm_data_2)
Togo_gadm_data_2_df<-as.data.frame(Togo_gadm_data_2)


Togo_cases_count_confirmed<- subset(Togo_cases_all_tests,Classification=='Confirmé') %>% group_by(district) %>%
  summarise(count=n())
Togo_cases_count_confirmed

Togo_cases_count_confirmed<-rbind(subset(Togo_cases_count_confirmed,!is.na(district)), data.frame(district = "Lomé", count = 881))


Togo_gadm_data_2_df<-left_join(Togo_gadm_data_2_df,Togo_cases_count_confirmed, by=c("NAME_2"="district"))
Togo_gadm_data_2_sf<-left_join(Togo_gadm_data_2_sf,Togo_cases_count_confirmed, by=c("NAME_2"="district"))


togo_sequence_count_district<-togo_sequence_metadata %>% group_by(District) %>%
  summarise(count=n())
togo_sequence_count_district

togo_sequence_count_district$District[togo_sequence_count_district$District=='Agoè']<-"Agoe-Nyive"
togo_sequence_count_district$District[togo_sequence_count_district$District=='Binah']<-"Bimah"
togo_sequence_count_district$District[togo_sequence_count_district$District=='Cinkasse']<-"Cinkassé"
togo_sequence_count_district$District[togo_sequence_count_district$District=='Golfe']<-"Lomé"
togo_sequence_count_district$District[togo_sequence_count_district$District=='Tone']<-"Tône"
#togo_sequence_count_district<-subset(Togo_cases_count_confirmed,!is.na(District))


Togo_gadm_data_2_sf2<-left_join(Togo_gadm_data_2_sf,togo_sequence_count_district, by=c("NAME_2"="District"))
Togo_gadm_data_2_sf2_centroids <- st_centroid(Togo_gadm_data_2_sf2)

Togo_epi_map<-ggplot() + 
  theme_void()+
  geom_sf(data=Togo_gadm_data_2_sf,aes(fill=log10(count))) +
  geom_sf(data=Togo_gadm_data_1_sf,colour='black',linewidth=0.5, fill=NA) +
  geom_point(data=Togo_gadm_data_2_sf2_centroids, aes(geometry=geometry,size=log2(count.x)),
             stat = "sf_coordinates",shape=21,colour='red3')+
  geom_text(data=Togo_gadm_data_2_sf2_centroids, aes(geometry=geometry,label=count.x),size=2,
             stat = "sf_coordinates",shape=21,colour='white')+
  scale_size_continuous(range = c(1, 10),breaks=c(1,5),labels=c(1,25),name='Genomic\nSampling')+  # Adjust min/max point sizes
scale_fill_distiller(palette = "PuBuGn", direction = 1,na.value = "white",
                       breaks = c(0, 1, 2), labels = c(1, 10, 100),
                       
                       name='Cases') +
  theme(legend.position = c(0.0,0.5))

Togo_epi_map



#Fig 1E-D

Togo_map1 <- ne_countries(type = "countries", country = "Togo",
                          scale = "medium", returnclass = "sf")
Togo_map2 <- rnaturalearth::ne_states(country = "Togo",
                                      returnclass = "sf")


#2023 transmission potential
TP_2023<-raster('data/transmisisonpotential_togo/Togo_TP_2023.tiff')
cols<-rev(brewer.pal(9,"RdBu"))
par()

plot(TP_2023,col=cols,axes=F,box=F,legend.args=list(text="Transmission Potential", cex=0.7, line=0.3, col="gray30"))
plot(Togo_map2$geometry,add=T,col='transparent')

TP_2023_monthly<-stack('data/transmisisonpotential_togo/Togo_monthly_TP_2023.tiff')

n <- nlayers(TP_2023_monthly)
ncol <- ceiling(sqrt(n))
nrow <- ceiling(n / ncol)

par(mfrow = c(nrow, ncol))
cols<-rev(brewer.pal(9,"RdBu"))
plot(TP_2023_monthly,col=cols,axes=F,box=F,smallplot=c(0.75,0.80,0.3,0.7),legend.args=list(text="Transmission\nPotential", cex=0.7, line=0.3, col="gray30"))
#plot(Togo_map2$geometry,add=T,col='transparent')

####Figure2

### Visualizing trees

#Fig2A
timetree <- read.newick("data/timetree.nwk")
timetree <- groupClade(timetree,.node=c(1359,1361))


p <- ggtree(timetree, mrsd="2024-09-27", as.Date=TRUE,aes(colour=group),size=0.3) + theme_tree2() +
  #geom_tiplab(size=0.5)+
  scale_colour_manual(values=c('cornsilk3','dodgerblue2','darkslateblue'))+
  scale_fill_manual(values=c('cornsilk3','dodgerblue2','darkslateblue'))+
  
  geom_tippoint(aes(fill=group),shape=21,stroke=0.04,size=2,colour='black')+
  #geom_tippoint(aes(
  #  subset=(grepl("Togo",label,fixed=TRUE)==TRUE)),size=2, align=F, fill='tomato1', colour='black', shape=21)+
  #geom_tippoint(aes(
  #  subset=(grepl('Tanzania',label,fixed=TRUE)==TRUE)),size=2, align=F, fill='dodgerblue3', colour='black', shape=21)+
  geom_tippoint(
    data = function(x) dplyr::filter(x, grepl("Togo", label, fixed = TRUE)),
    size = 2, align = FALSE, fill = 'red3', colour = 'black', shape = 21,stroke=0.05
  ) +
  geom_tippoint(
    data = function(x) dplyr::filter(x, grepl("Tanzania", label, fixed = TRUE)),
    size = 2, align = FALSE, fill = 'cornflowerblue', colour = 'black', shape = 21,stroke=0.05
  ) +
  theme(legend.position='none',
        axis.text = element_text(angle=90))+
  scale_x_date(date_labels = "%Y",date_breaks = "4 years")
  
  #geom_text(aes(label=node), hjust=1.5, size=2,colour='black')


p

#Fig2B


tree <- read.beast("data/beast_files/lineage_A.2_n483.relaxed.skygrid.500ml.Combined_MCC.tree")

p1 <- ggtree(tree, mrsd="2024-09-27", as.Date=TRUE,color='darkslateblue',size=0.5) + theme_tree2() +
  scale_color_manual(values=c('deeppink2','darkorange3','gold2','dodgerblue3'), labels=c("Omicron\nAncestor","BA.5","BA.4","BA.2"), name='')+
  #geom_tiplab(size=0.5)+
  theme(legend.position = "top")+
  theme(legend.direction = "horizontal")+
  geom_tippoint(
    data = function(x) dplyr::filter(x, grepl("Togo", label, fixed = TRUE)),
    size = 2, align = FALSE, fill = 'red3', colour = 'black', shape = 21,stroke=0.05
  ) +
  geom_tippoint(
    data = function(x) dplyr::filter(x, grepl("Tanzania", label, fixed = TRUE)),
    size = 2, align = FALSE, fill = 'cornflowerblue', colour = 'black', shape = 21,stroke=0.05
  ) +
  
  geom_tippoint(
    data = function(x) dplyr::filter(x, grepl("Cote_dIvoire", label, fixed = TRUE)),
    size = 2, align = FALSE, fill = 'hotpink3', colour = 'black', shape = 21,stroke=0.05
  ) +
  
  geom_tippoint(
    data = function(x) dplyr::filter(x, grepl("Benin", label, fixed = TRUE)),
    size = 2, align = FALSE, fill = 'darkgreen', colour = 'black', shape = 21,stroke=0.05
  ) +
  geom_tippoint(
    data = function(x) dplyr::filter(x, grepl("Burkina_Faso", label, fixed = TRUE)),
    size = 2, align = FALSE, fill = 'lightgreen', colour = 'black', shape = 21,stroke=0.05
  ) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    axis.text.x = element_text(angle=90)
    #legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    #legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+

  scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
  expand_limits(y=500)
#geom_text(aes(label=node), hjust=-.3, size=2)


p1


#Fig2C
pop_size<-read.table('data/beast_files/lineage_A.2_n483.relaxed.skygrid.500ml.Combined_pop_size.txt', header=T)
  
pop_size_plot<-ggplot(pop_size)+
  theme_bw()+
  geom_ribbon(aes(x=as.Date(date_decimal(time)),ymin=log10(lower),ymax=log10(upper)), fill='dodgerblue2',alpha=0.2, color=NA)+
  geom_line(aes(x=as.Date(date_decimal(time)),y=log10(median)),colour='dodgerblue2')+
  xlab("Date")+ylab("Effective Population Size (log)")

pop_size_plot



#Fig2D in other R code "Togo_Dengue_Seraphim.R"


# Supplementary Fig S2- lineage A.2 Tempest

tempest_df<-read_excel('data/lineage_A.2_IDs_n483_sequences_tempest.xlsx')
tempest_df$date2<-as.Date(date_decimal(as.numeric(tempest_df$date)))

ggplot(tempest_df,aes(date2,as.numeric(distance)))+
  theme_classic()+
  scale_fill_manual(values=c('darkgreen','lightgreen','grey','hotpink3','cornflowerblue','red3'), name='Country')+
  scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
  xlab('Sampling Date')+
  ylab("Average per site genetic\ndivergence from root")+
  geom_point(shape=21,size=3,aes(fill=country))+
  geom_smooth(method='lm', color='black',size=0.6)+
  annotate(geom = 'text',label='Correlation coefficient = 0.82', x=as.Date("2021/01/30"),y=0.022)+
  annotate(geom = 'text',label='R squared = 0.67', x=as.Date("2021/01/30"),y=0.02)









