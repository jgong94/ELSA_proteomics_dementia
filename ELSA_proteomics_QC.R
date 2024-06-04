rm(list = ls())

###LOADING LIBRARIES
#install.packages("OlinkAnalyze")
library(OlinkAnalyze)
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
#install.packages("readxl")
library(readxl)
#install.packages("ggplotify")
library(ggplotify)
#install.packages("pheatmap")
library(pheatmap)
#Downloading libraries
library(foreign)
#install.packages("haven")
library(haven)
#install.packages("stringi")
library(stringi)
library(magrittr)
library(readxl)
library(tidyr)

###LOADING DATASETS
setwd("~path")

###REMOVE CONTROL, QC_FAIL, OUTLIERS
###CVDII
data_CVDII <- read_NPX("Q-07545_Steptoe_NPX_CVDII.xlsx")
data_CVDII_remove_control <- data_CVDII[grepl("^BR", data_CVDII$SampleID), ] 
data_CVDII_pass_QC <- subset(data_CVDII_remove_control,QC_Warning=="Pass")
data_CVDII_clean=data_CVDII_pass_QC
variables_to_remove <- c("Index", "OlinkID", "UniProt", "MissingFreq", "Panel", "Panel_Version", "PlateID", "QC_Warning", "LOD", "Normalization")

data_CVDII_clean_1 <- data_CVDII_clean %>%
  dplyr::select(-variables_to_remove)

data_CVDII_clean_1$NPX<-as.numeric(data_CVDII_clean_1$NPX)

data_CVDII_clean_2 <- data_CVDII_clean_1 %>%
  filter(!(SampleID %in% c("BR973440G"))) #Remove outliers

length(unique(data_CVDII_clean_2$SampleID))

data_CVDII_wide <- data_CVDII_clean_2 %>%
  pivot_wider(names_from = Assay, values_from = NPX)

data_CVDII_wide <- data_CVDII_wide %>%
  mutate(across(where(is.double), as.character))

names(data_CVDII_wide)[names(data_CVDII_wide) == "SampleID"] <- "ELSA_ID"

#write.csv(data_CVDII_wide, "Q-07545_Steptoe_NPX_CVDII_wide.csv")

###NEUI
data_NEUI <- read_NPX("Q-07545_Steptoe_NPX_NEUI.xlsx")
data_NEUI_remove_control <- data_NEUI[grepl("^BR", data_NEUI$SampleID), ] 
data_NEUI_pass_QC <- subset(data_NEUI_remove_control,QC_Warning=="Pass")
data_NEUI_clean=data_NEUI_pass_QC

variables_to_remove <- c("Index", "OlinkID", "UniProt", "MissingFreq", "Panel", "Panel_Version", "PlateID", "QC_Warning", "LOD", "Normalization")

data_NEUI_clean_1 <- data_NEUI_clean %>%
  dplyr::select(-variables_to_remove)

data_NEUI_clean_1$NPX<-as.numeric(data_NEUI_clean_1$NPX)

data_NEUI_clean_2 <- data_NEUI_clean_1 %>%
  filter(!(SampleID %in% c("BR973440G", "BR984765M"))) #Rmove outliers

length(unique(data_NEUI_clean_2$SampleID))

data_NEUI_wide <- data_NEUI_clean_2 %>%
  pivot_wider(names_from = Assay, values_from = NPX)

data_NEUI_wide <- data_NEUI_wide %>%
  mutate(across(where(is.double), as.character))

names(data_NEUI_wide)[names(data_NEUI_wide) == "SampleID"] <- "ELSA_ID"

#write.csv(data_NEUI_wide, "Q-07545_Steptoe_NPX_NEUI_wide.csv")

###NEX
data_NEX_I <- read_NPX("Q-08033_Steptoe_NPX_NEX.xlsx")
data_NEX_II <- read_NPX("Q-08440_Steptoe_NPX_NEX.xlsx")
data_NEX <- read_NPX("Q_08033_Q_08440_Steptoe_NPX_NEX.csv")
data_NEX <- data_NEX %>% relocate(SampleID)

data_NEX_remove_bridge <- subset(data_NEX, Type=="NEX_DATA_1 Sample" | Type=="NEX_DATA_2 Sample" | Type=="NEX_DATA_1 Bridge") 

data_NEX_remove_control <- data_NEX_remove_bridge[grepl("^BR", data_NEX_remove_bridge$SampleID), ]

data_NEX_pass_QC <- subset(data_NEX_remove_control,QC_Warning=="Pass")

data_NEX_clean=data_NEX_pass_QC

data_NEX_below_LOD<-subset(data_NEX_clean,MissingFreq>0.5)

list(unique(data_NEX_below_LOD$Assay))

table(data_NEX_below_LOD$Assay, data_NEX_below_LOD$MissingFreq)

variables_to_remove <- c("Index", "OlinkID", "UniProt", "MissingFreq", "Panel", "Panel_Version", "PlateID", "QC_Warning", "LOD", "Normalization", "Type")

data_NEX_clean_1 <- data_NEX_clean %>%
  dplyr::select(-variables_to_remove)

data_NEX_clean_1$NPX <- as.numeric(data_NEX_clean_1$NPX)

data_NEX_aggregated <- data_NEX_clean_1 %>%
  group_by(SampleID, Assay) %>%
  summarise(NPX = mean(NPX, na.rm = TRUE)) %>%
  ungroup()

data_NEX_wide <- data_NEX_aggregated %>%
  pivot_wider(names_from = Assay, values_from = NPX)

data_NEX_wide <- data_NEX_wide %>%
  mutate(across(where(is.list), as.numeric))

names(data_NEX_wide)[names(data_NEX_wide) == "SampleID"] <- "ELSA_ID"

###CVDII
data_CVDII %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) %>% 
  olink_pca_plot(df = .,
                 color_g = "QC_Warning", byPanel = TRUE, outlierDefX = 5,
                 outlierDefY = 5, label_outliers = TRUE, outlierLines = TRUE)  

data_CVDII %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE'),
         Panel == 'Olink Cardiovascular II') %>% 
  olink_qc_plot(color_g = "QC_Warning", IQR_outlierDef = 5, median_outlierDef = 5, outlierLines = TRUE, label_outliers = TRUE)   

qc <- olink_qc_plot(data_CVDII, color_g = "QC_Warning", IQR_outlierDef = 5, median_outlierDef = 5)
qc$data %>% filter(Outlier == 1) %>% select(SampleID, Panel, IQR, sample_median, Outlier)

qc

###NEU I
data_NEUI %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) %>% 
  olink_pca_plot(df = .,
                 color_g = "QC_Warning", byPanel = TRUE, outlierDefX = 5,
                 outlierDefY = 5, label_outliers = TRUE, outlierLines = TRUE)  

data_NEUI %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE'),
         Panel == 'Olink Neurology') %>% 
  olink_qc_plot(color_g = "QC_Warning", IQR_outlierDef = 5, median_outlierDef = 5, label_outliers = TRUE, outlierLines = TRUE)   

qc <- olink_qc_plot(data_NEUI, color_g = "QC_Warning", IQR_outlierDef = 5, median_outlierDef = 5)
qc$data %>% filter(Outlier == 1) %>% select(SampleID, Panel, IQR, sample_median, Outlier)

qc

###NEX
data_NEX=data_NEX_remove_bridge

data_NEX %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE')) %>% 
  olink_pca_plot(df = .,
                 color_g = "QC_Warning", outlierDefX = 5,
                 outlierDefY = 5, label_outliers = TRUE, outlierLines = TRUE)  

data_NEX %>% 
  filter(!str_detect(SampleID, 'CONTROL_SAMPLE'),
         Panel == 'Olink Neuro Exploratory') %>% 
  olink_qc_plot(color_g = "QC_Warning", IQR_outlierDef = 5, median_outlierDef = 5, outlierLines = TRUE, label_outliers = TRUE)   

qc <- olink_qc_plot(data_NEX, color_g = "QC_Warning", IQR_outlierDef = 5, median_outlierDef = 5)
qc$data %>% filter(Outlier == 1) %>% select(SampleID, Panel, IQR, sample_median, Outlier)

qc

#write.csv(data_NEX_wide, "Q_08033_Q_08440_Steptoe_NPX_NEX_wide.csv", quote = TRUE, row.names = FALSE)

###Merge three panels
data_olink_all_1<-merge(data_NEUI_wide,data_CVDII_wide,by="ELSA_ID",all=T)
data_olink_all<-merge(data_olink_all_1,data_NEX_wide,by="ELSA_ID",all=T)
