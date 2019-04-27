## Script to clean up phi46.xlsx and save as parsable File
##  followed by stats analysis   2018-09-30
##
#--Environment----
#library(magrittr)    #dplyr automatically loads magrittr
#library(plyr)
#library(reshape)

library(dplyr)
library(data.table)
library(tidyr)
library(readxl) #part of tidyverse package, there is no write function yet. Thus us openxlsx to write to excel.
library(here)
#setwd("H:/GoogleDrive/VIMdb/Projects/PHIbase/")
#--MAIN---------------------------------------------------------
#phi=read.csv("phi46.csv", stringsAsFactors=FALSE, strip.white=TRUE, header=FALSE)
phi=read_excel(here("Phi46/phi46_0.xlsx"), sheet= "PHIbase_all")

#  ===TO TO ==== Beneath column removal not implemented in pipeline was doe manually in phi46  
#
#now deleted columns to allow PHI4SQL parsing later
#best now to change the header to MC header, then beneath the cleaning of columns can use the colheader names
#  the older headers need to be copied back to spreadsheet before parsing unless and empty row will suffice

#phi1=phi  %>% as.tbl #read dataset and trim empty Excel columns and header off
#remove additonal columns from Alayne not displayed in PHI4
#excludeColsPhi44=c(28,49,52,53,83,86:91)  #cols in spreadsheet but not in MC parser
#phi2=phi1 %>% select(-excludeColsPhi44)
#phi2=phi2 %>% select(1:64)

#write backup copy
#write.csv(phi2, file=paste("phi45ToParse",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")


# remover all leading+trailing+duplicated white space between words
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
replaceWhite <- function (x) gsub("\\s+", " ", x)
phi3=phi %>% mutate_all(funs(trim)) %>%mutate_all(funs(replaceWhite))
#%>% mutate_each(funs(tolower))
capitalFirstWord <- function(x) gsub("(?<=\\b)(^[a-z])", "\\U\\1", tolower(x), perl=TRUE)
colsForCapitals=c(20,23,26,52) #20=Pathogen_species, 23=Disease_name, 26=Experimental_Host_species, 52=Experimental_evidence
colsForLowerCase=c(35)    #35=mutant_phenotype

phi4=phi3 %>% mutate_at(.vars=vars(colsForCapitals), .funs=funs(capitalFirstWord))
phi5=phi4 %>% mutate_at(.vars=vars(colsForLowerCase), .funs=funs(tolower))
#delete col 1 and col 2, as records may be too long    #CURRENTLY NEED TO BE REINTRODUCED LATER FOR PARSING
#phi6=phi5 %>% select(3:64)
#
# ==To Do== save as Excel
#
write.csv(phi5, file=paste("phi46ToParseNow",Sys.Date(),".csv",sep=""), row.names=FALSE, na="") #save as .txt and import manually into excel to avoid ExcelDataMessing

#now add manually MC header file

#assign cleaned dataset
phi=phi5




str(phi)
#rename columns: phibase to R_phibase(more sensible for work in R)
#get new col names from Mapping table, take only col 4 and restrict to 86 cols
#newNames=read.table("colMap_PHI_R_MC_vers03.txt",header=T,sep="\t", stringsAsFactors=F)
newNames=fread("colMap_PHI_R_MC_vers05.csv",select=c(4),colClasses = list(character=4))
newNames=newNames[2:87]  #exclude header
newNames=unlist(newNames)
#phi=phi[-1,]

setnames(phi,colnames(phi),newNames)
#remove MC header row

phiSlim=phi %>%
  select(- contains("exclude"))

#write backup copy
write.csv(phiSlim, file=paste("phiSlim_All_",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
#phiSlim=read.csv("phiSlim.csv", stringsAsFactors=FALSE, strip.white=TRUE)


#remove plyr/reshape name space to allow dplyr to work with "summarise" function
#detach("package:reshape", unload=TRUE)
# tabulate tally: species x genes x interactions
data=phiSlim %>%
  group_by(PathSpecies) %>%
  summarise(genes=n_distinct(PhiAcc),
            PHI.base.accession=n())

data=arrange(data,desc(PHI.base.accession))
names(data)[names(data)=='PHI.base.accession']<-"interactions"
write.csv(data,file="speciesBYpathogen_geneBYinteraction_tally_PhiALLphi46.csv", row.names=FALSE)
# tabulate tally: genus x genes x interactions
#here split data "species name"
data2= data %>% separate(PathSpecies, c("Genus" ,"Species"), sep =" ", remove=FALSE, convert=FALSE)
data3=data2 %>%
  group_by(Genus ) %>%
  summarise(Species=n_distinct(Species))
data3=arrange(data3,desc(Species))

write.csv(data3, file=paste("genus_no_Species",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")


#===Global Statistics====
global_data=phiSlim %>%
  summarise(n_interactions=n(),
            n_pathogens=n_distinct(PathSpecies),
            n_phi_acc=n_distinct(PhiAcc),
            n_proteins_genes=n_distinct(ProtId),
            n_host=n_distinct(HostSpeciesTaxId), #contains typo
            n_LiteratureIds=n_distinct(Literature_Id),
            n_diseaseNames=n_distinct(DiseaseName))
write.csv(global_data, file=paste("phi46_ALL_globalStatsData",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
# remover all leading+trailing+duplicated white space between words; convert to LOWERCASE
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
replaceWhite <- function (x) gsub("\\s+", " ", x)
phiSlim=phiSlim %>% mutate_all(funs(trim)) %>% mutate_all(funs(replaceWhite))%>% mutate_all(funs(tolower))
write.csv(phiSlim, file=paste("phiSlim_bak_",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
phiSlim=read.csv("phiSlim.csv", stringsAsFactors=FALSE, strip.white=TRUE)

#
# To do ==> move path_species_all_list.csv to GITHUB  (PS: can be transformed to EXCEL!)
#
phiSpecies=read.csv("pathogen_species_all_list.csv", stringsAsFactors=FALSE, strip.white=TRUE)


####FOR Tab2   crosstabulation pathogen group --> PHI-Phenotypes
#count NAs in left_joined table for testing of succesful vlookup with lef_join in dplyr
#  Class is here the pathogen type!


##  Error produced by statement beneath, needs fixing: Error: becuase of integer/character incompatible data types
phiSlim %>%
  select(PhiAcc, PathSpecies,PathSpeciesTaxId,Phenotype)  %>%
  left_join(select(phiSpecies,PathSpeciesTaxId, Class), by=c("PathSpeciesTaxId"="PathSpeciesTaxId")) %>%
  summarise_each(funs(sum(is.na(.))))    #result should be "0" for all columns!

####################################BENEATH NEEDS CODE & FUNCTION REVIEW##############################################
# take all 'increased virulence (hypervirulence)' interations and subdivide the list into plant/animal/fungal infecting
# generate Excel spreadsheet

# IncrVirRecords=phiSlim %>%
#   filter(Phenotype=="increased virulence (hypervirulence)") 
# 
# IncrVirRecords %>% count #gives "514 records"
# 
# phiIncrVirByAniPlant=IncrVirRecords %>%
#   left_join(select(phiSpecies,PathSpeciesTaxId, ExperimentalHost), by=c("PathSpeciesTaxId"="PathSpeciesTaxId")) %>% 
#   group_by(ExperimentalHost) %>%
#   summarise(n_pathClassNum=n()) %>% arrange(desc(n_pathClassNum))
# 
# IncrVirRecords %>%
#   anti_join(phiSpecies, by=c("PathSpeciesTaxId"="PathSpeciesTaxId"))
# 
# ##now generate Excel sheet
# phiIncrVirByAniPlantEx=phi %>%
#   filter(Phenotype=="increased virulence (hypervirulence)") 
# 
# add_cols=phiSpecies %>%
#   select(PathogenTaxId, host1) %>% as.tbl
# add_cols$PathogenTaxId= as.character(add_cols$PathogenTaxId)
# exHypVir=left_join(phiIncrVirByAniPlantEx,add_cols, by=c("PathSpeciesTaxId"="PathogenTaxId")) 
# write.csv(exHypVir, file=paste("HyperVirByAnimalPlant",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
# 
# 
# #beneath left_join on INTEGER column PathogenTaxId
# #-anti-join pre-test to identify missing species in path_species_all_list
# phiSlim %>%
#   select(PathSpeciesTaxId) %>%
#   anti_join(phiSpecies, by=c("PathSpeciesTaxId"="PathogenTaxId"))
# 
# 
# tFus= phiSlim %>%
#   filter(PathSpecies=="fusarium graminearum") %>%
#   select(PhiAcc, PathSpecies, PathSpeciesTaxId,Phenotype)  %>%
#   group_by(Phenotype,PathSpecies) %>%
#   summarise(n_interactions=n())
# 
# # the Fusarium species divided out by phenotype 
# library(stringr)
# t2=phiSlim %>%
#   select(PhiAcc, PathSpecies, PathSpeciesTaxId,Phenotype) %>%
#   filter(str_detect(PathSpecies, regex("Fusarium",ignore_case=TRUE))) 
# 
# sFus=t2 %>% 
#   group_by(Phenotype,PathSpecies) %>%
#   summarise(n_interactions=n())
# write.csv(sFus, file=paste("FusSpecByPhenotyp",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
# 
# # the Cochliobulus/Bipolaris species divided out by phenotype 
# library(stringr)
# t3=phiSlim %>%
#   select(PhiAcc, PathSpecies, PathSpeciesTaxId,Phenotype) %>%
#   filter(str_detect(PathSpecies, regex("bipolaris|cochl",ignore_case=TRUE))) 
# 
# sCoch=t3 %>% 
#   group_by(Phenotype,PathSpecies) %>%
#   summarise(n_interactions=n())
# write.csv(sCoch, file=paste("CochBipSpecByPhenotyp",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
# 
# # the wheat entries species list ranked high to low divided out by phenotype. There are 1299.
# t4=phiSlim %>%
#   select(HostExpSpecies,PhiAcc, PathSpecies, PathSpeciesTaxId,Phenotype) %>%
#   filter(str_detect(HostExpSpecies, regex("wheat",ignore_case=TRUE)))
# 
# sWheat=t4 %>% 
#   group_by(PathSpecies,Phenotype) %>%
#   summarise(n_interactions=n())
# speciesRank=sWheat %>%
#   group_by(PathSpecies) %>%
#   summarise(TotInterActions=sum(n_interactions))
# ##now make lookup
# rankedWheatPath=sWheat %>%
#   left_join(speciesRank, by="PathSpecies") %>%
#   arrange(desc(TotInterActions), desc(n_interactions))
# write.csv(rankedWheatPath, file=paste("WheatPathRankedByInteratByPhenotype",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
# 
# 
# 
# 
# 
# sWheat=t4 %>% 
#   group_by(Phenotype,PathSpecies) %>%
#   summarise(n_interactions=n())
# 
# write.csv(sCoch, file=paste("CochBipSpecByPhenotyp",Sys.Date(),".csv",sep=""), row.names=FALSE, na="")
# 
# data=arrange(data,desc(PHI.base.accession))
# 
# 
# 
# 
# 
# #calculate Pathogen class numbers in phiSlim
# phiClassNum=phiSlim %>%
#   group_by(PathSpeciesTaxId) %>%
#   select(PathSpecies, PathSpeciesTaxId)  %>%
#   left_join(select(phiSpecies,PathogenTaxId, Class), by=c("PathSpeciesTaxId"="PathogenTaxId")) %>%
#   group_by(Class) %>%
#   summarise(n_pathClassNum=n_distinct(PathSpeciesTaxId)) %>% arrange(desc(n_pathClassNum))
# 
# #convert phenotype for bact,fung,nematode (PathogenClass) from narrow to wide data format
# #make sure the Phenotype  column is a factor
# t$Phenotype=factor(t$Phenotype)
# library(tidyr)
# t_wide=spread(t, Class, n_interactions)
# write.csv(t_wide,file="Tab2_phenotypByPathClass.csv", row.names=FALSE, na="0")
# #====
# ##Left join testing
# users <- data.frame(
#   user_id = c(1,2,4),
#   age = c(20,20,30)
# ) %>% as.tbl
# 
# items <- data.frame(
#   item_id = 1:3,
#   item = c("1+1","2*2","3/3")
# ) %>% as.tbl
# 
# log %>% left_join(users, "user_id") %>% left_join(items, "item_id")
# 
# ####FOR Tab1   crosstabulation: HOST --> PHI-Phenotypes
# #count NAs in left_joined table for testing of succesful vlookup with lef_join in dplyr
# #  Class is here the pathogen type!
# phiHosts=read.table("host_species_all.txt",header=T,sep="\t", stringsAsFactors=F) %>% as.tbl
# phiHosts=phiHosts %>% mutate_each(funs(trim)) %>%mutate_each(funs(replaceWhite))%>% mutate_each(funs(tolower))
# phiHosts$HostSpeciesTaxId=as.character(phiHosts$HostSpeciesTaxId)
# 
# #colnames(phiSlim)[15]="HostSpeciesTaxId"   #finally fix typo
# 
# 
# #calculate numbers
# phiHostNum=phiSlim %>%
#   group_by(HostSpeciesTaxId) %>%
#   select(HostSpeciesTaxId, Phenotype)  %>%
#   left_join(select(phiHosts,HostSpeciesTaxId, HostGroup), by=c("HostSpeciesTaxId"="HostSpeciesTaxId")) %>%
#   group_by(HostGroup) %>% na.omit() %>%
#   summarise(n_HostNum=n_distinct(HostSpeciesTaxId)) %>% arrange(desc(n_HostNum))
# 
# z= phiSlim %>%
#   select(PhiAcc,HostSpeciesTaxId,Phenotype)  %>%
#   left_join(select(phiHosts,HostSpeciesTaxId, HostGroup), by=c("HostSpeciesTaxId"="HostSpeciesTaxId")) %>%
#   group_by(Phenotype,HostGroup) %>%
#   summarise(n_interactions=n())
# z=na.omit(z)
# 
# #make sure the Phenotype  column is a factor
# z$Phenotype=factor(z$Phenotype)
# 
# library(tidyr)
# z_wide=spread(z,HostGroup, n_interactions)
# #t_wide=spread(t, Class, n_interactions)
# 
# write.csv(z_wide,file="Tab1_phenotypByHost.csv", row.names=FALSE, na="0")
# #====
# 
# 
# 
# #--hosts by pathogen type crosstabulation---------------------------------------------------------
# #phiSpecies=read.csv("PHI42_species_list.csv", stringsAsFactors=FALSE,header=TRUE,fill=TRUE)
# #write backup copy
# data2=phiSpecies %>%
#   group_by(HostExpSpecies_latin,pathogen.type) %>%
#   summarise(path_species=n_distinct(Pathogen.species))
# 
# 
# result=phiSpecies%>%
#   group_by(pathogen.type) %>%
#   arrange(desc(Pathogen.species))
# 
# write.csv(data2,file="hosts_by_organism_type.csv", row.names=FALSE)
# 
# #======================
# #==for Tabl 5====NAR===
# #=count interactions and unique genes per hosts
# #host_intact_genes=phiSlim%>%
# #  group_by(HostSeciesTaxId) %>%
# #  summarise(PHI.base.accessions=n(),genes=n_distinct(PhiAcc))
# 
# #write.csv(host_intact_genes,file="Tab5_hostsByInteractions_genes.csv", row.names=FALSE)
# 
# #==for Tab5 divide hosts by bacteria, fungi, oomycete, others===
# #=count interactions per hosts and dived by into path.typ (bacteria, fungi, oomycete)
# host.intact.ptype=phiSlim%>%
#   group_by(HostSpeciesTaxId, HostExpSpecies) %>%     
#   summarise(PhiAccessions=n())
# #  summarise(PHI.base.accessions=n(),genes=n_distinct(PHI.base.accession.no))
# 
# write.csv(host.intact.ptype,file="Tab5_hostsByIntactPatType.csv", row.names=FALSE)
# host.intact.ptype=read.csv("Tab5_hostsByIntactPatType.csv")
# 
# #example: colnames(select(phiSlim, matches("prot")))
# 
# 
# #=======beneath to next section may not be useful==#
# #convert host.intact.ptypet from narrow to wide data format
# #make sure the host  column is a factor
# host.intact.ptype$host=factor(host.intact.ptype$hosts)
# library(tidyr)
# host.intact.ptype_wide=spread(host.intact.ptype, path.type, PHI.base.accessions)
# write.csv(host.intact.ptype_wide,file="Tab5_hostsByIntactPatType_wide.csv", row.names=FALSE, na="0")
# #==================================================================
# ##tally 15 pathogens with most interactions
# freq_path_tally=phiSlim %>%
#   group_by(PathSpecies) %>% tally %>%
#   arrange(desc(n)) %>%head(n=20)
# write.csv(freq_path_tally,file="freq_path_tally20.csv", row.names=FALSE, na="0")
# 
# 

