#!/usr/bin/env Rscript


#base <- "~/Bureau/Prototype data analysis/Tools/Scripts/GINsim_heatmaps/"
#setwd(base)

## Get the required packages
ipkgs <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#Load existing packages, install missing packages.
required.packages <- c('here','data.table','ggplot2','tidyr','reshape','parallel','Matrix', 'gridExtra','reshape2','RColorBrewer','ggpubr','cowplot','grid','argparse','argparser')
ipkgs(required.packages)

setwd(here())
setwd("./GINsim_heatmaps")
print(getwd())





p <- arg_parser("A script making heatmaps based on GINsim results")


# specify our desired options

p <- add_argument(p, "-n", "--modelNodes", help = "Model nodes text file", flag = F)
p <- add_argument(p, "-i", "--inputs", help = "Input nodes text file (needs 2 inputs at least & at most)", flag = F)
p <- add_argument(p,"-b", "--booleanNodes", help = "Boolean nodes text file", flag = F)
p <- add_argument(p, "-M", "--multivaluedNodes", help = "Multivalued nodes text file", flag = F)
p <- add_argument(p, "-r", "--readouts", help = "Model readout nodes text file (the nodes we are interested in", flag = F)
p <- add_argument(p, "-o", "--outputs", help = "Model output nodes text file (the nodes we want to display in the pdf file)", flag = F)
p <- add_argument(p, "-ISo", "--IS_outputs", help = "Model output nodes, specific to Isc and Suf text file (displayed in the pdf file)", flag = F)


required_arguments <- c(4:10)

#Getting all the inputfiles in a parser
arguments <- commandArgs(FALSE)
true_args <- arguments[c(7:length(arguments))]
if (length(true_args)-6 != 7)
  {
    print(p)
    print("Usage: \n Rscript IsolatorMapper.R -n Model.txt -i Inputs.txt -b Boolean.txt -M Multivalued.txt -r Readouts.txt -o TheOutputs.txt -ISo IS_outputs.txt")
    quit(save = "no")
  } else 
  
  {
    argz <- parse_args(p)
    True_arguments <- argz[required_arguments]
    source("IsolatorHelper.R")
  }


ptm <- proc.time()










#GLOBAL VARIABLES

# Uncomment for manual diagnostic
#NAMES <- c("Fe_ext","O2","Fe_free","Fur","RyhB","ROS","OxyR","Hpx","Suf","IscR-A","IscR-H","IscSUA","ErpA","NfuA")
#INPUTS <- c("Fe_ext","O2")
#BOOLEAN <- c("ROS","Hpx","ErpA","NfuA","IscSUA","IscR-A","IscR-H","OxyR")
#MULTIVALUED <- c("Fe_ext","O2","Fe_free","Fur","RyhB","Suf")
#OUTPUTS <- c("Fe_free","ROS","IscSUA","ErpA") #For the Labeller function
#THEOUTPUTS <- c("Fe_free","ROS","IscSUA","Suf", "ErpA") #For the plotting part
#IS.OUTPUTS <- c("IscSUA","Suf")


#Automatic setting
NAMES <- readLines(True_arguments[[1]])
INPUTS <- readLines(True_arguments[[2]])
BOOLEAN <- readLines(True_arguments[[3]])
MULTIVALUED <- readLines(True_arguments[[4]])
OUTPUTS <- readLines(True_arguments[[5]])
THEOUTPUTS <- readLines(True_arguments[[6]])
IS.OUTPUTS <- readLines(True_arguments[[7]])


VALUESFUR <- c(1,2)
VALUESMULTI <- c(0,1,2)
VALUESBOOL <- c(0,1)

## Colorscales
COLORSCALE_MULTI <- c("#ffffff","#c8afff", "#cc00ff")
COLORSCALE_BOOL <- c("#ffffff","#cc00ff")


DARK_COLORSCALES <- list(COLORSCALE_BOOL,COLORSCALE_MULTI)




#### ISOLATOR MAPPER ####

# PHASE I: INTEGRATE EACH MUTANT ONE BY ONE IN A LIST OF DATAFRAMES.



#Recuperer les dossiers
Folderz <- list.dirs()
Folders <- Folderz[2:length(Folderz)]

#Recuperer les mutants
#Mu <- getmutants("./Report_01/")
#Hardtest <- Mu[[8]]

#Get all Mutants names
allmutNames <- getAllMutants(Folders)



#Recuperer tous les tableaux des dossiers
ALLMUT <- lapply(Folders, reader)

#Formatter tous les tableaux
ALLMUT <- G_Formatter(ALLMUT,NAMES)




#On sépare la liste de mutants en 3 (simple, 2x, 3x)
#TEST00 <- ALLMUT[[1]]
#TEST10 <- ALLMUT[[2]]
#Test20 <- ALLMUT[[3]]



#Prendre tous les allmutnames...

No_Isc_no_Suf <- lapply(ALLMUT, MultiDetection.No.IS ,out1 = "IscSUA", out2 = "Suf")

#Tous les IS mutants
IS <- mapply(function(partMutNames,is_mutant){return(partMutNames[which(is_mutant == T)])}, partMutNames = allmutNames, is_mutant = No_Isc_no_Suf)


ALL.IS.MUT <- mapply(reader2, chem1 = Folders, MutIS =  IS)
ALL.IS.MUT <- G_Formatter(ALL.IS.MUT,NAMES)



#For specific mutants
DarkUltimateIsolatorPlotter <- function(list2simulations, MutantISList)
{
  AllPlots4mutants <- lapply(list2simulations, DarkIsolatorPlotter, inp = INPUTS, o = OUTPUTS, t = THEOUTPUTS, v = VALUESMULTI)
  names(AllPlots4mutants) <- MutantISList
  return(AllPlots4mutants)
}

DarkUltimateIsolatorPlotter_IS <- function(list2simulations, MutantISList)
{
  AllPlots4mutants <- lapply(list2simulations, DarkIsolatorPlotter, inp = INPUTS, o = OUTPUTS, t = IS.OUTPUTS, v = VALUESMULTI)
  names(AllPlots4mutants) <- MutantISList
  return(AllPlots4mutants)
}




#Scalable depending of ram usage.
DarkMegaIsolatorPlotter <- function(allSimulations,allMutantList)
{
  Omega <- mapply(DarkUltimateIsolatorPlotter, list2simulations = allSimulations, MutantISList = allMutantList)
  return(Omega)
}


DarkMegaIsolatorPlotter_IS <- function(allSimulations,allMutantList)
{
  Omega <- mcmapply(DarkUltimateIsolatorPlotter_IS, list2simulations = allSimulations, MutantISList = allMutantList, mc.cores = 3)
  return(Omega)
}




Omegamon <- DarkMegaIsolatorPlotter_IS(ALL.IS.MUT,IS)

Omegamon.merciful <- DarkMegaIsolatorPlotter(ALLMUT,allmutNames)

#### METTRE EN PLACE LE PRINTER
ULTIMATE_PRINTER <- function(WholePlotList,WholeMutantList,folderlist)
{
  noms <- gsub("\\./", "", folderlist)
  #print(noms)
  N_A_M_E_S <- paste(noms,"_Readouts.pdf",sep = "")
  mapply(ThePlotPrinter, Megalist2plots = WholePlotList, list2mutants =WholeMutantList, plotname = N_A_M_E_S)
}
### For IS mutants only
ULTIMATE_PRINTER_IS <- function(WholePlotList,WholeMutantList,folderlist)
{
  noms <- gsub("\\./", "", folderlist)
  #print(noms)
  N_A_M_E_S <- paste(noms,"_IS_mutants.pdf",sep = "")
  mapply(ThePlotPrinter, Megalist2plots = WholePlotList, list2mutants =WholeMutantList, plotname = N_A_M_E_S)
}



ULTIMATE_PRINTER(Omegamon.merciful,allmutNames,Folders)

ULTIMATE_PRINTER_IS(Omegamon,IS,Folders)

proc.time() - ptm



