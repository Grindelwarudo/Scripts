#Load the functions needed...
# set working directory here


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
required.packages <- c('here','data.table','ggplot2','tidyr','reshape','parallel','Matrix', 'gridExtra','reshape2','RColorBrewer','ggpubr','cowplot','grid')
ipkgs(required.packages)

setwd(here())
setwd("./GINsim_heatmaps")
print(getwd())


source("IsolatorHelper.R")

#GLOBAL VARIABLES


NAMES <- c("Fe_ext","O2","Fe_free","Fur","RyhB","ROS","OxyR","Hpx","Suf","IscR-A","IscR-H","IscSUA","ErpA","NfuA")
BOOLEAN <- c("ROS","OxyR","Hpx","IscR-A","IscR-H","IscSUA","ErpA","NfuA")
MULTIVALUED <- c("Fe_ext","O2","Fe_free","Fur","RyhB","Suf")
INPUTS <- c("Fe_ext","O2")
OUTPUTS <- c("Fe_free","ROS","IscSUA","ErpA") #For the Labeller function
THEOUTPUTS <- c("Fe_free","ROS","IscSUA","Suf", "ErpA") #For the plotting part
IS.OUTPUTS <- c("IscSUA","Suf")

VALUESFUR <- c(1,2)
VALUESSUF <- c(0,1,2)
VALUESBOOL <- c(0,1)

## Colorscales
COLORSCALES_ROS <- c("#ffffff","#92b8ff")
COLORSCALES_FE <- c("#ffffff","#98ff76", "#98b776")
COLORSCALES_SUF <- c("#ffffff","#c8afff", "#cc00ff")
COLORSCALES_ISCSUA <- c("#ffffff","#ff5200")
COLORSCALES_ERPA <- c("#ffffff","#ffc66b")

COLORSCALE_MULTI <- c("#ffffff","#c8afff", "#cc00ff")
COLORSCALE_BOOL <- c("#ffffff","#cc00ff")
#Legends <- c("0","Osc 0-1","1","Osc 0-1-2","Osc 1-2","2")


DARK_COLORSCALES <- list(COLORSCALE_BOOL,COLORSCALE_MULTI)
COLORSCALES <- list(COLORSCALES_FE,COLORSCALES_ROS, COLORSCALES_ISCSUA, COLORSCALES_SUF, COLORSCALES_ERPA)




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


#ok <- ALLMUT[[1]]
#okk <- ok[[1]]

#oke <- ok[[8]]
#e <- allmutNames[[1]]

#On sépare la liste de mutants en 3 (simple, 2x, 3x)
TEST00 <- ALLMUT[[1]]
TEST10 <- ALLMUT[[2]]
Test20 <- ALLMUT[[3]]

TESTS <- lapply(ALLMUT, function(x){return(x)})

#Nom des mutants de chaque 1x,2x,3x
Helmut <- allmutNames[[1]]
Helmut2 <- allmutNames[[2]]
Helmut3 <- allmutNames[[3]]

#Detection des mutants no Isc no Suf
f <- MultiDetection.No.IS(TEST00,"IscSUA","Suf")
g <- MultiDetection.No.IS(TEST10,"IscSUA","Suf")
h <- MultiDetection.No.IS(Test20,"IscSUA","Suf")

#Nombre mutants (0) vs nb de mut IS (1)
f0 <- length(f)
g0 <- length(g)
h0 <- length(h)

f1 <- length(which(f ==T))
g1 <- length(which(g ==T))
h1 <- length(which(h ==T))

#Get IS mutnames
IS.1x <- Helmut[which(f ==T)]

IS.2x <- Helmut2[which(g == T)]
IS.3x <- Helmut3[which(h == T)]




#Assemblage de la liste des mutants
IS <- list(IS.1x,IS.2x,IS.3x)

ALL.IS.MUT <- mapply(reader2, chem1 = Folders, MutIS =  IS)
ALL.IS.MUT <- G_Formatter(ALL.IS.MUT,NAMES)


IS.ilie <- ALL.IS.MUT[[1]]
Test01 <- IS.ilie[[2]]


Test02 <- Labeller(Test01,OUTPUTS,VALUESSUF)

Test0X <- Labeller(Test01,OUTPUTS,VALUESSUF)
Test03 <- lapply(IS.OUTPUTS, DarkIsolator, i = INPUTS, df= Test02)
Testx3 <- DarkIsolator(Test02, INPUTS, "Suf")



IS.ilie.name <- IS[[1]]


Mightysilie <- DarkIsolatorPlotter(Test01,INPUTS,OUTPUTS,IS.OUTPUTS,VALUESSUF)
Mightysilie <- IsolatorPlotter(Test01,INPUTS,OUTPUTS,THEOUTPUTS,VALUESSUF)

TEST00 <- ALLMUT[[1]]
TEST10 <- ALLMUT[[2]]
Test20 <- ALLMUT[[3]]


Test01 <- TEST00[[11]]

### Suite des analyses


Test04 <- L_abeller(Test01,OUTPUTS,VALUESSUF)




n <- length(Might)
nCol <- floor(sqrt(n))
grid.arrange(Might[[1]],Might[[2]],Might[[3]],Might[[4]],Might[[5]], ncol = 5, left = "Cat E1")



Mighty <- DarkIsolatorPlotter(Test01,INPUTS,OUTPUTS,THEOUTPUTS,VALUESSUF)
Might <- IsolatorPlotter(Test01,INPUTS,OUTPUTS,THEOUTPUTS,VALUESSUF)
MightMKII <- UltimateIsolatorPlotter(TEST00, Folders[1])

MightMKII <- UltimateIsolatorPlotter(ALLMUT[[3]], Folders[3])



a <- grid.arrange(Might,ncol = 5, top = "WT")


#For all mutants
UltimateIsolatorPlotter <- function(list2simulations, Folder)
{
  AllPlots4mutants <- lapply(list2simulations, IsolatorPlotter, inp = INPUTS, o = OUTPUTS, t = THEOUTPUTS, v = VALUESSUF)
  allmutNames <- getmutants(Folder)
  names(AllPlots4mutants) <- mutnames
  return(AllPlots4mutants)
}

#For specific mutants
DarkUltimateIsolatorPlotter <- function(list2simulations, MutantISList)
{
  AllPlots4mutants <- lapply(list2simulations, DarkIsolatorPlotter, inp = INPUTS, o = OUTPUTS, t = THEOUTPUTS, v = VALUESSUF)
  allmutNames <- MutantISList
  names(AllPlots4mutants) <- allmutNames
  return(AllPlots4mutants)
}



YourWaifuIsShit <- DarkUltimateIsolatorPlotter(All.IS.MUT[[1]], IS[[1]] )


MegaIsolatorPlotter <- function(allSimulations,folders)
{
  Omega <- mapply(UltimateIsolatorPlotter, list2simulations = allSimulations, Folder = folders)
  return(Omega)
}

DarkMegaIsolatorPlotter <- function(allSimulations,allMutantList)
{
  Omega <- mcmapply(DarkUltimateIsolatorPlotter, list2simulations = allSimulations, MutantISList = allMutantList, mc.cores = 3)
  return(Omega)
}



Omegamon <- MegaIsolatorPlotter(ALLMUT,Folders)

Omegamon.merciful <- DarkMegaIsolatorPlotter(ALLMUT,allmutNames)

Alphamon <- DarkMegaIsolatorPlotter(ALL.IS.MUT,IS)




ULTIMATE_PRINTER <- function(WholePlotList,WholeMutantList,folder)
{
  noms <- gsub("\\./", "", folder)
  N_A_M_E_S <- paste(noms,"_Modules.pdf",sep = "")
  mcmapply(ThePlotPrinter, Megalist2plots = WholePlotList, list2mutants =WholeMutantList, plotname = N_A_M_E_S, mc.cores = 3)
}

ULTIMATE_PRINTER(ALLMUT,allmutNames,Folders)



ThePlotPrinter(ALLMUT[[3]], allmutNames[[3]],"Triple Mutants.pdf")





#Darkversion (IS ONLY)
ThePlotPrinter(Alphamon[[1]], IS[[1]],"IS Simple Mutants.pdf", IS.OUTPUTS)

ThePlotPrinter(Alphamon[[2]], IS[[2]],"IS Double Mutants.pdf", IS.OUTPUTS)
ThePlotPrinter(Alphamon[[3]], IS[[3]],"IS Triple Mutants.pdf", IS.OUTPUTS)

#Darkversion: All Nodes
ThePlotPrinter(Omegamon.merciful[[1]], allmutNames[[1]],"Simple Mutants.pdf", THEOUTPUTS)

ThePlotPrinter(Omegamon.merciful[[2]], allmutNames[[2]],"Double Mutants.pdf", THEOUTPUTS)
ThePlotPrinter(Omegamon.merciful[[3]], allmutNames[[3]],"Triple Mutants.pdf", THEOUTPUTS)



