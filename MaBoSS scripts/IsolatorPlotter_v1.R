#Ce script permet de pouvoir générer les plots des trajectoires et des entropies.


myArgs <- commandArgs(trailingOnly = TRUE)
#print(myArgs)
#print(typeof(myArgs))
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else {
  print("Let's analyze these files.")
}

print(myArgs[length(myArgs)])

#Manual setup (diagnostic)
base <- "~/Documents/MaBoSS/MaBoSS-env-2.0/Prototype"
#setwd(base)

#Les trajectoires sont enregistrées sous format pdf.
base <- myArgs[length(myArgs)]
setwd(base)


#Check if the packages are installed. if not: install it
ipkgs <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#Load existing packages, install missing packages.
required.packages <- c('data.table','ggplot2','tidyr','reshape','parallel','Matrix', 'gridExtra')
ipkgs(required.packages)



#Manual stuff...
#TableDSim <- read.csv("Table_Des_Simulations.csv")


#Après avoir récupéré la table des simulations via le script python...
TableDSim <- read.csv(myArgs[length(myArgs)-1])
TableDSim.t <- t(TableDSim[,2:length(TableDSim)])
nom2col <- TableDSim.t[1,]
TableDesSimulations <- as.data.frame(TableDSim.t[2:dim(TableDSim.t)[1],], colnames(nom2col))
colnames(TableDesSimulations) <- nom2col
TableDesSimulations.sortedc <- TableDesSimulations[,order(names(TableDesSimulations))]
TableDesSimulations.sortedL <- as.data.frame(lapply(TableDesSimulations.sortedc, sort))
#On commence à faire une mégaliste contenant tous les noms de dossier bien ordonnés
#base <- "~/Documents/MaBoSS/MaBoSS-env-2.0/Prototype/"
#setwd(base)
Listdesconditions <- as.character(TableDesSimulations.sortedL$WT)


if(length(Listdesconditions) > 9){
  stop("The Simulations Table is missing something... Check if the Mutation.txt file contains all the mutants listed in the directories \n
       If the analysis is for less than 9 conditions, please edit the IsolatorPlotter.R file at the line 36.")
}


#Et on fabrique la gigamatrice contenant toutes les dataframes <3
GigaDF <- list()
GigaErr <- list()
Thenamelist <- colnames(TableDesSimulations.sortedL)


i <- 1
for (i in 1:dim(TableDesSimulations.sortedL)[2]){
  TheList <- as.character(TableDesSimulations.sortedL[,i])
  
  
  ListofErrors <- lapply(X = TheList,FUN =  function(x){
    
    
    DataE <- as.data.frame(fread(paste(x, paste(x, "probtraj_table.csv", sep = "_"), sep = "/"),
                                 sep = "\t", header = TRUE))
    g <- grep("ErrProb", colnames(DataE))
    DataE <- DataE[,c(1,g)]
    if (dim(DataE)[1] == 999) {
      res <- DataE
    } else { #Ca veut dire 1000
      res <- DataE[1:(dim(DataE)[1]-1),]
    }
    res <- Matrix(as.matrix(res), sparse = TRUE) 
    return(res)
  })
  
  GigaErr[[i]] <- ListofErrors
  
  
  ListofConditions <- lapply(X = TheList, FUN =  function(x){
    print(x)
    DataJ <- as.data.frame(fread(paste(x, paste(x, "probtraj_table.csv", sep = "_"), sep = "/"),
                                 sep = "\t", header = TRUE))
    g <- grep("Prob", colnames(DataJ))
    f <- grep("Err", colnames(DataJ))
    DataJ <- DataJ[,c(1, g[!(g %in% f)])]
    if (dim(DataJ)[1] == 999) {
      res <- DataJ
    } else { #Ca veut dire 1000
      res <- DataJ[1:(dim(DataJ)[1]-1),]
    }
    res <- Matrix(as.matrix(res), sparse = TRUE) 
    return(res)
  })
  GigaDF[[i]] <-  ListofConditions
}

#Giga.P = OK
#Giga.P <- do.call(cbind,GigaDF)
# for (i in 1:dim(TableDesSimulations.sortedL)[1]){
#   X <- GigaDF[i]
#   for (j in 1:length(X)) {
#     Giga.P[j,i] <- X
#   }
# }


############################################################


###########################################################

#Maintenant, passons aux tâches à effectuer \o/




#Deja, on recupere toutes les dataframes d'une condition:
# List2cond <- Giga.P[,1]
# cond <- as.data.frame(List2cond[1])

purinharumaki <- list() #Mapper des différents outputs selon la condition étudiée. (pour tous les mutants)

#Thx to Reol, on va pouvoir avoir une bonne structure à gérer pour faire toutes les tâches :)
Reol <- function(Giga_DF, Giga_Err, TabSim) {
  
  TabSim <- TableDesSimulations.sortedL
  
  
  yugioh <- colnames(TabSim)
  i <- 1
  AlphaList <- list() #Accumulateur des matrices correspondant à chaque mutant.
  
  for (i in 1:dim(TabSim)[2] ) { #Pour chaque mutant
    #print(i)
    Yami <- list()
    GigaP.colonne <- Giga_DF[[i]]
    GigaE.colonne <- Giga_Err[[i]]
    OmegaMatrix <- list() #(nrow = nombre d'output, actuellement manuel mais possibilité de le rendre auto.)
    #print(paste("Je pose la carte", yugioh[i], "face cachée."))
    j <- 1
    for (j in 1:length(GigaP.colonne)) { #Pour chaque condition
      # /!\ En brut concernant le début.
      #print(j)
      
      length(GigaP.colonne)
      
      #Diag
      #cond <- deeF
      #End diag
      cond <- GigaP.colonne[[j]]
      E.cond <- GigaE.colonne[[j]]
      summ <- summary(cond)
      sume <- summary(E.cond)
      print(dim(cond))
      rownames(cond) <- seq_along(cond[,1])
      rownames(E.cond) <- seq_along(E.cond[,1])
      cond <- data.frame(Origin      = rownames(cond)[summ$i],
                         Destination = colnames(cond)[summ$j],
                         Weight      = summ$x)
      E.cond <- data.frame(Origin   = rownames(E.cond)[sume$i],
                         Destination = colnames(E.cond)[sume$j],
                         Weight      = sume$x)
      cond <- spread(data = cond, key = Destination, value = Weight)
      E.cond <- spread(data = E.cond, key = Destination, value = Weight)
      cond[is.na(cond)] <- 0
      E.cond[is.na(E.cond)] <- 0
      cond <- cond[order(cond$Time),]
      E.cond <- E.cond[order(cond$Time),]
      cond$Origin <- NULL
      E.cond$Origin <- NULL


      #cond.name <- substring(as.character(TabSim[j,i]), 18)
      
      #SUPER IMPORTANT POUR LE NOM DE L'OUTPUT!!!!!!!!!!!!
      nom_complet <- as.character(TabSim[j,i])
      
      cond.name <- substring(nom_complet, gregexpr(pattern = "Fe_[0-9]", nom_complet)[[1]][[1]])
      #colnames(cond)<- gsub(paste(gsub("-", ".", TabSim[j, i]), "[.]", sep = ""), "", colnames(cond))
      #colnames(cond)[1:10]
      #print(cond.name)
      #MANUEL
      conditio.n <- substring(cond.name,1,9)
      #print(conditio.n)
      #print(cond.name)
      #On récupère le nom de colonne le plus long afin de chercher chaque noeud 1 à 1.
      noms <- colnames(cond)
      #print(grep(cond.name,noms))
      nbcharnoms <-nchar(noms)
      Tailledunom <- nchar(noms)
      Taillemax <- unlist(max(nchar(noms)))
      Lepluslongnom <- noms[match(Taillemax,Tailledunom)]
      
      #On décompose ce nom pour avoir une liste contenant les différents noeuds.
      sep0 <- gsub("\\]", "", gsub(".*\\[", "", Lepluslongnom))
      separation1 <- unlist(strsplit(sep0, '[-]'))
      separation2 <- separation1[separation1 != ""]
      sep3 <- separation2[separation2 != "Prob"]
      sep3 <- sort(sep3)
      
      #print(sep3)
      Yami[[j]] <- sep3
      print(Yami)
      #Get the multivalued nodes for a different processing.
      multivalued_nodes <- grep("_b[0-9]$",sep3, value = T)
      multivalued_nodes.id <- unique(gsub("_b[0-9]$", "", multivalued_nodes))
      #print(multivalued_nodes.id)
      #Et ensuite on réactualise la liste des outputs d'une traite
      sep3.5 <- grep("_b[0-9]$",sep3, value = T, invert = T)
      sep4 <- sort(c(sep3.5, multivalued_nodes.id))
      
      #On sépare les dataframes en fonction des colonnes contenant le noeud en question
      List2output <- list()
      k = 1
      
      #print(dim(Di))
      
      print(cond.name)
      #W = cond or E.cond
      trieur <- function(x, w){
        
        #C'est ici qu'on va séparer les actions en fonctions la nature du noeud or not / S'il fait partie des mutants.
        if (any(grepl(x, cond.name)) && any(grepl(paste("d_",x,sep= ""),cond.name))) 
        {
          
          print(x)
          print("lol")
        } else {
          
          if (!x %in% multivalued_nodes) 
          {
            #si c'est pas multivalué
            noms.euh <- colnames(w)
            colonnesspe <- grep(x,noms.euh, value = T)
            Did <- w[,colonnesspe]
            
            return(Did)
            # print(paste("Je place l'output",sep3[k], "dans le mutant", cond.name))
            # print("OK, boolean")
            
          } else {
            #print(paste("Multivalued:",sep4[k]))
            # print(paste("Dimension Digimon:",dim(Di)[2]))
            
            nomsX <- colnames(w)
            # print(paste("Longueur des noms Digi... :",length(nomsX)))
            #On get l'identifiant du noeud multivalué
            multiN <- sub("_b[0-9]$", "",x)
            #Enlever les trajectoires interdites
            col.chel <- list()
            erds <- grep(paste(multiN,"_b2", sep = ""), nomsX)
            erd <- grep(paste(multiN,"_b1",sep = ""), nomsX)
            NOOOO <- erds[which(!erds %in% erd)]
            col.chel <- NOOOO
            
            
            #Faire une nouvelle DF en enlevant la merde
            col_biz <- unique(unlist(col.chel))
            numberz <- c(1:dim(w)[2])
            real.col_numb <- numberz[which(!numberz %in% col_biz)]
            #print(paste("Avant :",dim(Di)))
            Dx <- w[,real.col_numb]
            #print(paste("Après - Dimension DJ :",dim(Dx)))
            # print(paste("Dimension DJ:",dim(Dx)[2]))
            noms2n <- colnames(Dx)
            
            
            #Si j'ai bien séparé les noeuds c bon
            # print(sep3[k])#diagnostic
            # print(length(noms2n))#diagnostic
            
            
            ##### TEST SI NECESSAIRE
            
            
            #Separation obligatoire
            b2 <- grep(paste(multiN,"_b2", sep = ""), noms2n)
            b1_1 <- grep(paste(multiN,"_b1", sep = ""), noms2n)
            b1_2 <- b1_1[which(!b1_1 %in% b2)]
            
            
            #Si c'est le noeud B1 (0-1)
            if (grepl("_b1",x)){
              #On get le _b1 only
              return(Dx[,b1_2])
              # print(paste("Je place l'output",sep3[k], "dans le mutant", cond.name))
            }else { #Si c'est le noeud B2 (1-2)
              return(Dx[,b2])
              # print(paste("Je place l'output",sep3[k], "dans le mutant", cond.name))
            }
          }
        }
      }
      
      
      
      List2Err <- lapply(X = sep3, FUN =  trieur, w = E.cond)
        
      List2Output <- lapply(X = sep3, FUN =  trieur, w = cond)
            
            
      
      #Dy <- Dx[,c(b1_2,b2)]
            #List2output[[k]] <- Dy
            
            #Diagnostic
            # ora <- Dx[,b1_2]
            # orb <- Dx[,b2]
            # aaa <- rowSums(ora)
            # bbb <- rowSums(orb)
            # w <- aaa + bbb
            # zz <- rowSums(Dy)
            # tail(w)
            # tail(zz)
            # orc <- Dx[,b1_1]
            # ccc <- rowSums(orc)
            # wai <- ccc/2 + bbb / 2
           
      
      # print("Output ?")
      # print(length(List2output))
      Temps <- cond$Time
      
      #Pour chaque dataframe on fait la somme de toutes les proba
      Proba.4.all.statesinattractor <- list() #Ca c'est la liste des dataframes processées pour chaque output.
      Err_Proba.4.all.statesinattractor <- list()
      
      
      l <- 1
        #print(length(List2output))
      
      #x = List2Output/Err[i], y = digits, z = Err prob ou normal. Si normal: "", sinon: "Err_"
      separator <- function(x, y, z){
        
        node.names <- colnames(x)[1]
        node.names <- gsub("\\]", "", gsub(".*\\[", "", node.names))
        node.names <- unlist(strsplit(node.names, '[-]'))
        node.names <- node.names[node.names != ""]
        if (length(node.names) == 1 ) {
          true.node <- node.names
        } else {
          true.node <- node.names[length(node.names)]
        }
        
        multivalued_nodes_o <- grep("_b[0-9]$",true.node, value = T)
        multivalued_nodes_id <- unique(gsub("_b[0-9]$", "", multivalued_nodes_o))
        
        
        nom_complet <- as.character(TabSim[j,i])
        cond.name <- substring(nom_complet, gregexpr(pattern = "Fe_[0-9]", nom_complet)[[1]][[1]])
        conditio.n <- substring(cond.name,1,9)
        
        
        outputname <- conditio.n
        
        if (any(grepl(true.node,multivalued_nodes_id))) {
          #print(paste(sep_B[l],"est dans", multivalued_nodes_id))
          #print(any(grepl(sep_B[l],multivalued_nodes_id)))
          # print(Thenamelist[i])
          # print(conditio.n)
          # print(sepZ[l])
          Dja <- as.data.frame(x)
          #Dja <- Da
          
          Proba <- round((rowSums(Dja)/2),digits = y )
          Dka <- data.frame(Temps)
          Proba.name <- paste(z,outputname)
          Dka$Proba.name <- unlist(Proba)
          colnames(Dka)[2] <- Proba.name
          #print(tail(Dka))
          
          return(Dka)
        } else {
          Djb <- as.data.frame(x)
          #Djb <- Dy
          Proba <- round(rowSums(Djb),digits = y )
          Dkb <- data.frame(Temps)
          Proba.name <- paste(z,outputname)
          Dkb$Proba.name <- unlist(Proba)
          colnames(Dkb)[2] <- Proba.name
          #print(tail(Dkb))
          return(Dkb)
        }
        
        
        
        
      } #Fin du séparateur
      

      
      
      Proba.4.all.statesinattractor <- lapply(X = List2Output, FUN = separator, y = 3, z = "")
      Err_Proba.4.all.statesinattractor <- lapply(X = List2Err, FUN = separator, y = 5, z = "Err")
      The_Errors <- lapply(X = Err_Proba.4.all.statesinattractor, FUN = function(x){
        nom2col <- colnames(x)
        res <- as.data.frame(x[,2])
        colnames(res)[1] <- nom2col[2]
        
        
        
        return(res)})
      
      TheFinalProba <- mapply(function(x,y){
        Ko <- as.data.frame(x)
        truname <- colnames(y)
        Ko$Err <- unlist(y)
        colnames(Ko)[3] <- truname
        return(Ko)
        }, Proba.4.all.statesinattractor, The_Errors, SIMPLIFY = F)
      
      #print(length(Proba.4.all.statesinattractor))
      OmegaMatrix[[j]] <- TheFinalProba
      
      
    } #Fin de la boucle tournant sur les conditions
    print(Yami)
    purinharumaki[[i]] <- Yami
    OmegaMatrix <- do.call(cbind,OmegaMatrix)
    AlphaList[[i]] <- OmegaMatrix
    #On recupere le nom de la condition sur laquelle la boucle tourne. 
  }
  return(AlphaList)
}

Crusher.P <- Reol(GigaDF, GigaErr, TableDesSimulations.sortedL)
print(length(Crusher.P))




#diag
Julz <- Crusher.P[[1]]
Julzi <-Julz[[1]]
Julzz <- Julzi[[1]]



# #Preparing AllTasks :)
# GigaP.colonne <- Giga.P[,1]
# cond <- as.data.frame(GigaP.colonne[1])
# #On récupère le nom de colonne le plus long afin de chercher chaque noeud 1 à 1.
# noms <- colnames(cond)
# nbcharnoms <-nchar(noms)
# Tailledunom <- nchar(noms)
# Taillemax <- unlist(max(nchar(noms)))
# Lepluslongnom <- noms[match(Taillemax,Tailledunom)]
# 
# #On décompose ce nom pour avoir une liste contenant les différents noeuds.
# separation1 <- unlist(strsplit(Lepluslongnom, '[.]'))
# separation2 <- separation1[separation1 != ""]
# sep3 <- separation2[separation2 != "Prob"]
# sep3 <- sort(sep3)
# 
# #Get the multivalued nodes for a different processing.
# multivalued_nodes <- grep("_b[0-9]$",sep3, value = T)
# multivalued_nodes.id <- unique(gsub("_b[0-9]$", "", multivalued_nodes))
# print(multivalued_nodes.id)
# 
# #Et ensuite on réactualise la liste des outputs d'une traite
# sep3.5 <- grep("_b[0-9]$",sep3, value = T, invert = T)
# sep4 <- sort(c(sep3.5, multivalued_nodes.id))

#GO FOR TASK 1
#Task 1 : Faire des graphiques représentant l'évolution d'un output en fonction du temps, et ce pour les différentes conditions.
#Pour chaque output un fichier pdf sera généré.
#Pour chaque mutation un dossier sera généré, et il contiendra les différents PDF.



Thenamelist <- colnames(TableDesSimulations.sortedL)


################# READY FOR THE TASK 2 ##################

#Task2: Générer des barplot représentant la proportion à l'état stable de chaque output. Chaque plot = 1 output, 1 fichier = 1 mutant.

i <- 1
j <- 1
k <- 1
for (i in 1:length(Crusher.P)) {
  print(Thenamelist[i])
  #Pour chaque mutant, on doit d'abord get la liste des noeuds...
  GigaP.colonne <- GigaDF[[i]]
  cond <- GigaP.colonne[[1]]
  summ <- summary(cond)
  rownames(cond) <- seq_along(cond[,1])
  cond <- data.frame(Origin      = rownames(cond)[summ$i],
                     Destination = colnames(cond)[summ$j],
                     Weight      = summ$x)
  cond <- spread(data = cond, key = Destination, value = Weight)
  
  cond[is.na(cond)] <- 0
  cond$Origin <- NULL
  #On récupère le nom de colonne le plus long afin de chercher chaque noeud 1 à 1.
  noms <- colnames(cond)
  print(dim(cond))
  nbcharnoms <-nchar(noms)
  Tailledunom <- nchar(noms)
  Taillemax <- unlist(max(nchar(noms)))
  print(Taillemax)
  Lepluslongnom <- noms[match(Taillemax,Tailledunom)]
  
  #On décompose ce nom pour avoir une liste contenant les différents noeuds.
  sep0 <- gsub("\\]", "", gsub(".*\\[", "", Lepluslongnom))
  separation1 <- unlist(strsplit(sep0, '[-]'))
  separation2 <- separation1[separation1 != ""]
  sep3 <- separation2[separation2 != "Prob"]
  sep3 <- sort(sep3)
  
  #Get the multivalued nodes for a different processing.
  multivalued_nodes <- grep("_b[0-9]$",sep3, value = T)
  multivalued_nodes.id <- unique(gsub("_b[0-9]$", "", multivalued_nodes))
  print(multivalued_nodes.id)
  
  #Et ensuite on réactualise la liste des outputs d'une traite
  sep3.5 <- grep("_b[0-9]$",sep3, value = T, invert = T)
  sep4 <- sort(c(sep3.5, multivalued_nodes.id))
  
  namecol <- Thenamelist[i]
  print(namecol) #Nom du mutant 
  
  #Ensuite on commence les festivités
  gigaM <- Crusher.P[[i]]
  #Définissons des accumulateurs (List2df va contenir tous les outputs)
  list2df <- list()
  #Pour chaque output...
  for (j in 1:dim(gigaM)[1]) {
    
    #eldf contiendra 1 output & toutes les conditions à chaque fois.)
    edf <- as.data.frame(gigaM[j,1])
    print(sep3[j])
    naemu <- sub("X.", "", colnames(edf)[2])
    # print("Naemu-chan")
    # print(naemu)
    l1 <- as.numeric(edf[dim(edf)[1],c(2,3)])
    eldf <- data.frame(l1)
    colnames(eldf) <- naemu
    
    #Pour chaque condition...
    for (k in 2:dim(gigaM)[2]) {
      xdf <- as.data.frame(gigaM[j,k])
      #print(colnames(xdf))
      l2 <- as.numeric(xdf[dim(xdf)[1],c(2,3)])
      xddf <- data.frame(l2)
      eldf <- cbind(eldf,l2)
      naemux <- sub("X.", "", colnames(xdf)[2])
      #print("Watashi Wa...")
      #print(k)
      # print(naemux)
      colnames(eldf)[k] <- naemux
    }
    #print(dim(eldf))
    rownames(eldf) <- c("Proba_cumulee", "sd")
    list2df[[j]] <- eldf
    #print(list2df)
    print(tail(eldf))
    
    eldf <- data.frame()
  }
  #Faire les barplots représentant la proba en fonction des conditions.
  plotzu <- list()
  n <- 1
  for (n in 1:length(list2df)) {
    #print(list2df[n])
    #D'abord melt la df en question...
    tf2 <- as.data.frame(list2df[[n]])
    olol.f <- tf2[1,] 
    olol.t <- melt(olol.f)
    sx <- as.data.frame(t.data.frame(tf2))
    sd <- sx$sd
    sd <- replace(sd, sd ==0, NA)
    #Codé Manuellement, personnalisable
    titre2col <- c("Condition", "Proba_cumulee")
    colnames(olol.t) <- titre2col
    #On va chercher automatiquement l'output étudié pour l'inclure dans le titre (+ de lisibilité)
    futurtitreduplot <- as.character(olol.t[1,1])
    #print(futurtitreduplot)
    #Pour trouver l'output, on prend les 10 1res lettres de la condition, et on cherche dans sep3 (qui contient la liste des outputs).
    origin <- sub("_F.*", "", futurtitreduplot)
    p <-ggplot(data=as.data.frame(olol.t), aes(x=Condition, y=Proba_cumulee, fill = Proba_cumulee)) + ylim(0, 1.05)+ theme_minimal() + geom_bar(stat="identity") + geom_errorbar(aes(ymin = Proba_cumulee-sd, ymax = Proba_cumulee+sd), width = .2) + geom_text(aes(label=Proba_cumulee), vjust=2, color="white", size=2)  + theme(plot.title= element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=paste("Plot of Cumulated probabilites for",sep3[n]), 
                                                                    x="Condition", y = "Cumulated probability")

      plotzu[[n]] <- p
    
    
  }
  
  #geom_errorbarh(aes(xmin = upper, xmax = lower))
  
  
  if (!Thenamelist[i] %in% "WT") {
    output.name <- paste("Steady states - Cumulative probabilities in a",Thenamelist[i] ,"mutant.pdf") 
  } else {
    output.name <- paste("Steady states - Cumulative probabilities in a",Thenamelist[i], "strain.pdf")
  }
  pdf(output.name)
  bquiet = lapply(X = plotzu,FUN = print)
  dev.off()
}


# DYNAMIQUE SELON LE TEMPS
# i <- 1
# for (i in 1:length(Crusher.P)){
#   GigaP.colonne <- GigaDF[[i]]
#   cond <- GigaP.colonne[[1]]
#   summ <- summary(cond)
#   rownames(cond) <- seq_along(cond[,1])
#   cond <- data.frame(Origin      = rownames(cond)[summ$i],
#                      Destination = colnames(cond)[summ$j],
#                      Weight      = summ$x)
#   cond <- spread(data = cond, key = Destination, value = Weight)
#   
#   cond[is.na(cond)] <- 0
#   cond$Origin <- NULL
#   #On récupère le nom de colonne le plus long afin de chercher chaque noeud 1 à 1.
#   noms <- colnames(cond)
#   print(dim(cond))
#   nbcharnoms <-nchar(noms)
#   Tailledunom <- nchar(noms)
#   Taillemax <- unlist(max(nchar(noms)))
#   print(Taillemax)
#   Lepluslongnom <- noms[match(Taillemax,Tailledunom)]
#   
#   #On décompose ce nom pour avoir une liste contenant les différents noeuds.
#   sep0 <- gsub("\\]", "", gsub(".*\\[", "", Lepluslongnom))
#   separation1 <- unlist(strsplit(sep0, '[-]'))
#   separation2 <- separation1[separation1 != ""]
#   sep3 <- separation2[separation2 != "Prob"]
#   sep3 <- sort(sep3)
#   
#   
#   
#   #D'abord on prend le nom du mutant
#   nameofcol <- Thenamelist[i]
#   #Et la matrice contenant les dataframes (col = conditions, row = output)
#   gigaM <- Crusher.P[[i]]
#   
#   j <- 1
#   #On parcourt chaque output
#   
#   list2df <- list()
#   
#   for (j in 1:dim(gigaM)[1]) {
#     #On prend juste la 1re dataframe qu'on va ensuite bind aux autres DF =)
#     Arupha <- as.data.frame(gigaM[j,1])  
#     print(sep3[j])
#     naemu <- sub("X.", "", colnames(Arupha)[2])
#     # print("Naemu-chan")
#     # print(naemu)
#     l1 <- as.numeric(Arupha[,2])
#     eldf <- data.frame(cbind(Arupha[,1],l1))
#     colnames(eldf) <- c("Temps",naemu)
#     
#     k <- 2
#     for (k in 2:dim(gigaM)[2]) {
#       #Là on procède au bind, mais en faisant attention à prendre seulement la 2e colonne (la 1re = time)
#       
#       Beeta <- as.data.frame(gigaM[j,k])
#       naemu <- sub("X.", "", colnames(Beeta)[2])
#       
#       eldf <- cbind(eldf,Beeta[,2])
#       colnames(eldf)[k+1] <- naemu
#     }
#     #print(any(is.na(Arupha)))
#     #Des qu'Arupha = done, on procède à la création des plots ehwai!
#     #/!\ adapter le nombre (3) codé en brut au nb de plots qu'on veut faire!!
#     longueur <- length(eldf)-1
#     if (longueur%%3 == 0) {
#       a <- 1
#       d <- 1
#       plotz <- list()
#       while (a< (length(eldf))) {
#         b <- a +1
#         c <- a + 3  
#         pltName <- paste( 'a', a, sep = '' )
#         x <- data.frame(eldf[,c(1,b:c)])
#         #print(any(is.na(x)))
#         
#         #Je sais que les warnings sont déclenchés au MELT J'EN SUIS SUR!!
#         z <- ggplot(melt(x, id.vars = "Temps"), aes(x=Temps, y=value, group=variable, shape= variable,colour=variable)) + geom_point() + scale_y_continuous(limits=c(0,1))
#         plotz[[pltName]] <- z
#         a <- a + 3
#         d <- d + 1
#       }
#     } else {
#       a <- 1
#       d <- 1
#       last <- 0
#       plotz <- list()
#       while (a< (length(eldf) - (length(eldf)%%3))) {
#         b <- a +1
#         c <- a + 3  
#         pltName <- paste( 'a', a, sep = '' )
#         x <- data.frame(eldf[,c(1,b:c)])
#         #print(any(is.na(x)))
#         
#         #Je sais que les warnings sont déclenchés au MELT J'EN SUIS SUR!!
#         z <- ggplot(melt(x, id.vars = "Temps"), aes(x=Temps, y=value, group=variable, shape= variable,colour=variable)) + geom_point() + scale_y_continuous(limits=c(0,1))
#         plotz[[pltName]] <- z
#         a <- a + 3
#         d <- d + 1
#       }
#     }
#     last <- a
#     alpha <- data.frame(eldf[,c(1,last:length(eldf))])
#     #Je sais que les warnings sont déclenchés au MELT J'EN SUIS SUR!!
#     omega <- ggplot(melt(alpha, id.vars = "Temps"), aes(x=Temps, y=value, group=variable, shape= variable,colour=variable)) + geom_point() + scale_y_continuous(limits=c(0,1))
#     plotz[["Last"]] <- omega
#     
#     
#     #Necessite le package gridExtra chargé au préalable!
#     nom2fichier <- paste(nameofcol,sep3[j],".pdf",sep = "_")
#     pdf(nom2fichier)
#     bquiet = lapply(plotz, print)
#     dev.off()
#     
#   }
#   listofiles <- list.files(base)
#   nom2dossier <- paste(Thenamelist[i],"Traj_full_conditions",sep = "_")
#   dir.create(nom2dossier)
#   
# }











cat("Done")


###################################################################




  # Thenamelist2 <- Thenamelist2[which(!Thenamelist[i] %in% Thenamelist2)]



# #Partie du script permettant de montrer la répartition des états dans l'attracteur.
# Df_attractors <- D[nrow(D),]
# Df_attractors.simple <- subset(Df_attractors,,select= which(c(TRUE,(Df_attractors[which(Df_attractors$Time == "99.9"), ] >= 0.001)[2:ncol(Df_attractors)])))
# lbls <- paste(names(Df_attractors.simple2), sep="")
# Df_attractors.simple2 <- unlist(Df_attractors.simple)
# melt(Df_attractors.simple, id.vars = "Time")
# 
# #WARNING: NEED A PERCENTAGE FOR THE PIE TO WORK
# Attractors.prop <- as.numeric(Df_attractors.simple) * 100
# pie(Attractors.prop[2:length(Attractors.prop)], labels = lbls[-1],main = "Pie Chart of the states in the attractors")
# 
# pdf("Attractor repartition")
# bquiet = pie
# dev.off
# 
