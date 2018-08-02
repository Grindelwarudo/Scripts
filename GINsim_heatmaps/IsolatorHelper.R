#### All functions to make the analysis <3 ####




#Read content
#Fonction pour lire un fichier CSV selon le chemin
reader <- function(chem1) #chem1 = Dossier contenant un mutant donné
{
  #print(chem1)
  tp <- list.files(chem1)
  chemtp <- paste(chem1,"/",tp,sep = "")
  return(lapply(chemtp, read.table, sep = "," ,stringsAsFactors = FALSE, header = F))
}

#Lit les csv associés aux mutants IS, Specifique aux mutants IS
reader2 <- function(chem1, MutIS)
{
  print(chem1)
  tp <- list.files(chem1)
  tp2 <- get.IS.filenames(tp,MutIS)
  #print(tp2)
  chemtp <- paste(chem1,"/",tp2,sep = "")
  #print(chemtp)
  return(lapply(chemtp, read.table, sep = "," ,stringsAsFactors = FALSE, header = F))
}

#Extract the mutant name
extractor <- function(mot)
{
  xx <- substring(mot, gregexpr(pattern = "__", mot)[[1]][[1]]+2)
  xy <- substring(xx, 1, nchar(xx) - 4)
  return(xy)
}

#Extract all the mutants in a folder
getmutants <- function(chem1)
{
  tp <- list.files(chem1)
  return(lapply(tp, extractor))
}

#Extract all mutants in all folders
getAllMutants <- function(folders)
{
  allMut <- lapply(folders, getmutants)
  return(allMut)
}


#On formatte en enlevant la 1re colonne inutile & en mettant les noms des noeuds (hardcoded)
#Put the colnames and sort the inputs
Formatter <- function(list2DF,noms)
{
  a <- lapply(list2DF,function(x) {
    colnames(x) <- noms
    return(x)
  })

  b <- lapply(a, function(x) {
    x[order(x[,1]),]  })
  return(b)
}

#Apply the formatter to all folders
G_Formatter <- function(Bigstuff,nomsx)
{
  lapply(Bigstuff, Formatter, noms = nomsx)
}



#Vérifie si le nombre se trouve dans l'intervalle entre 2 nombres (exclusif)
checkInterval <- function(number, interval)
{
  #print(interval)
  #print(number > interval[1] && number < interval[2])
  #print(grepl("~", number) || grepl("SS",number) || grepl("X",number) || grepl("-",number))
  if (grepl("~", number) || grepl("SS",number) || grepl("X",number) || grepl("-",number)){
    return(FALSE)
  }
  else
  {
    numb <- as.numeric(number)
    return(numb > interval[1] && numb < interval[2])
  }
}

#Vérifie tous les intervalles pour voir si un nombre est dans l'un d'entre eux
AllIntervals <- function(number,intervals)
{
  return(lapply(intervals, checkInterval, number = number))
}


#Compares intervals... (to masterize!!!)
#False = Pas d'oscillation, True = oscillation
comparator <- function(number, values)
{
  intervalles <- unlist(AllIntervals(number,values))
  #print(intervalles)
  #print(number)
  return(any(intervalles))
}

splitWithOverlap <- function(vec, seg.length, overlap) {
  if (length(vec) > 2)
  {
    starts = seq(1, length(vec), by=seg.length-overlap)
    ends   = starts + seg.length - 1
    ends[ends > length(vec)] = length(vec)

    a <- lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
    return(a[1:length(a)-1])
  }
  else #Faudra que je prenne le cas où c'est moins de 2... mais normalement impossible vu que c'est censé être booléen à minima...
  {
    #print("\n")
    #print(vec)
    return(vec)
  }
}

#returns the nature of the attractor: (Cyclic attractor or Steady state)
isCA <- function(state, valor)
{
  valeurs <- splitWithOverlap(valor,2,1) #Hardcoded
  #print(valeurs)
  booltable <- unlist(lapply(state, comparator, values = valeurs))
  #print("\n")
  #print(any(booltable))
  if (any(booltable))
  {
    return("~")
  }
  else
  {
    return("SS")
  }
}


isOsc <- function(state, valor, output)
{
  #print(output)
  i <- grep(output,NAMES)
  nodevalue <- state[[i]]
  #print(nodevalue)
  valeurs <- splitWithOverlap(valor,2,1) #Hardcoded
  booltable <- comparator(nodevalue, values = valeurs)
  if (any(booltable))
  {
    #print(nodevalue)
    return("~")
  }
  else
  {
    return("-")
  }
}


#returns the table with a column representing the attractors. Valor = valeurs multivalué
CaSSetter <- function(table, value)
{
  Attractor_type <- apply(table,1, isCA, valor = value)
  #print(Attractor_type)
  #print(table)
  table$"Attracteur" <- Attractor_type
  return(table)
}

#Compares Suf value in order to determine if the oscillation range is 0-1/1-2 or 0-1-2
comp_suf <- function(state, valuesFur)
{
  valeurs <- splitWithOverlap(valuesFur,2,1) #Hardcoded
  #print(IsFurElevated(state["Fur"],valuesFur))
  #print((IsFurElevated(state["Fur"], valuesFur) && state["Suf"] == 1 && state["Attracteur"] == "~"))
  #print(state["Suf"])
  #print(as.numeric(state[["Suf"]]) == 1.0)
  #print(typeof(state[["Suf"]]))
  #print(identical(x = state[["Suf"]], y = 1.0, num.eq = T))
  #print(state["Attracteur"])
  Suf <- as.numeric(state[["Suf"]])
  #print(Suf)

  if (IsFurElevated(state["Fur"], valuesFur) && Suf == 1.0 && state["Attracteur"] == "~")
  {
    return("~")
  }
  else
  {
    if (state["Attracteur"] == "SS")
    {
      return("-")
    }
    else
    {

      if (comparator(Suf, valeurs))
      {
        return("~")
      }
      else
      {
        #print(state)
        #print(Suf)
        #print((comparator(Suf, valeurs)))
        return("-")
      }
    }
  }
}



IsFurElevated <- function(number, valuesFur)
{
  #print(valuesFur)
  #print(valuesFur[2])
  #print(valuesFur[3])
  if (number > valuesFur[2] && number < valuesFur[3])
  {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

SufSetter <- function(table,valeurs2Fur)
{
  TheSuf <- apply(table,1,comp_suf, valuesFur= valeurs2Fur)
  #print(TheSuf)
  table$"Suf_behavior" <- TheSuf
  return(table)
}

#Prend les outputs & décrit les oscillations.
OscillationDetector <- function(table, output) #value = boolean or multivalued
{
  #print(output)
  o <- paste(output,"behavior")
  #value <- Gateway(output)
  #print(value)
  Attractor_type <- apply(table,1, isOsc, valor = VALUESMULTI, output= output) #Value Hardcoded
  #print(Attractor_type)
  #print(table)
  table$o <- Attractor_type
  i <- grep("~",colnames(table))
  #print(colnames(table)[i])
  colnames(table)[i] <- o
  #print(colnames(table)[i])
  print(Attractor_type)
  #return(table)
}

Gateway <- function(x)
{
  if (x %in% BOOLEAN)
  {
    return(VALUESBOOL)
  }
  else if (x %in% MULTIVALUED)
  {
    return(VALUESMULTI)
  }
  else
  {
    print("lol")
  }
}

OutputBehaviors <- function(tableau, listoutputs)
{
  names <- paste(listoutputs, "behavior",sep = "_")
  Tables <- lapply(listoutputs, OscillationDetector, table = tableau)
  Table <- as.data.frame(lapply(Tables, cbind))
  colnames(Table) <- names
  Tablex <- cbind(tableau,Table)
  #print(Tablex)
  return(Tablex)

}

#3 étapes: 1) Définir si c'est un attracteur cyclique ou un état stable 2) Définir le cpt de Suf 3) Définir le cpt des autres outputs
L_abeller <- function(Tableau, OutputList, ValuesVector)
{
  a <- OutputBehaviors(Tableau, OutputList)
  #print("OK Output")
  b <- CaSSetter(a,ValuesVector)
  #print("OK CA vs SS")
  c <- SufSetter(b, ValuesVector)
  #print(typeof(c))
  return(c)
}

checkInterval2 <- function(number, interval) #comparaison entre 2 valeurs,
{
    numb <- as.numeric(number)
    return(numb > interval[1] && numb < interval[2])
}


#Subset data to get the desired output
Subsetter <- function(Tableau,inputlist,singleoutput)
{
  #print(inputlist)
  #print(singleoutput)
  #print(colnames(Tableau))
  O <- grep(singleoutput,colnames(Tableau))
  #print(O)
  I <- lapply(inputlist, grep, x = colnames(Tableau))
  #print(I)
  IO <- unlist(c(I,O))
  #print(IO)
  return(subset(Tableau, select = IO))
}

#Discretise les valeurs du tableau contenant l'output désiré
discretisator <- function(Table,outputname)
{
  beh <- grep("behavior",colnames(Table))
  bo <- grep(outputname,colnames(Table))
  OOF <- bo[which(!bo %in% beh)]
  OOJ <- bo[which(bo %in% beh)]
  e <- Table[,OOF]
  j <- Table[,OOJ]
  nomdetype <- paste(outputname,"type", sep = "_")
  if (outputname == "Suf")
  {
    ox <- lapply(e, discretisation_multi, valeursmulti = VALUESMULTI)
    oy <- as.factor(suf_discretisation(e,j))
    Table$X <- oy
    i <- grep("X",colnames(Table), fixed = T)
    colnames(Table)[i] <- nomdetype
    return(Table)
  }
  else if (any(grepl(outputname, x = MULTIVALUED)))
  {
    ox <- as.factor(unlist(lapply(e, discretisation_multi, valeursmulti = VALUESMULTI)))
    Table$X <- ox
    i <- grep("X",colnames(Table), fixed = T)
    colnames(Table)[i] <- nomdetype
    return(Table)
  }
  else
  {
    ox <- as.factor(unlist(lapply(e, discretisation_bool)))
    Table$X <- ox
    i <- grep("X",colnames(Table), fixed = T)
    colnames(Table)[i] <- nomdetype
    return(Table)
  }
}


suf_discretisation <-function(etat, typeOsc){
  a <- mapply( function(x,y) #x = valeur de Suf, y  = nature de l'attracteur, z = valeursmulti
  {
    #print(x)
    #print(discretisation_multi(x,z))
    if (discretisation_multi(x, VALUESMULTI) == "Low")
    {
      return("Low")
    }
    else if (discretisation_multi(x,VALUESMULTI) == "Medium" && y == "X")
    {
      return("High")
    }
    else if (discretisation_multi(x,VALUESMULTI) == "Medium" && y != "X")
    {
      return("Medium")
    }
    else
    {
      return("High")
    }
  }
  , x= etat, y = typeOsc)
  return(a)
}

discretisation_multi <- function(x, valeursmulti) #Valeursmulti = {0,1,2}, hardcoded
{
  n <- as.numeric(x)
  #print(valeursmulti)
  if (n == valeursmulti[1])
  {
    return("Low")
  }
  else if (n > valeursmulti [1] && n <= valeursmulti[2])
  {
    return("Medium")
  }
  else
  {
    return("High")
  }

}

discretisation_bool <- function(x)
{
  n <- as.numeric(x)
  if (n == 0)
  {
    return("OFF")
  }
  else
  {
    return("ON")
  }
}

Isolator <- function(df,i,o)
{
  S1 <- Subsetter(df,i,o)
  S2 <- discretisator(S1,o)
  return(S2)
}

#TestIsolator 1color scale continuous 4 plotting
DarkIsolator <- function(df,i,o)
{
    S1 <- Subsetter(df,i,o)
    return(S1)
}


plotter <- function(Table)
{

  y <- colnames(Table)[3]
  e <- grep(y,THEOUTPUTS)
  switch(e,
         {
           couleur <- COLORSCALES[[1]]
           behave <- grep("behavior", colnames(Table))
           type <- grep("type", colnames(Table))
           Arupha <- ggplot(data = Table, aes(x = O2, y = Fe_ext)) + theme_minimal() + geom_tile(aes(fill = Table[,dim(Table)[2]])) +  geom_text(aes(label = Table[,behave]), size = 10)+ scale_y_reverse() + scale_x_discrete(position = "top")
           Arupha <- Arupha + scale_fill_manual("",values = c("Low"=couleur[[1]],"Medium"=couleur[[2]], "High"=couleur[[3]] )) + theme(plot.caption= element_text(hjust = 0.5, size = 20)) + labs(caption =y)
         },
         {
           couleur <- COLORSCALES[[2]]
           #print(couleur[[1]])
           behave <- grep("behavior", colnames(Table))
           type <- grep("type", colnames(Table))
           Arupha <- ggplot(data = Table, aes(x = O2, y = Fe_ext)) + theme_minimal() + geom_tile(aes(fill = Table[,type])) +  geom_text(aes(label = Table[,behave]), size = 10)+ scale_y_reverse() + scale_x_discrete(position = "top")
           Arupha <- Arupha + scale_fill_manual("",values = c("OFF"=couleur[[1]],"ON"=couleur[[2]]))  + theme(plot.caption= element_text(hjust = 0.5, size = 20)) + labs(caption =y)
           return(Arupha)
         },
         {
           couleur <- COLORSCALES[[3]]
           #print(output)
           #print(couleur)
           behave <- grep("behavior", colnames(Table))
           type <- grep("type", colnames(Table))
           Arupha <- ggplot(data = Table, aes(x = O2, y = Fe_ext)) + theme_minimal() + geom_tile(aes(fill = Table[,type])) +  geom_text(aes(label = Table[,behave]), size = 10)+ scale_y_reverse() + scale_x_discrete(position = "top")
           Arupha <- Arupha + scale_fill_manual("",values = c("OFF"=couleur[[1]],"ON"=couleur[[2]]))  + theme(plot.caption= element_text(hjust = 0.5, size = 20)) + labs(caption =y)
         },
         {
           couleur <- COLORSCALES[[4]]
           behave <- grep("behavior", colnames(Table))
           type <- grep("type", colnames(Table))
           Arupha <- ggplot(data = Table, aes(x = O2, y = Fe_ext)) + theme_minimal() + geom_tile(aes(fill = Table[,dim(Table)[2]]))  +  geom_text(aes(label = Table[,behave]), size = 10)+ scale_y_reverse()  + scale_x_discrete(position = "top")
           Arupha <- Arupha + scale_fill_manual("",values = c("Low"=couleur[[1]],"Medium"=couleur[[2]], "High"=couleur[[3]] )) + theme(plot.caption= element_text(hjust = 0.5, size = 20)) + labs(caption =y)
         },
         {
           couleur <- COLORSCALES[[5]]
           behave <- grep("behavior", colnames(Table))
           type <- grep("type", colnames(Table))
           Arupha <- ggplot(data = Table, aes(x = O2, y = Fe_ext)) + theme_minimal() + geom_tile(aes(fill = Table[,type])) + geom_text(aes(label = Table[,behave]), size = 10)+ scale_y_reverse() + scale_x_discrete(position = "top")
           Arupha <- Arupha + scale_fill_manual("",values = c("OFF"=couleur[[1]],"ON"=couleur[[2]]))+ theme(plot.caption= element_text(hjust = 0.5, size = 20)) + labs(caption =y)
         }
  )
  return(Arupha)
}

#1 color, plusieurs teintes
Darkplotter <- function(Table)
{

  y <- colnames(Table)[3]
  #print(y)
  #e <- grep(y,IS.OUTPUTS)
  #print(grepl(y,BOOLEAN))
  if (any(grepl(y,BOOLEAN)))
         {
           couleur <- DARK_COLORSCALES[[1]]
           #print(couleur)
           behave <- grep("behavior", colnames(Table))
           #type <- grep("type", colnames(Table))
           Arupha <- ggplot(data = Table, aes(x = O2, y = Fe_ext)) + theme_minimal() + geom_tile(color="black",size=2,aes(fill = Table[,3])) +  geom_text(aes(label = Table[,behave]), size = 10)+ scale_y_reverse() + scale_x_discrete(position = "top")
           Arupha <- Arupha + scale_fill_continuous("value",low = COLORSCALE_BOOL[[1]],high= COLORSCALE_BOOL[[2]], limits = c(0,1), breaks = c(0,0.5,1)) + theme(plot.caption= element_text(hjust = 0.5, size = 20)) + labs(caption =y)

         }
  else
         {
           couleur <- DARK_COLORSCALES[[2]]
           #print(couleur[[1]])
           behave <- grep("behavior", colnames(Table))
#cheatcode : 3e colonne = la bonne
           Arupha <- ggplot(data = Table, aes(x = O2, y = Fe_ext)) + theme_minimal() + geom_tile(color="black",size=2,aes(fill = Table[,3]))+ scale_colour_gradient2() +  geom_text(aes(label = Table[,behave]), size = 10)+ scale_y_reverse() + scale_x_discrete(position = "top")
           Arupha <- Arupha + scale_fill_continuous("value",low = COLORSCALE_MULTI[[1]],high= COLORSCALE_MULTI[[3]], limits = c(0,2))+ theme(plot.caption= element_text(hjust = 0.5, size = 20)) + labs(caption =y)
           return(Arupha)
         }
  return(Arupha)
}



#A partir des résultats de simulation pour une situation, genère des plots pour les outputs étudiés ("THEOUTPUTS")
IsolatorPlotter <- function(SimulationResultsTable, inp, o,t, v) #i = INPUTS, o = OUTPUTS, t = THEOUTPUTS, v = Multivalued values
{
  #On indique les caractéristiques d'oscillation & des attracteurs pour nos outputs préférés
  a <- L_abeller(SimulationResultsTable,o,v)
  #On isole ensuite les outputs sous forme de tableaux différents, et on discrétise la valeur du noeud (double -> Qualitative)
  b <- lapply(t, Isolator, i = inp, df = a)
  #Et on plotte chaque output
  c <- lapply(b, plotter)
  return(c)
}

#A partir des résultats de simulation pour une situation, genère des plots pour les outputs étudiés ("IS.Outputs")
DarkIsolatorPlotter <- function(SimulationResultsTable, inp, o,t, v) #i = INPUTS, o = OUTPUTS, t = THEOUTPUTS, v = Multivalued values
{
  #On indique les caractéristiques d'oscillation & des attracteurs pour nos outputs préférés
  a <- L_abeller(SimulationResultsTable,o,v)
  #On isole ensuite les outputs sous forme de tableaux différents, et on discrétise la valeur du noeud (double -> Qualitative)
  b <- lapply(t, Isolator, i = inp, df = a)
  #Et on plotte chaque output
  c <- lapply(b, Darkplotter)
  return(c)
}

###### PRINT BIG SCALE PART #######

ThePlotPrinter <- function(Megalist2plots,list2mutants,plotname)
{
  # Evaluates the number of outputs in the plotlists
  n_cols = length(Megalist2plots[[1]])
  if (n_cols <= 2)
  {
    pdf(plotname, width = (4.5 * n_cols), height = (1.5 * (n_cols)))
    mapply(FUN = PlotPrinter, list2plots = Megalist2plots, mutant = list2mutants, n_col = n_cols)
    dev.off()  
  } else
  {
    pdf(plotname, width = (4.5 * n_cols), height = (1.5 * (n_cols/2)))
    mapply(FUN = PlotPrinter, list2plots = Megalist2plots, mutant = list2mutants, n_col = n_cols)
    
    dev.off() 
  }
  
}

PlotPrinter <- function(list2plots,mutant, n_col)
{
  #plot_grid(plotlist = Might, ncol = 5, )
  print(do.call("grid.arrange", list(grobs = list2plots, ncol= n_col, left = textGrob(mutant, gp=gpar(fontsize=15,font=8), rot = 90))))
}


#### Dark side of isolator plotter (1 color)

#### Detection No Isc No Suf
Detection.No.IS <- function(dataframe,output1, output2)
{
  #  dataframe <- okk
  #   output1 <- "IscSUA"
  #    output2 <- "Suf"
  o1 <- dataframe[,grep(output1, colnames(dataframe))]
  o2 <- dataframe[,grep(output2, colnames(dataframe))]
  bool <- mapply(comparison,o1,o2)
  return(any(bool))
}

comparison <- function(a,b){return(a == 0 && b == 0)}

#But = obtenir une liste de booléens où y'a ni isc ni suf, et l'index permettra d'avoir la liste des mutants où t'as aucun des 2
MultiDetection.No.IS <- function(listDF, out1,out2)
{
  X <- lapply(X = listDF, FUN = Detection.No.IS, output1 = out1, output2 = out2)
  return(X)
}

get.IS.filenames <- function(listOfFiles, MutIS)
{
  x <- lapply(MutIS, grep, x = listOfFiles, value = T)
  return(x)
}

MegaIsolatorPlotter <- function(allSimulations,folders)
{
  Omega <- mapply(UltimateIsolatorPlotter, list2simulations = allSimulations, Folder = folders)
  return(Omega)
}

UltimateIsolatorPlotter <- function(list2simulations, Folder)
{
  AllPlots4mutants <- lapply(list2simulations, IsolatorPlotter, inp = INPUTS, o = OUTPUTS, t = THEOUTPUTS, v = VALUESMULTI)
  allmutNames <- getmutants(Folder)
  names(AllPlots4mutants) <- mutnames
  return(AllPlots4mutants)
}
