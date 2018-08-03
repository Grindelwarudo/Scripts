import os
import re
import csv
import pandas
import subprocess

def GetDossiers(): #Recupere les dossiers
    
    cwd = os.getcwd()
    thelist = dirchecker(os.listdir(cwd))
    return(thelist)


def GetNodes(): #Recupere les noeuds du modèle, nécessite la spécification d'un fichier .cfg
    
    a = open(getCfg(),"r")
    line = a.readline()
    nodes = []
    while line != "":
        if "istate" in line:
            nodes = nodes + [line]
        line = a.readline()
    a.close()
    nodes = purgeuselessthings(nodes)
    return nodes


def getCfg(): #Recupere le fichier .cfg et s'il y en a plusieurs, renvoie une erreur
    a = []
    for file in os.listdir(os.getcwd()):
        if file.endswith(".cfg"):
            a = a + [file]
    if len(a) > 1 or len(a) == 0:
        sys.exit("Only one cfg file can be chosen. Pick it wisely...")
    else:
        print(a[0])
        return a[0]


def dirchecker(l): #Purge la liste contenant tous les fichiers des fichiers, ne laissant que des répertoires
    res = []
    for e in l:
        if os.path.isdir(e) == True:
            res = res + [e]
    return res




    
def lecteur2list(l): #Lit une liste (diagnostic)
    for e in l:
        print(e)
    return "ok"


def purgeuselessthings(l): #purge la liste des ".istate" j'saispasquoi.
    res = []
    for e in l:
        purgednode = re.sub('.istate=0;\n','',e)
        res = res + [purgednode]
    return res

def mutchecker(s,l1,l2):
    #Verifie les dossiers où le mutant est présent.
    #print("Nom du mutant à voir: ",s)
    #s = mutant à voir, l1 = liste de mutants, l2 = liste de fichiers
    res = []
    splitted = s.split(",")
    res = res + parkour(splitted, l2)
    return res

def getAllDico(d): #Recupere toutes les valeurs du dico.
    l = d.values()
    res = []
    for e in l:
        if e != []:
            res = res + e
    #print(len(res))
    return res

def parkour(splitted,l2): #Permet de parcourir tous les éléments du mutant
    #& de définir l'intersection entre toutes ces listes
    #print(len(l2))
    res = []
    helper = []
    i = 0
    while i < len(splitted):
        helper.append(parkour2(splitted[i],l2))
        #print(splitted[i], "\n")
        i = i + 1
    #print(helper[0])
    res = list(set(helper[0]).intersection(*helper))
    #print(set(helper[0]).intersection(*helper))
    #print(res)
    
    return res

def parkour2(et,l2):
    res = []
    for e in l2:
        if et in e:
            res = res + [e]
    return res


#Diag functions

def cptTailldico(d): #Compte les tailles de chanque élément du dico
    for e in d:
        print(e + ": %i" % len(d[e]))


#l2 = list2noeuds
def getdico(l1,l2,l3):
    dico = {}
    for e in l1:
        dico[e] = mutchecker(e,l1,l2)
    #cptTailldico(dico)
    return purj(dico,l2,l3)


def Wtchecker(d,l2): #Recupere les dossiers correspondant au WT
    #l2 = liste de noeuds...
    allV = getAllDico(d)
    WT = WTgetter(allV,l2)
    return WT

def WTgetter(l1,l2): #Recupere la difference entre 2 listes.
    #print(list(set(l1).symmetric_difference(set(l2))))
    return list(set(l1).symmetric_difference(set(l2)))

def MutationsGetter(): #Recupere tous les mutants, pas de discrimination de up/dwn
    a = open("Mutation.txt","r")
    line = a.readline()
    mutants = []
    while line != "":
        mutants = mutants + [line]
        line = a.readline()
    a.close()
    #print(mutants)
    mutants = purger(mutants)
    #print(mutants)
    mutants.sort(key = len, reverse = True)
    return mutants

def purger(liste): #purge la ligne spécifique du mutation.txt.
    res = []
    sub = {"u:":"u_","d:":"d_"}
    resi = [re.sub("u:",sub["u:"],e) for e in liste]
    res = [re.sub("d:",sub["d:"],e) for e in resi]
        
    res = [re.sub("\n","",f) for f in res]
    #print(res)
    return res #Retourne le résultat


def purj(dico, list2mutants,list2noeuds):
    dico2 = {}
    for e in dico:
        print(e)
        print(dico[e])
        dico2[e] = MutantVerif(dico[e],e,list2noeuds)
    dico2["WT"] = Wtchecker(dico,list2mutants)
    return dico2

def MutantVerif(dossiers,mutantetudie,list2noeuds):
#Verification du nombre de noeuds mutés. Si c'est identique, on vérifie que ce soit le bon
    #sinon... Tri.
    maxi = len(mutantetudie.split(","))
    print("\nNb max de noeuds du mutant: ", maxi)
    truelist = []
    for e in dossiers:
        i = checklength(e,list2noeuds) #verifie la taille des mutants
        print("Nb de noeuds du dossier: ", i)
        if i == maxi and checknames(e,mutantetudie,list2noeuds):
            truelist = truelist + [e]
    #print(truelist)
    #print(len(truelist))
    return truelist


            
def checklength(dossier,list2noeuds):
    res = 0
    for e in list2noeuds:
        if e in dossier:
            res = res + 1
    return res


def checknames(dossier,mutantetudie,list2noeuds):
    mutantNodes = mutantetudie.split(",")
    boo = [] #Accumulateur de mutant présent ou pas dans le dossier
    for e in mutantNodes:
        if e in dossier:
            boo = boo + [True]
        else:
            boo = boo + [False]
    return(all(boo))

def listmaker(d):
    res = []
    helper = []
    for e in d:
        helper = [e] + d[e]
        res = res + [helper]
        helper = []
    return res

def csvwriter(liste):
    pd = pandas.DataFrame(liste)
    pd.to_csv("Table_Des_Simulations.csv")
    return "Fichier crée, maintenant go sur R frère"


def CoreProgram():
    l1 = MutationsGetter()
    print(l1)
	
    l2 = GetDossiers()
    l3 = GetNodes()
    print("On créé le dictionnaire avec tous les dossiers correspondant aux mutants...")
    #C'est là qu'on assigne les mutants & les dossiers
    dico = getdico(l1,l2,l3)
    gigalist = listmaker(dico)
    cptTailldico(dico)
    print("... et on fout tout ça dans un fichier csv :)")
    csvwriter(gigalist)
    GraphMaker()




def GraphMaker(): #Travaille sur la table des simulations...
    current = os.getcwd()
    args = ["Table_Des_Simulations.csv", current]
    command = 'Rscript'
    #Path2script = script if cwd is the file where the R script is located.
    path2script = 'IsolatorPlotter-après.R'
    cmd = [command, path2script] + args
    #files = fileparser(".csv")
    #print(files)
    print(str(cmd))
    try:
        subprocess.check_output(cmd, universal_newlines=True, stderr = subprocess.STDOUT)
        print("Everything is okay now. The graphs are operationnal... (ヘ･_･)ヘ┳━┳")
    except subprocess.CalledProcessError as e: 
       print ("error>",e.output,'<') 
     
    


CoreProgram()


#l1 = MutationsGetter()
#l2 = GetDossiers()
#l3 = GetNodes()
#o = getdico(l1,l2,l3)
#test = "d_Fur_b1,d_Fur_b2,d_RyhB_b1,d_RyhB_b2"





#p = purj(o,l2,l3)
#cptTailldico(p)

#print(l1)

#dicoo = getdico2(l1,l2,l3)
#mut = ["d_Cat"]
#dico_keys = dicoo.keys()

#lecteur2list(getNamesInMutant(mut,dico_keys))

#Valeurs = dicoo["Cat"]
#splitted = ["Cat"]
#Ebola = WTgetter(splitted, l3)

#spaghetti = GetWTdir(l2,l1)
#print(spaghetti)
#print(len(spaghetti))
#print(Ebola)
##mut = ["d_Fur_b1"]
##doss = "Ox_R2_Fe-S-bool_mut_Fe_0-O2_0_d_Fur_b1"
##print(mutdsledoss(doss,mut))
