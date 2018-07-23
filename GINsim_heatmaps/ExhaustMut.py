import os
import itertools
from itertools import combinations
import re
import argparse

#Put the boolean nodes and multivalued nodes here
BoolNodes = ["H2O2","HO","OxyR","Hpx","IscSUA","IscR_A","IscR_H"]
MultiNodes = ["Fur","RyhB","Suf"]
NoIscRpls = ["IscR_H","Suf","Fur","RyhB"]

dicodor = {"Boolean": BoolNodes, "Multivalued": MultiNodes, "IscRHolo" : NoIscRpls}

mutNature = ["KO"]
#mutNature = ["KO","E1"]
boolNature = mutNature[:-1]
MultiNature = mutNature
MultiNatureX = ["KO"]
#MultiNatureX = ["KO",E2"]
def GetNodes():
    cwd = os.getcwd()
    nf = getnodefile(cwd)
    l = list2mut(nf)
    print(len(l))
    print(pow(2,len(l)))
    #print(l)
    l2 = listv2(l,dicodor)
    #print(l2)
    print(len(l2))
    print(pow(2,len(l2)))
    print(pow(2,len(l2)) / 12)
    #powa = list(powerset(l2))
    
    #power = powa[1:]
    #Embodiement = [list(row) for row in power]
    
    #print(Embodiement[1:25])
    Virtue = nest2notnest(Embodiement)
    #print(Virtue[1:25])
    #enumerer("Exhaustive_combinations", Virtue)
    return("Ok")


def GetNodes2(): #Fonction générant des fichiers 1 par 1, de combinaisons de n noeuds.
    cwd = os.getcwd()
    nf = getnodefile(cwd)
    l = list2mut(nf)
    #Only if we want ectopic mutants
    #l2 = listv2(l, dicodor)
    #print(l2)
    #print(len(l2))
    #combmaster(l2)
    #Only if we want KO mutants solo
    l3 = listv2_simple(l)
    combmaster(l3)
    return("done")

def combmaster(l):
    i = 1
    while i < len(l): #TO define manually, 5 is best for getting all the triple/quadruple mutants 
        c = []
        c.extend(combinations(l,i))
        print(len(c))
        if len(c) != 0:
            #print(c)
            filename = "Liste_mutants_%i" % i
            comb2 = [list(row) for row in c]
            comb3 = nest2notnest(comb2)
            #print(comb3)
            print("Printing the file: " + filename)
            enumerer(filename,comb3)
        i += 1
    return("Done")


def listv2_simple(listmutant): #Takes a list of mutants, and places the KO label (Knock-Out)
    res = []
    for e in listmutant:
        res = res + [e + " KO"]
    return res





#Takes into account the KO or E1 nature
def listv2(listmutant,d): #Prend une list de mutants, et greffe les natures des mutations dessus
    res = []
    for e in listmutant:
        if e in d["Boolean"] and not e in d["IscRHolo"]:
            a = e + " " + boolNature[0]
            b = e + " " + boolNature[1]
            res = res + [a] + [b]
        elif e in d["Boolean"] and e in d["IscRHolo"]:
            a = e + " " + boolNature[0]
            res = res + [a]
        elif e in d["Multivalued"] and not e in d["IscRHolo"]:
            a = e + " " + MultiNature[0]
            b = e + " " + MultiNature[1]
            c = e + " " + MultiNature[2]
            res = res + [a] + [b] + [c]
        elif e in d["Multivalued"] and e in d["IscRHolo"]:
            a = e + " " + MultiNatureX[0]
            b = e + " " + MultiNatureX[1]
            res = res + [a] + [b]
        else:
            return False
    return res


def getnodefile(cwd):
    for e in os.listdir(cwd):
        if e == 'Nodes':
            return 'Nodes'
    print('non')
    return(False)

def list2mut(file):
    a = open(file,"r")
    res = []
    line = a.readline()
    while line != "":
        line2 = line.strip()
        res = res + [line2]
        line = a.readline()
    a.close()
    return res[:-1]

def enumerer(file, l):
    a = open(file,"w")
    for e in l:
        a.write(e + "\n")
    a.close()
    return("Done")

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))


def it2list(iterable):
    return(list(iterable))

def nest2notnest(nes): #Transforme une nested list en une liste simple.
    res = []
    for e in nes:
        if len(e) == 1:
            res = res + e
        elif len(e) > 1:
            res = res + [','.join(e)] #concatene la liste en une string séparant les ets par une virgule
        else:
            return(False)
    return res

######### MAIN ##########

GetNodes2()



