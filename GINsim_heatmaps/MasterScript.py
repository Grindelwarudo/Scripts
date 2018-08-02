#MasterScript
#1) Launch ExhaustMut in order to get all mutants
#2) For each perturbation file file generated,
#   make the simulations using the Attractors_tables.py script and GINsim!
#3) Use the IsolatorMapper.R script in order to get plenty of heatmaps.

##############################################
#                                            #
#               Import Part                  #
#                                            #
##############################################

import os
import subprocess

##############################################
#                                            #
#    Script part for making all simulations  #
#                                            #
##############################################


#Detects if there is one or more models. If more, pick one from the input.
def detectModel():
    #Vérifie qu'il y a un seul modèle. S'il y en a deux, choisir.
    model = []
    for e in os.listdir(cwd):
        if ".zginml" in e:
            model += [e]            
    if len(model) == 1:
        return model[0]
    elif len(model) > 1:
        i = int(input(displayNumbers(model)))
        return(model[i-1])
    else:
        print("Warning: No model found in the directory!")

def displayNumbers(l):
    i = 0
    for e in l:
        print("%i : %s" % (i+1 , e))
        i += 1
    return()

def detectGINsimVersion():
    GINsim = []
    for e in os.listdir(cwd):
        if "GINsim" in e:
            GINsim += [e]                    
    if len(GINsim) == 1:
        return GINsim[0]
    elif len(GINsim) > 1:
        i = int(input(displayNumbers(GINsim)))
        return(GINsim[i-1])
    else:
        print("Warning: No GINsim software found in the directory!")



def detectAllPerturbations():
    all_perturbations = []
    try:
        subprocess.check_output("python3.5 ExhaustMut.py", shell=True, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
            print ("error>",e.output,'<')
            print ("Whoops... Looks like the perturbation files have not been generated!")
    for e in os.listdir(cwd):
        if "Liste_mutants_" in e:
            all_perturbations += [e]
    if len(all_perturbations) == 0:
        print("No perturbation file has been found...")
    else:
        return(sorted(all_perturbations))


def perturbationsSimulator():
    model = detectModel()
    GINsim = detectGINsimVersion()
    
    all_perturbations = detectAllPerturbations()
    script = 'Attractors_tables.py'
    for e in all_perturbations:
        args = [GINsim, script , model, e]
        #cmd = ['java -jar ',GINsim, ' -s ', script, ' ', model, ' -p ', e]
        cmd = 'java -jar %s -s %s %s -p %s' % (GINsim, script ,model, e)
        print(cmd)
        try:
            subprocess.check_output(cmd, shell=True, stderr = subprocess.STDOUT)
        except subprocess.CalledProcessError as e: 
            print ("error>",e.output,'<') 
    print("All simulations are done.")





cwd = os.getcwd()
perturbationsSimulator()

###################################################


##
##def GraphMaker(): #Travaille sur la table des simulations...
##    current = os.getcwd()
##    args = ["Table_Des_Simulations.csv", current]
##    command = 'Rscript'
##    #Path2script = script if cwd is the file where the R script is located.
##    path2script = 'IsolatorPlotter-après.R'
##    cmd = [command, path2script] + args
##    #files = fileparser(".csv")
##    #print(files)
##    print(str(cmd))
##    try:
##        subprocess.check_output(cmd, universal_newlines=True, stderr = subprocess.STDOUT)
##        print("Everything is okay now. The graphs are operationnal... (ヘ･_･)ヘ┳━┳")
##    except subprocess.CalledProcessError as e: 
##       print ("error>",e.output,'<') 

#java -jar GINsim-2.9.3-with-deps.jar -s Attractors_tables.py OxyR_Sufmulti_NoROS.zginml -p Liste_mutants_1
