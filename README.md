# GINsim_heatmap scripts

## Pre-requisites

This folder contians all the scripts used for analyzing a model with two inputs, with my GINsim model as an example.

### Hardware requirements:

* 8 Gb RAM
* 4 cores CPU (I could launch my analysis on my laptop, which has the "Intel® Core™ i5-4200H CPU @ 2.80GHz × 4" )

### Software requirements:

* python 3.5
* python 2.7
* R 3.4.4

Normally all R dependencies should be installed already (or R installs the prerequisite R packages otherwise).

### Files required

In order to use the scripts, make sure that:

* the model file is not configured to have predefined inputs
* there are 7 .txt files:
    + Model.txt
    + Boolean.txt
    + Multivalued.txt
    + TheInputs.txt
    + Readouts.txt
    + TheOutputs.txt
    + IS_outputs.txt
* The "Node" text file is present.

## How to use it:

```
# Make all the GINsim mutant perturbations
python3.5 Masterscript.py

# Analyse all the data made
Rscript IsolatorMapper.R -n Model.txt -i TheInputs.txt -b Boolean.txt -M Multivalued.txt -r Readouts.txt -o TheOutputs.txt -ISo IS_Outputs.txt

# -n = all model nodes
# -i = model inputs
# -b = boolean nodes
# -M = multivalued nodes
# -r = model readouts
# -o = model outputs (which will be printed in pdf files)
# -ISo = special model (for mutants where neither the outputs are present).
```

## Troubleshootings (and TO DO LIST)

> Oh noes the R script is too slow how do i make it faster

It is possible to launch it on a cluster, however the server stuff is not implemented yet. (TO DO)

In the meantime, it is possible to edit the IsolatorHelper.R script at the mc.cores line, and choose the number of cores used. (I should make an option for this...)


# MaBoSS analysis scripts

## Pre-requisites

Software:

* Python3.5
* R 3.4.4
* MaBoSS 2.0 and all its dependencies (please refer to [the MaBoSS installation page](https://maboss.curie.fr/)
* MaBoSS environment

For the hardware, please refer to the hardware requirements above.

Files:

* Booleanized MaBoSS model:

In order to get the booleanized MaBoSS model, I first exported my model in the SBML format. Then, I used bioLQM with the following command:

```
java -jar bioLQMlatestversion.jar SBMLfile -m booleanize ginmlfile
```

And then exported the ginml model to the MaBoSS format. 

NB: Make sure that no node name is included in the model filename, especially in the MaBoSS format!!

> I used this way in order to control the simulation inputs one by one. In the future I should be able to do the conversion in one step. - Firas

* The files containing MaBoSS generated data

To generate data with MaBoSS, please make sure that the model files, the SimulatorWT and Mutations.txt are in the same directory.

```
# Beforehand, activate the MaBoSS environment
source MaBoSS.env

#Then, prepare the model for the inputs. In the multivalued case, it should look like this:
MBSS_MutBndCfg Model.bnd Model.cfg "Input1_b1 Input1_b2 Input2_b1 Input2_b2"
# Note that it is the case for my model: MBSS_MutBndCfg Model.bnd Model.cfg "Fe_ext_b1 Fe_ext_b2 O2_b1 O2_b2" 
# You should obtain the same files with a "_mut" in the filename.

#Then, launch the simulations in the WT context and the mutants! It is necessary to have both for the next steps to work:

MBSS_MultipleSim.py Model_mut.cfg SimulatorWT

# You will have some stuff to prompt

#the first prompt input: 0
#the second prompt input (which situations you want to simulate): 1,2,3,4,5,6,7,8,9

# After making the simuations in the WT context, let's do the mutants:
MBSS_MultipleSim.py Model_mut.cfg SimulatorWT Mutation.txt

# Same prompt as above.

```

## How to use the analysis script:

```
# The analysis is made with the following cmd line, generating pdfs for each mutant. Each plot correspond to an output containing each condition:

python3.5 MaBoSS_Masterscript.py

```

## Troubleshooting (and TO DO LIST)

> Oh noes its too slow!

It is possible to replace all the "lapply" with mclapply, with mc.cores as the number of cores you need. Since I was limited by the RAM, I could not make some parallelization...

> I have put too much mutants, and I crashed...

This script uses plenty of memory since it makes a huge structure to store all the data coming from MaBoSS. The more mutants the bigger. I already optimized it using the SparseMatrix data type in R (it was even bigger before...). As for the GINsim scripts I should make an option for using it in a server...

