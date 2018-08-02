# GINsim_heatmap scripts

## Pre-requisites

This folder contians all the scripts used for analyzing a model with two inputs, with my GINsim model as an example.

Software requirements:

* python 3.5
* python 2.7
* R 3.4.4

Normally all R dependencies should be installed already (or R installs the prerequisite R packages otherwise).



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