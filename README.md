# SAXScreen
A set of scripts in Python and bash to analyse SAXS curves and extract binding patterns

# Introduction

Overall Procedure for SAXScreen

Terminology:
- Raw I(q)-file is the 1d curve directly output from the synchrotron software.
- Buffer-subtracted I(q) is the initial processed one after correction for combinations between sample and ligand buffer.
- Smoothed I(q) is further processed in ATSAS.

Requirements:
- Recent Python with scipy and numpy. On EMBL machines, install a python virtal-environment and use it with this program.

# Step-wise Instructions

## 0) Prepare files and scripts.

You need to sort and figure out four things:
- List of titration curves and their identities (titration_dictionary)
    - This file contains four columns
    - The name of the titration series (e.g. name of the ligand) to identify and track all data.
    - The sample-ligand ratio indicating the titration point.
    - The volume-wise fraction of ligand aliquot added to the sample during measuement (e.g. 20uL+20uL = 0.5 volume fraction)
    - The location of the raw I(q) file corresponding to this titration point.

- List of ligand curves and their identities    (ligand_dictionary)
    - This file contains two columns. The name of the titration series linked to the ligand SAXS-curves
- List of Sample buffer                         (sample_buffer_list)
    - This file is a simple list of buffer SAXS curves.
- List of Ligand buffer                         (ligand_buffer_list)
    - This file is a simple list of buffer SAXS curves.

Examine the example files for a better sense of what each file is meant to contain.
We will be calling these dictionary files.
The script step0-make-simple-saxs-filelist.bash is here to help
automatically generate the file listing based on the given title/code of the synchrotron output data.

## 1) Then run the scripts in order. Once every parameter is set, then they should all work OK.

The only file to copy over are general-settings.txt, and optionally all the step\*.bash scripts.
- general-settings.txt attempts to keep all the variable that you're expected to modify in a single place.
- bash scripts are the top-level interface to shell scripts and python commands.

