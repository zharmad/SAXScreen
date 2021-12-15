#!/bin/bash

scrfold=$HOME/scripts/SAXScreen/scripts

# = = Generate the average apo-sample and buffer curves for fitting later.
$scrfold/step0B-generate-average-curves.bash

# = = Copy over template dictionary files. For the tutorial, all dictionary files have been prepared for you
# $scrfold/step0-make-simple-saxs-filelist.bash

# = = Run manual subtraction with optional corrections for additional buffer
$scrfold/step1-run-manual-subtraction.bash

# = = Dedicated script to run chi_lin calculations. Takes some time.
$scrfold/step2A-run-linear-fitting.bash

# = = Automation of AUTOGNOM
$scrfold/step2B-run-autognom-atsas.bash

# = = Dedicated script to generate similarity matrices for clustering
$scrfold/step3B-construct-matrix.bash

# = = Collate or compute structural quantities. Options include: VR, chi, Vc, Rg, PV, chi_free, ... 
$scrfold/step3-collate-measurements.bash chi
$scrfold/step3-collate-measurements.bash V_R

# = = Run final fitting for ranking and affinity estimation. Takes some time.
$scrfold/step4-fit-affinity-curves.bash chi
$scrfold/step4-fit-affinity-curves.bash V_R
