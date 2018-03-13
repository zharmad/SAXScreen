#!/bin/bash

scrfold="../scripts"

# = = Generate the average apo-sample and buffer curves for fitting later.
$scrfold/step0-generate-average-curves.bash

# = = Copy over template dictionary files
# $scrfold/step0-make-simple-saxs-filelist.bash

# = = Run manual subtraction with optional corrections for additional buffer
$scrfold/step1-run-manual-subtraction.bash

# = = Dedicated script to run chi_lin calculations. Takes some time.
$scrfold/step2A-run-linear-fitting.bash

# = = Automation of AUTOGNOM
$scrfold/step2B-run-autognom-atsas.bash

# = = Dedicated script to generate similarity matrices for clustering
$scrfold/step3B-construct-chimatrix.bash

for quant in chi Vc PV Rg ; do
    # = = Collate or compute structural quantities
    $scrfold/step3-collate-measurements.bash $quant

    # = = Run final fitting for ranking and affinity estimation. Takes some time.
    $scrfold/step4-fit-curves.bash $quant
done
