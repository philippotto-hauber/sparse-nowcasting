#!/bin/bash
module load matlab/2020b
mcc -m batch_eval_hpc.m

mv batch_eval_hpc ./../../eval-hpc
rm mccExcludedFiles.log requiredMCRProducts.txt readme.txt run_batch_eval_hpc.sh
