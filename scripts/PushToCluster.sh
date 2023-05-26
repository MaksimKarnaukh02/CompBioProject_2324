#!/bin/sh
# This script will push the files to the cluster
rsync -av -r -P --exclude 'scripts' --exclude 'win_venv' --exclude '.idea' --exclude '.git' ./ vsc20886@login.hpc.uantwerpen.be:/scratch/antwerpen/208/vsc20886/CompBioProject_2324