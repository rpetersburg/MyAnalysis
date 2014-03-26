#!/bin/bash

cd ../../../

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
echo "Local Atlas Setup"
setupATLAS
echo "Athena Setup"
asetup 17.7.0,here,64,slc5

cd MyAnalysis/Run/RunScript

#echo "Other Stuff Setup"
#source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh
#source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh
