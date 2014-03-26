#!/bin/bash

cd ../../../

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

voms-proxy-init -valid 48:00 --voms atlas:/atlas/phys-higgs/Role=production -pwstdin < ~/.pas

localSetupDQ2Client --skipConfirm 
#source /afs/cern.ch/atlas/offline/external/GRID/ddm/DQ2Clients/setup.sh
source /afs/cern.ch/atlas/offline/external/GRID/DA/panda-client/latest/etc/panda/panda_setup.sh 

cd MyAnalysis/Run/RunScript

#echo "Other Stuff Setup"
#source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh
#source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh
