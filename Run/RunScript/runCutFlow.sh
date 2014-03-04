#!/bin/bash

cd ../../../

$ROOTCOREDIR/scripts/find_packages.sh
$ROOTCOREDIR/scripts/compile.sh

cd MyAnalysis/Run

root -l -b 'runAnalysisRoot.cxx'
cd RunScript
