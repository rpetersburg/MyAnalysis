#!/bin/bash
cd ../../../

$ROOTCOREDIR/scripts/find_packages.sh
$ROOTCOREDIR/scripts/compile.sh

cd MyAnalysis/Run

rm Output.txt
root -l 'runAnalysis.cxx' >> Output.txt

cd RunScript
