#!/bin/bash

cd ../../../

$ROOTCOREDIR/scripts/find_packages.sh
$ROOTCOREDIR/scripts/compile.sh

cd MyAnalysis/Run

make
for i in {1..16}
do
	./runAnalysis $i
done	   

cd Output
hadd -f mc12a_Mix_Combined.root mc12a_Mix_1.root mc12a_Mix_2.root mc12a_Mix_3.root
cd ..
cd RunScript
