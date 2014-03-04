#!/bin/bash

echo "Going into the folder"

inputSTR=$1

echo "Input Container Name: $inputSTR"

if [ "$1" == "" ];
  then
    echo "No Input argument"
	inputSTR="MC12"
fi

useNewGeo=false

if [ "$2" != "" ];
then
	echo "Use this new geo"
	useNewGeo=true
fi

cd ../ 

echo $inputSTR > inputName.txt
root -l -q -b 'runAnalysisGrid.cxx('$useNewGeo')'

cd RunScript
