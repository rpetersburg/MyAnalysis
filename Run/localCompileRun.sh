#!/bin/bash
cd ../cmt
make -f Makefile.RootCore
cd ../Run
root -l -q -b runAnalysisRoot.cxx
