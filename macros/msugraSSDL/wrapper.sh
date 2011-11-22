#!/bin/bash

source /swshare/ROOT/root_v5.30.00_40807_slc5_amd64//bin/thisroot.sh

cd /shome/mdunser/workspace/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/msugraSSDL

#root -l -b -q /shome/mdunser/workspace/CMSSW_4_2_8/src/DiLeptonAnalysis/NTupleProducer/macros/msugraScan/run_scan.c\(${1},${2}\)
root -l -b -q run_scan.c\(${1}\)
