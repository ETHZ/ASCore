#!/bin/bash

LOG="cmssw"
scram setup lhapdffull
touch $CMSSW_BASE/src/DiLeptonAnalysis/NTupleProducer/BuildFile.xml
eval `scram ru -sh`
scram b

cmsRun -j crab_fjr.xml ../../ntupleproducer_cfg.py runon=MC ModelScan=True