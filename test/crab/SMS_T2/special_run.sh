#!/bin/bash

LOG="cmssw"
scram setup lhapdffull
eval `scram ru -sh`
scram b

cmsRun -j  $RUNTIME_AREA/crab_fjr_$NJob.xml pset.py
