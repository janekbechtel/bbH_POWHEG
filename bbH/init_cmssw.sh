#!/bin/bash

# Export CMSSW
export SCRAM_ARCH=slc6_amd64_gcc491
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scram project CMSSW CMSSW_8_0_22
cd CMSSW_8_0_22/src


# Build
scram b -j 24
scram b python
