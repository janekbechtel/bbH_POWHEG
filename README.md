git clone https://github.com/janekbechtel/bbH_POWHEG .

cd bbH

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch

source $VO_CMS_SW_DIR/cmsset_default.sh

scram project CMSSW_8_0_22

cd CMSSW_8_0_22/src

cmsenv

cd -

make pwhg_main

