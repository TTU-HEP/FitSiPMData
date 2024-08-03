# Instructions


## Setup code

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_14_0_13
cd CMSSW_14_0_13
cmsenv
```

```
cd $CMSSW_BASE/src
git clone git@github.com:TTU-HEP/FitSiPMData.git
cd FitSiPMData/test/
autoconf
./configure
make clean
make -j 4
```

```
cd $CMSSW_BASE/src/FitSiPMData/test/
cp /uscms_data/d3/cmadrid/DREAM/run742_V22000.root .
./mipFitsSiPM
```

