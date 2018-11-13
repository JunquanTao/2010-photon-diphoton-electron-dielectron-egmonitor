# CMS OpenData Analysis with photon(s)/electron(s) in the final status

A brief introduction for this repo which contains different examples:

a) Firstly two simple codes for the photon/diphoton analysis ("DecaysToPhotons") and electron/dielectron analysis ("DecaysToElectrons"). You can run the examples according to the README in each work directory. These codes were developed for quick check on the validation plots of photon/diphoton/electron/dielectron objects. Only histograms are stored and you can run the codes with Python based command. You can run the codes locally but cann't submit the CRAB jobs, so much suitable for the check of small part of the datasets.

b) Then new codes are developed to validate 2010B datasets : Photon, Electron and EGMonitor. The codes "PhotonElectronAnalyzer" include Standard EDAnalyzer and store five minitrees for general event information ("eventTree"), single photon ("photonTree"), diphoton ("diphotonTree"), single electron ("electronTree") and di-electron ("dielectronTree") after some selections. You can run "cmsRun PhotonElectronAnalyzer.py" in the sub-directory test/ to test the codes, "crabConfig_\*.py" in test/ directory for crab job submission, and some ROOT and C++ based scripts ("Draw\*.C") in test to draw the plots, based on the five minitrees.


The following recipe are based on case b) with "PhotonElectronAnalyzer".

From now on it is assumed that you will work on a VM properly contextualized for CMS, the current version "CMS-OpenData-1.1.2" for the validation of 2010B dataset. You can use SLC5 shell for code compilation and test. You need to change to SLC7 for CRAB job submission. Please follow the recipe https://twiki.cern.ch/twiki/bin/viewauth/CMS/DPOALegacyDataJobSubmission .


## Creating the Working Area

This step is only needed the first time.

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_4_2_8
cd CMSSW_4_2_8/src
git init
git clone git://github.com/JunquanTao/CMSOpenDataAnalysis.git
scram b
```

## Sourcing the environment 

This step is needed each time you want to run the exercise.

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSW_4_2_8/src
cmsenv
```

## Compile the "PhotonElectronAnalyzer" codes with SLC5 as mentioned above

```
cd  PhotonElectronAnalyzer/
scramv1 b
``` 

## Running the example in the test to check if the codes can be ran sucessfully

```
cd CMSOpenDataAnalysis/PhotonElectronAnalyzer/test/
cmsRun  PhotonElectronAnalyzer.py
```
After the running of this example, you will get an output root file, PhotonElectron_Test.root, containing the five minitrees as mention above.


## Running CRAB jobs to analysis the full dataset found in DAS, with SLC7 as mentioned above

You need to source the CRAB env with SLC7 firstly:

```
cd CMSSW_4_2_8/src/
source /cvmfs/grid.cern.ch/umd-c7ui-latest/etc/profile.d/setup-c7-ui-example.sh
eval `scramv1 runtime -sh`
source /cvmfs/cms.cern.ch/crab3/crab_standalone.sh
```

You can use the "crabConfig_\*.py" in the test directory to create the CRAB job work directory for different dataset, "crabConfig_Photon.py" for 2010B Photon dataset, "crabConfig_Electron.py" for 2010B Electron dataset and "crabConfig_Egamma.py" for 2010B EGMonitor dataset. By dafault, in these cfg file, the output will be stored to my CERN T2 area, you need to change the followng two lines according to the T2 site you can access:

```
config.Data.outLFNDirBase = "/store/user/jtao/CMSOpenData2010/"
config.Site.storageSite = 'T2_CH_CERN'  # you might need to change this to a site you have acces too
```

Please follow the CRAB commands for Grid job submission, 
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial .


## Running the ROOT and C++ based scripts ("Draw\*.C") in test/ to draw the plots, based on the five minitrees.

Please modify the codes to take the outputs (rootuples) of your CRAB jobs as the input root files of these scripts. Then simply using for example 
```
cd PhotonElectronAnalyzer/test/
root -l -b -q DrawEvtTreeVar.C 
```
to obtain the final validation plots.



Please enjoy!
