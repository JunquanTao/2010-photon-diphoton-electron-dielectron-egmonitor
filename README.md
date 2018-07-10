# CMS OpenData Analysis with photon(s) in the final status

A brief introduction for this repo which contains two examples for:

a) the single photon analysis with the producer "SinglePhotonAnalyzer"

b) the double photon analysis with the producer "TwoPhotonAnalyzer"

Both examples exercises running on Run2010B photon data in the AOD format, http://opendata.cern.ch/record/12.

From now on it is assumed that you will work on a VM properly contextualized for CMS.

## Creating the Working Area

This step is only needed the first time.

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_4_2_8
cd CMSSW_4_2_8/src
git init
git clone git://github.com/JunquanTao/CMSOpenDataAnalysis.git
```

## Sourcing the environment 

This step is needed each time you want to run the exercise.

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSW_4_2_8/src
cmsenv
```

## Running the Exercise

The path to the pattuples is giving in the file:
```
CMSOpenDataAnalysis/DecaysToPhotons/python/sources.py
``` 

To run the code you have to move to the directory:

```
CMSOpenDataAnalysis/DecaysToPhotons/run/
```

You must specify in your run.py code from the path above which analysis you want to run:

```python
# Import the Analyzer you want to run:
# SinglePhotonAnalyzer or TwoPhotonAnalyzer
# by uncommenting the appropiate line below.
#from CMSOpenDataAnalysis.DecaysToPhotons.TwoPhotonAnalyzer import TwoPhotonAnalyzer as MyAnalyzer
from CMSOpenDataAnalysis.DecaysToPhotons.SinglePhotonAnalyzer import SinglePhotonAnalyzer as MyAnalyzer
``` 

The number of events to be analyzed can be modified in the run.py file.

```python
for sample in sources:
    # maxEv defines the maximum number of events to analyze
    # set it to -1 to analyze all available events; 
    analyzer.processSample(sample, maxEv=100)
```

To get enough events in the plots, you would need to run over all available samples. That takes time.

```
python run.py 
```

or in a non-interactive mode:

```
python run.py 
```
At the beginning you will get a message like: 

```python
Processing Files
['root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/004DDBA5-7471-E011-A381-0017A4770C08.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/00C1B689-F670-E011-B054-1CC1DE1CF1BA.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/063C79CD-3271-E011-934D-0025B3E022C2.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/065409FD-0071-E011-B185-001F296B758E.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/06C1851C-8E71-E011-83BB-00237DA16C42.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/08118EBE-9B71-E011-915E-001F296B758E.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0A6EBDBB-8371-E011-A8EE-0017A477003C.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0C59D16B-F670-E011-A3D4-1CC1DE1CDDBC.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0C9AE5C4-3371-E011-85B7-1CC1DE1CEDB2.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0CC571FD-7571-E011-98CB-1CC1DE046F78.root', 'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Photon/AOD/Apr21ReReco-v1/0005/0E1ED66B-3971-E011-9C33-00237DA1EDE0.root']
90090 events available for processing
1 events processed in 0.000134944915771 seconds
90090 events processed in 80.6395988464 seconds
```

By default all the histrograms defined in your analyzer will be plotted and saved in png format (with the name of the histogram). In the interactive mode, you can also plot them typing (depending on the analysis your are running):

```python
analyzer.makeAllPlots()
```
to obtain all the plots by default from the codes or for example 
```
analyzer.makePlot("diphoton_mass")
```
from the TwoPhotonAnalyzer to get only the mass distribution.

You can exit the python session by typing exit() or ctrl+d.

Events selection can be modified in the SinglePhotonAnalyzer.py and TwoPhotonAnalyzer.py codes:
```
CMSOpenDataAnalysis/DecaysToPhotons/python/SinglePhotonAnalyzer.py
CMSOpenDataAnalysis/DecaysToPhotons/python/TwoPhotonAnalyzer.py
```
For single photon, the loose photon id selections in PAS-EGM-10-006 (http://cdsweb.cern.ch/record/1324545) together with photon pT > 21 GeV are used. No electron-veto requirement is asked by default in the codes. The output plots include the photon pT and eta, and some variables related to shower shape and photon isolation (11 variables in EB and 12 variables in EE at present) used for photon identification.

For diphoton, I applied the event selections based on the cross section measurement with 2010 data (JHEP01(2012)133) but more simple and direct way. No electron-veto selection is required. The output plots include diphoton mass and pT, delta_Phi and cos_theta_star between two photons, which were used in the cross section measurement, and also additional the pt and eta of leading and subleading photons, with 8 variables in total. From the mass plot, you can see the Z peak and maybe also the Higgs :)



