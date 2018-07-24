###JTao test####
import ROOT
ROOT.gROOT.ProcessLine(".x ./hggPaperStyle.C")

from CMSOpenDataAnalysis.DecaysToElectrons.sources import sources

# Import the Analyzer you want to run:
# SingleElectronAnalyzer or TwoElectronAnalyzer
# by uncommenting the appropiate line below.
from CMSOpenDataAnalysis.DecaysToElectrons.TwoElectronAnalyzer import TwoElectronAnalyzer as MyAnalyzer
#from CMSOpenDataAnalysis.DecaysToElectrons.SingleElectronAnalyzer import SingleElectronAnalyzer as MyAnalyzer

analyzer = MyAnalyzer()

analyzer.declareHistos()

for sample in sources:
    # maxEv defines the maximum number of events to analyze
    # set it to -1 to analyze all available events; 
    analyzer.processSample(sample, maxEv=-1)
    #analyzer.processSample(sample, maxEv=100000)


analyzer.makeAllPlots()

# uncommet line below to export selected data to a json file
#analyzer.exportData()
