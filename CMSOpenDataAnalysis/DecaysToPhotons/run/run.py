###JTao test####
import ROOT
ROOT.gROOT.ProcessLine(".x hggPaperStyle.C")

from CMSOpenDataAnalysis.DecaysToPhotons.sources import sources

# Import the Analyzer you want to run:
# SinglePhotonAnalyzer or TwoPhotonAnalyzer
# by uncommenting the appropiate line below.
#from CMSOpenDataAnalysis.DecaysToPhotons.TwoPhotonAnalyzer import TwoPhotonAnalyzer as MyAnalyzer
from CMSOpenDataAnalysis.DecaysToPhotons.SinglePhotonAnalyzer import SinglePhotonAnalyzer as MyAnalyzer

analyzer = MyAnalyzer()

analyzer.declareHistos()

for sample in sources:
    # maxEv defines the maximum number of events to analyze
    # set it to -1 to analyze all available events; 
    analyzer.processSample(sample, maxEv=-1)


analyzer.makeAllPlots()

# uncommet line below to export selected data to a json file
#analyzer.exportData()
