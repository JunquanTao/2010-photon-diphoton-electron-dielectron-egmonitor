######JTao#########
import itertools
import math
import ROOT

from CMSOpenDataAnalysis.DecaysToElectrons.Analyzer import Analyzer, Object

class SingleElectronAnalyzer(Analyzer):
    """
    XXX
    Some comment here.
    """
    def __init__(self):
        super(SingleElectronAnalyzer, self).__init__()

    #####CHANGE THIS METHOD TO CHANGE PHOTON ID######
    def electronID(self, electron):
        if electron.pt() < 20 or abs(electron.superCluster().position().Eta()) > 2.5:
            return False
         # excluded the EB/EE gap region
        if abs(electron.superCluster().position().Eta()) > 1.4442 and abs(electron.superCluster().position().Eta()) < 1.566:
            return False

        #Employ loose selections of PAS-EGM-10-004   http://cds.cern.ch/record/1299116
        #The electron identification variables that have been found to be the most powerful, and are used in the selection, are: the energy-momentum match between the seed cluster and the track E_seed/p_in, the variables measuring spatial matching between the track and the supercluster, DeltaEta_in and DeltaPhi_in, the supercluster eta_width, sigma_ietaieta (as taken from the covariance matrix using logarithmic weights), and the hadronic leakage variable H/E.
        #Isolation variables are computed in three sub-detectors: the tracker, the ECAL, and the HCAL. Transverse energy/momentum sums are evaluated in regions of DR < 0.3.
        # "WP85" used  (WP95 and WP80, Cut-in-category Cic Loose (95% eff) and Cic Loose (80% eff) from W/Z)
        # https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID

	#trackIso = electron.dr03TkSumPt()           # calculated with   electronTrackIsolationScone_cfi.py
	#ecalIso  = electron.dr03EcalRecHitSumEt()   # calculated with   electronEcalRecHitIsolationScone_cfi.py
	#hcalIso  = electron.dr03HcalTowerSumEt()    # calculated with   electronHcalTowerIsolationScone_cfi.py
           
	trackIsoRel = electron.dr03TkSumPt()/electron.pt()              
	ecalIsoRel  = electron.dr03EcalRecHitSumEt()/electron.pt()
	hcalIsoRel  = electron.dr03HcalTowerSumEt()/electron.pt()  

	#Combined Isolation as is defined with the following way:
	RelCombinedIsoEB = ( electron.dr03TkSumPt() + max(0., electron.dr03EcalRecHitSumEt() - 1.) + electron.dr03HcalTowerSumEt() ) / electron.pt() 
	RelCombinedIsoEE = ( electron.dr03TkSumPt() + electron.dr03EcalRecHitSumEt() + electron.dr03HcalTowerSumEt() ) / electron.pt() 
	#where the -1 is the pedestal subtraction and appears only in the barrel.

	#Shower shape:
	HoE = electron.hadronicOverEm()
	SigIeIe = electron.sigmaIetaIeta()
	#Track-cluster matching:
	DeltaPhi = electron.deltaPhiSuperClusterTrackAtVtx()
	DeltaEta = electron.deltaEtaSuperClusterTrackAtVtx()
	#Conversion rejection: Number of missing hits (the number of missing expected hits in front of the innermost valid hit) - If NumberOfExpectedInnerHits is greater than 1, then the electron is vetoed as from a converted photon and should be rejected in an analysis looking for prompt photons.   
	Nmisshit = electron.gsfTrack().trackerExpectedHitsInner().numberOfHits()
	#Minimum distance between conversion tracks:
	DisConv = electron.convDist()
	#$ \Delta cot \theta $ between conversion tracks at conversion vertex:
	DeltaCotTheta = electron.convDcot()

        if Nmisshit > 1:
	    return False
        if DisConv > 0.02:
            return False
        if DeltaCotTheta > 0.02:
            return False

        #EB      
        if abs(electron.superCluster().position().Eta()) < 1.4442:
            if RelCombinedIsoEB > 0.09:
                return False
            if trackIsoRel > 0.09:
               return False
	    if ecalIsoRel > 0.08:
               return False
	    if hcalIsoRel > 0.10:
              return False
            if SigIeIe > 0.01:
              return False
            if DeltaPhi > 0.06:
              return False
            if DeltaEta > 0.006:
              return False
            if HoE > 0.04:
              return False

        # EE:
        if abs(electron.superCluster().position().Eta()) > 1.566:
            if RelCombinedIsoEB > 0.06:
                return False
            if trackIsoRel > 0.05:
               return False
            if ecalIsoRel > 0.05:
               return False
            if hcalIsoRel > 0.025:
              return False
            if SigIeIe > 0.03:
              return False
            if DeltaPhi > 0.04:
              return False
            if DeltaEta > 0.007:
              return False
            if HoE > 0.025:
              return False

        return True


    #####ANALYSIS######
    def analyze(self, box):

        # Now check if there are at least 1 preselcted electron:
        if len(box.selectedElectrons) < 1:
            return False

        #first select electrons
        box.FinalSelectedElectrons = []
        for ele in box.selectedElectrons:
            try:
#                if self.electronID(ele, box.vertex):
                if self.electronID(ele):
                    box.FinalSelectedElectrons.append(ele)
            except ZeroDivisionError:
                continue

        if len(box.FinalSelectedElectrons) < 1:
            return False

        return True


    def declareHistos(self):
        super(SingleElectronAnalyzer, self).declareHistos()

        ###ADD YOUR HISTOGRAMS AFTER THIS LINE AS AbOVE#####
        self.declareHisto('electron_n', 10, -0.5, 9.5, "Electron multiplicity")
        self.declareHisto('electron_pt', 100, 15, 115, "p_{T,e} [GeV]")
        self.declareHisto('electron_eta', 52, -2.6, 2.6, "#eta_{e}")
        self.declareHisto('electron_SCeta', 52, -2.6, 2.6, "#eta of electron supercluster")
        ###EB
        #self.declareHisto('electron_R9_EB', 110, 0, 1.1, "e R_{9}")
        self.declareHisto('electron_SigIeIe_EB', 100, 0.005, 0.015, "Electron #sigma_{i#etai#eta}")
        self.declareHisto('electron_HoE_EB', 100, 0, 0.2, "Electron H/E")
        
        self.declareHisto('electron_RelCombinedIso03_EB', 100, 0, 0.1, "Relative combined ISO (#DeltaR=0.3)")
        self.declareHisto('electron_trackIsoRel03_EB', 100, 0, 0.1, "Relative track ISO (#DeltaR=0.3)")
        self.declareHisto('electron_ecalIsoRel03_EB', 100, 0, 0.1, "Relative ECAL ISO (#DeltaR=0.3)")
        self.declareHisto('electron_hcalIsoRel03_EB', 100, 0, 0.1, "Relative HCAL ISO (#DeltaR=0.3)")

        self.declareHisto('electron_DeltaEta_EB', 100, 0, 0.01, "Track-cluster matching : #Delta#eta")
        self.declareHisto('electron_DeltaPhi_EB', 100, 0, 0.1, "Track-cluster matching : #Delta#phi")
        
        self.declareHisto('electron_Nmisshit_EB', 5, -0.5, 4.5, "Number of missing expected hits")
        self.declareHisto('electron_DistanceConv_EB', 50, 0., 0.05, "Minimum distance between conversion tracks")
        self.declareHisto('electron_DeltaCotTheta_EB', 50, 0., 0.05, "#Deltacot#theta between conversion tracks at conversion vertex")

        ###EE
        #self.declareHisto('electron_R9_EE', 110, 0, 1.1, "e R_{9}")
        self.declareHisto('electron_SigIeIe_EE', 100, 0.005, 0.035, "Electron #sigma_{i#etai#eta}")
        self.declareHisto('electron_HoE_EE', 100, 0, 0.2, "Electron H/E")
        
        self.declareHisto('electron_RelCombinedIso03_EE', 100, 0, 0.1, "Relative combined ISO (#DeltaR=0.3)")
        self.declareHisto('electron_trackIsoRel03_EE', 100, 0, 0.1, "Relative track ISO (#DeltaR=0.3)")
        self.declareHisto('electron_ecalIsoRel03_EE', 100, 0, 0.1, "Relative ECAL ISO (#DeltaR=0.3)")
        self.declareHisto('electron_hcalIsoRel03_EE', 100, 0, 0.1, "Relative HCAL ISO (#DeltaR=0.3)")

        self.declareHisto('electron_DeltaEta_EE', 100, 0, 0.01, "Track-cluster matching : #Delta#eta")
        self.declareHisto('electron_DeltaPhi_EE', 100, 0, 0.1, "Track-cluster matching : #Delta#phi")
        
        self.declareHisto('electron_Nmisshit_EE', 5, -0.5, 4.5, "Number of missing expected hits")
        self.declareHisto('electron_DistanceConv_EE', 50, 0., 0.05, "Minimum distance between conversion tracks")
        self.declareHisto('electron_DeltaCotTheta_EE', 50, 0., 0.05, "#Deltacot#theta between conversion tracks at conversion vertex")
       #ES
        self.declareHisto('electron_EesOEsc_EE', 100, 0, 0.4, "Electron E_{es}/E_{SC}")


    def fillHistos(self, box, sample, weight=1):
        super(SingleElectronAnalyzer, self).fillHistos(box, sample, weight)
 
        self.fillHisto('electron_n', sample, len(box.FinalSelectedElectrons), weight)

        for electron in box.FinalSelectedElectrons:
            self.fillHisto('electron_pt', sample, electron.pt(), weight)
            self.fillHisto('electron_eta', sample, electron.eta(), weight)
            self.fillHisto('electron_SCeta', sample, electron.superCluster().position().Eta(), weight)

            trackIsoRel = electron.dr03TkSumPt()/electron.pt()
            ecalIsoRel  = electron.dr03EcalRecHitSumEt()/electron.pt()
            hcalIsoRel  = electron.dr03HcalTowerSumEt()/electron.pt()
            RelCombinedIsoEB = ( electron.dr03TkSumPt() + max(0., electron.dr03EcalRecHitSumEt() - 1.) + electron.dr03HcalTowerSumEt() ) / electron.pt()
            RelCombinedIsoEE = ( electron.dr03TkSumPt() + electron.dr03EcalRecHitSumEt() + electron.dr03HcalTowerSumEt() ) / electron.pt()
            DeltaPhi = electron.deltaPhiSuperClusterTrackAtVtx()
            DeltaEta = electron.deltaEtaSuperClusterTrackAtVtx()
            Nmisshit = electron.gsfTrack().trackerExpectedHitsInner().numberOfHits()
            DisConv = electron.convDist()
            DeltaCotTheta = electron.convDcot()
         
            ##EB
            if (abs(electron.superCluster().position().Eta())<1.45):
                #self.fillHisto('electron_R9_EB', sample, electron.superCluster().r9(), weight)
                self.fillHisto('electron_SigIeIe_EB', sample, electron.sigmaIetaIeta(), weight)
                self.fillHisto('electron_HoE_EB', sample, electron.hadronicOverEm(), weight)
                
                self.fillHisto('electron_RelCombinedIso03_EB', sample, RelCombinedIsoEB, weight)
                self.fillHisto('electron_trackIsoRel03_EB', sample, trackIsoRel, weight)
                self.fillHisto('electron_ecalIsoRel03_EB', sample, ecalIsoRel, weight)
                self.fillHisto('electron_hcalIsoRel03_EB', sample, hcalIsoRel, weight)
                
                self.fillHisto('electron_DeltaEta_EB', sample, DeltaEta, weight)
                self.fillHisto('electron_DeltaPhi_EB', sample, DeltaPhi, weight)
                
                self.fillHisto('electron_Nmisshit_EB', sample, Nmisshit, weight)
                self.fillHisto('electron_DistanceConv_EB', sample, DisConv, weight)
                self.fillHisto('electron_DeltaCotTheta_EB', sample, DeltaCotTheta, weight)
                #self.fillHisto('', sample, , weight)
            ##EE   
            if (abs(electron.superCluster().position().Eta())>1.56):
                #self.fillHisto('electron_R9_EE', sample, electron.superCluster().r9(), weight)
                self.fillHisto('electron_SigIeIe_EE', sample, electron.sigmaIetaIeta(), weight)
                self.fillHisto('electron_HoE_EE', sample, electron.hadronicOverEm(), weight)
                
                self.fillHisto('electron_RelCombinedIso03_EE', sample, RelCombinedIsoEE, weight)
                self.fillHisto('electron_trackIsoRel03_EE', sample, trackIsoRel, weight)
                self.fillHisto('electron_ecalIsoRel03_EE', sample, ecalIsoRel, weight)
                self.fillHisto('electron_hcalIsoRel03_EE', sample, hcalIsoRel, weight)
                
                self.fillHisto('electron_DeltaEta_EE', sample, DeltaEta, weight)
                self.fillHisto('electron_DeltaPhi_EE', sample, DeltaPhi, weight)
                
                self.fillHisto('electron_Nmisshit_EE', sample, Nmisshit, weight)
                self.fillHisto('electron_DistanceConv_EE', sample, DisConv, weight)
                self.fillHisto('electron_DeltaCotTheta_EE', sample, DeltaCotTheta, weight)

                self.fillHisto('electron_EesOEsc_EE', sample, electron.superCluster().preshowerEnergy()/electron.superCluster().rawEnergy(), weight)


    def addEvent(self, box):
        self.data.append(len(box.FinalSelectedElectrons))
