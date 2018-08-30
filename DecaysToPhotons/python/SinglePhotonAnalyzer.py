######JTao#########
import itertools
import math
import ROOT

from CMSOpenDataAnalysis.DecaysToPhotons.Analyzer import Analyzer, Object

class SinglePhotonAnalyzer(Analyzer):
    """
    XXX
    Some comment here.
    """
    def __init__(self):
        super(SinglePhotonAnalyzer, self).__init__()

    #####CHANGE THIS METHOD TO CHANGE PHOTON ID######
    def photonID(self, photon):
        if photon.pt() < 21 or abs(photon.eta()) > 2.5:
            return False
         # excluded the EB/EE gap region
        if abs(photon.eta()) > 1.4442 and abs(photon.eta()) < 1.566:
            return False

        #Employ loose selections of PAS-EGM-10-006 http://cdsweb.cern.ch/record/1324545
        #sigmaIetaIeta < 0.01/0.03 in EB/EE 
        #HoE (hadronicOverEm) < 0.05 
        #Track ISO (cone DeltaR=0.4) trkSumPtHollowConeDR04 < 2.0GeV for both EB/EE
        #ECAL ISO (cone DeltaR=0.4) ecalRecHitSumEtConeDR03 < 4.2GeV for both EB/EE
        #HCAL ISO (cone DeltaR=0.4) hcalTowerSumEtConeDR04 < 2.2GeV  for both EB/EE

        #common selections for both EB and EE
        if photon.hadronicOverEm() > 0.05:
            return False
        if photon.trkSumPtHollowConeDR04() > 2.:
            return False
        if photon.ecalRecHitSumEtConeDR04() > 4.2:
            return False
        if photon.hcalTowerSumEtConeDR04() > 2.2:
            return False
        # pixel-seed e-veto 
        if photon.hasPixelSeed():
            return False

        # EB: 
        #if abs(photon.eta()) < 1.4442:
        if abs(photon.superCluster().position().Eta()) < 1.4442:
            if photon.sigmaIetaIeta() > 0.01:
                return False

        # EE:
        #if abs(photon.eta()) > 1.566:
        if abs(photon.superCluster().position().Eta()) > 1.566:
            if photon.sigmaIetaIeta() > 0.03:
                return False

        return True


    #####ANALYSIS######
    def analyze(self, box):

        # Now check if there are at least 1 preselcted photon:
        if len(box.selectedPhotons) < 1:
            return False

        #first select photons
        box.FinalSelectedPhotons = []
        for pho in box.selectedPhotons:
            try:
#                if self.photonID(pho, box.vertex):
                if self.photonID(pho):
                    box.FinalSelectedPhotons.append(pho)
            except ZeroDivisionError:
                continue

        if len(box.FinalSelectedPhotons) < 1:
            return False

        return True


    def declareHistos(self):
        super(SinglePhotonAnalyzer, self).declareHistos()

        ###ADD YOUR HISTOGRAMS AFTER THIS LINE AS AbOVE#####
        self.declareHisto('photon_n', 10, -0.5, 9.5, "#gamma multiplicity")
        self.declareHisto('photon_pt', 100, 15, 115, "p_{T,#gamma} [GeV]")
        self.declareHisto('photon_eta', 52, -2.6, 2.6, "#eta_{#gamma}")
        self.declareHisto('photon_SCeta', 52, -2.6, 2.6, "#eta of #gamma supercluster")
        ###EB
        self.declareHisto('photon_R9_EB', 110, 0, 1.1, "#gamma R_{9}")
        self.declareHisto('photon_SigIeIe_EB', 100, 0.005, 0.015, "#gamma #sigma_{i#etai#eta}")
        self.declareHisto('photon_etaWidth_EB', 100, 0, 0.025, "#gamma #eta_{width}")
        self.declareHisto('photon_phiWidth_EB', 100, 0, 0.04, "#gamma #phi_{width}")
        self.declareHisto('photon_HoE_EB', 100, 0, 0.1, "#gamma H/E")
        self.declareHisto('photon_TrkIsoHollow04_EB', 100, 0, 5.0, "Track ISO with Hollow cone 0.4 [GeV]")
        self.declareHisto('photon_EcalIso04_EB', 100, 0, 5.0, "ECAL ISO (#DeltaR=0.4) [GeV]")
        self.declareHisto('photon_HcalIso04_EB', 100, 0, 5.0, "HCAL ISO (#DeltaR=0.4) [GeV]")
        #pf iso
        self.declareHisto('photon_pfChargedIso_EB', 100, 0, 10.0, "pf charged hardon ISO [GeV]")
        self.declareHisto('photon_pfNeutralIso_EB', 100, 0, 10.0, "pf neutral hardon ISO [GeV]")
        self.declareHisto('photon_pfphotonIso_EB', 100, 0, 10.0, "pf photon ISO [GeV]")

        ###EE
        self.declareHisto('photon_R9_EE', 110, 0, 1.1, "#gamma R_{9}")
        self.declareHisto('photon_SigIeIe_EE', 100, 0.015, 0.035, "#gamma #sigma_{i#etai#eta}")
        self.declareHisto('photon_etaWidth_EE', 100, 0, 0.05, "#gamma #eta_{width}")
        self.declareHisto('photon_phiWidth_EE', 100, 0, 0.1, "#gamma #phi_{width}")
        self.declareHisto('photon_HoE_EE', 100, 0, 0.1, "#gamma H/E")
        self.declareHisto('photon_TrkIsoHollow04_EE', 100, 0, 5.0, "Track ISO with Hollow cone 0.4 [GeV]")
        self.declareHisto('photon_EcalIso04_EE', 100, 0, 5.0, "ECAL ISO (#DeltaR=0.4) [GeV]")
        self.declareHisto('photon_HcalIso04_EE', 100, 0, 5.0, "HCAL ISO (#DeltaR=0.4) [GeV]")
        #pf iso
        self.declareHisto('photon_pfChargedIso_EE', 100, 0, 10.0, "pf charged hardon ISO [GeV]")
        self.declareHisto('photon_pfNeutralIso_EE', 100, 0, 10.0, "pf neutral hardon ISO [GeV]")
        self.declareHisto('photon_pfphotonIso_EE', 100, 0, 10.0, "pf photon ISO [GeV]")
        #ES
        self.declareHisto('photon_EesOEsc_EE', 100, 0, 0.4, "#gamma E_{es}/E_{SC}")


    def fillHistos(self, box, sample, weight=1):
        super(SinglePhotonAnalyzer, self).fillHistos(box, sample, weight)
 
        self.fillHisto('photon_n', sample, len(box.FinalSelectedPhotons), weight)

        for pho in box.FinalSelectedPhotons:
            self.fillHisto('photon_pt', sample, pho.pt(), weight)
            self.fillHisto('photon_eta', sample, pho.eta(), weight)
            self.fillHisto('photon_SCeta', sample, pho.superCluster().position().Eta(), weight)
            ##EB
            if (abs(pho.superCluster().position().Eta())<1.45):
                self.fillHisto('photon_R9_EB', sample, pho.r9(), weight)
                self.fillHisto('photon_SigIeIe_EB', sample, pho.sigmaIetaIeta(), weight)
                self.fillHisto('photon_etaWidth_EB', sample, pho.superCluster().etaWidth(), weight)
                self.fillHisto('photon_phiWidth_EB', sample, pho.superCluster().phiWidth(), weight)
                self.fillHisto('photon_HoE_EB', sample, pho.hadronicOverEm(), weight)
                self.fillHisto('photon_TrkIsoHollow04_EB', sample, pho.trkSumPtHollowConeDR04(), weight)
                self.fillHisto('photon_EcalIso04_EB', sample, pho.ecalRecHitSumEtConeDR04(), weight)
                self.fillHisto('photon_HcalIso04_EB', sample, pho.hcalTowerSumEtConeDR04(), weight)
                self.fillHisto('photon_pfChargedIso_EB', sample, pho.chargedHadronIso(), weight)
                self.fillHisto('photon_pfNeutralIso_EB', sample, pho.neutralHadronIso(), weight)
                self.fillHisto('photon_pfphotonIso_EB', sample, pho.photonIso(), weight)
            if (abs(pho.superCluster().position().Eta())>1.56):
                self.fillHisto('photon_R9_EE', sample, pho.r9(), weight)
                self.fillHisto('photon_SigIeIe_EE', sample, pho.sigmaIetaIeta(), weight)
                self.fillHisto('photon_etaWidth_EE', sample, pho.superCluster().etaWidth(), weight)
                self.fillHisto('photon_phiWidth_EE', sample, pho.superCluster().phiWidth(), weight)
                self.fillHisto('photon_HoE_EE', sample, pho.hadronicOverEm(), weight)
                self.fillHisto('photon_TrkIsoHollow04_EE', sample, pho.trkSumPtHollowConeDR04(), weight)
                self.fillHisto('photon_EcalIso04_EE', sample, pho.ecalRecHitSumEtConeDR04(), weight)
                self.fillHisto('photon_HcalIso04_EE', sample, pho.hcalTowerSumEtConeDR04(), weight)
                self.fillHisto('photon_pfChargedIso_EE', sample, pho.chargedHadronIso(), weight)
                self.fillHisto('photon_pfNeutralIso_EE', sample, pho.neutralHadronIso(), weight)
                self.fillHisto('photon_pfphotonIso_EE', sample, pho.photonIso(), weight)
                self.fillHisto('photon_EesOEsc_EE', sample, pho.superCluster().preshowerEnergy()/pho.superCluster().rawEnergy(), weight)


    def addEvent(self, box):
        self.data.append(len(box.FinalSelectedPhotons))
