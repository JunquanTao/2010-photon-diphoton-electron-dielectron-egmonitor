######JTao#########
import itertools
import math
import ROOT

from CMSOpenDataAnalysis.DecaysToElectrons.Analyzer import Analyzer, Object

##def deltaPhi
def twoobj_deltaPhi(phi1, phi2):
	deltaPhi = phi1 - phi2
        if (deltaPhi >= math.pi):
            deltaPhi -= 2.0*math.pi
        if (deltaPhi < -1.0*math.pi):
            deltaPhi += 2.0*math.pi
        return deltaPhi;

##def deltaR
def twoobj_deltaR(eta1, phi1, eta2, phi2):
	deltaEta = eta1-eta2
        deltaPhi = twoobj_deltaPhi(phi1, phi2)
        deltaR = math.sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi)
        return deltaR

##cosThetaStar
def getCosThetaCS(g1, g2, sqrtS=7):
        #print "JTao: pT1 = ",g1.Pt(),", pT2 = ",g2.Pt()," and sqrt = ",sqrtS
	beamE = 500.*sqrtS	
	b1 = ROOT.TLorentzVector(0., 0., beamE, beamE)
        b2 = ROOT.TLorentzVector(0., 0., -1.0*beamE, beamE)
        dielectron=g1+g2;
        boostToDipartFrame = -dielectron.BoostVector();
	#boost to dielectron frame
	refDiParticle_g1 = g1
	refDiParticle_g1.Boost(boostToDipartFrame)
	refDiParticle_b1 = b1
	refDiParticle_b1.Boost(boostToDipartFrame)
	refDiParticle_b2 = b2
	refDiParticle_b2.Boost(boostToDipartFrame)
	#Getting beam 3-vector from 4-vectors
	refDiParticle_vb1_direction = refDiParticle_b1.Vect().Unit()
	refDiParticle_vb2_direction = refDiParticle_b2.Vect().Unit()
	#Definition of zz directions
	direction_cs = (refDiParticle_vb1_direction - refDiParticle_vb2_direction).Unit() # CS direction
	
	return math.cos(direction_cs.Angle(refDiParticle_g1.Vect()))
  

class TwoElectronAnalyzer(Analyzer):
    """
    XXX
    Some comment here.
    """
    def __init__(self):
        super(TwoElectronAnalyzer, self).__init__()

    #####CHANGE THIS METHOD TO CHANGE PHOTON ID######
#    def electronID(self, electron, vertex):
    def electronID(self, electron):
        if electron.pt() < 10 or abs(electron.superCluster().position().Eta()) > 2.5:
            return False
         # excluded the EB/EE gap region
        if abs(electron.superCluster().position().Eta()) > 1.4442 and abs(electron.superCluster().position().Eta()) < 1.566:
            return False

        #"WP85" used : https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID

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


    def select_dielectroncandidates(self, box):
        dielectroncandidates = []
        ###dielectron XS measurement with 2010 data : JHEP01(2012)133
        ### PT>(23,20)GeV, DeltaR>0.45
        for ele1, ele2 in itertools.combinations(box.electrons, 2):
            if max(ele1.pt(), ele2.pt())<20. or min(ele1.pt(), ele2.pt())<10.:
                continue
            #if twoobj_deltaR(ele1.eta(), ele1.phi(), ele2.eta(), ele2.phi()) < 0.45:
            #    continue
            # now create a dielectron object and check the variable - mass
            dpcand = Object(ele1, ele2)
            # mass range?
            #if not (hcand.mass() > 100 and hcand.mass() < 150):
            #    continue
            dielectroncandidates.append(dpcand)
        return dielectroncandidates

    #####ANALYSIS######
    def analyze(self, box):

        #####START FROM A bOX CONTAINING SELECTED PHOTONS and MAKE DiParticleTON CANDIDATES

        # Now check if there are at least two electrons:
	if len(box.selectedElectrons) < 2:
            return False

        # Now create dielectron candidates and apply cuts:
        box.dielectroncandidates = self.select_dielectroncandidates(box)
        if len(box.dielectroncandidates) == 0:
            return False

        # OK if there are more than one dielectron candidate
        # pick the one with the hightest Sum PT of two electrons
        sortedDPs = sorted(box.dielectroncandidates,
                          key=lambda x: x.ele1.pt() + x.ele2.pt(), reverse=True)
        box.DP = sortedDPs[0]

        # create the selected dielectron candiate
        #box.DP = Object(box.DP) #### ele1=box.DP, ele2 = null

        return True

    def declareHistos(self):
        super(TwoElectronAnalyzer, self).declareHistos()

        ###ADD YOUR HISTOGRAMS AFTER THIS LINE AS AbOVE#####
        self.declareHisto('dielectron_mass', 100, 0, 200, "m_{ee} [GeV]")
        self.declareHisto('dielectron_mass0to6', 60, 0, 6, "m_{ee} [GeV]")
        self.declareHisto('dielectron_mass6to20', 70, 6, 20, "m_{ee} [GeV]")
        self.declareHisto('dielectron_mass20to60', 80, 20, 60, "m_{ee} [GeV]")
        self.declareHisto('dielectron_mass60to120', 60, 60, 120, "m_{ee} [GeV]")
        self.declareHisto('dielectron_mass120to150', 60, 120, 150, "m_{ee} [GeV]")
        self.declareHisto('dielectron_mass150to300', 75, 150, 300, "m_{ee} [GeV]")
	
        self.declareHisto('dielectron_pt', 100, 0, 200, "p_{T,ee} [GeV]")
        self.declareHisto('dielectron_deltaPhi', 100, 0, 3.1416, "#Delta#phi_{ee}")
        self.declareHisto('dielectron_cosThetaStar', 100, 0, 1, "cos#theta^{*}_{ee}")
	self.declareHisto('dielectron_lead_pt', 100, 20, 120, "p_{T,e}^{leading} [GeV]")
        self.declareHisto('dielectron_lead_eta', 52, -2.6, 2.6, "#eta_{e}^{leading}")
        self.declareHisto('dielectron_sublead_pt', 100, 15, 115, "p_{T,e}^{subleading} [GeV]")
        self.declareHisto('dielectron_sublead_eta', 52, -2.6, 2.6, "#eta_{e}^{subleading}")
	

    def fillHistos(self, box, sample, weight=1):
        super(TwoElectronAnalyzer, self).fillHistos(box, sample, weight)
	#print "JTao: mass - ",box.DP.mass()," and pT - ",box.DP.pt()
        self.fillHisto('dielectron_mass', sample, box.DP.mass(), weight)
        self.fillHisto('dielectron_mass0to6', sample, box.DP.mass(), weight)
        self.fillHisto('dielectron_mass6to20', sample, box.DP.mass(), weight)
        self.fillHisto('dielectron_mass20to60', sample, box.DP.mass(), weight)
        self.fillHisto('dielectron_mass60to120', sample, box.DP.mass(), weight)
        self.fillHisto('dielectron_mass120to150', sample, box.DP.mass(), weight)
        self.fillHisto('dielectron_mass150to300', sample, box.DP.mass(), weight)
	
        self.fillHisto('dielectron_pt', sample, box.DP.pt(), weight)
	#print "JTao: dielectron_phi -",box.DP.phi()
	#print "JTao : ele1_phi - ",box.DP.ele1.phi()," and ele2_phi - ",box.DP.ele2.phi()
        self.fillHisto('dielectron_deltaPhi', sample, abs(twoobj_deltaPhi(box.DP.ele1.phi(), box.DP.ele2.phi())), weight)
        box.ele1Obj = Object(box.DP.ele1)
        box.ele2Obj = Object(box.DP.ele2)
        #print "JTao: ele1Obj_pT = ",box.ele1Obj.pt()," and ele2Obj_pT = ",box.ele2Obj.pt()
	TL_ele1 = ROOT.TLorentzVector(box.ele1Obj.px(), box.ele1Obj.py(), box.ele1Obj.pz(), box.ele1Obj.energy())
	TL_ele2 = ROOT.TLorentzVector(box.ele2Obj.px(), box.ele2Obj.py(), box.ele2Obj.pz(), box.ele2Obj.energy())
        #costhetastar = getCosThetaCS(box.ele1Obj.p4(), box.ele2Obj.p4(), 7)
	costhetastar = getCosThetaCS(TL_ele1, TL_ele2, 7)
        #print "JTao : costhetastar - ",costhetastar
        self.fillHisto('dielectron_cosThetaStar', sample, abs(costhetastar), weight)
	leadpt = box.DP.ele1.pt()
	leadeta = box.DP.ele1.eta()
	subleadpt = box.DP.ele2.pt()
        subleadeta = box.DP.ele2.eta() 
	if (box.DP.ele2.pt()>box.DP.ele1.pt()):
		leadpt = box.DP.ele2.pt()
		leadeta = box.DP.ele2.eta() 
		subleadpt = box.DP.ele1.pt()
		subleadeta = box.DP.ele1.eta()
	self.fillHisto('dielectron_lead_pt', sample, leadpt, weight)
	self.fillHisto('dielectron_lead_eta', sample, leadeta, weight)
        self.fillHisto('dielectron_sublead_pt', sample, subleadpt, weight)
        self.fillHisto('dielectron_sublead_eta', sample, subleadeta, weight)


    def addEvent(self, box):
        self.data.append(box.DP.mass())
        
