######JTao#########
import itertools
import math
import ROOT

from CMSOpenDataAnalysis.DecaysToPhotons.Analyzer import Analyzer, Object

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
        diphoton=g1+g2;
        boostToDiphotonFrame = -diphoton.BoostVector();
	#boost to diphoton frame
	refDIPHO_g1 = g1
	refDIPHO_g1.Boost(boostToDiphotonFrame)
	refDIPHO_b1 = b1
	refDIPHO_b1.Boost(boostToDiphotonFrame)
	refDIPHO_b2 = b2
	refDIPHO_b2.Boost(boostToDiphotonFrame)
	#Getting beam 3-vector from 4-vectors
	refDIPHO_vb1_direction = refDIPHO_b1.Vect().Unit()
	refDIPHO_vb2_direction = refDIPHO_b2.Vect().Unit()
	#Definition of zz directions
	direction_cs = (refDIPHO_vb1_direction - refDIPHO_vb2_direction).Unit() # CS direction
	
	return math.cos(direction_cs.Angle(refDIPHO_g1.Vect()))
  

class TwoPhotonAnalyzer(Analyzer):
    """
    XXX
    Some comment here.
    """
    def __init__(self):
        super(TwoPhotonAnalyzer, self).__init__()

    #####CHANGE THIS METHOD TO CHANGE PHOTON ID######
#    def photonID(self, photon, vertex):
    def photonID(self, photon):
        if photon.pt() < 20 or abs(photon.eta()) > 2.5:
            return False
         # excluded the EB/EE gap region
        if abs(photon.eta()) > 1.4442 and abs(photon.eta()) < 1.566:
            return False

        #To make it simple, I select (ref. but a little diff. from diphoton XS measurement with 2010 data - JHEP01(2012)133):
        #sigmaIetaIeta < 0.01/0.03 in EB/EE and HoE (hadronicOverEm) < 0.05 
        #(loose selections of PAS-EGM-10-006 http://cdsweb.cern.ch/record/1324545)
        #Track ISO (cone DeltaR=0.4) trkSumPtHollowConeDR04 < 2./4. GeV in EB/EE
        #HCAL ISO (cone DeltaR=0.4) hcalTowerSumEtConeDR04 < 2./4. GeV in EB/EE
        #ECAL ISO (cone DeltaR=0.3) ecalRecHitSumEtConeDR03 < 0.3 * photon PT
 
        #common selections for both EB and EE:
        if photon.hadronicOverEm() > 0.05:
            return False
        if photon.ecalRecHitSumEtConeDR03() / photon.pt() > 0.3:
            return False
 	# No electron-veto required : bool hasPixelSeed() to keep more/conversion photons	
        # EB: 
 	#if abs(photon.eta()) < 1.4442:
        if abs(photon.superCluster().position().Eta()) < 1.4442:
            if photon.trkSumPtHollowConeDR04() > 2.:
                return False
            if photon.hcalTowerSumEtConeDR04() > 2.:
                return False
            if photon.sigmaIetaIeta() > 0.01:
                return False
        # EE:
	#if abs(photon.eta()) > 1.566:
        if abs(photon.superCluster().position().Eta()) > 1.566:
            if photon.trkSumPtHollowConeDR04() > 4.:
                return False
            if photon.hcalTowerSumEtConeDR04() > 4.:
                return False
            if photon.sigmaIetaIeta() > 0.03:
                return False

        return True


    def select_diphotoncandidates(self, box):
        diphotoncandidates = []
        ###diphoton XS measurement with 2010 data : JHEP01(2012)133
        ### PT>(23,20)GeV, DeltaR>0.45
        for pho1, pho2 in itertools.combinations(box.photons, 2):
            if max(pho1.pt(), pho2.pt())<23. or min(pho1.pt(), pho2.pt())<20.:
                continue
            if twoobj_deltaR(pho1.eta(), pho1.phi(), pho2.eta(), pho2.phi()) < 0.45:
                continue
            # now create a diphoton object and check the variable - mass
            dpcand = Object(pho1, pho2)
            # mass range?
            #if not (hcand.mass() > 100 and hcand.mass() < 150):
            #    continue
            diphotoncandidates.append(dpcand)
        return diphotoncandidates

    #####ANALYSIS######
    def analyze(self, box):

        #####START FROM A bOX CONTAINING SELECTED PHOTONS and MAKE DIPHOTON CANDIDATES

        # Now check if there are at least two photons:
	if len(box.selectedPhotons) < 2:
            return False

        # Now create diphoton candidates and apply cuts:
        box.diphotoncandidates = self.select_diphotoncandidates(box)
        if len(box.diphotoncandidates) == 0:
            return False

        # OK if there are more than one diphoton candidate
        # pick the one with the hightest Sum PT of two photons
        sortedDPs = sorted(box.diphotoncandidates,
                          key=lambda x: x.pho1.pt() + x.pho2.pt(), reverse=True)
        box.DP = sortedDPs[0]

        # create the selected diphoton candiate
        #box.DP = Object(box.DP) #### pho1=box.DP, pho2 = null

        return True

    def declareHistos(self):
        super(TwoPhotonAnalyzer, self).declareHistos()

        ###ADD YOUR HISTOGRAMS AFTER THIS LINE AS AbOVE#####
        self.declareHisto('diphoton_mass', 100, 0, 300, "m_{#gamma#gamma} [GeV]")
        self.declareHisto('diphoton_pt', 100, 0, 200, "p_{T,#gamma#gamma} [GeV]")
        self.declareHisto('diphoton_deltaPhi', 100, 0, 3.1416, "#Delta#phi_{#gamma#gamma}")
        self.declareHisto('diphoton_cosThetaStar', 100, 0, 1, "cos#theta^{*}_{#gamma#gamma}")
	self.declareHisto('diphoton_lead_pt', 100, 20, 120, "p_{T,#gamma}^{leading} [GeV]")
        self.declareHisto('diphoton_lead_eta', 52, -2.6, 2.6, "#eta_{#gamma}^{leading}")
        self.declareHisto('diphoton_sublead_pt', 100, 15, 115, "p_{T,#gamma}^{subleading} [GeV]")
        self.declareHisto('diphoton_sublead_eta', 52, -2.6, 2.6, "#eta_{#gamma}^{subleading}")
	

    def fillHistos(self, box, sample, weight=1):
        super(TwoPhotonAnalyzer, self).fillHistos(box, sample, weight)
	#print "JTao: mass - ",box.DP.mass()," and pT - ",box.DP.pt()
        self.fillHisto('diphoton_mass', sample, box.DP.mass(), weight)
        self.fillHisto('diphoton_pt', sample, box.DP.pt(), weight)
	#print "JTao: diphoton_phi -",box.DP.phi()
	#print "JTao : pho1_phi - ",box.DP.pho1.phi()," and pho2_phi - ",box.DP.pho2.phi()
        self.fillHisto('diphoton_deltaPhi', sample, abs(twoobj_deltaPhi(box.DP.pho1.phi(), box.DP.pho2.phi())), weight)
        box.pho1Obj = Object(box.DP.pho1)
        box.pho2Obj = Object(box.DP.pho2)
        #print "JTao: pho1Obj_pT = ",box.pho1Obj.pt()," and pho2Obj_pT = ",box.pho2Obj.pt()
	TL_pho1 = ROOT.TLorentzVector(box.pho1Obj.px(), box.pho1Obj.py(), box.pho1Obj.pz(), box.pho1Obj.energy())
	TL_pho2 = ROOT.TLorentzVector(box.pho2Obj.px(), box.pho2Obj.py(), box.pho2Obj.pz(), box.pho2Obj.energy())
        #costhetastar = getCosThetaCS(box.pho1Obj.p4(), box.pho2Obj.p4(), 7)
	costhetastar = getCosThetaCS(TL_pho1, TL_pho2, 7)
        #print "JTao : costhetastar - ",costhetastar
        self.fillHisto('diphoton_cosThetaStar', sample, abs(costhetastar), weight)
	leadpt = box.DP.pho1.pt()
	leadeta = box.DP.pho1.eta()
	subleadpt = box.DP.pho2.pt()
        subleadeta = box.DP.pho2.eta() 
	if (box.DP.pho2.pt()>box.DP.pho1.pt()):
		leadpt = box.DP.pho2.pt()
		leadeta = box.DP.pho2.eta() 
		subleadpt = box.DP.pho1.pt()
		subleadeta = box.DP.pho1.eta()
	self.fillHisto('diphoton_lead_pt', sample, leadpt, weight)
	self.fillHisto('diphoton_lead_eta', sample, leadeta, weight)
        self.fillHisto('diphoton_sublead_pt', sample, subleadpt, weight)
        self.fillHisto('diphoton_sublead_eta', sample, subleadeta, weight)


    def addEvent(self, box):
        self.data.append(box.DP.mass())
        
