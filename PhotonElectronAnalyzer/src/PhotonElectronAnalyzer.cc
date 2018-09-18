// -*- C++ -*-
//
// Package:    PhotonElectronAnalyzer
// Class:      PhotonElectronAnalyzer
// 
/**\class PhotonElectronAnalyzer PhotonElectronAnalyzer.cc CMSOpenDataAnalysis/PhotonElectronAnalyzer/src/PhotonElectronAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Junquan Tao
//         Created:  Mon Jul 16 09:15:26 CEST 2018
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "CLHEP/Vector/LorentzVector.h"
#include <HepMC/WeightContainer.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//HLT
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//Photon
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

//Electron
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//
// class decleration
//

#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;
using namespace reco;
using namespace edm;

#define TAODEBUG 0
#define TAODEBUGHLT 0

class PhotonElectronAnalyzer : public edm::EDAnalyzer {
public:
  explicit PhotonElectronAnalyzer(const edm::ParameterSet&);
  ~PhotonElectronAnalyzer();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  double twoobj_deltaPhi(double phi1, double phi2);
  double twoobj_deltaR(double eta1, double phi1, double eta2, double phi2);

  // ----------member data ---------------------------
  //edm::InputTag puInforProducer_;
  //edm::InputTag rhoCollection_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag vertexProducer_;
  edm::InputTag photonProducer_;
  edm::InputTag electronProducer_;

  int nEvt;// used to count the number of events

  // to be used for root output tree
  std::string outputFile; // output file
  TFile *outputfile;
  //====================
  TTree *eventTree;
  TTree *photonTree;
  TTree *diphotonTree;
  TTree *electronTree;
  TTree *dielectronTree;

  std::vector<std::string> HLTPaths_photon;
  std::vector<std::string> HLTPaths_diphoton;

  // Variables to fill
  Int_t run, lumis, event, IsData;
  int nvtx;
  double  gen_weight;
  int pu_n;

  int event_hlt_n;
  int event_PassHLT_SinglePhoton;
  int event_PassHLT_DiPhoton;


  int n_photon;

  //single photon
  double photon_pt, photon_eta, photon_SCeta, photon_phi, photon_energy, photon_Escraw;
  double photon_r9, photon_hoe, photon_sieie; 
  int photon_hasPixelSeed;
  double photon_etaWidth, photon_phiWidth;
  double photon_TrkIsoHollow04, photon_EcalIso04, photon_HcalIso04;
 

  //diphoton
  int n_diphoton_sel;
  int dipho_pho1_ind;  int dipho_pho2_ind;

  double dipho_pt,  dipho_eta,  dipho_phi, dipho_energy, dipho_mass;
  double dipho_CosThetaStar, dipho_DeltaPhi, dipho_DeltaR;

  double dipho_pho1_pt, dipho_pho1_eta, dipho_pho1_SCeta, dipho_pho1_phi, dipho_pho1_energy;
  double dipho_pho1_r9, dipho_pho1_hoe, dipho_pho1_sieie;
  int dipho_pho1_hasPixelSeed;
  double dipho_pho1_TrkIsoHollow04, dipho_pho1_EcalIso04, dipho_pho1_HcalIso04;
  double dipho_pho1_EcalIso03OverPT ;

  double dipho_pho2_pt, dipho_pho2_eta, dipho_pho2_SCeta, dipho_pho2_phi, dipho_pho2_energy;
  double dipho_pho2_r9, dipho_pho2_hoe, dipho_pho2_sieie;
  int dipho_pho2_hasPixelSeed;
  double dipho_pho2_TrkIsoHollow04, dipho_pho2_EcalIso04, dipho_pho2_HcalIso04;
  double dipho_pho2_EcalIso03OverPT ;

  //=========================
  int n_electron;
  //single electron
  double electron_pt, electron_eta, electron_phi, electron_SCeta;
  int electron_q;
  double electron_trckIso03Rel, electron_ecalIso03Rel, electron_hcalIso03Rel, electron_RelCombinedIso03;
  double electron_hoe, electron_sigieie;
  double electron_DeltaPhiSCTrk, electron_DeltaEtaSCTrk;
  int electron_Nmisshit;
  double electron_DisConv, electron_DeltaCotTheta;
  double electron_fbrem,  electron_ESCoPin, electron_AbsInvEmInvPin;

  //dielectron
  int n_dielectron_sel;
  int diele_ele1_ind;  int diele_ele2_ind;

  double diele_ele1_pt, diele_ele1_eta, diele_ele1_phi, diele_ele1_e, diele_ele1_SCeta, diele_ele1_Escraw;
  int diele_ele1_q;
  double diele_ele1_trckIso03Rel, diele_ele1_ecalIso03Rel, diele_ele1_hcalIso03Rel, diele_ele1_RelCombinedIso03;
  double diele_ele1_hoe, diele_ele1_sigieie, diele_ele1_DeltaPhiSCTrk, diele_ele1_DeltaEtaSCTrk;
  int diele_ele1_Nmisshit;
  double diele_ele1_DisConv, diele_ele1_DeltaCotTheta;
  double diele_ele1_fbrem, diele_ele1_ESCoPin, diele_ele1_AbsInvEmInvPin;

  double diele_ele2_pt, diele_ele2_eta, diele_ele2_phi, diele_ele2_e, diele_ele2_SCeta, diele_ele2_Escraw;
  int diele_ele2_q;
  double diele_ele2_trckIso03Rel, diele_ele2_ecalIso03Rel, diele_ele2_hcalIso03Rel, diele_ele2_RelCombinedIso03;
  double diele_ele2_hoe, diele_ele2_sigieie, diele_ele2_DeltaPhiSCTrk, diele_ele2_DeltaEtaSCTrk;
  int diele_ele2_Nmisshit;
  double diele_ele2_DisConv, diele_ele2_DeltaCotTheta;
  double diele_ele2_fbrem, diele_ele2_ESCoPin, diele_ele2_AbsInvEmInvPin;

  double diele_pt,  diele_eta,  diele_phi, diele_energy, diele_mass;
  double diele_DeltaR;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PhotonElectronAnalyzer::PhotonElectronAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  nEvt=0;

  //CMS AN AN-10-268
  string hlt_pho1="HLT_Photon10_Cleaned_L1R"; 
  string hlt_pho2="HLT_Photon15_Cleaned_L1R"; 
  string hlt_pho3="HLT_Photon20_Cleaned_L1R";   
  string hlt_pho4="HLT_Photon30_Cleaned_L1R";
  string hlt_pho5="HLT_Photon50_NoHE_Cleaned_L1R";
  HLTPaths_photon.push_back(hlt_pho1);
  HLTPaths_photon.push_back(hlt_pho2);
  HLTPaths_photon.push_back(hlt_pho3);
  HLTPaths_photon.push_back(hlt_pho4);
  HLTPaths_photon.push_back(hlt_pho5);

  //CMS AN -2011/032
  string hlt_dipho1="HLT_DoublePhoton15_L1R";  //141956 - 144116  2.9pb-1
  string hlt_dipho2="HLT_DoublePhoton17_L1R";  //146428 - 148058  14.52pb-1
  string hlt_dipho3="HLT_Photon17_Isol_SC17HE_L1R";  //148822 - 149294 18.44pb-1
  HLTPaths_diphoton.push_back(hlt_dipho1);
  HLTPaths_diphoton.push_back(hlt_dipho2);
  HLTPaths_diphoton.push_back(hlt_dipho3);

  outputFile= iConfig.getUntrackedParameter<std::string>("outputfileName");
  outputfile= TFile::Open(outputFile.c_str(),"RECREATE"); // open output file
  triggerResultsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag");
  vertexProducer_ = iConfig.getUntrackedParameter<edm::InputTag>("vertexProducer");
  photonProducer_ = iConfig.getUntrackedParameter<edm::InputTag>("photonProducer");
  electronProducer_ = iConfig.getUntrackedParameter<edm::InputTag>("electronProducer");

  outputfile->cd();
  eventTree = new TTree("eventTree", "per-event tree");
  photonTree = new TTree("photonTree", "single photon tree");
  diphotonTree = new TTree("diphotonTree", "di-photon tree");
  electronTree = new TTree("electronTree", "single electron tree");
  dielectronTree  = new TTree("dielectronTree", "di-electron tree");

  //----event tree----
  eventTree->Branch( "IsData", &IsData, "IsData/I" );
  eventTree->Branch( "run", &run, "run/I" );
  eventTree->Branch( "event", &event, "event/I" );
  eventTree->Branch( "lumis", &lumis, "lumis/I" );
  eventTree->Branch( "pu_n", &pu_n, "pu_n/I" );
  eventTree->Branch( "nvtx", &nvtx, "nvtx/I" );
  eventTree->Branch( "gen_weight", &gen_weight, "gen_weight/D" );
  eventTree->Branch( "event_hlt_n", &event_hlt_n, "event_hlt_n/I" );
  eventTree->Branch( "event_PassHLT_SinglePhoton", &event_PassHLT_SinglePhoton, "event_PassHLT_SinglePhoton/I" );
  eventTree->Branch( "event_PassHLT_DiPhoton", &event_PassHLT_DiPhoton, "event_PassHLT_DiPhoton/I" );

  //----single photon------
  photonTree->Branch( "IsData", &IsData, "IsData/I" );
  photonTree->Branch( "run", &run, "run/I" );
  photonTree->Branch( "event", &event, "event/I" );
  photonTree->Branch( "lumis", &lumis, "lumis/I" );
  photonTree->Branch( "pu_n", &pu_n, "pu_n/I" );
  photonTree->Branch( "nvtx", &nvtx, "nvtx/I" );
  photonTree->Branch( "gen_weight", &gen_weight, "gen_weight/D" );
  photonTree->Branch( "event_hlt_n", &event_hlt_n, "event_hlt_n/I" );
 
  photonTree->Branch( "event_PassHLT_SinglePhoton", &event_PassHLT_SinglePhoton, "event_PassHLT_SinglePhoton/I" );
 
  photonTree->Branch("n_photon", &n_photon,"n_photon/I");
  photonTree->Branch("photon_pt", &photon_pt,"photon_pt/D");
  photonTree->Branch("photon_eta", &photon_eta,"photon_eta/D");
  photonTree->Branch("photon_SCeta", &photon_SCeta,"photon_SCeta/D");
  photonTree->Branch("photon_phi", &photon_phi,"photon_phi/D");
  photonTree->Branch("photon_energy", &photon_energy,"photon_energy/D");
  photonTree->Branch("photon_Escraw", &photon_Escraw,"photon_Escraw/D");
  photonTree->Branch("photon_r9", &photon_r9,"photon_r9/D");
  photonTree->Branch("photon_hoe", &photon_hoe,"photon_hoe/D");
  photonTree->Branch("photon_sieie", &photon_sieie,"photon_sieie/D");
  photonTree->Branch("photon_hasPixelSeed", &photon_hasPixelSeed,"photon_hasPixelSeed/I");
  photonTree->Branch("photon_etaWidth", &photon_etaWidth,"photon_etaWidth/D");
  photonTree->Branch("photon_phiWidth", &photon_phiWidth,"photon_phiWidth/D");

  photonTree->Branch("photon_TrkIsoHollow04", &photon_TrkIsoHollow04,"photon_TrkIsoHollow04/D");
  photonTree->Branch("photon_EcalIso04", &photon_EcalIso04,"photon_EcalIso04/D");
  photonTree->Branch("photon_HcalIso04", &photon_HcalIso04,"photon_HcalIso04/D");

  //---diphoton----
  diphotonTree->Branch( "IsData", &IsData, "IsData/I" );
  diphotonTree->Branch( "run", &run, "run/I" );
  diphotonTree->Branch( "event", &event, "event/I" );
  diphotonTree->Branch( "lumis", &lumis, "lumis/I" );
  diphotonTree->Branch( "pu_n", &pu_n, "pu_n/I" );
  diphotonTree->Branch( "nvtx", &nvtx, "nvtx/I" );
  diphotonTree->Branch( "gen_weight", &gen_weight, "gen_weight/D" );
  diphotonTree->Branch( "event_hlt_n", &event_hlt_n, "event_hlt_n/I" );
 
  diphotonTree->Branch("event_PassHLT_DiPhoton", &event_PassHLT_DiPhoton, "event_PassHLT_DiPhoton/I" );
  diphotonTree->Branch("n_diphoton_sel", &n_diphoton_sel,"n_diphoton_sel/I");
  diphotonTree->Branch("dipho_pho1_ind", &dipho_pho1_ind,"dipho_pho1_ind/I");
  diphotonTree->Branch("dipho_pho2_ind", &dipho_pho2_ind,"dipho_pho2_ind/I");

  diphotonTree->Branch("dipho_pt", &dipho_pt,"dipho_pt/D");
  diphotonTree->Branch("dipho_eta", &dipho_eta,"dipho_eta/D");
  diphotonTree->Branch("dipho_phi", &dipho_phi,"dipho_phi/D");
  diphotonTree->Branch("dipho_energy", &dipho_energy,"dipho_energy/D");
  diphotonTree->Branch("dipho_mass", &dipho_mass,"dipho_mass/D");
  diphotonTree->Branch("dipho_CosThetaStar", &dipho_CosThetaStar,"dipho_CosThetaStar/D");
  diphotonTree->Branch("dipho_DeltaPhi", &dipho_DeltaPhi,"dipho_DeltaPhi/D");
  diphotonTree->Branch("dipho_DeltaR", &dipho_DeltaR,"dipho_DeltaR/D");

  diphotonTree->Branch("dipho_pho1_pt", &dipho_pho1_pt,"dipho_pho1_pt/D");
  diphotonTree->Branch("dipho_pho1_eta", &dipho_pho1_eta,"dipho_pho1_eta/D");
  diphotonTree->Branch("dipho_pho1_SCeta", &dipho_pho1_SCeta,"dipho_pho1_SCeta/D");
  diphotonTree->Branch("dipho_pho1_phi", &dipho_pho1_phi,"dipho_pho1_phi/D");
  diphotonTree->Branch("dipho_pho1_energy", &dipho_pho1_energy,"dipho_pho1_energy/D");
  diphotonTree->Branch("dipho_pho1_r9", &dipho_pho1_r9,"dipho_pho1_r9/D");
  diphotonTree->Branch("dipho_pho1_hoe", &dipho_pho1_hoe,"dipho_pho1_hoe/D");
  diphotonTree->Branch("dipho_pho1_sieie", &dipho_pho1_sieie,"dipho_pho1_sieie/D");
  diphotonTree->Branch("dipho_pho1_hasPixelSeed", &dipho_pho1_hasPixelSeed,"dipho_pho1_hasPixelSeed/I");
  diphotonTree->Branch("dipho_pho1_TrkIsoHollow04", &dipho_pho1_TrkIsoHollow04,"dipho_pho1_TrkIsoHollow04/D");
  diphotonTree->Branch("dipho_pho1_EcalIso04", &dipho_pho1_EcalIso04,"dipho_pho1_EcalIso04/D");
  diphotonTree->Branch("dipho_pho1_HcalIso04", &dipho_pho1_HcalIso04,"dipho_pho1_HcalIso04/D");
  diphotonTree->Branch("dipho_pho1_EcalIso03OverPT", &dipho_pho1_EcalIso03OverPT,"dipho_pho1_EcalIso03OverPT/D");

  diphotonTree->Branch("dipho_pho2_pt", &dipho_pho2_pt,"dipho_pho2_pt/D");
  diphotonTree->Branch("dipho_pho2_eta", &dipho_pho2_eta,"dipho_pho2_eta/D");
  diphotonTree->Branch("dipho_pho2_SCeta", &dipho_pho2_SCeta,"dipho_pho2_SCeta/D");
  diphotonTree->Branch("dipho_pho2_phi", &dipho_pho2_phi,"dipho_pho2_phi/D");
  diphotonTree->Branch("dipho_pho2_energy", &dipho_pho2_energy,"dipho_pho2_energy/D");
  diphotonTree->Branch("dipho_pho2_r9", &dipho_pho2_r9,"dipho_pho2_r9/D");
  diphotonTree->Branch("dipho_pho2_hoe", &dipho_pho2_hoe,"dipho_pho2_hoe/D");
  diphotonTree->Branch("dipho_pho2_sieie", &dipho_pho2_sieie,"dipho_pho2_sieie/D");
  diphotonTree->Branch("dipho_pho2_hasPixelSeed", &dipho_pho2_hasPixelSeed,"dipho_pho2_hasPixelSeed/I");
  diphotonTree->Branch("dipho_pho2_TrkIsoHollow04", &dipho_pho2_TrkIsoHollow04,"dipho_pho2_TrkIsoHollow04/D");
  diphotonTree->Branch("dipho_pho2_EcalIso04", &dipho_pho2_EcalIso04,"dipho_pho2_EcalIso04/D");
  diphotonTree->Branch("dipho_pho2_HcalIso04", &dipho_pho2_HcalIso04,"dipho_pho2_HcalIso04/D");
  diphotonTree->Branch("dipho_pho2_EcalIso03OverPT", &dipho_pho2_EcalIso03OverPT,"dipho_pho2_EcalIso03OverPT/D");

  //diphotonTree->Branch("", &,"/D");

  //----------
  electronTree->Branch( "IsData", &IsData, "IsData/I" );
  electronTree->Branch( "run", &run, "run/I" );
  electronTree->Branch( "event", &event, "event/I" );
  electronTree->Branch( "lumis", &lumis, "lumis/I" );
  electronTree->Branch( "pu_n", &pu_n, "pu_n/I" );
  electronTree->Branch( "nvtx", &nvtx, "nvtx/I" );
  electronTree->Branch( "gen_weight", &gen_weight, "gen_weight/D" );

  electronTree->Branch("n_electron", &n_electron,"n_electron/I");

  electronTree->Branch("electron_pt", &electron_pt,"electron_pt/D");
  electronTree->Branch("electron_eta", &electron_eta,"electron_eta/D");
  electronTree->Branch("electron_phi", &electron_phi,"electron_phi/D");
  electronTree->Branch("electron_SCeta", &electron_SCeta,"electron_SCeta/D");
  electronTree->Branch("electron_q", &electron_q,"electron_q/I");

  electronTree->Branch("electron_trckIso03Rel", &electron_trckIso03Rel,"electron_trckIso03Rel/D");
  electronTree->Branch("electron_ecalIso03Rel", &electron_ecalIso03Rel,"electron_ecalIso03Rel/D");
  electronTree->Branch("electron_hcalIso03Rel", &electron_hcalIso03Rel,"electron_hcalIso03Rel/D");
  electronTree->Branch("electron_RelCombinedIso03", &electron_RelCombinedIso03,"electron_RelCombinedIso03/D");

  electronTree->Branch("electron_hoe", &electron_hoe,"electron_hoe/D");
  electronTree->Branch("electron_sigieie", &electron_sigieie,"electron_sigieie/D");
  electronTree->Branch("electron_DeltaPhiSCTrk", &electron_DeltaPhiSCTrk,"electron_DeltaPhiSCTrk/D");
  electronTree->Branch("electron_DeltaEtaSCTrk", &electron_DeltaEtaSCTrk,"electron_DeltaEtaSCTrk/D");
  electronTree->Branch("electron_Nmisshit", &electron_Nmisshit,"electron_Nmisshit/I");
  electronTree->Branch("electron_DisConv", &electron_DisConv,"electron_DisConv/D");
  electronTree->Branch("electron_DeltaCotTheta", &electron_DeltaCotTheta,"electron_DeltaCotTheta/D");

  electronTree->Branch("electron_fbrem", &electron_fbrem,"electron_fbrem/D");
  electronTree->Branch("electron_ESCoPin", &electron_ESCoPin,"electron_ESCoPin/D");
  electronTree->Branch("electron_AbsInvEmInvPin", &electron_AbsInvEmInvPin,"electron_AbsInvEmInvPin/D");
  //electronTree->Branch("", &,"/D");

  //di-electron
  dielectronTree->Branch( "IsData", &IsData, "IsData/I" );
  dielectronTree->Branch( "run", &run, "run/I" );
  dielectronTree->Branch( "event", &event, "event/I" );
  dielectronTree->Branch( "lumis", &lumis, "lumis/I" );
  dielectronTree->Branch( "pu_n", &pu_n, "pu_n/I" );
  dielectronTree->Branch( "nvtx", &nvtx, "nvtx/I" );
  dielectronTree->Branch( "gen_weight", &gen_weight, "gen_weight/D" );

  dielectronTree->Branch("n_dielectron_sel", &n_dielectron_sel,"n_dielectron_sel/I");
  dielectronTree->Branch("diele_ele1_ind", &diele_ele1_ind,"diele_ele1_ind/I");
  dielectronTree->Branch("diele_ele2_ind", &diele_ele2_ind,"diele_ele2_ind/I");

  dielectronTree->Branch("diele_ele1_pt", &diele_ele1_pt,"diele_ele1_pt/D");
  dielectronTree->Branch("diele_ele1_eta", &diele_ele1_eta,"diele_ele1_eta/D");
  dielectronTree->Branch("diele_ele1_phi", &diele_ele1_phi,"diele_ele1_phi/D");
  dielectronTree->Branch("diele_ele1_e", &diele_ele1_e,"diele_ele1_e/D");
  dielectronTree->Branch("diele_ele1_SCeta", &diele_ele1_SCeta,"diele_ele1_SCeta/D");
  dielectronTree->Branch("diele_ele1_Escraw", &diele_ele1_Escraw,"diele_ele1_Escraw/D");
  dielectronTree->Branch("diele_ele1_q", &diele_ele1_q,"diele_ele1_q/I");
  dielectronTree->Branch("diele_ele1_trckIso03Rel", &diele_ele1_trckIso03Rel,"diele_ele1_trckIso03Rel/D");
  dielectronTree->Branch("diele_ele1_ecalIso03Rel", &diele_ele1_ecalIso03Rel,"diele_ele1_ecalIso03Rel/D");
  dielectronTree->Branch("diele_ele1_hcalIso03Rel", &diele_ele1_hcalIso03Rel,"diele_ele1_hcalIso03Rel/D");
  dielectronTree->Branch("diele_ele1_RelCombinedIso03", &diele_ele1_RelCombinedIso03,"diele_ele1_RelCombinedIso03/D");
  dielectronTree->Branch("diele_ele1_hoe", &diele_ele1_hoe,"diele_ele1_hoe/D");
  dielectronTree->Branch("diele_ele1_sigieie", &diele_ele1_sigieie,"diele_ele1_sigieie/D");
  dielectronTree->Branch("diele_ele1_DeltaPhiSCTrk", &diele_ele1_DeltaPhiSCTrk,"diele_ele1_DeltaPhiSCTrk/D");
  dielectronTree->Branch("diele_ele1_DeltaEtaSCTrk", &diele_ele1_DeltaEtaSCTrk,"diele_ele1_DeltaEtaSCTrk/D");
  dielectronTree->Branch("diele_ele1_Nmisshit", &diele_ele1_Nmisshit,"diele_ele1_Nmisshit/I");
  dielectronTree->Branch("diele_ele1_DisConv", &diele_ele1_DisConv,"diele_ele1_DisConv/D");
  dielectronTree->Branch("diele_ele1_DeltaCotTheta", &diele_ele1_DeltaCotTheta,"diele_ele1_DeltaCotTheta/D");
  dielectronTree->Branch("diele_ele1_fbrem", &diele_ele1_fbrem,"diele_ele1_fbrem/D");
  dielectronTree->Branch("diele_ele1_ESCoPin", &diele_ele1_ESCoPin,"diele_ele1_ESCoPin/D");
  dielectronTree->Branch("diele_ele1_AbsInvEmInvPin", &diele_ele1_AbsInvEmInvPin,"diele_ele1_AbsInvEmInvPin/D");
 
  dielectronTree->Branch("diele_ele2_pt", &diele_ele2_pt,"diele_ele2_pt/D");
  dielectronTree->Branch("diele_ele2_eta", &diele_ele2_eta,"diele_ele2_eta/D");
  dielectronTree->Branch("diele_ele2_phi", &diele_ele2_phi,"diele_ele2_phi/D");
  dielectronTree->Branch("diele_ele2_e", &diele_ele2_e,"diele_ele2_e/D");
  dielectronTree->Branch("diele_ele2_SCeta", &diele_ele2_SCeta,"diele_ele2_SCeta/D");
  dielectronTree->Branch("diele_ele2_Escraw", &diele_ele2_Escraw,"diele_ele2_Escraw/D");
  dielectronTree->Branch("diele_ele2_q", &diele_ele2_q,"diele_ele2_q/I");
  dielectronTree->Branch("diele_ele2_trckIso03Rel", &diele_ele2_trckIso03Rel,"diele_ele2_trckIso03Rel/D");
  dielectronTree->Branch("diele_ele2_ecalIso03Rel", &diele_ele2_ecalIso03Rel,"diele_ele2_ecalIso03Rel/D");
  dielectronTree->Branch("diele_ele2_hcalIso03Rel", &diele_ele2_hcalIso03Rel,"diele_ele2_hcalIso03Rel/D");
  dielectronTree->Branch("diele_ele2_RelCombinedIso03", &diele_ele2_RelCombinedIso03,"diele_ele2_RelCombinedIso03/D");
  dielectronTree->Branch("diele_ele2_hoe", &diele_ele2_hoe,"diele_ele2_hoe/D");
  dielectronTree->Branch("diele_ele2_sigieie", &diele_ele2_sigieie,"diele_ele2_sigieie/D");
  dielectronTree->Branch("diele_ele2_DeltaPhiSCTrk", &diele_ele2_DeltaPhiSCTrk,"diele_ele2_DeltaPhiSCTrk/D");
  dielectronTree->Branch("diele_ele2_DeltaEtaSCTrk", &diele_ele2_DeltaEtaSCTrk,"diele_ele2_DeltaEtaSCTrk/D");
  dielectronTree->Branch("diele_ele2_Nmisshit", &diele_ele2_Nmisshit,"diele_ele2_Nmisshit/I");
  dielectronTree->Branch("diele_ele2_DisConv", &diele_ele2_DisConv,"diele_ele2_DisConv/D");
  dielectronTree->Branch("diele_ele2_DeltaCotTheta", &diele_ele2_DeltaCotTheta,"diele_ele2_DeltaCotTheta/D");
  dielectronTree->Branch("diele_ele2_fbrem", &diele_ele2_fbrem,"diele_ele2_fbrem/D");
  dielectronTree->Branch("diele_ele2_ESCoPin", &diele_ele2_ESCoPin,"diele_ele2_ESCoPin/D");
  dielectronTree->Branch("diele_ele2_AbsInvEmInvPin", &diele_ele2_AbsInvEmInvPin,"diele_ele2_AbsInvEmInvPin/D");

  dielectronTree->Branch("diele_pt", &diele_pt,"diele_pt/D");
  dielectronTree->Branch("diele_eta", &diele_eta,"diele_eta/D");
  dielectronTree->Branch("diele_phi", &diele_phi,"diele_phi/D");
  dielectronTree->Branch("diele_energy", &diele_energy,"diele_energy/D");
  dielectronTree->Branch("diele_mass", &diele_mass,"diele_mass/D");
  dielectronTree->Branch("diele_DeltaR", &diele_DeltaR,"diele_DeltaR/D");

  //dielectronTree->Branch("", &,"/I");
  //tree->Branch("", &,"/I");

}


PhotonElectronAnalyzer::~PhotonElectronAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


double PhotonElectronAnalyzer::twoobj_deltaPhi(double phi1, double phi2){
  double deltaPhi = phi1 - phi2;
  if(deltaPhi >= TMath::Pi()) deltaPhi -= 2.0*TMath::Pi();
  if(deltaPhi < -1.0*TMath::Pi()) deltaPhi += 2.0*TMath::Pi();
  return deltaPhi;
}

double PhotonElectronAnalyzer::twoobj_deltaR(double eta1, double phi1, double eta2, double phi2){
  double deltaEta = eta1 - eta2;
  double deltaPhi = twoobj_deltaPhi(phi1, phi2);
  double deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
  return deltaR;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PhotonElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;

  nEvt++;

  if(TAODEBUG==1) cout<<"JTao: event info"<<endl;
  run = iEvent.id().run();
  event = iEvent.id().event();
  lumis = iEvent.eventAuxiliary().luminosityBlock();
  IsData = iEvent.isRealData();

  gen_weight = 1.0;
  pu_n = -1;
  if( ! IsData ) { //MC
    edm::Handle<GenEventInfoProduct> genEvent;
    try{
      iEvent.getByLabel("generator", genEvent);
      if( !genEvent.isValid() ){
	cout <<  "   ===> No Valid GenEventInfoProduct, skip PDF Infos" << endl;
      }else{
	//const auto &weights = genEvent->weights();
	//gen_weight = weights[0];
	if(genEvent->weights().size()>0 ) gen_weight=genEvent->weight();
      }
    }catch ( cms::Exception& ex ) { 
      edm::LogError("TaoAnalyzer") <<"Error! can't get collection with label : generator"; 
    }
    //PU infor

    edm::Handle<std::vector< PileupSummaryInfo> > PupInfo;
    try {
      iEvent.getByLabel("addPileupInfo", PupInfo);
    } catch ( cms::Exception& ex ) { 
      edm::LogError("TaoAnalyzer") <<"Error! can't get collection with label : addPileupInfo"; 
    }
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int PU_Bunchcrossing = PVI->getBunchCrossing();
      if(PU_Bunchcrossing == 0) {
	pu_n = PVI->getPU_NumInteractions();
	break;
      }
    }
  }
 

  //RecVertices
  nvtx = 0;   
  Handle<reco::VertexCollection> recVtxs;
  //edm::Handle<std::vector<reco::Vertex> > recVtxs;
  try {
    iEvent.getByLabel(vertexProducer_, recVtxs);
  } catch ( cms::Exception& ex ) {
    edm::LogError("TaoAnalyzer") <<"Error! can't get collection with label : Vertices"; }
  reco::VertexCollection::const_iterator myvtx;
  //std::vector<reco::Vertex>::const_iterator myvtx;
  //nvtx = recVtxs->size();
  for(myvtx=recVtxs->begin(); myvtx!=recVtxs->end(); myvtx++){
    nvtx++; 
  }

  //HLT
  event_hlt_n = 0;
  event_PassHLT_SinglePhoton = 0;
  edm::Handle<edm::TriggerResults> trigResults;
  try {
    iEvent.getByLabel(triggerResultsTag_,trigResults);
  } catch ( cms::Exception& ex ) {
    edm::LogError("TaoAnalyzer") <<"Error! can't get collection with label : HLT trigger";
  }
  if (trigResults.isValid()) {
    int ntrigs = int(trigResults->size());
    event_hlt_n = ntrigs;
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*trigResults);
    for (int itrig = 0; itrig != ntrigs; ++itrig){
      const std::string trigName = triggerNames.triggerName(itrig);
      if(TAODEBUGHLT == 1 && nEvt == 1) cout<<"JTao: HLT name - "<<trigName<<" for itrig="<<itrig<<endl;
      //single photon HLT
      for(unsigned int i=0; i<HLTPaths_photon.size(); i++){
	string MyHLTName = HLTPaths_photon[i];
	size_t iffind = trigName.find(MyHLTName);
	if ((iffind != std::string::npos)) {
	  if(trigResults->accept(itrig)) event_PassHLT_SinglePhoton = 1;
	}
      }
      //diphoton HLT     
      for(unsigned int i=0; i<HLTPaths_diphoton.size(); i++){
	string MyHLTName = HLTPaths_diphoton[i];
	size_t iffind = trigName.find(MyHLTName);
	if ((iffind != std::string::npos)) {
	  if(trigResults->accept(itrig)) event_PassHLT_DiPhoton = 1;
	}
      }
      //
    }
  }
  // //////////////////////////////////Fill event tree
  eventTree->Fill();

  ///////////////////////////////Photons/////////////////////////
  n_photon = 0;
  edm::Handle<reco::PhotonCollection> photonHandle;
  iEvent.getByLabel(photonProducer_,photonHandle);
  if( photonHandle.isValid() ) {
    const reco::PhotonCollection* phoCollection = photonHandle.product();  
    reco::PhotonCollection::const_iterator photons;
    
    //======Single photon analyzer QCD-10-019 :   Phys. Rev. Lett. 106, 082001 
    int SelectedPhotonIndex = -1; double maxPhotonPT = 0;
    n_photon =  0; //phoCollection->size();
    int Diphoton_LeadInd = -1, Diphoton_SubLeadInd = -1; 
    n_diphoton_sel =  0;
    dipho_pho1_ind = -1;  dipho_pho2_ind = -1; double maxSumPT =0.;
    int photonIndex = 0;
    
    
    for ( photons = phoCollection->begin(); photons != phoCollection->end(); ++photons ) {
      if(TAODEBUG==1) cout<<"JTao : in the photon loops!"<<endl;
      n_photon ++;


      //Single photon analyzer QCD-10-019 :   Phys. Rev. Lett. 106, 082001  
      //Kinematics
      double et_p = photons->et();
      double scEta_p = photons->superCluster()->position().Eta();
      double phi_p = photons->phi();
      //pixel-seed e-veto
      //sigmaIetaIeta < 0.011/0.03 in EB/EE
      //H/E (hadronicOverEm) < 0.05 
      //IsoTRK < 2 GeV, 0.04-0.4 hollow cone
      //IsoECAL < 4.2 GeV, 0.06-0.4 cone
      //IsoHCAL < 2.2 GeV, 0.15-0.4 cone
      //if ( et_p > 21. && fabs(scEta_p) < 1.45   && !photons->hasPixelSeed() && photons->hadronicOverEm() < 0.05 && photons->trkSumPtHollowConeDR04() < 2. && photons->ecalRecHitSumEtConeDR04() < 4.2 && photons->hcalTowerSumEtConeDR04() < 2.2 ) {
      if ( et_p > 21. && ( fabs(scEta_p) < 1.4442 || abs(scEta_p)>1.566 && fabs(scEta_p)<2.5 )  && !photons->hasPixelSeed() && photons->hadronicOverEm() < 0.05 && photons->trkSumPtHollowConeDR04() < 2. && photons->ecalRecHitSumEtConeDR04() < 4.2 && photons->hcalTowerSumEtConeDR04() < 2.2 ) {  //2010B
        if(TAODEBUG==1) cout<<"JTao : passed the single photon selection!"<<endl;
	if( et_p > maxPhotonPT ){
	  SelectedPhotonIndex = photonIndex;
          maxPhotonPT = et_p;
	}
      }


      //diphoton analyzer  QCD-10-035 : JHEP01(2012)133)
      //To make it simple, I select (ref. but a little diff. from diphoton XS measurement with 2010 data - JHEP01(2012)133)
      //(loose selections of PAS-EGM-10-006 http://cdsweb.cern.ch/record/1324545)
      //pixel-seed e-veto
      //sigmaIetaIeta < 0.011/0.03 in EB/EE and HoE (hadronicOverEm) < 0.05 
      //Track ISO (cone DeltaR=0.4) trkSumPtHollowConeDR04 < 2./4. GeV in EB/EE
      //HCAL ISO (cone DeltaR=0.4) hcalTowerSumEtConeDR04 < 2./4. GeV in EB/EE
      //ECAL ISO (cone DeltaR=0.3) ecalRecHitSumEtConeDR03 < 0.2 * photon PT     
      if(et_p > 20. && (fabs(scEta_p)<1.44442 || fabs(scEta_p)>1.566 && fabs(scEta_p)<2.5) ){ //kinematics
	if(TAODEBUG==1) cout<<"JTao : passed the first photon kinematics in diphoton selection!"<<endl;
	if(photons->hadronicOverEm()<0.05 && !photons->hasPixelSeed() && photons->ecalRecHitSumEtConeDR03()/et_p < 0.2 ){//common selections
	  if( (fabs(scEta_p)<1.4442 && photons->trkSumPtHollowConeDR04() < 2. && photons->hcalTowerSumEtConeDR04() < 2. && photons->sigmaIetaIeta() < 0.01) || (fabs(scEta_p)>1.566 && fabs(scEta_p)<2.5 && photons->trkSumPtHollowConeDR04() < 4. && photons->hcalTowerSumEtConeDR04() < 4. && photons->sigmaIetaIeta() < 0.03) ){ //EB or EE isolation
	    if(TAODEBUG==1) cout<<"JTao : passed the first photon iso in diphoton selection with leading index = "<<photonIndex<<endl;
	    //photon2
            int photon2Ind = 0;
            for ( reco::PhotonCollection::const_iterator photons2 = phoCollection->begin(); photons2 != phoCollection->end(); ++photons2 ) {
	      if(photon2Ind>photonIndex){//not photon1
		if(TAODEBUG==1) cout<<"JTao : starting find the second photon index in diphoton selection!"<<endl;
		double et_p2 = photons2->et();
		double scEta_p2 = photons2->superCluster()->position().Eta();
		double phi_p2 = photons2->phi();
		if(et_p2 > 20. && (fabs(scEta_p2)<1.44442 || fabs(scEta_p2)>1.566 && fabs(scEta_p2)<2.5) ){ //kinematics
                  if(TAODEBUG==1) cout<<"JTao : passed the second photon kinematics in diphoton selection!"<<endl;        
		  if(photons2->hadronicOverEm()<0.05 && !photons2->hasPixelSeed() && photons2->ecalRecHitSumEtConeDR03()/et_p2 < 0.2 ){//common selections
		    if( (fabs(scEta_p2)<1.4442 && photons->trkSumPtHollowConeDR04() < 2. && photons2->hcalTowerSumEtConeDR04() < 2. && photons2->sigmaIetaIeta() < 0.01) || (fabs(scEta_p2)>1.566 && fabs(scEta_p2)<2.5 && photons2->trkSumPtHollowConeDR04() < 4. && photons2->hcalTowerSumEtConeDR04() < 4. && photons2->sigmaIetaIeta() < 0.03) ){ //EB or EE isolation
		      if(TAODEBUG==1) cout<<"JTao : passed the second photon iso in diphoton selection!"<<endl;
		      //if(TAODEBUG==1) cout<<"JTao : inputs of twoobj_deltaR scEta_p= "<<scEta_p<<" phi_p= "<<phi_p<<" scEta_p2= "<<scEta_p2<<" and phi_p2= "<<phi_p2<<endl;
		      double DelataRgg = twoobj_deltaR(scEta_p, phi_p, scEta_p2, phi_p2);
		      if(TAODEBUG==1) cout<<"JTao : in diphoton selection DelataRgg = "<<DelataRgg<<endl;
		      //PT>(23,20)GeV, DeltaR>0.45
		      if(max(et_p, et_p2)>23. && DelataRgg>0.45){
			if(TAODEBUG==1) cout<<"JTao : passed diphoton selection with pho1_ind="<<photonIndex<<" and pho2_ind="<<photon2Ind<<endl;
			n_diphoton_sel++;
			Diphoton_LeadInd = photonIndex; Diphoton_SubLeadInd = photon2Ind;
			if(et_p2 > et_p){
			  Diphoton_LeadInd = photon2Ind; Diphoton_SubLeadInd = photonIndex; 
			}
			if( et_p2 + et_p > maxSumPT ) { // if more than 2 selected diphotons, select the pair with the largest sum PT
			  maxSumPT = et_p2 + et_p;
			  dipho_pho1_ind = Diphoton_LeadInd;  dipho_pho2_ind = Diphoton_SubLeadInd;  
			}
		      }
		    }
		  }
		}
	      }//not photon1
	      photon2Ind++;
	    }//photons2 loop

	  }
	}
      }
      
      photonIndex++;
    }//photons loop

    //=======Single photon tree===========
    photonIndex = 0;
    for ( photons = phoCollection->begin(); photons != phoCollection->end(); ++photons ) {
      if(photonIndex == SelectedPhotonIndex){
        if(TAODEBUG==1) cout<<"JTao : trying to fill the single photon tree!"<<endl;
	photon_pt = photons->et(); photon_eta = photons->eta();  photon_SCeta = photons->superCluster()->position().Eta();  photon_phi = photons->phi();  photon_energy = photons->energy();  photon_Escraw = photons->superCluster()->rawEnergy();;
	photon_r9 = photons->r9();  photon_hoe = photons->hadronicOverEm();  photon_sieie = photons->sigmaIetaIeta();  
	photon_hasPixelSeed = photons->hasPixelSeed();
        photon_etaWidth = photons->superCluster()->etaWidth(); photon_phiWidth = photons->superCluster()->phiWidth();
	photon_TrkIsoHollow04 = photons->trkSumPtHollowConeDR04();  photon_EcalIso04 =  photons->ecalRecHitSumEtConeDR04(); photon_HcalIso04 = photons->hcalTowerSumEtConeDR04();
	//Fill Single photon tree
        photonTree->Fill();
      }
      photonIndex++;
    }

    //=======Di-photon tree===========
    int pho1Ind = 0;
    for ( photons = phoCollection->begin(); photons != phoCollection->end(); ++photons ) {
      if(pho1Ind == dipho_pho1_ind){ // find leading 
        if(TAODEBUG==1) cout<<"JTao : find the leading photon!"<<endl;
	int pho2Ind = 0;
	for ( reco::PhotonCollection::const_iterator photons2 = phoCollection->begin(); photons2 != phoCollection->end(); ++photons2 ) {
	  if(pho2Ind == dipho_pho2_ind){// find subleading
            if(TAODEBUG==1) cout<<"JTao : find the subleading photon!"<<endl;

	    dipho_pho1_pt = photons->et(); dipho_pho1_eta = photons->eta(); dipho_pho1_SCeta = photons->superCluster()->position().Eta();  dipho_pho1_phi = photons->phi(); dipho_pho1_energy= photons->energy();
	    dipho_pho1_r9 = photons->r9(); dipho_pho1_hoe = photons->hadronicOverEm(); dipho_pho1_sieie = photons->sigmaIetaIeta();
            dipho_pho1_hasPixelSeed = photons->hasPixelSeed();  dipho_pho1_TrkIsoHollow04 = photons->trkSumPtHollowConeDR04(); dipho_pho1_EcalIso04 =  photons->ecalRecHitSumEtConeDR04(); dipho_pho1_HcalIso04 = photons->hcalTowerSumEtConeDR04();
            dipho_pho1_EcalIso03OverPT = photons->ecalRecHitSumEtConeDR03()/photons->et();

	    dipho_pho2_pt = photons2->et(); dipho_pho2_eta = photons2->eta(); dipho_pho2_SCeta = photons2->superCluster()->position().Eta();  dipho_pho2_phi = photons2->phi(); dipho_pho2_energy= photons2->energy();
	    dipho_pho2_r9 = photons2->r9(); dipho_pho2_hoe = photons2->hadronicOverEm(); dipho_pho2_sieie = photons2->sigmaIetaIeta();
            dipho_pho2_hasPixelSeed = photons2->hasPixelSeed();  dipho_pho2_TrkIsoHollow04 = photons2->trkSumPtHollowConeDR04(); dipho_pho2_EcalIso04 =  photons2->ecalRecHitSumEtConeDR04(); dipho_pho2_HcalIso04 = photons2->hcalTowerSumEtConeDR04();
            dipho_pho2_EcalIso03OverPT = photons2->ecalRecHitSumEtConeDR03()/photons2->et();

	    TLorentzVector lead_p4; lead_p4.SetPtEtaPhiM(dipho_pho1_pt, dipho_pho1_eta, dipho_pho1_phi, 0);
	    TLorentzVector sublead_p4; sublead_p4.SetPtEtaPhiM(dipho_pho2_pt, dipho_pho2_eta, dipho_pho2_phi, 0);
	    TLorentzVector diphoton_p4 = lead_p4 + sublead_p4;
	    dipho_pt = diphoton_p4.Pt();  dipho_eta = diphoton_p4.Eta();  dipho_phi = diphoton_p4.Phi(); dipho_energy = diphoton_p4.E();  dipho_mass = diphoton_p4.M();
	    dipho_CosThetaStar = fabs(TMath::TanH((lead_p4.Rapidity() - sublead_p4.Rapidity())/2.0)); //2010 definition
	    dipho_DeltaPhi = twoobj_deltaPhi(dipho_pho1_phi, dipho_pho2_phi);  dipho_DeltaR = twoobj_deltaR(dipho_pho1_eta, dipho_pho1_phi, dipho_pho2_eta, dipho_pho2_phi);
	    //Fill Single photon tree
	    diphotonTree->Fill();
	  }
       
	  pho2Ind++;
	}
      }
      pho1Ind++;
    }


  }//photonHandle.isValid()

  ///////////////////////////////Electrons/////////////////////////
  edm::Handle<reco::GsfElectronCollection> pElectrons;
  iEvent.getByLabel(electronProducer_, pElectrons);
  if( pElectrons.isValid() ) {

    const reco::GsfElectronCollection* eleCollection = pElectrons.product();  
    reco::GsfElectronCollection::const_iterator electrons;

    int SelectedEleIndex = -1; double maxElePT = 0;
    n_electron =  0; 
    int electronIndex = 0;

    int Dielectron_LeadInd = -1, Dielectron_SubLeadInd = -1; 
    n_dielectron_sel =  0;
    diele_ele1_ind = -1;  diele_ele2_ind = -1; double DeltaclosestZMass = 9999999.;
 
    for ( electrons = eleCollection->begin(); electrons != eleCollection->end(); ++electrons ) {
      if(TAODEBUG==1) cout<<"JTao : in the electron loops!"<<endl;
      n_electron++;

      //Employ loose selections of PAS-EGM-10-004   http://cds.cern.ch/record/1299116
      //The electron identification variables that have been found to be the most powerful, and are used in the selection, are: the energy-momentum match between the seed cluster and the track E_seed/p_in, the variables measuring spatial matching between the track and the supercluster, DeltaEta_in and DeltaPhi_in, the supercluster eta_width, sigma_ietaieta (as taken from the covariance matrix using logarithmic weights), and the hadronic leakage variable H/E.
      //Isolation variables are computed in three sub-detectors: the tracker, the ECAL, and the HCAL. Transverse energy/momentum sums are evaluated in regions of DR < 0.3.
      //"WP85" used  (WP95 and WP80, Cut-in-category Cic Loose (95% eff) and Cic Loose (80% eff) from W/Z)
      //https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
      double pt_ele = electrons->pt();
      double scEta_ele = electrons->superCluster()->position().Eta();
      if(  pt_ele > 20 &&  fabs(scEta_ele) < 2.5 && !(fabs(scEta_ele)>1.4442 && fabs(scEta_ele)<1.566) ){ //kinematic
 
        //double trackIso = electrons->dr03TkSumPt()           # calculated with   electronTrackIsolationScone_cfi.py
        //double ecalIso  = electrons->dr03EcalRecHitSumEt()   # calculated with   electronEcalRecHitIsolationScone_cfi.py
        //double hcalIso  = electrons->dr03HcalTowerSumEt()    # calculated with   electronHcalTowerIsolationScone_cfi.py
        double trackIsoRel = electrons->dr03TkSumPt()/electrons->pt();
	double ecalIsoRel  = electrons->dr03EcalRecHitSumEt()/electrons->pt();
	double hcalIsoRel  = electrons->dr03HcalTowerSumEt()/electrons->pt();
	//Combined Isolation as is defined with the following way:
	double RelCombinedIsoEB = ( electrons->dr03TkSumPt() + max(0., electrons->dr03EcalRecHitSumEt() - 1.) + electrons->dr03HcalTowerSumEt() ) / electrons->pt();
	double RelCombinedIsoEE = ( electrons->dr03TkSumPt() + electrons->dr03EcalRecHitSumEt() + electrons->dr03HcalTowerSumEt() ) / electrons->pt();
	//where the -1 is the pedestal subtraction and appears only in the barrel.
	//Shower shape:
	double HoE = electrons->hadronicOverEm();
	double SigIeIe = electrons->sigmaIetaIeta();
	//Track-cluster matching:
	double DeltaPhi = electrons->deltaPhiSuperClusterTrackAtVtx();
	double DeltaEta = electrons->deltaEtaSuperClusterTrackAtVtx();
	//Conversion rejection: Number of missing hits (the number of missing expected hits in front of the innermost valid hit) - If NumberOfExpectedInnerHits is greater than 1, then the electron is vetoed as from a converted photon and should be rejected in an analysis looking for prompt photons. 
         int Nmisshit = electrons->gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 
	//int Nmisshit = electrons->gsfTrack().trackerExpectedHitsInner().numberOfHits();
         // default value is -9999 if conversion partner not found
	//Minimum distance between conversion tracks:
	double DisConv = electrons->convDist();
	// \Delta cot \theta $ between conversion tracks at conversion vertex:
	double  DeltaCotTheta = electrons->convDcot();
       
        int ifpassedselection = 1; 
        if( Nmisshit > 1 ) ifpassedselection = 0;
	//if( fabs(DisConv) < 0.02 ) ifpassedselection = 0; 
        //if( fabs(DeltaCotTheta) < 0.02 ) ifpassedselection = 0;
        Bool_t isConv = fabs(DisConv) < 0.02 && fabs(DeltaCotTheta) < 0.02;
        if(isConv) ifpassedselection = 0;

        if ( fabs(scEta_ele) < 1.4442){ //EB
	  if ( RelCombinedIsoEB > 0.09 ) ifpassedselection = 0;
	  if ( trackIsoRel > 0.09 )      ifpassedselection = 0;
	  if ( ecalIsoRel > 0.08 )       ifpassedselection = 0;
	  if ( hcalIsoRel > 0.10 )       ifpassedselection = 0;
	  if ( SigIeIe > 0.01 )          ifpassedselection = 0;
	  if ( DeltaPhi > 0.06 )         ifpassedselection = 0;
	  if ( DeltaEta > 0.006 )        ifpassedselection = 0;
	  if ( HoE > 0.04 )              ifpassedselection = 0;
	} else{ //EE
	  if ( RelCombinedIsoEE > 0.06 ) ifpassedselection = 0;
	  if ( trackIsoRel > 0.05 )      ifpassedselection = 0;
	  if ( ecalIsoRel > 0.05 )       ifpassedselection = 0;
	  if ( hcalIsoRel > 0.025 )      ifpassedselection = 0;
	  if ( SigIeIe > 0.03 )          ifpassedselection = 0;
	  if ( DeltaPhi > 0.04 )         ifpassedselection = 0;
	  if ( DeltaEta > 0.007 )        ifpassedselection = 0;
	  if ( HoE > 0.025 )             ifpassedselection = 0;
	}

	if( ifpassedselection == 1 && pt_ele > maxElePT ){
	  SelectedEleIndex = electronIndex;
	  maxElePT = pt_ele;
	}

	//di-electron : same selections of PAS-EGM-10-004   http://cds.cern.ch/record/1299116
	if(ifpassedselection == 1){
          int electron2Ind = 0;
	  for ( reco::GsfElectronCollection::const_iterator electrons2 = eleCollection->begin(); electrons2 != eleCollection->end(); ++electrons2 ) {
	    if(electron2Ind != electronIndex){
	      double pt_ele2 = electrons2->pt();
	      double scEta_ele2 = electrons2->superCluster()->position().Eta();
              int QProd = electrons->charge()*electrons2->charge();
	      if(  pt_ele2 > 20 &&  fabs(scEta_ele2) < 2.5 && !(fabs(scEta_ele2)>1.4442 && fabs(scEta_ele2)<1.566) && QProd<0){ //kinematic
                n_dielectron_sel ++;
		TLorentzVector lead_p4; lead_p4.SetPtEtaPhiM(pt_ele, electrons->eta(), electrons->phi(), electrons->mass());
		TLorentzVector sublead_p4; sublead_p4.SetPtEtaPhiM(pt_ele2, electrons2->eta(), electrons2->phi(), electrons2->mass());
		TLorentzVector dielectron_p4 = lead_p4 + sublead_p4;
		double DeltaMass = fabs(dielectron_p4.M()-91.1876);
		if( DeltaMass < DeltaclosestZMass ){
		  DeltaclosestZMass = DeltaMass;
                  diele_ele1_ind = electronIndex;  diele_ele2_ind = electron2Ind;
		}
	      }
	    }
	    electron2Ind++;

	  }//ele2 loop
	}//ele1 passed the selection
      }//fulfill the kinematics
      electronIndex++;
    }// electron loop

    //=========single electron tree==============
    electronIndex = 0;
    for ( electrons = eleCollection->begin(); electrons != eleCollection->end(); ++electrons ) {
      //single electron
      if(TAODEBUG==1) cout<<"JTao : in the single electron tree filling!"<<endl;
      if(electronIndex == SelectedEleIndex){
        electron_pt = electrons->pt();
        electron_eta = electrons->eta();
        electron_phi = electrons->phi();
        electron_SCeta = electrons->superCluster()->position().Eta();
        electron_q = electrons->charge();
        electron_trckIso03Rel = electrons->dr03TkSumPt()/electrons->pt();
        electron_ecalIso03Rel  = electrons->dr03EcalRecHitSumEt()/electrons->pt();
	electron_hcalIso03Rel  = electrons->dr03HcalTowerSumEt()/electrons->pt();
	double RelCombinedIsoEB = ( electrons->dr03TkSumPt() + max(0., electrons->dr03EcalRecHitSumEt() - 1.) + electrons->dr03HcalTowerSumEt() ) / electrons->pt();
	double RelCombinedIsoEE = ( electrons->dr03TkSumPt() + electrons->dr03EcalRecHitSumEt() + electrons->dr03HcalTowerSumEt() ) / electrons->pt();
        electron_RelCombinedIso03 = RelCombinedIsoEB;
        if(fabs(electron_SCeta)>1.5)  electron_RelCombinedIso03 = RelCombinedIsoEE;
        
   	electron_hoe = electrons->hadronicOverEm();
	electron_sigieie = electrons->sigmaIetaIeta();
     
 	electron_DeltaPhiSCTrk = electrons->deltaPhiSuperClusterTrackAtVtx();
	electron_DeltaEtaSCTrk  = electrons->deltaEtaSuperClusterTrackAtVtx();
        electron_Nmisshit = electrons->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	electron_DisConv = fabs(electrons->convDist());
	electron_DeltaCotTheta = fabs(electrons->convDcot());
	electron_fbrem = electrons->fbrem();
        electron_ESCoPin = electrons->eSuperClusterOverP();
        electron_AbsInvEmInvPin = fabs(1.0/electrons->ecalEnergy() - 1.0/electrons->trackMomentumAtVtx().R());
	//Fill Single electron tree
        electronTree->Fill();

     }

      //di-electron
      if(electronIndex == diele_ele1_ind){       
	int electron2Ind = 0;
	for ( reco::GsfElectronCollection::const_iterator electrons2 = eleCollection->begin(); electrons2 != eleCollection->end(); ++electrons2 ) {
	  if( electron2Ind == diele_ele2_ind){

	    //ele 1
	    diele_ele1_pt = electrons->pt();   diele_ele1_eta = electrons->eta();  diele_ele1_phi = electrons->phi();  diele_ele1_e = electrons->energy();
            diele_ele1_SCeta = electrons->superCluster()->position().Eta();        diele_ele1_Escraw = electrons->superCluster()->rawEnergy();
            diele_ele1_q = electrons->charge();   
	    diele_ele1_trckIso03Rel = electrons->dr03TkSumPt()/electrons->pt();
	    diele_ele1_ecalIso03Rel  = electrons->dr03EcalRecHitSumEt()/electrons->pt();
	    diele_ele1_hcalIso03Rel  = electrons->dr03HcalTowerSumEt()/electrons->pt();
	    double RelCombinedIsoEB = ( electrons->dr03TkSumPt() + max(0., electrons->dr03EcalRecHitSumEt() - 1.) + electrons->dr03HcalTowerSumEt() ) / electrons->pt();
	    double RelCombinedIsoEE = ( electrons->dr03TkSumPt() + electrons->dr03EcalRecHitSumEt() + electrons->dr03HcalTowerSumEt() ) / electrons->pt();
	    diele_ele1_RelCombinedIso03 = RelCombinedIsoEB;
	    if(fabs(diele_ele1_SCeta)>1.5)  diele_ele1_RelCombinedIso03 = RelCombinedIsoEE;
	    diele_ele1_hoe = electrons->hadronicOverEm();
	    diele_ele1_sigieie = electrons->sigmaIetaIeta();

	    diele_ele1_DeltaPhiSCTrk = electrons->deltaPhiSuperClusterTrackAtVtx();
	    diele_ele1_DeltaEtaSCTrk  = electrons->deltaEtaSuperClusterTrackAtVtx();
	    diele_ele1_Nmisshit = electrons->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	    diele_ele1_DisConv = fabs(electrons->convDist());
	    diele_ele1_DeltaCotTheta = fabs(electrons->convDcot());
	    diele_ele1_fbrem = electrons->fbrem();
	    diele_ele1_ESCoPin = electrons->eSuperClusterOverP();
	    diele_ele1_AbsInvEmInvPin = fabs(1.0/electrons->ecalEnergy() - 1.0/electrons->trackMomentumAtVtx().R());
            
            //ele 2
 	    diele_ele2_pt = electrons2->pt();   diele_ele2_eta = electrons2->eta();  diele_ele2_phi = electrons2->phi();  diele_ele2_e = electrons2->energy();
            diele_ele2_SCeta = electrons2->superCluster()->position().Eta();        diele_ele2_Escraw = electrons2->superCluster()->rawEnergy();
            diele_ele2_q = electrons2->charge();   
	    diele_ele2_trckIso03Rel = electrons2->dr03TkSumPt()/electrons2->pt();
	    diele_ele2_ecalIso03Rel  = electrons2->dr03EcalRecHitSumEt()/electrons2->pt();
	    diele_ele2_hcalIso03Rel  = electrons2->dr03HcalTowerSumEt()/electrons2->pt();
	    double RelCombinedIso2EB = ( electrons2->dr03TkSumPt() + max(0., electrons2->dr03EcalRecHitSumEt() - 1.) + electrons2->dr03HcalTowerSumEt() ) / electrons2->pt();
	    double RelCombinedIso2EE = ( electrons2->dr03TkSumPt() + electrons2->dr03EcalRecHitSumEt() + electrons2->dr03HcalTowerSumEt() ) / electrons2->pt();
	    diele_ele2_RelCombinedIso03 = RelCombinedIso2EB;
	    if(fabs(diele_ele2_SCeta)>1.5)  diele_ele2_RelCombinedIso03 = RelCombinedIso2EE;
	    diele_ele2_hoe = electrons2->hadronicOverEm();
	    diele_ele2_sigieie = electrons2->sigmaIetaIeta();

	    diele_ele2_DeltaPhiSCTrk = electrons2->deltaPhiSuperClusterTrackAtVtx();
	    diele_ele2_DeltaEtaSCTrk  = electrons2->deltaEtaSuperClusterTrackAtVtx();
	    diele_ele2_Nmisshit = electrons2->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
	    diele_ele2_DisConv = fabs(electrons2->convDist());
	    diele_ele2_DeltaCotTheta = fabs(electrons2->convDcot());
	    diele_ele2_fbrem = electrons2->fbrem();
	    diele_ele2_ESCoPin = electrons2->eSuperClusterOverP();
	    diele_ele2_AbsInvEmInvPin = fabs(1.0/electrons2->ecalEnergy() - 1.0/electrons2->trackMomentumAtVtx().R());

	    //di-electron
	    TLorentzVector lead_p4; lead_p4.SetPtEtaPhiM(diele_ele1_pt, diele_ele1_eta, diele_ele1_phi, electrons->mass());
	    TLorentzVector sublead_p4; sublead_p4.SetPtEtaPhiM(diele_ele2_pt, diele_ele2_eta, diele_ele2_phi, electrons2->mass());
	    TLorentzVector dielectron_p4 = lead_p4 + sublead_p4;
	    diele_pt = dielectron_p4.Pt();  diele_eta = dielectron_p4.Eta();   diele_phi = dielectron_p4.Phi();  diele_energy = dielectron_p4.E();   diele_mass = dielectron_p4.M();         
            diele_DeltaR = twoobj_deltaR(diele_ele1_eta, diele_ele1_phi, diele_ele2_eta, diele_ele2_phi);
	    //=============           
	    dielectronTree->Fill();
	  }
	  electron2Ind++;
	}//find ele2

      }//find ele1

      electronIndex++;
    }

  }  //electron valid

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonElectronAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonElectronAnalyzer::endJob() {
  std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "analyzed " << nEvt << " events " << std::endl;
  std::cout << "writing information into file: " << outputfile->GetName() << std::endl;
  outputfile->Write();
  outputfile->Close();
}



//define this as a plug-in
DEFINE_FWK_MODULE(PhotonElectronAnalyzer);
