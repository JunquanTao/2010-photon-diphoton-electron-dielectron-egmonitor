//#include "setTDRStyle.C"
#include "TMath.h"
//#include "TROOT.h"
#include <TH1D.h>
#include <TH1D.h>
#include <TH1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLeafF.h>
#include <TChain.h>
#include <TFile.h>
#include "TSystem.h"
#include <TChain.h>
#include "TSystem.h"
#include <TString.h>
#include <iostream>
#include <vector>
#include <TPostScript.h>
#include <iostream>
#include <iomanip>  //for precision

//==============
const int debug=1;
const string TreeName = "photonTree";
const string PlotPreName = "Photon2010B_";
const string PlotLabel = "Photon 2010B";
const string InputFile = "/eos/cms/store/user/jtao/CMSOpenData2010/Photon_2010B-Apr21ReReco.root";

const TString PrintInfor1="#bf{CMS} #it{} #it{Preliminary}";
//const TString PrintInfor2="2.9 pb^{-1} (7 TeV)";
const TString PrintInfor2="2010 (7 TeV)";

const TString PrintEtaEB="|#eta|<1.4442";
const TString PrintEtaEE="|#eta|>1.566";


void DrawMyPlots(string Object, string Selections,  string XTitle, string YUnit, string PlotName, int BinTotal, double BinXLow, double BinXHig, int LegendLR=0, int IfLogY=0, int IfEB=-1){

  TCanvas *c1 = new TCanvas("reconstruction1","reconstruction1");
  c1->SetFillColor(0);

  TLegend *legend;
  if(LegendLR==0) legend = new TLegend(0.65,0.75,0.9,0.85);
  else legend = new TLegend(0.24,0.75,0.49,0.85);
  legend->SetFillColor(0);

  //====add root file=============
  TChain *Data_Tree=new TChain(TreeName.c_str());
  // Data_Tree->Add(Form("%s/*.root", InputDir.c_str()));
  Data_Tree->Add(InputFile.c_str());

  //=========entries================
  int entries_Data = Data_Tree->GetEntries();
  if(debug==1) cout <<"JTao: nEntries_Data = "<<entries_Data<<endl;

  //=========
  c1->cd();

  //=================
  char *myLimits= new char[100];
  sprintf(myLimits,"(%d,%f,%f)",BinTotal,BinXLow,BinXHig);
  TString Taolimits(myLimits);

  //====data=======
  TString variable_Data = Object + ">>Histo_Data_temp" + Taolimits;
  Data_Tree->Draw(variable_Data, Selections.c_str());
  TH1D *h_data = (TH1D*)gDirectory->Get("Histo_Data_temp");
  h_data->SetTitle("");
  c1->Clear();


  double Ntot_Data=h_data->Integral();
  if( debug==1 ) cout<<"JTao: N_Data= "<<Ntot_Data<<endl;

  h_data->SetLineColor(1);
  h_data->SetFillColor(0);
  h_data->SetLineStyle(1);
  h_data->SetLineWidth(2);
  double WidthBin=(BinXHig-BinXLow)/BinTotal;
  string PreTitleY( Form("Events / %.2g ",WidthBin) );
  string TitleY =  PreTitleY + YUnit;
  h_data->GetYaxis()->SetTitle(TitleY.c_str());
  h_data->GetXaxis()->SetTitle(XTitle.c_str());

  h_data->SetTitleSize(0.05,"X");
  h_data->SetTitleSize(0.05,"Y");
  h_data->SetTitleOffset(1.6, "Y");
  //h_data->SetTitleOffset(1.1, "Y");

  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerSize(1.0);
  h_data->SetMarkerStyle(20);

  legend->AddEntry(h_data,PlotLabel.c_str(),"pe");
  //h_data->GetXaxis()->SetLabelColor(0);
  h_data->SetNdivisions(510 ,"X");

  double maxY=h_data->GetMaximum();
  double minY=h_data->GetMinimum();
  if(IfLogY==1 && minY<1) minY = 1.;
  h_data->GetYaxis()->SetRangeUser(0.9*minY, 1.05*maxY);
  if(IfLogY==1) h_data->GetYaxis()->SetRangeUser(0.1*minY, 5.0*maxY);

  gPad->SetLogy(IfLogY);
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  h_data->Draw("PE1");
  legend->Draw("same");
  h_data->Draw("samePE1");
  h_data->Draw("Axissame");

  TLatex *tex1 = new TLatex(0.16,0.94, PrintInfor1);
  tex1->SetNDC();
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.045);
  tex1->SetLineWidth(2);
  tex1->Draw();

  TLatex *tex2 = new TLatex(0.62,0.933, PrintInfor2);
  tex2 = new TLatex(0.72,0.94, PrintInfor2);
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.045);
  tex2->SetLineWidth(2);
  tex2->Draw();

  if(IfEB>=0){
    TLatex *texEta = new TLatex(0.22,0.65, PrintEtaEB);
    if(IfEB==0) texEta = new TLatex(0.22,0.65, PrintEtaEE);
    if(LegendLR==0) {
      texEta = new TLatex(0.72,0.65, PrintEtaEB);
      if(IfEB==0) texEta = new TLatex(0.72,0.65, PrintEtaEE);
    }
    if(LegendLR==-1) {
      texEta = new TLatex(0.46,0.65, PrintEtaEB);
      if(IfEB==0) texEta = new TLatex(0.46,0.65, PrintEtaEE);
    }

    texEta->SetNDC();
    texEta->SetTextFont(42);
    texEta->SetTextSize(0.045);
    texEta->SetLineWidth(2);
    if(IfEB>-1) texEta->Draw();
  }

  //
  string nameplots=PlotPreName+PlotName+".png";
  c1->Print(nameplots.c_str());

  string nameplotspdf=PlotPreName+PlotName+".pdf";
  c1->Print(nameplotspdf.c_str());

}


void DrawSinglePhoton(){

  gROOT->ProcessLine(".x hggPaperStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);

  //string BasicSelections="event_PassHLT_SinglePhoton == 1 && run <= 144116"; //2.9pb-1
  string BasicSelections="event_PassHLT_SinglePhoton == 1"; //2.9pb-1

  DrawMyPlots("photon_pt", BasicSelections, "p_{T} of #gamma (GeV)", "GeV",  "photon_pt", 60, 18, 78);
  DrawMyPlots("photon_SCeta", BasicSelections,"#eta_{SC}", "", "photon_SCeta",60, -3, 3);
  DrawMyPlots("photon_r9", BasicSelections, "R_{9} of #gamma", "", "photon_r9", 81, 0.2, 1.01, 1);
  DrawMyPlots("photon_hoe", BasicSelections, "H/E of #gamma", "", "photon_hoe", 120, 0.0, 0.06, 0, 1);
  DrawMyPlots("photon_sieie", BasicSelections,"#sigma_{i#etai#eta}", "", "photon_sieie", 60, 0.005, 0.035);
  DrawMyPlots("photon_TrkIsoHollow04", BasicSelections, "Track Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_TrkIsoHollow04", 110, -0.1, 2.1, 0, 1);
  DrawMyPlots("photon_EcalIso04", BasicSelections, "ECAL Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_EcalIso04", 90, -0.1, 4.4);
  DrawMyPlots("photon_HcalIso04", BasicSelections, "HCAL Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_HcalIso04", 50, -0.1, 2.4, 0, 1);
    
  string BasicSelectionsEB="event_PassHLT_SinglePhoton == 1 && fabs(photon_SCeta)<1.5"; //2.9pb-1
  DrawMyPlots("photon_pt", BasicSelectionsEB, "p_{T} of #gamma (GeV)", "GeV",  "photon_pt_EB", 60, 18, 78, 0, 0, 1);
  DrawMyPlots("photon_r9", BasicSelectionsEB, "R_{9} of #gamma", "", "photon_r9_EB", 81, 0.2, 1.01, 1, 0, 1);
  DrawMyPlots("photon_hoe", BasicSelectionsEB, "H/E of #gamma", "", "photon_hoe_EB", 120, 0.0, 0.06, 0, 1, 1);
  DrawMyPlots("photon_sieie", BasicSelectionsEB,"#sigma_{i#etai#eta}", "", "photon_sieie_EB", 60, 0.0, 0.03, 0, 0, 1);
  DrawMyPlots("photon_TrkIsoHollow04", BasicSelectionsEB, "Track Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_TrkIsoHollow04_EB", 110, -0.1, 2.1, 0, 1, 1);
  DrawMyPlots("photon_EcalIso04", BasicSelectionsEB, "ECAL Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_EcalIso04_EB", 90, -0.1, 4.4, 0, 0, 1);
  DrawMyPlots("photon_HcalIso04", BasicSelectionsEB, "HCAL Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_HcalIso04_EB", 50, -0.1, 2.4, 0, 1, 1);

  string BasicSelectionsEE="event_PassHLT_SinglePhoton == 1 && fabs(photon_SCeta)>1.5"; //2.9pb-1
  DrawMyPlots("photon_pt", BasicSelectionsEE, "p_{T} of #gamma (GeV)", "GeV",  "photon_pt_EE", 60, 18, 78, 0, 0, 0);
  DrawMyPlots("photon_r9", BasicSelectionsEE, "R_{9} of #gamma", "", "photon_r9_EE", 81, 0.2, 1.01, 1, 0, 0);
  DrawMyPlots("photon_hoe", BasicSelectionsEE, "H/E of #gamma", "", "photon_hoe_EE", 120, 0.0, 0.06, 0, 1, 0);
  DrawMyPlots("photon_sieie", BasicSelectionsEE,"#sigma_{i#etai#eta}", "", "photon_sieie_EE", 100, 0.01, 0.06, 0, 0, 0);
  DrawMyPlots("photon_TrkIsoHollow04", BasicSelectionsEE, "Track Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_TrkIsoHollow04_EE", 110, -0.1, 2.1, 0, 1, 0);
  DrawMyPlots("photon_EcalIso04", BasicSelectionsEE, "ECAL Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_EcalIso04_EE", 90, -0.1, 4.4, 0, 0, 0);
  DrawMyPlots("photon_HcalIso04", BasicSelectionsEE, "HCAL Iso (#DeltaR=0.4) (GeV)", "GeV", "photon_HcalIso04_EE", 50, -0.1, 2.4, 0, 1, 0);
}
