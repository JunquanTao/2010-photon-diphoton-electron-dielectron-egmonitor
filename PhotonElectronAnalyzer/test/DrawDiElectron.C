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
const string TreeName = "dielectronTree";
const string PlotPreName = "Electron2010B_"; //"Photon2010B_";
const string PlotLabel = "Electron 2010B"; //"Photon 2010B";
//const string InputFile = "/eos/cms/store/user/jtao/CMSOpenData2010/Photon_2010B-Apr21ReReco.root";
const string InputFile = "/eos/cms/store/user/jtao/CMSOpenData2010/Electron_Run2010B-Apr21ReReco.root";

const TString PrintInfor1="#bf{CMS} #it{} #it{Preliminary}";
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
  //h_data->SetNdivisions(510 ,"X");
  h_data->SetNdivisions(505 ,"X");

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


void DrawDiElectron(){

  gROOT->ProcessLine(".x hggPaperStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);

  string BasicSelections="diele_pt>0.";  //no effect

  DrawMyPlots("diele_pt", BasicSelections, "p_{T} (GeV)", "GeV",  "diele_pt", 100, 20, 120);
  DrawMyPlots("diele_mass", BasicSelections, "M_{ee} (GeV)", "GeV",  "diele_mass", 100, 20, 120);
  DrawMyPlots("diele_DeltaR", BasicSelections, "#DeltaR_{e,e}", "",  "diele_DeltaR", 70, 0,7.0);

  DrawMyPlots("diele_ele1_pt", BasicSelections, "p_{T} (GeV)", "GeV",  "diele_ele1_pt", 60, 18, 78);
  DrawMyPlots("diele_ele1_SCeta", BasicSelections,"#eta_{SC}", "", "diele_ele1_SCeta",60, -3, 3);
  DrawMyPlots("diele_ele1_hoe", BasicSelections, "H/E", "", "diele_ele1_hoe", 100, 0.0, 0.05, 0, 1);
  DrawMyPlots("diele_ele1_sigieie", BasicSelections,"#sigma_{i#etai#eta}", "", "diele_ele1_sigieie", 60, 0.005, 0.035);
 
  DrawMyPlots("diele_ele1_trckIso03Rel", BasicSelections, "Track Iso (#DeltaR=0.3) / p_{T}^{e}", "", "diele_ele1_trckIso03Rel", 100, 0., 0.1, 0, 1);
  DrawMyPlots("diele_ele1_ecalIso03Rel", BasicSelections, "ECAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "diele_ele1_ecalIso03Rel", 120, -0.02, 0.1);
  DrawMyPlots("diele_ele1_hcalIso03Rel", BasicSelections, "HCAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "diele_ele1_hcalIso03Rel", 100,  0., 0.1, 0, 1);
  DrawMyPlots("diele_ele1_RelCombinedIso03", BasicSelections, "Combined Iso (#DeltaR=0.3) / p_{T}^{e}", "", "diele_ele1_RelCombinedIso03", 100,  0., 0.1, 0, 1);
  DrawMyPlots("diele_ele1_DeltaEtaSCTrk", BasicSelections, "#Delta#eta(SC, trk at vtx)", "", "diele_ele1_DeltaEtaSCTrk", 100,  0., 0.01);
  DrawMyPlots("diele_ele1_DeltaPhiSCTrk", BasicSelections, "#Delta#phi(SC, trk at vtx)", "", "diele_ele1_DeltaPhiSCTrk", 100,  0., 0.08);
 
  DrawMyPlots("diele_ele1_Nmisshit", BasicSelections, "Number of missing hits", "", "diele_ele1_Nmisshit", 3,  -0.5, 2.5);
  DrawMyPlots("diele_ele1_DisConv", BasicSelections, "Minimum distance between conversion tracks", "", "diele_ele1_DisConv", 60,  -0.03, 0.03);
  DrawMyPlots("diele_ele1_DeltaCotTheta", BasicSelections, "#Deltacot#theta between conversion tracks at conversion vertex:", "", "diele_ele1_DeltaCotTheta", 80,  -0.02, 0.02);

 
}
