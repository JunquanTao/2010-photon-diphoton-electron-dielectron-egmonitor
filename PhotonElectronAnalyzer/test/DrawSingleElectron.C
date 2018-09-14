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
const string TreeName = "electronTree";
const string PlotPreName = "Photon2010B_";
const string PlotLabel = "Photon 2010B";
const string InputFile = "/eos/cms/store/user/jtao/CMSOpenData2010/Photon_2010B-Apr21ReReco.root";

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


void DrawSingleElectron(){

  gROOT->ProcessLine(".x hggPaperStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);

  string BasicSelections="electron_pt>10.";  //no effect

  DrawMyPlots("electron_pt", BasicSelections, "p_{T} (GeV)", "GeV",  "electron_pt", 60, 18, 78);
  DrawMyPlots("electron_SCeta", BasicSelections,"#eta_{SC}", "", "electron_SCeta",60, -3, 3);
  DrawMyPlots("electron_hoe", BasicSelections, "H/E", "", "electron_hoe", 100, 0.0, 0.05, 0, 1);
  DrawMyPlots("electron_sigieie", BasicSelections,"#sigma_{i#etai#eta}", "", "electron_sigieie", 60, 0.005, 0.035);
 
  DrawMyPlots("electron_trckIso03Rel", BasicSelections, "Track Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_trckIso03Rel", 100, 0., 0.1, 0, 1);
  DrawMyPlots("electron_ecalIso03Rel", BasicSelections, "ECAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_ecalIso03Rel", 120, -0.02, 0.1);
  DrawMyPlots("electron_hcalIso03Rel", BasicSelections, "HCAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_hcalIso03Rel", 100,  0., 0.1, 0, 1);
  DrawMyPlots("electron_RelCombinedIso03", BasicSelections, "Combined Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_RelCombinedIso03", 100,  0., 0.1, 0, 1);
  DrawMyPlots("electron_DeltaEtaSCTrk", BasicSelections, "#Delta#eta(SC, trk at vtx)", "", "electron_DeltaEtaSCTrk", 100,  0., 0.01);
  DrawMyPlots("electron_DeltaPhiSCTrk", BasicSelections, "#Delta#phi(SC, trk at vtx)", "", "electron_DeltaPhiSCTrk", 100,  0., 0.08);
 
  DrawMyPlots("electron_Nmisshit", BasicSelections, "Number of missing hits", "", "electron_Nmisshit", 3,  -0.5, 2.5);
  DrawMyPlots("electron_DisConv", BasicSelections, "Minimum distance between conversion tracks", "", "electron_DisConv", 60,  -0.03, 0.03);
  DrawMyPlots("electron_DeltaCotTheta", BasicSelections, "#Deltacot#theta between conversion tracks at conversion vertex:", "", "electron_DeltaCotTheta", 80,  -0.02, 0.02);

  string BasicSelectionsEB="electron_pt>10. && fabs(electron_SCeta)<1.5"; //EB
  DrawMyPlots("electron_pt", BasicSelectionsEB, "p_{T} (GeV)", "GeV",  "electron_pt_EB", 60, 18, 78);
  DrawMyPlots("electron_SCeta", BasicSelectionsEB,"#eta_{SC}", "", "electron_SCeta_EB",60, -3, 3);
  DrawMyPlots("electron_hoe", BasicSelectionsEB, "H/E", "", "electron_hoe_EB", 100, 0.0, 0.05, 0, 1);
  DrawMyPlots("electron_sigieie", BasicSelectionsEB,"#sigma_{i#etai#eta}", "", "electron_sigieie_EB", 60, 0.005, 0.035);
 
  DrawMyPlots("electron_trckIso03Rel", BasicSelectionsEB, "Track Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_trckIso03Rel_EB", 100, 0., 0.1, 0, 1);
  DrawMyPlots("electron_ecalIso03Rel", BasicSelectionsEB, "ECAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_ecalIso03Rel_EB", 120, -0.02, 0.1);
  DrawMyPlots("electron_hcalIso03Rel", BasicSelectionsEB, "HCAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_hcalIso03Rel_EB", 100,  0., 0.1, 0, 1);
  DrawMyPlots("electron_RelCombinedIso03", BasicSelectionsEB, "Combined Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_RelCombinedIso03_EB", 100,  0., 0.1, 0, 1);
  DrawMyPlots("electron_DeltaEtaSCTrk", BasicSelectionsEB, "#Delta#eta(SC, trk at vtx)", "", "electron_DeltaEtaSCTrk_EB", 100,  0., 0.01);
  DrawMyPlots("electron_DeltaPhiSCTrk", BasicSelectionsEB, "#Delta#phi(SC, trk at vtx)", "", "electron_DeltaPhiSCTrk_EB", 100,  0., 0.08);
 
  DrawMyPlots("electron_Nmisshit", BasicSelectionsEB, "Number of missing hits", "", "electron_Nmisshit_EB", 3,  -0.5, 2.5);


  string BasicSelectionsEE="electron_pt>10. && fabs(electron_SCeta)>1.5"; //EE
  DrawMyPlots("electron_pt", BasicSelectionsEE, "p_{T} (GeV)", "GeV",  "electron_pt_EE", 60, 18, 78);
  DrawMyPlots("electron_SCeta", BasicSelectionsEE,"#eta_{SC}", "", "electron_SCeta_EE",60, -3, 3);
  DrawMyPlots("electron_hoe", BasicSelectionsEE, "H/E", "", "electron_hoe_EE", 100, 0.0, 0.05, 0, 1);
  DrawMyPlots("electron_sigieie", BasicSelectionsEE,"#sigma_{i#etai#eta}", "", "electron_sigieie_EE", 60, 0.005, 0.035);
 
  DrawMyPlots("electron_trckIso03Rel", BasicSelectionsEE, "Track Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_trckIso03Rel_EE", 100, 0., 0.1, 0, 1);
  DrawMyPlots("electron_ecalIso03Rel", BasicSelectionsEE, "ECAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_ecalIso03Rel_EE", 120, -0.02, 0.1);
  DrawMyPlots("electron_hcalIso03Rel", BasicSelectionsEE, "HCAL Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_hcalIso03Rel_EE", 100,  0., 0.1, 0, 1);
  DrawMyPlots("electron_RelCombinedIso03", BasicSelectionsEE, "Combined Iso (#DeltaR=0.3) / p_{T}^{e}", "", "electron_RelCombinedIso03_EE", 100,  0., 0.1, 0, 1);
  DrawMyPlots("electron_DeltaEtaSCTrk", BasicSelectionsEE, "#Delta#eta(SC, trk at vtx)", "", "electron_DeltaEtaSCTrk_EE", 100,  0., 0.01);
  DrawMyPlots("electron_DeltaPhiSCTrk", BasicSelectionsEE, "#Delta#phi(SC, trk at vtx)", "", "electron_DeltaPhiSCTrk_EE", 100,  0., 0.08);
 
  DrawMyPlots("electron_Nmisshit", BasicSelectionsEE, "Number of missing hits", "", "electron_Nmisshit_EE", 3,  -0.5, 2.5);
}
