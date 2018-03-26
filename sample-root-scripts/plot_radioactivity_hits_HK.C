#include <iostream>
#include <TH1F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TStyle.h"


int main(){

  TString RadioactivityDirec = Form("outputs/radioactivity/");
  //  TString RadioactivityDirec = Form("outputs/SuperKamiokandeRadioactivity_3.5/");
  //TString RadioactivityDirec = Form("outputs/radioactivity_BiWater1event/");
  TString RadioactivityFile= Form("output.root");
  TString ElectronsDirec = Form("outputs/eminus_3MeV_noise_random_4pi_000_40per/");
  TString ElectronsFile= Form("output.root");

  TFile output("radioactivity_hits.root","RECREATE");

  TFile Radioactivityinput((RadioactivityDirec + RadioactivityFile).Data(),"READ");
  TFile Electronsinput((ElectronsDirec + ElectronsFile).Data(),"READ");

  std::vector<Int_t> * trigger_number_digitized_hits = 0;
  std::vector<std::vector<Float_t> > * digitized_hit_Q = 0;

  double maxQ;
  double nhits_min=36000;
  double nhits_max=41000;

  TTree * Radioactivity_primary_events_tree = (TTree*)Radioactivityinput.Get("primary_events_tree");
  Radioactivity_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits); 
  Radioactivity_primary_events_tree->SetBranchAddress("digitized_hit_Q",&digitized_hit_Q);
  TH1F *h_Radioactivity_trigger_number_digitized_hits = new TH1F("h_Radioactivity_trigger_number_digitized_hits","radioactivity",100,nhits_min,nhits_max);
  h_Radioactivity_trigger_number_digitized_hits->SetLineColor(kRed);
  h_Radioactivity_trigger_number_digitized_hits->SetLineStyle(1);
  h_Radioactivity_trigger_number_digitized_hits->SetLineWidth(2);
  h_Radioactivity_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
  TH1F *h_Radioactivity_charge_hit = new TH1F("h_Radioactivity_charge_hit","radioactivity",100,0,250);
  h_Radioactivity_charge_hit->SetLineColor(kRed);
  h_Radioactivity_charge_hit->SetLineStyle(1);
  h_Radioactivity_charge_hit->SetLineWidth(2);
  h_Radioactivity_charge_hit->GetXaxis()->SetTitle("charge per hit");
  TH1F *h_Radioactivity_charge_max_hit = new TH1F("h_Radioactivity_charge_max_hit","radioactivity",100,0,250);
  h_Radioactivity_charge_max_hit->SetLineColor(kRed);
  h_Radioactivity_charge_max_hit->SetLineStyle(1);
  h_Radioactivity_charge_max_hit->SetLineWidth(2);
  h_Radioactivity_charge_max_hit->GetXaxis()->SetTitle("maximum charge per hit");
  for(int ievent=0; ievent<Radioactivity_primary_events_tree->GetEntries(); ievent++){
    Radioactivity_primary_events_tree->GetEvent(ievent); 
    for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
      h_Radioactivity_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
      maxQ = 0.;
      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_Q->at(itrigger)).size(); idigitizedhit++){
	h_Radioactivity_charge_hit->Fill((digitized_hit_Q->at(itrigger)).at(idigitizedhit));
      	if( (digitized_hit_Q->at(itrigger)).at(idigitizedhit) > maxQ ) maxQ = (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
      }
      h_Radioactivity_charge_max_hit->Fill(maxQ);
    }
  }
  h_Radioactivity_trigger_number_digitized_hits->Scale(1./h_Radioactivity_trigger_number_digitized_hits->Integral());
  h_Radioactivity_charge_max_hit->Scale(1./h_Radioactivity_charge_max_hit->Integral());
  h_Radioactivity_charge_hit->Scale(1./h_Radioactivity_charge_hit->Integral());


  TTree * Electrons_primary_events_tree = (TTree*)Electronsinput.Get("primary_events_tree");
  Electrons_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits);
  Electrons_primary_events_tree->SetBranchAddress("digitized_hit_Q",&digitized_hit_Q);
  TH1F *h_Electrons_trigger_number_digitized_hits = new TH1F("h_Electrons_trigger_number_digitized_hits","3 MeV electrons",100,nhits_min,nhits_max);
  h_Electrons_trigger_number_digitized_hits->SetLineColor(kBlue);
  h_Electrons_trigger_number_digitized_hits->SetLineStyle(1);
  h_Electrons_trigger_number_digitized_hits->SetLineWidth(2);
  h_Electrons_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
  TH1F *h_Electrons_charge_hit = new TH1F("h_Electrons_charge_hit","3 MeV electrons",38,0,76);
  h_Electrons_charge_hit->SetLineColor(kBlue);
  h_Electrons_charge_hit->SetLineStyle(1);
  h_Electrons_charge_hit->SetLineWidth(2);
  h_Electrons_charge_hit->GetXaxis()->SetTitle("charge per hit");
  TH1F *h_Electrons_charge_max_hit = new TH1F("h_Electrons_charge_max_hit","3 MeV electrons",38,0,76);
  h_Electrons_charge_max_hit->SetLineColor(kBlue);
  h_Electrons_charge_max_hit->SetLineStyle(1);
  h_Electrons_charge_max_hit->SetLineWidth(2);
  h_Electrons_charge_max_hit->GetXaxis()->SetTitle("maximum charge per hit");
  for(int ievent=0; ievent<Electrons_primary_events_tree->GetEntries(); ievent++){
    Electrons_primary_events_tree->GetEvent(ievent); 
    for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
      h_Electrons_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
      maxQ = 0.;
      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_Q->at(itrigger)).size(); idigitizedhit++){
	h_Electrons_charge_hit->Fill((digitized_hit_Q->at(itrigger)).at(idigitizedhit));
      	if( (digitized_hit_Q->at(itrigger)).at(idigitizedhit) > maxQ ) maxQ = (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
      }
      h_Electrons_charge_max_hit->Fill(maxQ);
    }
  }
  h_Electrons_trigger_number_digitized_hits->Scale(1./h_Electrons_trigger_number_digitized_hits->Integral());
  h_Electrons_charge_max_hit->Scale(1./h_Electrons_charge_max_hit->Integral());
  h_Electrons_charge_hit->Scale(1./h_Electrons_charge_hit->Integral());



  output.cd();
  h_Radioactivity_trigger_number_digitized_hits->Write();
  h_Electrons_trigger_number_digitized_hits->Write();
  h_Radioactivity_charge_max_hit->Write();
  h_Radioactivity_charge_hit->Write();
  h_Electrons_charge_max_hit->Write();
  h_Electrons_charge_hit->Write();
  output.Write();

  TCanvas canvas("canvas", "canvas", 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend legend(0.5,0.5,0.85,0.85);
  legend.SetFillColor(0);

  h_Electrons_trigger_number_digitized_hits->Draw();
  h_Radioactivity_trigger_number_digitized_hits->Draw("same");
  legend.AddEntry(h_Radioactivity_trigger_number_digitized_hits, h_Radioactivity_trigger_number_digitized_hits->GetTitle(), "l");
  legend.AddEntry(h_Electrons_trigger_number_digitized_hits, h_Electrons_trigger_number_digitized_hits->GetTitle(), "l");

  legend.Draw("same");

  canvas.Print("n-hits-radioactivity-electrons.png");

  std::clog << " n digitized hits from radioactivity " << h_Radioactivity_trigger_number_digitized_hits->GetMean() << std::endl;
  std::clog << " n digitized hits from 3 MeV electrons " << h_Electrons_trigger_number_digitized_hits->GetMean() << std::endl;


  TLegend legend2(0.5,0.5,0.85,0.85);
  legend2.SetFillColor(0);

  h_Radioactivity_charge_max_hit->Draw();
  h_Radioactivity_charge_max_hit->GetYaxis()->SetRangeUser(1.e-2,1.1);
  canvas.SetLogy(1);
  legend2.AddEntry(h_Radioactivity_charge_max_hit, h_Radioactivity_charge_max_hit->GetTitle(), "l");
  h_Electrons_charge_max_hit->Draw("same");
  legend2.AddEntry(h_Electrons_charge_max_hit, h_Electrons_charge_max_hit->GetTitle(), "l");

  legend2.Draw("same");

  canvas.Print("max_charge_radioactivity.png");


  exit(0);
  return 1;
}



