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

  TString direc = Form("radioactivity_e61_detector_3inch/");
  TString BiPMTFile= Form("output_Bi214_PMT_1event_e61_detector.root");
  TString TlPMTFile= Form("output_Tl208_PMT_1event_e61_detector.root");
  TString KPMTFile= Form("output_K40_PMT_1event_e61_detector.root");
  TString BiWaterFile= Form("output_Bi214_Water_1event_e61_detector.root");
  TString NeutronFile= Form("output_neutron.root");

  TFile output("radioactivity_hits.root","RECREATE");

  TFile BiPMTinput((direc + BiPMTFile).Data(),"READ");
  TFile TlPMTinput((direc + TlPMTFile).Data(),"READ");
  TFile KPMTinput((direc + KPMTFile).Data(),"READ");
  TFile BiWaterinput((direc + BiWaterFile).Data(),"READ");
  TFile Neutroninput((direc + NeutronFile).Data(),"READ");

  std::vector<Int_t> * trigger_number_digitized_hits = 0;
  std::vector<std::vector<Int_t> > * digitized_hit_tube_id = 0;
  int nentries = 14000;
  int nhits_max_pmt[14000];

  int nbins = 20;

  TTree * BiPMT_primary_events_tree = (TTree*)BiPMTinput.Get("primary_events_tree");
  BiPMT_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits); 
  BiPMT_primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  TH1F *h_BiPMT_trigger_number_digitized_hits = new TH1F("h_BiPMT_trigger_number_digitized_hits","Bi^{214} on PMT",nbins,1,-1);
  h_BiPMT_trigger_number_digitized_hits->SetLineColor(kRed);
  h_BiPMT_trigger_number_digitized_hits->SetLineStyle(1);
  h_BiPMT_trigger_number_digitized_hits->SetLineWidth(2);
  h_BiPMT_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
  TH1F *h_BiPMT_nhits_max_pmt = new TH1F("h_BiPMT_nhits_max_pmt","Bi^{214} on PMT",nbins,1,-1);
  h_BiPMT_nhits_max_pmt->SetLineColor(kRed);
  h_BiPMT_nhits_max_pmt->SetLineStyle(1);
  h_BiPMT_nhits_max_pmt->SetLineWidth(2);
  h_BiPMT_nhits_max_pmt->GetXaxis()->SetTitle("maximum n hits per PMT");
  for(int ievent=0; ievent<BiPMT_primary_events_tree->GetEntries(); ievent++){
    BiPMT_primary_events_tree->GetEvent(ievent); 
    for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
      h_BiPMT_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
      for(int i = 0; i < nentries; ++i) nhits_max_pmt[i] = 0;
      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_tube_id->at(itrigger)).size(); idigitizedhit++){
      	nhits_max_pmt[(digitized_hit_tube_id->at(itrigger)).at(idigitizedhit)]++;
      }
      for(int i = 1; i < nentries; ++i) if(nhits_max_pmt[0] < nhits_max_pmt[i]) nhits_max_pmt[0] = nhits_max_pmt[i];
      h_BiPMT_nhits_max_pmt->Fill(nhits_max_pmt[0]);
    }
  }
  h_BiPMT_trigger_number_digitized_hits->Scale(1./h_BiPMT_trigger_number_digitized_hits->Integral(3,100));


  TTree * TlPMT_primary_events_tree = (TTree*)TlPMTinput.Get("primary_events_tree");
  TlPMT_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits);
  TlPMT_primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  TH1F *h_TlPMT_trigger_number_digitized_hits = new TH1F("h_TlPMT_trigger_number_digitized_hits","Tl^{208} on PMT",nbins,1,-1);
  h_TlPMT_trigger_number_digitized_hits->SetLineColor(kBlue);
  h_TlPMT_trigger_number_digitized_hits->SetLineStyle(1);
  h_TlPMT_trigger_number_digitized_hits->SetLineWidth(2);
  h_TlPMT_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
  TH1F *h_TlPMT_nhits_max_pmt = new TH1F("h_TlPMT_nhits_max_pmt","Tl^{208} on PMT",nbins,1,-1);
  h_TlPMT_nhits_max_pmt->SetLineColor(kRed);
  h_TlPMT_nhits_max_pmt->SetLineStyle(1);
  h_TlPMT_nhits_max_pmt->SetLineWidth(2);
  h_TlPMT_nhits_max_pmt->GetXaxis()->SetTitle("maximum n hits per PMT");
  for(int ievent=0; ievent<TlPMT_primary_events_tree->GetEntries(); ievent++){
    TlPMT_primary_events_tree->GetEvent(ievent); 
    for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
      h_TlPMT_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
      for(int i = 0; i < nentries; ++i) nhits_max_pmt[i] = 0;
      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_tube_id->at(itrigger)).size(); idigitizedhit++){
      	nhits_max_pmt[(digitized_hit_tube_id->at(itrigger)).at(idigitizedhit)]++;
      }
      for(int i = 1; i < nentries; ++i) if(nhits_max_pmt[0] < nhits_max_pmt[i]) nhits_max_pmt[0] = nhits_max_pmt[i];
      h_TlPMT_nhits_max_pmt->Fill(nhits_max_pmt[0]);
    }
  }
  h_TlPMT_trigger_number_digitized_hits->Scale(1./h_TlPMT_trigger_number_digitized_hits->Integral(3,100));


  TTree * KPMT_primary_events_tree = (TTree*)KPMTinput.Get("primary_events_tree");
  KPMT_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits);
  KPMT_primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  TH1F *h_KPMT_trigger_number_digitized_hits = new TH1F("h_KPMT_trigger_number_digitized_hits","K^{40} on PMT",nbins,1,-1);
  h_KPMT_trigger_number_digitized_hits->SetLineColor(kGreen);
  h_KPMT_trigger_number_digitized_hits->SetLineStyle(1);
  h_KPMT_trigger_number_digitized_hits->SetLineWidth(2);
  h_KPMT_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
  TH1F *h_KPMT_nhits_max_pmt = new TH1F("h_KPMT_nhits_max_pmt","K^{40} on PMT",nbins,1,-1);
  h_KPMT_nhits_max_pmt->SetLineColor(kRed);
  h_KPMT_nhits_max_pmt->SetLineStyle(1);
  h_KPMT_nhits_max_pmt->SetLineWidth(2);
  h_KPMT_nhits_max_pmt->GetXaxis()->SetTitle("maximum n hits per PMT");
  for(int ievent=0; ievent<KPMT_primary_events_tree->GetEntries(); ievent++){
    KPMT_primary_events_tree->GetEvent(ievent); 
    for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
      h_KPMT_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
      for(int i = 0; i < nentries; ++i) nhits_max_pmt[i] = 0;
      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_tube_id->at(itrigger)).size(); idigitizedhit++){
      	nhits_max_pmt[(digitized_hit_tube_id->at(itrigger)).at(idigitizedhit)]++;
      }
      for(int i = 1; i < nentries; ++i) if(nhits_max_pmt[0] < nhits_max_pmt[i]) nhits_max_pmt[0] = nhits_max_pmt[i];
      h_KPMT_nhits_max_pmt->Fill(nhits_max_pmt[0]);
    }
  }
  h_KPMT_trigger_number_digitized_hits->Scale(1./h_KPMT_trigger_number_digitized_hits->Integral(3,100));


  TTree * BiWater_primary_events_tree = (TTree*)BiWaterinput.Get("primary_events_tree");
  BiWater_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits);
  BiWater_primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  TH1F *h_BiWater_trigger_number_digitized_hits = new TH1F("h_BiWater_trigger_number_digitized_hits","Bi^{214} in water",nbins,1,-1);
  h_BiWater_trigger_number_digitized_hits->SetLineColor(kRed);
  h_BiWater_trigger_number_digitized_hits->SetLineStyle(2);
  h_BiWater_trigger_number_digitized_hits->SetLineWidth(2);
  h_BiWater_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
  TH1F *h_BiWater_nhits_max_pmt = new TH1F("h_BiWater_nhits_max_pmt","Bi^{214} in water",nbins,1,-1);
  h_BiWater_nhits_max_pmt->SetLineColor(kRed);
  h_BiWater_nhits_max_pmt->SetLineStyle(1);
  h_BiWater_nhits_max_pmt->SetLineWidth(2);
  h_BiWater_nhits_max_pmt->GetXaxis()->SetTitle("maximum n hits per PMT");
  for(int ievent=0; ievent<BiWater_primary_events_tree->GetEntries(); ievent++){
    BiWater_primary_events_tree->GetEvent(ievent); 
    for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
      h_BiWater_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
      for(int i = 0; i < nentries; ++i) nhits_max_pmt[i] = 0;
      for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_tube_id->at(itrigger)).size(); idigitizedhit++){
      	nhits_max_pmt[(digitized_hit_tube_id->at(itrigger)).at(idigitizedhit)]++;
      }
      for(int i = 1; i < nentries; ++i) if(nhits_max_pmt[0] < nhits_max_pmt[i]) nhits_max_pmt[0] = nhits_max_pmt[i];
      h_BiWater_nhits_max_pmt->Fill(nhits_max_pmt[0]);
    }
  }
  h_BiWater_trigger_number_digitized_hits->Scale(1./h_BiWater_trigger_number_digitized_hits->Integral(3,100));
 
  TTree * Neutron_primary_events_tree = (TTree*)Neutroninput.Get("primary_events_tree");
  Neutron_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits);
  Neutron_primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);
  TH1F *h_Neutron_trigger_number_digitized_hits = new TH1F("h_Neutron_trigger_number_digitized_hits","1 MeV neutron",nbins,1,-1);
  h_Neutron_trigger_number_digitized_hits->SetLineColor(kBlack);
  h_Neutron_trigger_number_digitized_hits->SetLineStyle(1);
  h_Neutron_trigger_number_digitized_hits->SetLineWidth(2);
  h_Neutron_trigger_number_digitized_hits->SetFillColor(kBlack);
  h_Neutron_trigger_number_digitized_hits->SetFillStyle(3003);
  h_Neutron_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
  for(int ievent=0; ievent<Neutron_primary_events_tree->GetEntries(); ievent++){
    Neutron_primary_events_tree->GetEvent(ievent); 
    for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
      h_Neutron_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
    }
  }
  h_Neutron_trigger_number_digitized_hits->Scale(1./h_Neutron_trigger_number_digitized_hits->Integral(3,100));

  output.cd();
  h_BiPMT_trigger_number_digitized_hits->Write();
  h_TlPMT_trigger_number_digitized_hits->Write();
  h_KPMT_trigger_number_digitized_hits->Write();
  h_BiWater_trigger_number_digitized_hits->Write();
  h_Neutron_trigger_number_digitized_hits->Write();
  h_BiPMT_nhits_max_pmt->Write();
  h_TlPMT_nhits_max_pmt->Write();
  h_KPMT_nhits_max_pmt->Write();
  h_BiWater_nhits_max_pmt->Write();
  output.Write();

  TCanvas canvas("canvas", "canvas", 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend legend(0.5,0.5,0.85,0.85);
  legend.SetFillColor(0);

  h_Neutron_trigger_number_digitized_hits->Draw();
  h_Neutron_trigger_number_digitized_hits->GetYaxis()->SetRangeUser(0.,h_BiPMT_trigger_number_digitized_hits->GetMaximum()*1.1);
  legend.AddEntry(h_Neutron_trigger_number_digitized_hits, h_Neutron_trigger_number_digitized_hits->GetTitle(), "f");
  h_BiPMT_trigger_number_digitized_hits->Draw("same");
  legend.AddEntry(h_BiPMT_trigger_number_digitized_hits, h_BiPMT_trigger_number_digitized_hits->GetTitle(), "l");
  h_TlPMT_trigger_number_digitized_hits->Draw("same");
  legend.AddEntry(h_TlPMT_trigger_number_digitized_hits, h_TlPMT_trigger_number_digitized_hits->GetTitle(), "l");
  h_KPMT_trigger_number_digitized_hits->Draw("same");
  legend.AddEntry(h_KPMT_trigger_number_digitized_hits, h_KPMT_trigger_number_digitized_hits->GetTitle(), "l");
  h_BiWater_trigger_number_digitized_hits->Draw("same");
  legend.AddEntry(h_BiWater_trigger_number_digitized_hits, h_BiWater_trigger_number_digitized_hits->GetTitle(), "l");

  legend.Draw("same");

  canvas.Print("n-hits-bi214-tl208.png");

  std::clog << " n digitized hits from 1 MeV neutron " << h_Neutron_trigger_number_digitized_hits->GetMean() << std::endl;
  std::clog << " n digitized hits from Bi214 on PMT " << h_BiPMT_trigger_number_digitized_hits->GetMean() << std::endl;
  std::clog << " n digitized hits from Tl208 on PMT " << h_TlPMT_trigger_number_digitized_hits->GetMean() << std::endl;
  std::clog << " n digitized hits from K214 on PMT " << h_KPMT_trigger_number_digitized_hits->GetMean() << std::endl;
  std::clog << " n digitized hits from Bi214 in water " << h_BiWater_trigger_number_digitized_hits->GetMean() << std::endl;

  exit(0);
  return 1;
}



