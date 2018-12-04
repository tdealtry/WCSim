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


  TString Isotopes[16];


  // alpha
  // Isotopes[0] = Form("Rn220");
  // Isotopes[1] = Form("Po216");
  // Isotopes[2] = Form("Bi212");
  // Isotopes[3] = Form("Po212");
  // Isotopes[4] = Form("Rn222");
  // Isotopes[5] = Form("Po218");
  // Isotopes[6] = Form("At218");
  // Isotopes[7] = Form("Bi214");
  // Isotopes[8] = Form("Po214");
  // Isotopes[9] = Form("Pb210");
  // Isotopes[10] = Form("Bi210");
  // Isotopes[11] = Form("Po210");
  // Isotopes[12] = Form("Rn219");
  // Isotopes[13] = Form("Po215");
  // Isotopes[14] = Form("At215");
  // Isotopes[15] = Form("Bi211");
  // Isotopes[16] = Form("Po211");


  // beta
  Isotopes[0] = Form("Pb212");
  Isotopes[1] = Form("Bi212");
  Isotopes[2] = Form("Tl208");
  Isotopes[3] = Form("Po218");
  Isotopes[4] = Form("Pb214");
  Isotopes[5] = Form("Bi214");
  Isotopes[6] = Form("Tl210");
  Isotopes[7] = Form("Pb210");
  Isotopes[8] = Form("Bi210");
  Isotopes[9] = Form("Hg206");
  Isotopes[10] = Form("Tl206");
  Isotopes[11] = Form("Po215");
  Isotopes[12] = Form("Pb211");
  Isotopes[13] = Form("Bi211");
  Isotopes[14] = Form("Tl207");
  Isotopes[15] = Form("K40");


  TFile output("radioactivity_comparison.root","RECREATE");

  std::vector<Int_t> * trigger_number_digitized_hits = 0;
  std::vector<std::vector<Float_t> > * digitized_hit_Q = 0;

  double maxQ, totQ;
  double nhits_min=-0.5;//36000;
  double nhits_max=99.5;//41000;

  for(int i=0; i<16; i++){

    TString RadioactivityDirec = Form("outputs/radioactivity_")+Isotopes[i]+Form("Pmt1event/");
    TString RadioactivityFile= Form("output.root");
    

    TFile Radioactivityinput((RadioactivityDirec + RadioactivityFile).Data(),"READ");

    TTree * Radioactivity_primary_events_tree = (TTree*)Radioactivityinput.Get("primary_events_tree");
    Radioactivity_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits); 
    Radioactivity_primary_events_tree->SetBranchAddress("digitized_hit_Q",&digitized_hit_Q);
    TH1F *h_Radioactivity_trigger_number_digitized_hits = new TH1F("h_Radioactivity_trigger_number_digitized_hits",(Isotopes[i]).Data(),100,nhits_min,nhits_max);
    h_Radioactivity_trigger_number_digitized_hits->SetLineColor(1+i);
    h_Radioactivity_trigger_number_digitized_hits->SetLineStyle(1);
    h_Radioactivity_trigger_number_digitized_hits->SetLineWidth(2);
    h_Radioactivity_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
    TH1F *h_Radioactivity_charge_hit = new TH1F("h_Radioactivity_charge_hit",(Isotopes[i]).Data(),350,-0.5,349.5);
    h_Radioactivity_charge_hit->SetLineColor(1+i);
    h_Radioactivity_charge_hit->SetLineStyle(1);
    h_Radioactivity_charge_hit->SetLineWidth(2);
    h_Radioactivity_charge_hit->GetXaxis()->SetTitle("charge per hit [p.e.]");
    TH1F *h_Radioactivity_charge_max_hit = new TH1F("h_Radioactivity_charge_max_hit",(Isotopes[i]).Data(),350,-0.5,349.5);
    h_Radioactivity_charge_max_hit->SetLineColor(1+i);
    h_Radioactivity_charge_max_hit->SetLineStyle(1);
    h_Radioactivity_charge_max_hit->SetLineWidth(2);
    h_Radioactivity_charge_max_hit->GetXaxis()->SetTitle("charge of hit with maximum charge [p.e.]");
    TH1F *h_Radioactivity_charge_all_hits = new TH1F("h_Radioactivity_charge_all_hits",(Isotopes[i]).Data(),350,-0.5,349.5);
    h_Radioactivity_charge_all_hits->SetLineColor(1+i);
    h_Radioactivity_charge_all_hits->SetLineStyle(1);
    h_Radioactivity_charge_all_hits->SetLineWidth(2);
    h_Radioactivity_charge_all_hits->GetXaxis()->SetTitle("charge of all hits in the event [p.e.]");
    
    std::clog << " isotope " << (Isotopes[i]).Data() << " events " << Radioactivity_primary_events_tree->GetEntries() << std::endl;
    for(int ievent=0; ievent<Radioactivity_primary_events_tree->GetEntries(); ievent++){
      
      Radioactivity_primary_events_tree->GetEvent(ievent); 
      for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){
	h_Radioactivity_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
	maxQ = 0.;
	totQ = 0.;
	for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_Q->at(itrigger)).size(); idigitizedhit++){
	  h_Radioactivity_charge_hit->Fill((digitized_hit_Q->at(itrigger)).at(idigitizedhit));
	  totQ += (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
	  if( (digitized_hit_Q->at(itrigger)).at(idigitizedhit) > maxQ ) maxQ = (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
	}
	h_Radioactivity_charge_max_hit->Fill(maxQ);
	h_Radioactivity_charge_all_hits->Fill(totQ);
      }
    }
    h_Radioactivity_trigger_number_digitized_hits->Scale(1./h_Radioactivity_trigger_number_digitized_hits->Integral());
    h_Radioactivity_charge_max_hit->Scale(1./h_Radioactivity_charge_max_hit->Integral());
    h_Radioactivity_charge_all_hits->Scale(1./h_Radioactivity_charge_all_hits->Integral());
    h_Radioactivity_charge_hit->Scale(1./h_Radioactivity_charge_hit->Integral());

    output.cd();
    h_Radioactivity_trigger_number_digitized_hits->Write();
    h_Radioactivity_charge_max_hit->Write();
    h_Radioactivity_charge_all_hits->Write();
    h_Radioactivity_charge_hit->Write();
    output.Write();

  }

  exit(0);
  return 1;
}



