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


  TString Isotopes[46];


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
  // Isotopes[0] = Form("Pb212");
  // Isotopes[1] = Form("Bi212");
  // Isotopes[2] = Form("Tl208");
  // Isotopes[3] = Form("Po218");
  // Isotopes[4] = Form("Pb214");
  // Isotopes[5] = Form("Bi214");
  // Isotopes[6] = Form("Tl210");
  // Isotopes[7] = Form("Pb210");
  // Isotopes[8] = Form("Bi210");
  // Isotopes[9] = Form("Hg206");
  // Isotopes[10] = Form("Tl206");
  // Isotopes[11] = Form("Po215");
  // Isotopes[12] = Form("Pb211");
  // Isotopes[13] = Form("Bi211");
  // Isotopes[14] = Form("Tl207");
  // Isotopes[15] = Form("K40");


  // all
  Isotopes[0] = Form("Tl210");// 5.4
  Isotopes[1] = Form("Tl208"); // 5.0
  Isotopes[2] = Form("Bi214"); //3.27
  Isotopes[3] = Form("Pa234"); // 2.3
  Isotopes[4] = Form("Bi215"); // 2.2
  Isotopes[5] = Form("Ac228"); // 2.1
  Isotopes[6] = Form("Tl206"); // 1.5
  Isotopes[7] = Form("K40"); // 1.3
  Isotopes[8] = Form("Tl207"); // 1.4
  Isotopes[9] = Form("Pb211"); // 0.57
  Isotopes[10] = Form("Hg206"); // 1.3
  Isotopes[11] = Form("Bi210"); // 1.2
  Isotopes[12] = Form("Bi212"); // 2.3
  Isotopes[13] = Form("Rn220");
  Isotopes[14] = Form("Po216");
  Isotopes[15] = Form("Po212");
  Isotopes[16] = Form("Rn222");
  Isotopes[17] = Form("Po218");
  Isotopes[18] = Form("Pb210");
  Isotopes[19] = Form("Po215");
  Isotopes[20] = Form("Bi211");
  Isotopes[21] = Form("Pb212");
  Isotopes[22] = Form("At218");
  Isotopes[23] = Form("Po214");
  Isotopes[24] = Form("Po210");
  Isotopes[25] = Form("Rn219");
  Isotopes[26] = Form("At215");
  Isotopes[27] = Form("Po211");
  Isotopes[28] = Form("Pb214");
  Isotopes[29] = Form("Th232");
  Isotopes[30] = Form("Ra228");
  Isotopes[31] = Form("Th228");
  Isotopes[32] = Form("Ra224");
  Isotopes[33] = Form("U238");
  Isotopes[34] = Form("Th234");
  Isotopes[35] = Form("U234");
  Isotopes[36] = Form("Th230");
  Isotopes[37] = Form("Ra226");
  Isotopes[38] = Form("U235");
  Isotopes[39] = Form("Th231");
  Isotopes[40] = Form("Pa231");
  Isotopes[41] = Form("Ac227");
  Isotopes[42] = Form("Th227");
  Isotopes[43] = Form("Fr223");
  Isotopes[44] = Form("Ra223");
  Isotopes[45] = Form("At219");


  TFile output("radioactivity_comparison.root","RECREATE");

  std::vector<Int_t> * trigger_number_digitized_hits = 0;
  std::vector<std::vector<Float_t> > * digitized_hit_Q = 0;
  std::vector<std::vector<Float_t> > * trigger_vtx_x = 0; std::vector<std::vector<Float_t> > * trigger_vtx_y = 0; std::vector<std::vector<Float_t> > * trigger_vtx_z = 0;
  std::vector<std::vector<Int_t> > * digitized_hit_tube_id = 0;

  double maxQ, totQ, OnlyInitialPMT_maxQ, OnlyInitialPMT_totQ;
  double nhits_min=-0.5;//36000;
  double nhits_max=99.5;//41000;

  TString InitialRadioactivityDirec = Form("outputs/radioactivity_")+Isotopes[0]+Form("Pmt1event/");
  //  TString InitialRadioactivityDirec = Form("outputs/radioactivity_SK_")+Isotopes[0]+Form("Pmt1event/");
  TString InitialRadioactivityFile= Form("output.root");
  TFile InitialRadioactivityinput((InitialRadioactivityDirec + InitialRadioactivityFile).Data(),"READ");

  TTree * all_pmts_tree = (TTree*)InitialRadioactivityinput.Get("all_pmts_tree");
  Int_t pmt_number, pmt_location;
  Float_t pmt_ux, pmt_uy, pmt_uz, pmt_x, pmt_y, pmt_z;
  all_pmts_tree->SetBranchAddress("pmt_number",&pmt_number);
  all_pmts_tree->SetBranchAddress("pmt_location",&pmt_location);
  all_pmts_tree->SetBranchAddress("pmt_ux",&pmt_ux);
  all_pmts_tree->SetBranchAddress("pmt_uy",&pmt_uy);
  all_pmts_tree->SetBranchAddress("pmt_uz",&pmt_uz);
  all_pmts_tree->SetBranchAddress("pmt_x",&pmt_x);
  all_pmts_tree->SetBranchAddress("pmt_y",&pmt_y);
  all_pmts_tree->SetBranchAddress("pmt_z",&pmt_z);

  TTree * geom_tree = (TTree*)InitialRadioactivityinput.Get("geom_tree");
  Float_t detector_length, detector_radius, pmt_radius;
  Int_t number_of_pmts;
  geom_tree->SetBranchAddress("detector_length",&detector_length);
  geom_tree->SetBranchAddress("detector_radius",&detector_radius);
  geom_tree->SetBranchAddress("pmt_radius",&pmt_radius);
  geom_tree->SetBranchAddress("number_of_pmts",&number_of_pmts);
  geom_tree->GetEntry(0);

  double vtx_x, vtx_y, vtx_z;
  double dist;
  bool found_pmt;
  int found_ipmt;
  size_t n_found_hits;


  for(int i=0; i<46; i++){

    TString RadioactivityDirec = Form("outputs/radioactivity_")+Isotopes[i]+Form("Pmt1event/");
    //    TString RadioactivityDirec = Form("outputs/radioactivity_SK_")+Isotopes[i]+Form("Pmt1event/");
    TString RadioactivityFile= Form("output.root");
    

    TFile Radioactivityinput((RadioactivityDirec + RadioactivityFile).Data(),"READ");

    TTree * Radioactivity_primary_events_tree = (TTree*)Radioactivityinput.Get("primary_events_tree");
    Radioactivity_primary_events_tree->SetBranchAddress("trigger_number_digitized_hits",&trigger_number_digitized_hits); 
    Radioactivity_primary_events_tree->SetBranchAddress("digitized_hit_Q",&digitized_hit_Q);
    Radioactivity_primary_events_tree->SetBranchAddress("trigger_vtx_x",&trigger_vtx_x); 
    Radioactivity_primary_events_tree->SetBranchAddress("trigger_vtx_y",&trigger_vtx_y);
    Radioactivity_primary_events_tree->SetBranchAddress("trigger_vtx_z",&trigger_vtx_z);
    Radioactivity_primary_events_tree->SetBranchAddress("digitized_hit_tube_id",&digitized_hit_tube_id);

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
    
    TH1F *h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits = new TH1F("h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits",(Isotopes[i]).Data(),100,nhits_min,nhits_max);
    h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->SetLineColor(1+i);
    h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->SetLineStyle(1);
    h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->SetLineWidth(2);
    h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->GetXaxis()->SetTitle("number of digitized hits");
    TH1F *h_OnlyInitialPMT_Radioactivity_charge_hit = new TH1F("h_OnlyInitialPMT_Radioactivity_charge_hit",(Isotopes[i]).Data(),350,-0.5,349.5);
    h_OnlyInitialPMT_Radioactivity_charge_hit->SetLineColor(1+i);
    h_OnlyInitialPMT_Radioactivity_charge_hit->SetLineStyle(1);
    h_OnlyInitialPMT_Radioactivity_charge_hit->SetLineWidth(2);
    h_OnlyInitialPMT_Radioactivity_charge_hit->GetXaxis()->SetTitle("charge per hit [p.e.]");
    TH1F *h_OnlyInitialPMT_Radioactivity_charge_max_hit = new TH1F("h_OnlyInitialPMT_Radioactivity_charge_max_hit",(Isotopes[i]).Data(),350,-0.5,349.5);
    h_OnlyInitialPMT_Radioactivity_charge_max_hit->SetLineColor(1+i);
    h_OnlyInitialPMT_Radioactivity_charge_max_hit->SetLineStyle(1);
    h_OnlyInitialPMT_Radioactivity_charge_max_hit->SetLineWidth(2);
    h_OnlyInitialPMT_Radioactivity_charge_max_hit->GetXaxis()->SetTitle("charge of hit with maximum charge [p.e.]");
    TH1F *h_OnlyInitialPMT_Radioactivity_charge_all_hits = new TH1F("h_OnlyInitialPMT_Radioactivity_charge_all_hits",(Isotopes[i]).Data(),350,-0.5,349.5);
    h_OnlyInitialPMT_Radioactivity_charge_all_hits->SetLineColor(1+i);
    h_OnlyInitialPMT_Radioactivity_charge_all_hits->SetLineStyle(1);
    h_OnlyInitialPMT_Radioactivity_charge_all_hits->SetLineWidth(2);
    h_OnlyInitialPMT_Radioactivity_charge_all_hits->GetXaxis()->SetTitle("charge of all hits in the event [p.e.]");
    
    std::clog << " isotope " << (Isotopes[i]).Data() << " events " << Radioactivity_primary_events_tree->GetEntries() << std::endl;
    for(int ievent=0; ievent<Radioactivity_primary_events_tree->GetEntries(); ievent++){
      
      Radioactivity_primary_events_tree->GetEvent(ievent); 
      for(size_t itrigger=0; itrigger<trigger_number_digitized_hits->size(); itrigger++){

	vtx_x = (trigger_vtx_x->at(itrigger)).at(0);
	vtx_y = (trigger_vtx_y->at(itrigger)).at(0);
	vtx_z = (trigger_vtx_z->at(itrigger)).at(0);
	
	found_pmt = false;
	for(int ipmt = 0; ipmt < all_pmts_tree->GetEntries(); ipmt++){
	  all_pmts_tree->GetEntry(ipmt);
	  dist = std::sqrt(std::pow(vtx_x - pmt_x,2) + std::pow(vtx_y - pmt_y,2) + std::pow(vtx_z - pmt_z,2));
	  if( dist < 2*pmt_radius ){
	    found_pmt = true;
	    found_ipmt = ipmt;
	    break;
	  }
	}
	if( ! found_pmt ){
	  std::clog << " problem: cannot find pmt for vertex (" << vtx_x << ", " << vtx_y << ", " << vtx_z << ") " << std::endl;
	  break;
	}

	h_Radioactivity_trigger_number_digitized_hits->Fill(trigger_number_digitized_hits->at(itrigger));
	maxQ = 0.;
	totQ = 0.;
	OnlyInitialPMT_maxQ=0.;
	OnlyInitialPMT_totQ=0.;
	n_found_hits = 0;
	for(size_t idigitizedhit=0; idigitizedhit<(digitized_hit_Q->at(itrigger)).size(); idigitizedhit++){
	  h_Radioactivity_charge_hit->Fill((digitized_hit_Q->at(itrigger)).at(idigitizedhit));
	  totQ += (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
	  if( (digitized_hit_Q->at(itrigger)).at(idigitizedhit) > maxQ ) maxQ = (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
	  if( (digitized_hit_tube_id->at(itrigger)).at(idigitizedhit) == found_ipmt + 1 ){
	    n_found_hits ++;
	    h_OnlyInitialPMT_Radioactivity_charge_hit->Fill((digitized_hit_Q->at(itrigger)).at(idigitizedhit));
	    OnlyInitialPMT_totQ += (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
	    if( (digitized_hit_Q->at(itrigger)).at(idigitizedhit) > OnlyInitialPMT_maxQ ) OnlyInitialPMT_maxQ = (digitized_hit_Q->at(itrigger)).at(idigitizedhit);
	  }
	}
	h_Radioactivity_charge_max_hit->Fill(maxQ);
	h_Radioactivity_charge_all_hits->Fill(totQ);
	h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->Fill(n_found_hits);
	h_OnlyInitialPMT_Radioactivity_charge_max_hit->Fill(OnlyInitialPMT_maxQ);
	h_OnlyInitialPMT_Radioactivity_charge_all_hits->Fill(OnlyInitialPMT_totQ);
      }
    }
    h_Radioactivity_trigger_number_digitized_hits->Scale(1./h_Radioactivity_trigger_number_digitized_hits->Integral());
    h_Radioactivity_charge_max_hit->Scale(1./h_Radioactivity_charge_max_hit->Integral());
    h_Radioactivity_charge_all_hits->Scale(1./h_Radioactivity_charge_all_hits->Integral());
    h_Radioactivity_charge_hit->Scale(1./h_Radioactivity_charge_hit->Integral());
    h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->Scale(1./h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->Integral());
    h_OnlyInitialPMT_Radioactivity_charge_max_hit->Scale(1./h_OnlyInitialPMT_Radioactivity_charge_max_hit->Integral());
    h_OnlyInitialPMT_Radioactivity_charge_all_hits->Scale(1./h_OnlyInitialPMT_Radioactivity_charge_all_hits->Integral());
    h_OnlyInitialPMT_Radioactivity_charge_hit->Scale(1./h_OnlyInitialPMT_Radioactivity_charge_hit->Integral());

    output.cd();
    h_Radioactivity_trigger_number_digitized_hits->Write();
    h_Radioactivity_charge_max_hit->Write();
    h_Radioactivity_charge_all_hits->Write();
    h_Radioactivity_charge_hit->Write();
    h_OnlyInitialPMT_Radioactivity_trigger_number_digitized_hits->Write();
    h_OnlyInitialPMT_Radioactivity_charge_max_hit->Write();
    h_OnlyInitialPMT_Radioactivity_charge_all_hits->Write();
    h_OnlyInitialPMT_Radioactivity_charge_hit->Write();
    output.Write();

  }

  exit(0);
  return 1;
}



