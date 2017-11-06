#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TFile.h>
#include <TTree.h>

using namespace std;


#define NEVENTS 10000
#define NHITS 200
#define NPMTS 50000

int main(){

  TFile input("angles.root","READ");
  TFile * output = new TFile("angles_plot.root","RECREATE");

  TTree * angles_tree = (TTree*)input.Get("angles_tree");
  Float_t vtx_x, vtx_y, vtx_z;
  Float_t dir_x, dir_y, dir_z;
  Float_t kinetic;
  Int_t n_hits;
  Float_t hit_x[NHITS];
  Float_t hit_y[NHITS];
  Float_t hit_z[NHITS];
  Float_t hit_t[NHITS];
  angles_tree->SetBranchAddress("vtx_x",&vtx_x);
  angles_tree->SetBranchAddress("vtx_y",&vtx_y);
  angles_tree->SetBranchAddress("vtx_z",&vtx_z);
  angles_tree->SetBranchAddress("dir_x",&dir_x);
  angles_tree->SetBranchAddress("dir_y",&dir_y);
  angles_tree->SetBranchAddress("dir_z",&dir_z);
  angles_tree->SetBranchAddress("kinetic",&kinetic);
  angles_tree->SetBranchAddress("n_hits",&n_hits);
  angles_tree->SetBranchAddress("hit_x",hit_x);
  angles_tree->SetBranchAddress("hit_y",hit_y);
  angles_tree->SetBranchAddress("hit_z",hit_z);
  angles_tree->SetBranchAddress("hit_t",hit_t);

  double cos_theta;
  double dist_x, dist_y, dist_z, dist;
  double tof, tcorr;

  double theta;
  double theta0 = 41.4;
  double pi = acos(-1.);
  double c = 30/1.3; // cm / ns
  double tmin = 950.;
  double tmax = 970.;

  TH1F h_t("h_t","t;time [ns]",200,1,-1);
  TH1F h_tcorr("h_tcorr","tcorr;corrected time [ns]",200,1,-1);

  TH1F h_cos_theta_1("h_cos_theta_1","1 MeV;cos#theta",200,-1,1);
  h_cos_theta_1.SetLineColor(kGray);
  h_cos_theta_1.SetLineWidth(2);
  TH1F h_theta_1("h_theta_1","1 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_1("h_theta_from_theta0_1","1 MeV;#theta - #theta_{0}",200,-360.,360.);
  TH1F h_cos_theta_2("h_cos_theta_2","2 MeV;cos#theta",200,-1,1);
  h_cos_theta_2.SetLineColor(kBlack);
  h_cos_theta_2.SetLineWidth(2);
  TH1F h_theta_2("h_theta_2","2 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_2("h_theta_from_theta0_2","2 MeV;#theta - #theta_{0}",200,-360.,360.);
  TH1F h_cos_theta_3("h_cos_theta_3","3 MeV;cos#theta",200,-1,1);
  h_cos_theta_3.SetLineColor(kRed);
  h_cos_theta_3.SetLineWidth(2);
  TH1F h_theta_3("h_theta_3","3 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_3("h_theta_from_theta0_3","3 MeV;#theta - #theta_{0}",200,-360.,360.);
  TH1F h_cos_theta_4("h_cos_theta_4","4 MeV;cos#theta",200,-1,1);
  h_cos_theta_4.SetLineColor(kBlue);
  h_cos_theta_4.SetLineWidth(2);
  TH1F h_theta_4("h_theta_4","4 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_4("h_theta_from_theta0_4","4 MeV;#theta - #theta_{0}",200,-360.,360.);
  TH1F h_cos_theta_5("h_cos_theta_5","5 MeV;cos#theta",200,-1,1);
  h_cos_theta_5.SetLineColor(kGreen);
  h_cos_theta_5.SetLineWidth(2);
  TH1F h_theta_5("h_theta_5","5 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_5("h_theta_from_theta0_5","5 MeV;#theta - #theta_{0}",200,-360.,360.);

  TH1F h_hits_fraction_1("h_hits_fraction_1","1 MeV;hits_fraction",50,0,1);
  h_hits_fraction_1.SetLineColor(kGray);
  h_hits_fraction_1.SetLineWidth(2);
  TH1F h_hits_fraction_2("h_hits_fraction_2","2 MeV;hits_fraction",50,0,1);
  h_hits_fraction_2.SetLineColor(kBlack);
  h_hits_fraction_2.SetLineWidth(2);
  TH1F h_hits_fraction_3("h_hits_fraction_3","3 MeV;hits_fraction",50,0,1);
  h_hits_fraction_3.SetLineColor(kRed);
  h_hits_fraction_3.SetLineWidth(2);
  TH1F h_hits_fraction_4("h_hits_fraction_4","4 MeV;hits_fraction",50,0,1);
  h_hits_fraction_4.SetLineColor(kBlue);
  h_hits_fraction_4.SetLineWidth(2);
  TH1F h_hits_fraction_5("h_hits_fraction_5","5 MeV;hits_fraction",50,0,1);
  h_hits_fraction_5.SetLineColor(kGreen);
  h_hits_fraction_5.SetLineWidth(2);

  double cos_theta_min = 0.72;
  double cos_theta_max = 0.78;

  for(int ientry = 0; ientry < angles_tree->GetEntries(); ientry++){
    angles_tree->GetEntry(ientry);

    double n_hits_in_range=0;

    for(int ihit=0; ihit < n_hits; ihit++){
      dist_x = hit_x[ihit] - vtx_x;
      dist_y = hit_y[ihit] - vtx_y;
      dist_z = hit_z[ihit] - vtx_z;
      dist = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));
      tof = dist/c;
      tcorr = hit_t[ihit] - tof;
      h_t.Fill(hit_t[ihit]);
      h_tcorr.Fill(tcorr);

      //   if( tcorr < tmin || tcorr > tmax ) continue;

      cos_theta = (dist_x*dir_x + dist_y*dir_y + dist_z*dir_z)/dist;
      theta = acos(cos_theta)*180./pi;

      if( cos_theta >= cos_theta_min && cos_theta <= cos_theta_max ){
	n_hits_in_range++;
      }

    
      if( kinetic == 1. ){
	h_cos_theta_1.Fill(cos_theta);
	h_theta_1.Fill(theta);
	h_theta_from_theta0_1.Fill(theta - theta0);
      }else if( kinetic == 2. ){
	h_cos_theta_2.Fill(cos_theta);
	h_theta_2.Fill(theta);
	h_theta_from_theta0_2.Fill(theta - theta0);
      }else if( kinetic == 3. ){
	h_cos_theta_3.Fill(cos_theta);
	h_theta_3.Fill(theta);
	h_theta_from_theta0_3.Fill(theta - theta0);
      }else if( kinetic == 4. ){
	h_cos_theta_4.Fill(cos_theta);
	h_theta_4.Fill(theta);
	h_theta_from_theta0_4.Fill(theta - theta0);
      }else if( kinetic == 5. ){
	h_cos_theta_5.Fill(cos_theta);
	h_theta_5.Fill(theta);
	h_theta_from_theta0_5.Fill(theta - theta0);
      }else{
	clog << " unknown kinetic " << kinetic << endl;
      }


    }
    
    
    if( kinetic == 1. )
      h_hits_fraction_1.Fill(n_hits_in_range/(double)n_hits);
    else if( kinetic == 2. )
      h_hits_fraction_2.Fill(n_hits_in_range/(double)n_hits);
    else if( kinetic == 3. )
      h_hits_fraction_3.Fill(n_hits_in_range/(double)n_hits);
    else if( kinetic == 4. )
      h_hits_fraction_4.Fill(n_hits_in_range/(double)n_hits);
    else if( kinetic == 5. )
      h_hits_fraction_5.Fill(n_hits_in_range/(double)n_hits);

  }


  output->Write();

  return 0;

}

