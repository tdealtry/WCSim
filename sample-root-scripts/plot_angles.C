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

  double theta;
  double theta0 = 41.4;
  double pi = acos(-1.);

  TH1F h_cos_theta("h_cos_theta","cos#theta",100,-1,1);
  TH1F h_theta_from_theta0("h_theta_from_theta0","#theta - #theta_{0}",100,-360.,360.);

  for(int ientry = 0; ientry < angles_tree->GetEntries(); ientry++){
    angles_tree->GetEntry(ientry);

    for(int ihit=0; ihit < n_hits; ihit++){
      dist_x = hit_x[ihit] - vtx_x;
      dist_y = hit_y[ihit] - vtx_y;
      dist_z = hit_z[ihit] - vtx_z;
      dist = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));
      cos_theta = (dist_x*dir_x + dist_y*dir_y + dist_z*dir_z)/dist;
      h_cos_theta.Fill(cos_theta);
      theta = acos(cos_theta)*180./pi;
      h_theta_from_theta0.Fill(theta - theta0);
    }

  }

  output->Write();

  return 0;

}

