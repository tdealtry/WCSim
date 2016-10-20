#include <iostream>
#include <TH1F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TSystem.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"
#endif

#ifdef __MAKECINT__
#pragma link C++ class std::vector<std::vector<Int_t> >+;
#pragma link C++ class std::vector<std::vector<Float_t> >+;
#pragma link C++ class std::vector<std::vector<std::vector<Int_t> > >+;
#pragma link C++ class std::vector<std::vector<std::vector<Float_t> > >+;
#endif

void load_library();
TFile * get_input_file(char *filename);
TString create_filename(const char * prefix, TString& filename_string);

void create_pmt_list(char *filename="../wcsim.root", bool verbose = false)
{
  // load library
  load_library();

  // open input file
  TFile *f = get_input_file(filename);
  
  // load geometry tree
  WCSimRootGeom* geom = new WCSimRootGeom();
  TTree *input_geom_tree = (TTree*)f->Get("wcsimGeoT");
  TBranch *gb = input_geom_tree->GetBranch("wcsimrootgeom");
  gb->SetAddress(&geom);
  input_geom_tree->GetEntry(0);

  //open output file
  TString ofilename = geom->GetDetectorName().c_str();
  ofilename += "_pmts.list";
  std::ofstream output(ofilename, std::ofstream::out);

  // all pmts list
  const int number_of_pmts = geom->GetWCNumPMT();
  Int_t pmt_number, pmt_location;
  Float_t pmt_ux, pmt_uy, pmt_uz, pmt_x, pmt_y, pmt_z;
  for(int i=0; i<number_of_pmts; i++){
    WCSimRootPMT pmt = geom->GetPMT(i);
    pmt_number = pmt.GetTubeNo();
    //pmt_location = pmt.GetCylLoc();
    //pmt_ux = pmt.GetOrientation(0);
    //pmt_uy = pmt.GetOrientation(1);
    //pmt_uz = pmt.GetOrientation(2);
    pmt_x = pmt.GetPosition(0);
    pmt_y = pmt.GetPosition(1);
    pmt_z = pmt.GetPosition(2);
    output << pmt_number
	   << "\t" << pmt_x
	   << "\t" << pmt_y
	   << "\t" << pmt_z
	   << endl;
  }

  output.close();
}


void load_library(){

#if !defined(__MAKECINT__)
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  }else{
    gSystem->Load("../libWCSimRoot.so");
  }
#endif
  return;

}

TString create_filename(const char * prefix, TString& filename_string)
{
  //std::cout << "Creating filename from prefix " << prefix << " and filename_string " << filename_string << std::endl;                                                                                                                                                  
  TString prefix_string(prefix);
  TString outfilename = prefix_string + filename_string;
  return outfilename;
}

TFile * get_input_file(char *filename){

  TFile *f;
  // Open the file                                                                                                            
  f = new TFile(filename,"read");

  if (!f->IsOpen()){
    std::cerr << "Error, could not open input file: " << filename << std::endl;
    exit(0);
  }
  
  return f;

}



