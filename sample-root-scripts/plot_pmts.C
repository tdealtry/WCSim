#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TExec.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"
#endif

void FillHistogram(int ipmt, WCSimRootGeom * geo, const double z_limit,
		   TH2I * htop, TH2I * hbottom, TH2I * hside, bool pmtid_as_weight)
{
  //get PMT info
  WCSimRootPMT pmt = geo->GetPMT(ipmt);
  const int pmtno = pmt.GetTubeNo();
  assert(pmtno == ipmt + 1);
  const double x0 = pmt.GetPosition(0);
  const double y0 = pmt.GetPosition(1);
  const double z0 = pmt.GetPosition(2);
  const double phi0 = atan2(y0, x0) * 180. / acos(-1.);
  const int weight = pmtid_as_weight ? pmtno : 1;
  if(z0 > z_limit){ // top
    if(htop->GetBinContent(htop->FindBin(x0, y0))) {
      cerr << "PMT with top position " << x0 << " " << y0 << " already exists" << endl;
      exit(-1);
    }
    htop->Fill(x0, y0, weight);
  }
  else if(z0 < - z_limit){ // bottom
    if(hbottom->GetBinContent(hbottom->FindBin(x0, y0))) {
      cerr << "PMT with bottom position " << x0 << " " << y0 << " already exists" << endl;
      exit(-1);
    }
    hbottom->Fill(x0, y0, weight);
  }
  else{ // side
    if(hside->GetBinContent(hside->FindBin(phi0, z0))) {
      cerr << "PMT with side position " << phi0 << " " << z0 << " already exists" << endl;
      exit(-1);
    }
    hside->Fill(phi0, z0, weight);
  }
}

void SetStyle(TH2I * h, Color_t color, int size, int style, bool exec)
{
  h->SetFillColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerSize(size);
  h->SetMarkerStyle(style);
  h->GetXaxis()->SetTitleSize(0.02);
  h->GetXaxis()->SetLabelSize(0.02);
  h->GetYaxis()->SetTitleSize(0.02);
  h->GetYaxis()->SetLabelSize(0.02);
  h->GetXaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetTickLength(0.01);
  if(exec) {
    TExec * ex = new TExec("ex", "drawtext()");
    h->GetListOfFunctions()->Add(ex);
  }
}

string gHist;
double gTextSize;

void drawtext()
{
  Int_t p;
  Double_t x,y;
  TLatex *l;
  TH2I *h = (TH2I*)gPad->GetListOfPrimitives()->FindObject(gHist.c_str());
  for(int ix = 0; ix <= h->GetNbinsX(); ix++) {
    for(int iy = 0; iy <= h->GetNbinsY(); iy++) {
      int p = h->GetBinContent(ix, iy);
      if(p) {
	x = h->GetXaxis()->GetBinCenter(ix);
	y = h->GetYaxis()->GetBinCenter(iy);
	l = new TLatex(x, y, Form("%d", p));
	l->SetTextColor(h->GetMarkerColor());
	l->SetTextSize(gTextSize);
	l->SetTextFont(42);
 	l->SetTextAlign(22);
	l->SetTextAngle(45);
	l->Paint();
      }
    }//iy
  }//ix
}


int plot_pmts(const char * filename = "../wcsim.root", TString deadpmt_filename = "",
	      double textsize = 0.005, bool verbose = true)
{
  gTextSize = textsize;

  // open input file
  TFile * f = new TFile(filename, "READ");

  // Geometry tree - only need 1 "event"
  TTree * geotree = 0;
  f->GetObject("wcsimGeoT", geotree);
  WCSimRootGeom * geo = 0; 
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose)
    std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
    exit(9);
  }
  geotree->GetEntry(0);

  //get detector sizes, pmt sizes, npmts, ...
  const Float_t detector_radius = geo->GetWCCylRadius();
  const Float_t detector_length = geo->GetWCCylLength(); 
  const Float_t r_limit = detector_radius;
  const Float_t z_limit = detector_length / 2.;
  const Int_t   npmts = geo->GetWCNumPMT();
  const Float_t pmt_radius = geo->GetWCPMTRadius();
  const TString detname = geo->GetDetectorName();

  // define space binning so that different pmts are sure to end up in different bins
  int nbins_x   = (int)(3. * detector_radius / pmt_radius);
  int nbins_y   = nbins_x;
  int nbins_z   = (int)(3. * detector_length / (2. * pmt_radius));
  int nbins_phi = (int)(3. * acos(-1.) * detector_radius / pmt_radius);
  const double expansion_factor = 1.05;

  //open output file
  TString ofilenamebase = TString::Format("%s", detname.Data());
  bool deadpmtlist_exists = !string(deadpmt_filename).empty();
  if(deadpmtlist_exists) {
    TObjArray * tokens = deadpmt_filename.Tokenize("/");
    tokens = ((TObjString*)(tokens->Last()))->String().Tokenize(".");
    TString deadfilestub = ((TObjString*)(tokens->First()))->String();
    ofilenamebase += TString::Format("_%s", deadfilestub.Data());
  }
  TFile * of = new TFile(TString::Format("%s.root", ofilenamebase.Data()), "RECREATE");

  //setup all-pmts histograms
  TH2I *  h_y0_vs_x0_top = new TH2I("h_y0_vs_x0_top","y0_vs_x0_top;x0 (mm);y0 (mm)",
				    nbins_x, -r_limit * expansion_factor, r_limit * expansion_factor,
				    nbins_y, -r_limit * expansion_factor, r_limit * expansion_factor);
  SetStyle(h_y0_vs_x0_top, kBlack, 0, 0, true);
  TH2I *  h_y0_vs_x0_bottom = new TH2I("h_y0_vs_x0_bottom","y0_vs_x0_bottom;x0 (mm);y0 (mm)",
				       nbins_x, -r_limit * expansion_factor, r_limit * expansion_factor,
				       nbins_y, -r_limit * expansion_factor, r_limit * expansion_factor);
  SetStyle(h_y0_vs_x0_bottom, kBlack, 0, 0, true);
  TH2I *  h_z0_vs_phi = new TH2I("h_z0_vs_phi","z0_vs_phi;phi (degrees);z0 (mm)",
				 nbins_phi, -180. * expansion_factor, 180. * expansion_factor,
				 nbins_z, -z_limit * expansion_factor, z_limit * expansion_factor);
  SetStyle(h_z0_vs_phi, kBlack, 0, 0, true);

  //and the dead ones
  TH2I *  h_y0_vs_x0_top_dead = new TH2I("h_y0_vs_x0_top_dead","y0_vs_x0_top_dead;x0 (mm);y0 (mm)",
					 nbins_x, -r_limit * expansion_factor, r_limit * expansion_factor,
					 nbins_y, -r_limit * expansion_factor, r_limit * expansion_factor);
  SetStyle(h_y0_vs_x0_top_dead, kRed, 0, 0, false);
  TH2I *  h_y0_vs_x0_bottom_dead = new TH2I("h_y0_vs_x0_bottom_dead","y0_vs_x0_bottom_dead;x0 (mm);y0 (mm)",
					    nbins_x, -r_limit * expansion_factor, r_limit * expansion_factor,
					    nbins_y, -r_limit * expansion_factor, r_limit * expansion_factor);
  SetStyle(h_y0_vs_x0_bottom_dead, kRed, 0, 0, false);
  TH2I *  h_z0_vs_phi_dead = new TH2I("h_z0_vs_phi_dead","z0_vs_phi_dead;phi (degrees);z0 (mm)",
				      nbins_phi, -180. * expansion_factor, 180. * expansion_factor,
				      nbins_z, -z_limit * expansion_factor, z_limit * expansion_factor);
  SetStyle(h_z0_vs_phi_dead, kRed, 0, 0, false);

  //fill the all-pmts histograms
  for(int ipmt = 0; ipmt < npmts; ipmt++) {
    FillHistogram(ipmt, geo, z_limit - 1E-6, h_y0_vs_x0_top, h_y0_vs_x0_bottom, h_z0_vs_phi, true);
  }//ipmt

  //fill the dead-pmts histograms
  if(deadpmtlist_exists) {
    ifstream ifs(deadpmt_filename.Data());
    int pmtid;
    vector<int> pmts;
    if(ifs.good()) {
      while(!ifs.eof()) {
	ifs >> pmtid;
	pmts.push_back(pmtid);
      }
    }
    std::sort(pmts.begin(), pmts.end());
    std::vector<int>::iterator it = std::unique(pmts.begin(), pmts.end());
    pmts.resize(std::distance(pmts.begin(), it));
    for(size_t i = 0; i < pmts.size(); i++)
      FillHistogram(pmts[i] - 1, geo, z_limit - 1E-6, h_y0_vs_x0_top_dead, h_y0_vs_x0_bottom_dead, h_z0_vs_phi_dead, false);
    cout << pmts.size() << " total dead PMTs" << endl;
  }

  f->Close();

  TCanvas * c_display = new TCanvas("c_display","c_display",600,600);
  double BottomMargin = .05;   // fraction of pad height
  double LeftMargin   = .05;
  double RightMargin  = .02;
  double TopMargin    = .02;
  c_display->SetBottomMargin(BottomMargin);
  c_display->SetLeftMargin(LeftMargin);
  c_display->SetRightMargin(RightMargin);
  c_display->SetTopMargin(TopMargin);

  
  c_display->Print(TString::Format("%s.pdf[", ofilenamebase.Data()));

  gHist = "h_y0_vs_x0_top";
  h_y0_vs_x0_top->Draw("AXIS");
  h_y0_vs_x0_top_dead->Draw("SAME COL");
  c_display->Modified();
  c_display->Update();
  c_display->Print(TString::Format("%s.pdf", ofilenamebase.Data()));

  gHist = "h_z0_vs_phi";
  h_z0_vs_phi->Draw("AXIS");
  h_z0_vs_phi_dead->Draw("SAME COL");
  c_display->Modified();
  c_display->Update();
  c_display->Print(TString::Format("%s.pdf", ofilenamebase.Data()));

  gHist = "h_y0_vs_x0_bottom";
  h_y0_vs_x0_bottom->Draw("AXIS");
  h_y0_vs_x0_bottom_dead->Draw("SAME COL");
  c_display->Modified();
  c_display->Update();
  c_display->Print(TString::Format("%s.pdf", ofilenamebase.Data()));

  c_display->Print(TString::Format("%s.pdf]", ofilenamebase.Data()));

  of->cd();
  h_y0_vs_x0_top->Write();
  h_y0_vs_x0_bottom->Write();
  h_z0_vs_phi->Write();

  return 1;

}

