#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "TGaxis.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TTree.h"
#include "TStyle.h"
#include <cassert>

using namespace std;

Color_t cols  [3] = {kBlack, kRed, kBlue};
Style_t styles[3] = {1,2,3};
int energies[3] = {0, 4, 10};

double efficiencies[3] = {0.50, 0.90, 0.99};

TTree * trees[3];
double var_max_ndigits[3];
double var_max_ndigits_time[3];
double var_max_ndigits_1slice[3];
double var_max_itc[3];
double var_max_itc_time[3];
double var_max_itc_1slice[3];
TH1D * h1itc_fine_sig = new TH1D("h1itc_fine_sig","",200000,0,2);
TH1D * h1itc_fine_bkg = new TH1D("h1itc_fine_bkg","",200000,0,2);

void SetMargins(TVirtualPad & c, double left = 0.15, double right = 0.05, double top = 0.05, double bottom = 0.15)
{
  c.SetLeftMargin(left);
  c.SetRightMargin(right);
  c.SetTopMargin(top);
  c.SetBottomMargin(bottom);
}

double smallfracs[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9};
int  largewindows[] = {200, 300, 400, 500, 600, 800, 1000, 1250, 1500, 1750, 2000};
const int nsmallfracs   = sizeof(smallfracs)   / sizeof(double);
const int nlargewindows = sizeof(largewindows) / sizeof(int);

const double leftmargin   = 0.10;
const double rightmargin  = 0.05;
const double topmargin    = 0.05;
const double splitfrac    = 0.40;
const double legfrac      = 0.80;
const double bottommargin = 0.15;

TLegend * leffs[2][3];
TCanvas * ceffs[2][3];
TPad    * peffs[2][3][4];
TH1D * heffs[nsmallfracs + 1][2][3];
TH1D * hcuts[nsmallfracs + 1][2][3];
Color_t effcols  [nsmallfracs + 1] = {kBlack, kMagenta, kBlue, kCyan, kOrange-1, kGreen+2, kGray+2, kRed};
Color_t effstyles[nsmallfracs + 1] = {20, 21, 22, 23, 24, 25, 26, 27};
void InitEffHists() {
  for(int ismall = 0; ismall < nsmallfracs + 1; ismall++) {
    for(int ien = 0; ien < 2; ien++) {
      for(int ieff = 0; ieff < 3; ieff++) {
	heffs[ismall][ien][ieff] = new TH1D(TString::Format("heff_%.1f_%d_%d", smallfracs[(ismall < nsmallfracs ? ismall : ismall+999)], ien, ieff), ";Large window size;Noise rejection efficiency (%)", 4001, 0, 2001);
	heffs[ismall][ien][ieff]->GetYaxis()->SetRangeUser(0,100);
	heffs[ismall][ien][ieff]->SetMarkerColor(effcols[ismall]);
	heffs[ismall][ien][ieff]->SetMarkerStyle(effstyles[ismall]);
	heffs[ismall][ien][ieff]->GetYaxis()->SetTitleOffset(0.8);
	heffs[ismall][ien][ieff]->GetXaxis()->SetLabelOffset(0.8);
	hcuts[ismall][ien][ieff] = new TH1D(TString::Format("hcut_%.1f_%d_%d", smallfracs[(ismall < nsmallfracs ? ismall : ismall+999)], ien, ieff), ";Large window size;Cut value", 4001, 0, 2001);
	hcuts[ismall][ien][ieff]->GetYaxis()->SetRangeUser(0,1);
	hcuts[ismall][ien][ieff]->SetMarkerColor(effcols[ismall]);
	hcuts[ismall][ien][ieff]->SetMarkerStyle(effstyles[ismall]);
	hcuts[ismall][ien][ieff]->GetYaxis()->SetTitleOffset(0.8);
	if(ismall == 0) {
	  ceffs[ien][ieff] = new TCanvas();
	  ceffs[ien][ieff]->SetCanvasSize(700,1000);
	  peffs[ien][ieff][0] = new TPad(TString::Format("peffs_0_%d_%d", ien, ieff), "", 0, legfrac,   1, 1);
	  peffs[ien][ieff][1] = new TPad(TString::Format("peffs_1_%d_%d", ien, ieff), "", 0, splitfrac, 1, legfrac);
	  peffs[ien][ieff][2] = new TPad(TString::Format("peffs_2_%d_%d", ien, ieff), "", 0, 0,         1, splitfrac);
	  peffs[ien][ieff][3] = new TPad(TString::Format("peffs_2_%d_%d", ien, ieff), "", 0, 0,         1, splitfrac);
	  peffs[ien][ieff][3]->SetFillStyle(4000); //will be transparent
	  peffs[ien][ieff][2]->SetTicks(1,0);
	  for(int ipad = 0; ipad < 4; ipad++)
	    peffs[ien][ieff][ipad]->Draw();
	  SetMargins(*(peffs[ien][ieff][1]), leftmargin, rightmargin, topmargin, 0.025);
	  SetMargins(*(peffs[ien][ieff][2]), leftmargin, rightmargin, 0.025,     bottommargin);
	  SetMargins(*(peffs[ien][ieff][3]), leftmargin, rightmargin, 0.025,     bottommargin);
	  leffs[ien][ieff] = new TLegend(0, 0, 1, 1);
	  leffs[ien][ieff]->SetNColumns(3);
	  leffs[ien][ieff]->SetHeader(TString::Format("Noise rejection efficiency for cut allowing %d%% of %d MeV signal", (int)(100 * efficiencies[ieff]), energies[ien+1]));
	  cout << TString::Format("Noise rejection efficiency for cut allowing %d%% of %d MeV signal", (int)(100 * efficiencies[ieff]), energies[ien+1]) << endl;
	}
	if(ismall == 1 || ismall == 3)
	  continue;
	else if(ismall < nsmallfracs)
	  leffs[ien][ieff]->AddEntry(heffs[ismall][ien][ieff], TString::Format("Small = large * %.1f", smallfracs[ismall]), "p");
	else
	  leffs[ien][ieff]->AddEntry(heffs[ismall][ien][ieff], "NHITS", "p");
      }//ieff
    }//ien
  }//ismall
}

void DrawEffHists(bool oneslice) {
  TString plotname = (oneslice ? "efficiency_oneslice.pdf" : "efficiency.pdf");
  ceffs[0][0]->SaveAs(TString::Format("%s[", plotname.Data()));
  for(int ien = 0; ien < 2; ien++) {
    for(int ieff = 0; ieff < 3; ieff++) {
      /*
      cout << "Drawing ccuts" << endl;
      ccuts[ien][ieff]->Draw();
      cout << "Drawn ccuts" << endl;
      */
      for(int ismall = 0; ismall < nsmallfracs + 1; ismall++) {
	if(ismall == 1 || ismall == 3)
	  continue;
	peffs[ien][ieff][1]->cd();
	if(ismall == 0)
	  heffs[ismall][ien][ieff]->Draw("P");
	else
	  heffs[ismall][ien][ieff]->Draw("P SAME");
	/*
	cout << "drawn heffs" << endl;
	ccuts[ien][ieff]->cd(0);
	cout << "cded to cuts pad" << endl;
	*/
	peffs[ien][ieff][2]->cd();
	if(ismall == 0)
	  hcuts[ismall][ien][ieff]->Draw("P");
	else if(ismall == nsmallfracs) {
	  /*
	    peffs[ien][ieff][3]->cd();
	    double ymin = hcuts[0][ien][ieff]->GetMinimum() - 5;
	    double ymax = hcuts[0][ien][ieff]->GetMaximum() + 5;
	    double dy   = (ymax-ymin) / (1 - (0.025 + bottommargin));
	    double xmin = hcuts[0][ien][ieff]->GetXaxis()->GetXmin();
	    double xmax = hcuts[0][ien][ieff]->GetXaxis()->GetXmax();
	    double dx   = (ymax-ymin) / (1 - (rightmargin + leftmargin));
	    peffs[ien][ieff][3]->Range(xmin - leftmargin  * dx, ymin - bottommargin * dy,
	    xmax + rightmargin * dx, ymax + 0.025        * dy);
	    hcuts[0][ien][ieff]->Draw("][sames");
	    TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,0,100,50510,"+L");
	    axis->Draw();
	  */
	  peffs[ien][ieff][2]->cd();
	  peffs[ien][ieff][2]->Update();
	  double max = 1.1 * hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMaximumBin());
	  double min = 0.9 * hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMinimumBin());
	  for(int ix = 0; ix <= hcuts[ismall][ien][ieff]->GetNbinsX() + 1; ix++)
	    if(hcuts[ismall][ien][ieff]->GetBinContent(ix))
	      cout << "Bin " << ix << " content " << hcuts[ismall][ien][ieff]->GetBinContent(ix) 
		   << " min " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMinimumBin()) 
		   << " max " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMaximumBin()) 
		   << endl;
	  cout << "HIST min " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMinimumBin()) 
	       << " max " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMaximumBin()) 
	       << " scaling by " << gPad->GetUymax() / max << endl;
	  hcuts[ismall][ien][ieff]->Scale(gPad->GetUymax() / max);
	  cout << "HIST min " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMinimumBin()) 
	       << " max " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMaximumBin()) << endl;
	  hcuts[ismall][ien][ieff]->Draw("P SAME");
	  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), min, max, 510, "+L");
	  cout << "Ux " << gPad->GetUxmin() << ", " << gPad->GetUxmax() << endl
	       << "Uy " << gPad->GetUymin() << ", " << gPad->GetUymax() << endl
	       << "min,max " << min << ", " << max << endl;
	  axis->SetLabelFont(132);
	  axis->SetLabelSize(0.05);
	  axis->Draw();

	  for(int ix = 0; ix <= hcuts[ismall][ien][ieff]->GetNbinsX() + 1; ix++)
	    if(hcuts[ismall][ien][ieff]->GetBinContent(ix))
	      cout << "Bin " << ix << " content " << hcuts[ismall][ien][ieff]->GetBinContent(ix) 
		   << " min " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMinimumBin()) 
		   << " max " << hcuts[ismall][ien][ieff]->GetBinContent(hcuts[ismall][ien][ieff]->GetMaximumBin()) 
		   << endl;
	}
	else
	  hcuts[ismall][ien][ieff]->Draw("P SAME");
	//cout << heffs[ismall][ien][ieff]->GetEntries() << "\t" << heffs[ismall][ien][ieff]->GetMean() << endl;
      }//ismall
      peffs[ien][ieff][0]->cd();
      leffs[ien][ieff]->Draw();
      ceffs[ien][ieff]->SaveAs(TString::Format("%s", plotname.Data()));
    }//ieff
  }//ien
  ceffs[0][0]->SaveAs(TString::Format("%s]", plotname.Data()));
}



double GetBkgRejection(TH1D * sig, TH1D * bkg, double sigeff)
{
  //first, find the bin corresponding to the cut value that gives the given signal efficiency
  // need the right hand side ("allowed") to be above the given sigeff
  int cutbin = -1;
  for(int i = 0; i <= sig->GetNbinsX(); i++) {
    double integral = sig->Integral(0, i);
    double efficiency = integral / (double)sig->Integral();//GetEntries();
    cout << "Efficiency from bin 0 to bin " << i << " is " << efficiency
	 << "\t looking for " << (1 - sigeff) << "\t" << sigeff << endl;
    if(efficiency > (1. - sigeff)) {
      cout << "SIGNAL ";
      cout << "Efficiency from bin 0 to bin " << i << " is " << efficiency << endl;
      cout << "Setting cut bin to " << i - 1 << endl;
      cutbin = i - 1;
      break;
    }
  }//i
  if(cutbin < 0) {
    cerr << "Could not find a good cut value" << endl;
    exit(-1);
  }
  //now, see what background rejection we get
  double integral = bkg->Integral(0, cutbin);
  double rejection = integral / (double)bkg->GetEntries();
  return rejection;
}

double GetBkgRejectionTree(int isig, double sigeff, bool oneslice)
{
  assert(isig == 1 || isig == 2);
  const char * one = (oneslice ? "_1slice" : "");
  trees[0]   ->Draw(TString::Format("var_max_itc%s>>h1itc_fine_bkg(200000,0,2)", one), "", "GOFF");
  trees[isig]->Draw(TString::Format("var_max_itc%s>>h1itc_fine_sig(200000,0,2)", one), "", "GOFF");

  if(h1itc_fine_bkg->GetEntries() <= 0 || h1itc_fine_sig->GetEntries() <= 0)
    return -99;
  return GetBkgRejection(h1itc_fine_sig, h1itc_fine_bkg, sigeff);
}

void GetBkgRejectionTree2(ofstream & ss, bool oneslice, int ismall, double largewindow)
{
  vector<double> bkg, sig1, sig2;
  if(oneslice) {
    for(int i = 0; i < trees[0]->GetEntries(); i++) {
      trees[0]->GetEntry(i);
      bkg.push_back(var_max_itc_1slice[0]);
    }//i
    for(int i = 0; i < trees[1]->GetEntries(); i++) {
      trees[1]->GetEntry(i);
      sig1.push_back(var_max_itc_1slice[1]);
    }//i
    for(int i = 0; i < trees[2]->GetEntries(); i++) {
      trees[2]->GetEntry(i);
      sig2.push_back(var_max_itc_1slice[2]);
    }//i
  }
  else {
    for(int i = 0; i < trees[0]->GetEntries(); i++) {
      trees[0]->GetEntry(i);
      bkg.push_back(var_max_itc[0]);
    }//i
    for(int i = 0; i < trees[1]->GetEntries(); i++) {
      trees[1]->GetEntry(i);
      sig1.push_back(var_max_itc[1]);
    }//i
    for(int i = 0; i < trees[2]->GetEntries(); i++) {
      trees[2]->GetEntry(i);
      sig2.push_back(var_max_itc[2]);
    }//i
  }
  sort(sig1.begin(), sig1.end());
  sort(sig2.begin(), sig2.end());
  sort(bkg.begin(), bkg.end());

  for(int ien = 0; ien < 2; ien++) {
    for(int ieff = 0; ieff < 3; ieff++) {
      double sigeff = efficiencies[ieff];
      //get the cut value
      double cut = (!ien ? sig1[(1 - sigeff) * sig1.size()] : sig2[(1 - sigeff) * sig2.size()]);

      //find the position in the list of sig where this cut is
      vector<double>::iterator sig_cut_it = std::lower_bound((!ien ? sig1.begin() : sig2.begin()), (!ien ? sig1.end() : sig2.end()), cut);
      //and the number of sig events below this value
      int nsig_below = sig_cut_it - (!ien ? sig1.begin() : sig2.begin());
      //and the sacrifice
      double sacrifice = (nsig_below / (double)(!ien ? sig1.size() : sig2.size())) * 100;
      
      //find the position in the list of background where this cut is
      vector<double>::iterator bkg_cut_it = std::lower_bound(bkg.begin(), bkg.end(), cut);
      //and the number of bkg events below this value
      int nbkg_below = bkg_cut_it - bkg.begin();
      //and the rejection efficiency
      double rejection = (nbkg_below / (double)bkg.size()) * 100;

      ss << " & " << rejection;
      heffs[ismall][ien][ieff]->SetBinContent(heffs[ismall][ien][ieff]->FindBin(largewindow), rejection);
      hcuts[ismall][ien][ieff]->SetBinContent(hcuts[ismall][ien][ieff]->FindBin(largewindow), cut);

      bool verbose = 1;
      if(verbose)
	cout << "Noise rejection for " << energies[ien+1] << " MeV @ " << 100 * efficiencies[ieff] << "% with " 
	     << largewindow << " ns * " << smallfracs[ismall]
	     << " is " << rejection << "%" << endl
	     << "Cut at " << cut << endl
	     << "Bkg below " << nbkg_below << " out of " << bkg.size() << endl
	     << "Sig below " << nsig_below << " out of " << (!ien ? sig1.size() : sig2.size()) << " sacrificing " << sacrifice << "%"
	     << endl;
    }//ieff
  }//ien
}

void GetBkgRejectionTree2NDigits(bool oneslice, double largewindow)
{
  const int ismall = nsmallfracs;
  vector<double> bkg, sig1, sig2;
  if(oneslice) {
    for(int i = 0; i < trees[0]->GetEntries(); i++) {
      trees[0]->GetEntry(i);
      bkg.push_back(var_max_ndigits_1slice[0]);
    }//i
    for(int i = 0; i < trees[1]->GetEntries(); i++) {
      trees[1]->GetEntry(i);
      sig1.push_back(var_max_ndigits_1slice[1]);
    }//i
    for(int i = 0; i < trees[2]->GetEntries(); i++) {
      trees[2]->GetEntry(i);
      sig2.push_back(var_max_ndigits_1slice[2]);
    }//i
  }
  else {
    for(int i = 0; i < trees[0]->GetEntries(); i++) {
      trees[0]->GetEntry(i);
      bkg.push_back(var_max_ndigits[0]);
    }//i
    for(int i = 0; i < trees[1]->GetEntries(); i++) {
      trees[1]->GetEntry(i);
      sig1.push_back(var_max_ndigits[1]);
    }//i
    for(int i = 0; i < trees[2]->GetEntries(); i++) {
      trees[2]->GetEntry(i);
      sig2.push_back(var_max_ndigits[2]);
    }//i
  }
  sort(sig1.begin(), sig1.end());
  sort(sig2.begin(), sig2.end());
  sort(bkg.begin(), bkg.end());

  for(int ien = 0; ien < 2; ien++) {
    for(int ieff = 0; ieff < 3; ieff++) {
      double sigeff = efficiencies[ieff];
      //get the cut value
      double cut = (!ien ? sig1[(1 - sigeff) * sig1.size()] : sig2[(1 - sigeff) * sig2.size()]);

      //find the position in the list of sig where this cut is
      vector<double>::iterator sig_cut_it = std::lower_bound((!ien ? sig1.begin() : sig2.begin()), (!ien ? sig1.end() : sig2.end()), cut);
      //and the number of sig events below this value
      int nsig_below = sig_cut_it - (!ien ? sig1.begin() : sig2.begin());
      //and the sacrifice
      double sacrifice = (nsig_below / (double)(!ien ? sig1.size() : sig2.size())) * 100;

      //find the position in the list of background where this cut is
      vector<double>::iterator bkg_cut_it = std::lower_bound(bkg.begin(), bkg.end(), cut);
      //and the number of bkg events below this value
      int nbkg_below = bkg_cut_it - bkg.begin();
      //and the rejection efficiency
      double rejection = ((nbkg_below / (double)bkg.size())) * 100;
      //ss << " & " << rejection;
      heffs[ismall][ien][ieff]->SetBinContent(heffs[ismall][ien][ieff]->FindBin(largewindow), rejection);
      hcuts[ismall][ien][ieff]->SetBinContent(hcuts[ismall][ien][ieff]->FindBin(largewindow), cut);

      cout << "Noise rejection for " << energies[ien+1] << " MeV @ " << 100 * efficiencies[ieff] << "% with " 
	   << largewindow << " ns "
	   << " is " << rejection << "%" << endl
	   << "Cut at " << cut << endl
	   << "Bkg below " << nbkg_below << " out of " << bkg.size() << endl
	   << "Sig below " << nsig_below << " out of " << (!ien ? sig1.size() : sig2.size()) << " sacrificing " << sacrifice << "%"
	   << endl;
    }//ieff
  }//ien
}

void get_trigger_efficiencies(TString filenames, bool oneslice = false)
{
  gStyle->SetOptStat(0);
  InitEffHists();

  TString fname;
  Ssiz_t from = 0;
  vector<TFile*> files;
  while(filenames.Tokenize(fname, from, "[$]")) {
    cout << fname << endl;
    files.push_back(new TFile(fname.Data(), "READ"));
  }//loop over filenames
  const int nfiles = files.size();
  assert(nfiles == 3);
  TH1D * hists[nfiles];
  TH1D * hists_ndigits[nfiles];

  TCanvas c, cndigits;
  SetMargins(c);
  SetMargins(cndigits);
  TString plotname = (oneslice ? "plots_itc_oneslice.pdf" : "plots_itc.pdf");
  TString plotname_ndigits = (oneslice ? "plots_ndigits_oneslice.pdf" : "plots_ndigits.pdf");
  c.SaveAs(TString::Format("%s[", plotname.Data()));
  cndigits.SaveAs(TString::Format("%s[", plotname_ndigits.Data()));
  TLegend l(0.15,0.75,0.45,0.95);
  TLegend lndigits(0.15,0.75,0.35,0.95);

  ofstream ss((oneslice ? "itc_table_oneslice.tex" : "itc_table.tex"));
  ss << "\\tiny"
     << "\\begin{tabular}{| l l | c c c | c c c |}" << endl
     << "\\hline" << endl
     << "Big window & Small window & \\multicolumn{6}{c|}{Rejection @ efficiency} \\\\" << endl
     << " & & \\multicolumn{3}{c|}{4 MeV} & \\multicolumn{3}{c|}{10 MeV} \\\\" << endl
     << " & & 50\\% & 90\\% & 100\\% & 50\\% & 90\\% & 100\\% \\\\" << endl
     << "\\hline" << endl;

  for(int ilarge = 0; ilarge < nlargewindows; ilarge++) {
    int largewindow = largewindows[ilarge];
    //get the ndigits histogram name & setup plotting
    TString hname_ndigits = TString::Format("h1_max_ndigits%d", largewindow);
    if(oneslice)
      hname_ndigits += "_1slice";
    lndigits.Clear();
    lndigits.SetHeader(TString::Format("NHITS %d ns", largewindow));

    for(int ismall = 0; ismall < nsmallfracs; ismall++) {
      int smallwindow = smallfracs[ismall]*largewindows[ilarge];

      //get the itc histogram name & setup plotting
      TString hname = TString::Format("h1_max_itc%d_%d_0", smallwindow, largewindow);
      TString tname = TString::Format("t_itc%d_%d_0", smallwindow, largewindow);
      if(oneslice)
	hname += "_1slice";
      l.Clear();
      l.SetHeader(TString::Format("ITC %d ns %d ns", smallwindow, largewindow));

      //loop over files
      for(int ifile = 0; ifile < nfiles; ifile++) {
	//get tree
	trees[ifile] = 0;
	files[ifile]->GetObject(tname.Data(), trees[ifile]);
	if(!trees[ifile]) {
	  cout << tname << " not found" << endl;
	  continue;
	}
	trees[ifile]->SetBranchAddress("var_max_ndigits", &(var_max_ndigits[ifile]));
	trees[ifile]->SetBranchAddress("var_max_ndigits_time", &(var_max_ndigits_time[ifile]));
	trees[ifile]->SetBranchAddress("var_max_ndigits_1slice", &(var_max_ndigits_1slice[ifile]));
	trees[ifile]->SetBranchAddress("var_max_itc", &(var_max_itc[ifile]));
	trees[ifile]->SetBranchAddress("var_max_itc_time", &(var_max_itc_time[ifile]));
	trees[ifile]->SetBranchAddress("var_max_itc_1slice", &(var_max_itc_1slice[ifile]));

	//get itc histograms
	hists[ifile] = 0;
	files[ifile]->GetObject(hname.Data(), hists[ifile]);
	if(!hists[ifile]) {
	  cout << hname << " not found" << endl;
	  continue;
	}
	//plot itc
	c.cd(0);
	hists[ifile]->SetLineColor(cols[ifile]);
	hists[ifile]->SetLineStyle(styles[ifile]);
	if(ifile == 0)
	  hists[ifile]->Draw();
	else
	  hists[ifile]->Draw("SAME");
	l.AddEntry(hists[ifile], TString::Format("%d MeV", energies[ifile]), "l");

	//get ndigits histograms
	if(ismall == 0) {
	  hists_ndigits[ifile] = 0;
	  files[ifile]->GetObject(hname_ndigits.Data(), hists_ndigits[ifile]);
	  if(!hists_ndigits[ifile]) {
	    cout << hname << " not found" << endl;
	    continue;
	  }
	  //plot ndigits
	  cndigits.cd(0);
	  hists_ndigits[ifile]->SetLineColor(cols[ifile]);
	  hists_ndigits[ifile]->SetLineStyle(styles[ifile]);
	  if(ifile == 0)
	    hists_ndigits[ifile]->Draw();
	  else
	    hists_ndigits[ifile]->Draw("SAME");
	  lndigits.AddEntry(hists_ndigits[ifile], TString::Format("%d MeV", energies[ifile]), "l");
	}//ismall == 0 (ndigits stuff)
      }//ifile
      c.cd(0);
      l.Draw();
      c.SaveAs(TString::Format("%s", plotname.Data()));

      if(ismall == 0) {
	cndigits.cd(0);
	lndigits.Draw();
	cndigits.SaveAs(TString::Format("%s", plotname_ndigits.Data()));
	GetBkgRejectionTree2NDigits(oneslice, largewindow);
      }

      //now all the histograms are here, work out efficiencies
      ss << "\\hline" << endl
	 << largewindow << " & " << smallwindow;
      GetBkgRejectionTree2(ss, oneslice, ismall, largewindow);
      /*
      for(int isig = 1; isig <= 2; isig++) //4 MeV, 10 MeV
	for(int ieff = 0; ieff < 3; ieff++) //50%, 90%, 99%
	  //ss << " & " << GetBkgRejection(hists[isig], hists[0], efficiencies[ieff]);
	  ss << " & " << GetBkgRejectionTree(isig, efficiencies[ieff], oneslice);
      */
      ss << " \\\\" << endl;
    }//ismall
  }//ilarge

  c.SaveAs(TString::Format("%s]", plotname.Data()));
  cndigits.SaveAs(TString::Format("%s]", plotname_ndigits.Data()));
  ss << "\\hline" << endl;
  ss.close();

  DrawEffHists(oneslice);
}
