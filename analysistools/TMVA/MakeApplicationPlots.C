#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TEventList.h"
#include "TMath.h"

void DrawStacks(THStack * hs, TString opt, TString pdfname, TCanvas * c);
TH1D * MakeHistogram(const int mode, const int ie,
		     const TString hname, const TString htitle, const int nbins, const double xmax);
TEventList * MakeEventList(TTree * t, const TString elistname, const TString ecut);

void MakeApplicationPlots(TString variables = "Bdt,BdtB,BdtG,BdtD,Cuts", int verbose = 1, const char * filename = "TMVApp.root")
{
  TFile f(filename);
  TTree * t = 0;
  f.GetObject("tmva_results", t);
  if(!t) {
    cerr << "Could not get tree" << endl;
    return;
  }
  cout << "Tree contains " << t->GetEntries() << " events" << endl;

  TFile fout("TMVAPlots.root", "RECREATE");
  TCanvas * c = new TCanvas();
  c->SetTopMargin(0.2);
  c->SaveAs("TMVAPlots.pdf[");

  const int nvariables = 10;
  const int nenergies  = 10;
  //book histograms
  TH1D * h_value[nvariables][nenergies];
  THStack * hs_value[nvariables];
  TH1D * h_efficiency[nvariables][nenergies];
  THStack * hs_efficiency[nvariables];
  TH1D * h_purity[nvariables][nenergies];
  THStack * hs_purity[nvariables];
  TH1D * h_ep[nvariables][nenergies];
  THStack * hs_ep[nvariables];
  //book everything else
  double max_axis[nvariables];
  TEventList * eventlists[nenergies], * eventlists_sig[nenergies], * eventlists_bkg[nenergies];
  Long64_t nentries[nenergies];
  //loop variables
  int iv = 0;
  Ssiz_t from = 0;  //dummy
  TString variable;
  while(variables.Tokenize(variable, from, ",")) {
    cout << variable << endl;
    hs_value[iv] = new THStack(TString::Format("hs_value_%s", variable.Data()),
			       TString::Format("%s;%s value;Number of events", variable.Data(), variable.Data()));
    hs_efficiency[iv] = new THStack(TString::Format("hs_efficiency_%s", variable.Data()),
			       TString::Format("%s;%s value;Efficiency", variable.Data(), variable.Data()));
    hs_purity[iv] = new THStack(TString::Format("hs_purity_%s", variable.Data()),
				TString::Format("%s;%s value;Purity", variable.Data(), variable.Data()));
    hs_ep[iv] = new THStack(TString::Format("hs_ep_%s", variable.Data()),
				TString::Format("%s;%s value;Efficiency #times Purity", variable.Data(), variable.Data()));
    double min = t->GetMinimum(variable);
    double max = t->GetMaximum(variable);
    max_axis[iv] = TMath::Ceil(TMath::Max(TMath::Abs(min), TMath::Abs(max)));

    for(int ie = 0; ie < nenergies; ie++) {
      //apply energy cut
      if(iv == 0) {
	eventlists[ie] = MakeEventList(t, TString::Format("elist%d", ie),
				       TString::Format("TRUE_e_energy>=%d&&TRUE_e_energy<%d", ie, ie+1));
	nentries[ie] = eventlists[ie]->GetN();

	//and for the purities...
	eventlists_sig[ie] = MakeEventList(t, TString::Format("elist_sig%d", ie),
					   TString::Format("TRUE_e_energy>=%d", ie));
	eventlists_bkg[ie] = MakeEventList(t, TString::Format("elist_bkg%d", ie),
					   TString::Format("TRUE_e_energy<%d", ie));
      }

      //make the histograms & make them pretty
      TString htitle_ecut = TString::Format("E = %d-%d MeV", ie, ie+1);
      TString htitle_ecut_purity = TString::Format("E >= %d MeV", ie);
      //value
      h_value[iv][ie] = MakeHistogram(0, ie, TString::Format("h_value_%s_e%d", variable.Data(), ie),
				      htitle_ecut, 100, max_axis[iv]);
      h_efficiency[iv][ie] = MakeHistogram(1, ie, TString::Format("h_efficiency_%s_e%d", variable.Data(), ie),
					   htitle_ecut, 100, max_axis[iv]);
      h_purity[iv][ie] = MakeHistogram(1, ie, TString::Format("h_purity_%s_e%d", variable.Data(), ie),
				       htitle_ecut_purity, 100, max_axis[iv]);
      h_ep[iv][ie] = MakeHistogram(1, ie, TString::Format("h_ep_%s_e%d", variable.Data(), ie),
				   htitle_ecut, 100, max_axis[iv]);

      //fill the value histogram
      cout << "Event list " << eventlists[ie]->GetTitle() << " has " << nentries[ie] << " events" << endl;
      t->SetEventList(0);
      t->SetEventList(eventlists[ie]);
      t->Draw(TString::Format("%s>>%s", variable.Data(), h_value[iv][ie]->GetName()),
	      "", "goff");

      //fill the efficiency & purity & ep histograms
      for(int ib = 1; ib <= h_purity[iv][ie]->GetNbinsX(); ib++) {
	t->SetEventList(eventlists[ie]);
	TString varcut = TString::Format("%s>=%f", variable.Data(), h_purity[iv][ie]->GetBinLowEdge(ib));
	//efficiency
	Long64_t nentries_passed = t->GetEntries(varcut);
	double efficiency = (double)nentries_passed / (double)nentries[ie];
	h_efficiency[iv][ie]->SetBinContent(ib, efficiency);
	//purity
	t->SetEventList(eventlists_sig[ie]);
	Long64_t nentries_passed_sig = t->GetEntries(varcut);
	t->SetEventList(eventlists_bkg[ie]);
	Long64_t nentries_passed_bkg = t->GetEntries(varcut);
	double purity = ((nentries_passed_sig + nentries_passed_bkg == 0) ? 0 
			 : ((double)nentries_passed_sig / (double)(nentries_passed_sig + nentries_passed_bkg)));
	h_purity[iv][ie]->SetBinContent(ib, purity);
	//ep
	h_ep[iv][ie]->SetBinContent(ib, efficiency * purity);

	if(verbose > 2)
	  cout << h_purity[iv][ie]->GetBinLowEdge(ib) << endl
	       << "\tefficiency " << efficiency << "\t=\t" << nentries_passed << "\t/\t" << nentries[ie] << endl
	       << "\tpurity     " << purity     << "\t=\t" << nentries_passed_sig << "\t/\t" 
	       << nentries_passed_sig + nentries_passed_bkg << endl;
      }//ib

      //populate the stacks
      hs_value     [iv]->Add(h_value[iv][ie]);
      hs_efficiency[iv]->Add(h_efficiency[iv][ie]);
      hs_purity    [iv]->Add(h_purity[iv][ie]);
      hs_ep        [iv]->Add(h_ep[iv][ie]);

      //save the histograms
      h_value     [iv][ie]->Write();
      h_efficiency[iv][ie]->Write();
      h_purity    [iv][ie]->Write();
      h_ep        [iv][ie]->Write();

      //cout << "finishing energy loop " << ie << endl;
    }//ie
    t->SetEventList(0);


    //draw the histograms & save them
    DrawStacks(hs_value[iv], "", "TMVAPlots.pdf", c);
    DrawStacks(hs_efficiency[iv], "NOSTACK", "TMVAPlots.pdf", c);
    DrawStacks(hs_purity[iv], "NOSTACK", "TMVAPlots.pdf", c);
    DrawStacks(hs_ep[iv], "NOSTACK", "TMVAPlots.pdf", c);

    //loop
    iv++;
  }//loop over variables

  c->SaveAs("TMVAPlots.pdf]");
  fout.Close();
  f.Close();
}

void DrawStacks(THStack * hs, TString opt, TString pdfname, TCanvas * c)
{
  hs->Write();
  hs->Draw(opt);
  hs->GetYaxis()->SetTitleOffset(1.1);
  TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), 1 - gPad->GetRightMargin(), 1.0, "");//, "F");
  l->SetNColumns(5);
  gPad->Update();
  c->SaveAs(pdfname);
}


TH1D * MakeHistogram(const int mode, const int ie,
		     const TString hname, const TString htitle, const int nbins, const double xmax)
{
  TH1D * h = new TH1D(hname, htitle, nbins, -xmax, +xmax);

  if(mode == 0) {
    h->SetFillColor(1 + ie);
    h->SetMarkerColor(h->GetLineColor());
    h->SetMarkerStyle(1);
  }
  else if(mode == 1) {
    h->SetLineColor(ie == 9 ? kGray : 1 + ie);
    h->SetMarkerColor(h->GetLineColor());
    h->SetMarkerStyle(1);
  }
  return h;
}

TEventList * MakeEventList(TTree * t, const TString elistname, const TString ecut)
{
  TEventList * elist = new TEventList(elistname, ecut);
  t->SetEventList(0);
  t->Draw(TString::Format(">>%s", elistname.Data()), ecut);
  return elist;
}
