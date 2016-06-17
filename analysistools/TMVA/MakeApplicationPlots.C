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
  TCanvas c;
  c.SetTopMargin(0.2);
  c.SaveAs("TMVAPlots.pdf[");

  const int nvariables = 1;
  const int nenergies  = 10;
  TH1D * h_value[nvariables][nenergies];
  THStack * hs_value[nvariables];
  TH1D * h_efficiency[nvariables][nenergies];
  THStack * hs_efficiency[nvariables];
  TH1D * h_purity[nvariables][nenergies];
  THStack * hs_purity[nvariables];
  TH1D * h_ep[nvariables][nenergies];
  THStack * hs_ep[nvariables];
  TEventList * eventlists[nenergies], * eventlists_sig[nenergies], * eventlists_bkg[nenergies];
  Long64_t nentries[nenergies];
  double max_axis[nvariables];
  int iv = 0;
  Ssiz_t from = 0;
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
	t->SetEventList(0);
	TString elistname = TString::Format("elist%d", ie);
	eventlists[ie] = new TEventList(elistname);
	t->Draw(TString::Format(">>%s", elistname.Data()),
		TString::Format("TRUE_e_energy>=%d&&TRUE_e_energy<%d", ie, ie+1));
	nentries[ie] = eventlists[ie]->GetN();
	cout << "Event list " << elistname << " has " << nentries[ie] << " events" << endl;

	//and for the purities...
	t->SetEventList(0);
	TString elistname_sig = TString::Format("elist_sig%d", ie);
	eventlists_sig[ie] = new TEventList(elistname_sig);
	t->Draw(TString::Format(">>%s", elistname_sig.Data()), TString::Format("TRUE_e_energy>=%d", ie));
	t->SetEventList(0);
	TString elistname_bkg = TString::Format("elist_bkg%d", ie);
	eventlists_bkg[ie] = new TEventList(elistname_bkg);
	t->Draw(TString::Format(">>%s", elistname_bkg.Data()), TString::Format("TRUE_e_energy<%d", ie));
      }
      t->SetEventList(eventlists[ie]);

      //make the histograms pretty
      TString htitle_ecut = TString::Format("E = %d-%d MeV", ie, ie+1);
      TString htitle_ecut_purity = TString::Format("E >= %d MeV", ie);
      //value
      TString hname = TString::Format("h_value_%s_e%d", variable.Data(), ie);
      h_value[iv][ie] = new TH1D(hname, htitle_ecut, 100, -max_axis[iv], +max_axis[iv]);
      h_value[iv][ie]->SetFillColor(1 + ie);
      h_value[iv][ie]->SetMarkerColor(h_value[iv][ie]->GetLineColor());
      h_value[iv][ie]->SetMarkerStyle(1);
      //efficiency
      hname = TString::Format("h_efficiency_%s_e%d", variable.Data(), ie);
      h_efficiency[iv][ie] = new TH1D(hname, htitle_ecut, 100, -max_axis[iv], +max_axis[iv]);
      h_efficiency[iv][ie]->SetLineColor(ie == 9 ? kGray : 1 + ie);
      h_efficiency[iv][ie]->SetMarkerColor(h_efficiency[iv][ie]->GetLineColor());
      h_efficiency[iv][ie]->SetMarkerStyle(1);
      //purity
      hname = TString::Format("h_purity_%s_e%d", variable.Data(), ie);
      h_purity[iv][ie] = new TH1D(hname, htitle_ecut_purity, 100, -max_axis[iv], +max_axis[iv]);
      h_purity[iv][ie]->SetLineColor(ie == 9 ? kGray : 1 + ie);
      h_purity[iv][ie]->SetMarkerColor(h_purity[iv][ie]->GetLineColor());
      h_purity[iv][ie]->SetMarkerStyle(1);
      //ep
      hname = TString::Format("h_ep_%s_e%d", variable.Data(), ie);
      h_ep[iv][ie] = new TH1D(hname, htitle_ecut, 100, -max_axis[iv], +max_axis[iv]);
      h_ep[iv][ie]->SetLineColor(ie == 9 ? kGray : 1 + ie);
      h_ep[iv][ie]->SetMarkerColor(h_ep[iv][ie]->GetLineColor());
      h_ep[iv][ie]->SetMarkerStyle(1);

      //fill the value histogram
      t->Draw(TString::Format("%s>>%s", variable.Data(), hname.Data()),
	      "", "goff");

      //fill the efficiency histogram
      for(int ib = 1; ib <= h_efficiency[iv][ie]->GetNbinsX(); ib++) {
	//cout << h_efficiency[iv][ie]->GetBinLowEdge(ib) << "\t"
	//   << nentries_passed << "\t" << nentries[ie] << "\t" << efficiency << endl;
      }//ib

      //fill the efficiency & purity & ep histograms
      for(int ib = 1; ib <= h_purity[iv][ie]->GetNbinsX(); ib++) {
	TString varcut = TString::Format("%s>=%f", variable.Data(), h_purity[iv][ie]->GetBinLowEdge(ib));
	//efficiency
	t->SetEventList(eventlists[ie]);
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

	h_ep[iv][ie]->SetBinContent(ib, efficiency * purity);

	if(verbose > 2)
	  cout << h_purity[iv][ie]->GetBinLowEdge(ib) << endl
	       << "\tefficiency " << efficiency << "\t=\t" << nentries_passed << "\t/\t" << nentries[ie] << endl
	       << "\tpurity     " << purity     << "\t=\t" << nentries_passed_sig << "\t/\t" 
	       << nentries_passed_sig + nentries_passed_bkg << endl;
      }//ib

      //populate the stacks
      hs_value[iv]->Add(h_value[iv][ie]);
      hs_efficiency[iv]->Add(h_efficiency[iv][ie]);
      hs_purity[iv]->Add(h_purity[iv][ie]);
      hs_ep[iv]->Add(h_ep[iv][ie]);

      //cout << "finishing energy loop " << ie << endl;
    }//ie
    t->SetEventList(0);


    //draw the histograms & save them
    hs_value[iv]->Write();
    hs_efficiency[iv]->Write();
    hs_purity[iv]->Write();

    hs_value[iv]->Draw();
    hs_value[iv]->GetYaxis()->SetTitleOffset(1.1);
    TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), 1 - gPad->GetRightMargin(), 1.0, "");//, "F");
    l->SetNColumns(5);
    gPad->Update();
    c.SaveAs("TMVAPlots.pdf");

    hs_efficiency[iv]->Draw("NOSTACK");
    hs_efficiency[iv]->GetYaxis()->SetTitleOffset(1.1);
    l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), 1 - gPad->GetRightMargin(), 1.0, "");//, "F");
    l->SetNColumns(5);
    gPad->Update();
    c.SaveAs("TMVAPlots.pdf");

    hs_purity[iv]->Draw("NOSTACK");
    hs_purity[iv]->GetYaxis()->SetTitleOffset(1.1);
    l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), 1 - gPad->GetRightMargin(), 1.0, "");//, "F");
    l->SetNColumns(5);
    gPad->Update();
    c.SaveAs("TMVAPlots.pdf");

    hs_ep[iv]->Draw("NOSTACK");
    hs_ep[iv]->GetYaxis()->SetTitleOffset(1.1);
    l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), 1 - gPad->GetRightMargin(), 1.0, "");//, "F");
    l->SetNColumns(5);
    gPad->Update();
    c.SaveAs("TMVAPlots.pdf");

    //loop
    iv++;
  }//loop over variables

  c.SaveAs("TMVAPlots.pdf]");
  fout.Close();
  f.Close();
}
