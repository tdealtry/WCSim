#include "TF1.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1D.h"
#include <iostream>

using namespace std;

int i = 0;
const int n20 = 21784;
const int n3  = 18604;
const int nM  = 10894;

void FitResult(TTree * t, TString axes, const char * selection, TCanvas * c2, const char * plotname)
{
  cout << axes << endl;
  c2->cd();
  t->Draw(TString::Format("%s>>hfun%d", axes.Data(), i), selection, "COLZ");
  TH1 * graph = (TH1*)gDirectory->Get(TString::Format("hfun%d", i));
  if(graph) {
    TObjArray * tokens = axes.Tokenize(":");
    graph->GetXaxis()->SetTitle(((TObjString*)tokens->At(1))->String());
    graph->GetYaxis()->SetTitle(((TObjString*)tokens->At(0))->String());
    TF1 * f = new TF1("f", "[0]*x + [1]");
    graph->Fit("f");
    TF1 * fitinfo = graph->GetFunction("f");
    cout << "y = " << fitinfo->GetParameter(0) << " x + " << fitinfo->GetParameter(1) << endl;
  }
  c2->SaveAs(plotname);
  i++;
}

void PlotHitId(TTree * t, TString axes, const char * selection, int nbins, double min, double max, TCanvas * c, int nev, const char * plotname)
{
  cout << axes << endl;
  c->cd();
  cout << selection << endl;
  t->Draw(TString::Format("%s>>hfun%d(%d,%d,%d)", axes.Data(), i, nbins, min, max), selection);
  TH1 * graph = (TH1*)gDirectory->Get(TString::Format("hfun%d", i));
  if(graph) {
    graph->GetXaxis()->SetTitle(axes.Data());
    graph->GetYaxis()->SetTitle(TString::Format("N hits per PMT across %d events", nev));
  }
  c->SaveAs(plotname);
  i++;
}

void PlotHitDist(TTree * t, const char * axes, TCanvas * c, const char * plotname)
{
  t->Draw(axes);
  TH1F * htemp = (TH1F*)gPad->GetPrimitive("htemp");
  cout << axes << endl
       << "\tMean  \t" << htemp->GetMean() << endl
       << "\tMedian\t" << htemp->GetBinCenter(htemp->GetMaximumBin()) << endl
       << "\tRMS   \t" << htemp->GetRMS() << endl;
  c->SaveAs(plotname);
}

void analyse_plots(const char * filename = "plots.root", const char * tag = "", bool mode=0)
{
  TFile fin(filename, "READ");
  TTree * event = 0;
  fin.GetObject("event", event);
  TTree * hit = 0;
  fin.GetObject("hit", hit);
  TTree * truehit = 0;
  fin.GetObject("truehit", truehit);
  TTree * global = 0;
  fin.GetObject("global", global);
  const char * plotname = Form("plots%s.pdf", tag);
  
  gStyle->SetOptFit(1);
  TCanvas * c = new TCanvas();
  TCanvas * c2 = new TCanvas();
  c2->SetRightMargin(0.15);
  c->SaveAs(Form("%s[", plotname));
  c->cd();
  event->Draw("X"); c->SaveAs(plotname);
  event->Draw("Y"); c->SaveAs(plotname);
  c2->cd();
  event->Draw("Y:X","","COLZ"); c2->SaveAs(plotname);  
  c->cd();
  event->Draw("Z"); c->SaveAs(plotname);
  event->Draw("E"); c->SaveAs(plotname);
  PlotHitDist(event, "NDigi20", c, plotname);
  PlotHitDist(event, "NDigiM", c, plotname);
  event->Draw("NTrue20"); c->SaveAs(plotname);
  event->Draw("NTrueM"); c->SaveAs(plotname);
  hit->Draw("Q"); c->SaveAs(plotname);
  truehit->Draw("PMTId_True"); c->SaveAs(plotname);
  truehit->Draw("mPMTId_True"); c->SaveAs(plotname);

  //do linear fits here
  c2->cd();
  FitResult(event, "NDigi20:E", "" ,c2, plotname);
  FitResult(event, "NDigiM:E", "" ,c2, plotname);
  FitResult(event, "NDigiM:NDigi20", "" ,c2, plotname);

  //calculate total number of events
  if(!mode)
    global->Draw("NEv>>hnew","NEvM>0||NEv20>0","goff");
  else
    global->Draw("NEv>>hnew","","goff");
  TH1 * hnew = (TH1*)gDirectory->Get("hnew");
  long nev = hnew->GetEntries() * hnew->GetMean();

  //hit plots
  //PlotHitId(hit, "PMTId",  TString::Format("(1/%f)*(mPMTId<0)",  (double)nev), n20, 0, n20, c);
  //PlotHitId(hit, "PMTId",  TString::Format("(1/%f)*(mPMTId>=0)", (double)nev), nM, n20, n20+nM, c);
  //PlotHitId(hit, "mPMTId", TString::Format("(1/%f)*(mPMTId>=0)", (double)nev), nM*19, 0, nM*19, c);
  if(!mode) {
    PlotHitId(hit, "PMTId",  "mPMTId<0",  n20, 0, n20,     c, nev, plotname);
    PlotHitId(hit, "PMTId",  "mPMTId>=0", nM*19, n20, n20+(nM*19), c, nev, plotname);
    PlotHitId(hit, "mPMTId", "mPMTId>=0", nM, 0, nM, c, nev, plotname);
  }
  else{
    PlotHitId(hit, "PMTId",  "",  n3, 0, n3,     c, nev, plotname);
    PlotHitId(hit, "PMTId",  "", nM*19, n20, n20+(nM*19), c, nev, plotname);
    PlotHitId(hit, "mPMTId", "", nM, 0, nM, c, nev, plotname);
  }

  //time since last truehit plot
  int thEv;
  int thPMTId, thmPMTId, thTrueTime;
  truehit->SetBranchAddress("Ev", &thEv);
  truehit->SetBranchAddress("TrueTime", &thTrueTime);
  truehit->SetBranchAddress("PMTId_True", &thPMTId);
  truehit->SetBranchAddress("mPMTId_True", &thmPMTId);
  float last_ev;
  int last_id, last_mid, last_time;
  int count_by_id = 1, count_by_mid = 1;
  const int ntimebins = mode ? 140 : 140;
  const int nhitsbins = mode ? 10 : 30;
  TH1D * hTimeSinceLast20 = new TH1D("hTimeSinceLast20", ";Time since last hit on this 20\" PMT (ns);N", ntimebins, 0, +ntimebins);
  TH1D * hTimeSinceLastM  = new TH1D("hTimeSinceLastM",  ";Time since last hit on this mPMT module (ns);N", ntimebins, 0, +ntimebins);
  TH1D * hTimeSinceLast3  = new TH1D("hTimeSinceLast3",  ";Time since last hit on this 3\" mPMT (ns);N", ntimebins, 0, ntimebins);
  TH1D * hNHitsPerPMT20 = new TH1D("hNHitsPerPMT20", ";Number of hits per event on this 20\" PMT;N", nhitsbins - 2, 2, nhitsbins);
  TH1D * hNHitsPerPMTM  = new TH1D("hNHitsPerPMTM",  ";Number of hits per event on this mPMT module;N", nhitsbins - 2, 2, nhitsbins);
  TH1D * hNHitsPerPMT3  = new TH1D("hNHitsPerPMT3",  ";Number of hits per event on this 3\" PMT;N", nhitsbins - 2, 2, nhitsbins);
  if(mode) {
    hTimeSinceLast20->GetXaxis()->SetTitle("Time since last hit on this OD PMT (ns)");
    hNHitsPerPMT20->GetXaxis()->SetTitle("Number of hits per event on this OD PMT");
  }
  for(long i = 0; i < truehit->GetEntries(); i++) {
    truehit->GetEntry(i);
    if(i) {
      //hit on same pmt (or mpmt module)
      if(last_ev == thEv && last_id == thPMTId) {
	//20 inch
	if(mode || thmPMTId < 0)
	  hTimeSinceLast20->Fill(TMath::Abs(thTrueTime - last_time));
	//mPMT
	else {
	  hTimeSinceLastM->Fill(thTrueTime - last_time);
	  if(last_mid == thmPMTId)
	    hTimeSinceLast3->Fill(TMath::Abs(thTrueTime - last_time));
	}//tube type
	count_by_id++;
      }//hit on same pmt (or mpmt module)
      else {
	//20 inch
	if(mode || thmPMTId < 0)
	  hNHitsPerPMT20->Fill(count_by_id);
	else
	  hNHitsPerPMTM->Fill(count_by_id);
	if(count_by_id > 20)
	  cout << count_by_id << endl;
	count_by_id = 1;
      }
      if(!mode) {
	if(last_ev == thEv && last_id == thPMTId && last_mid == thmPMTId && thmPMTId >= 0)
	  count_by_mid++;
	else {
	  hNHitsPerPMT3->Fill(count_by_mid);
	  count_by_mid = 1;
	}
      }
    }//i != 0
    if(count_by_mid > 1)
      cout << "ev num " << thEv << " pmt " << thPMTId << " mpmt " << thmPMTId << " count " << count_by_id << " mcount " << count_by_mid << endl;
    last_ev = thEv;
    last_time = thTrueTime;
    last_id = thPMTId;
    last_mid = thmPMTId;
  }//i
  hTimeSinceLast20->Draw(); c->SaveAs(plotname);
  //hTimeSinceLastM->Draw(); c->SaveAs(plotname); //this doesn't work properly. wrong logic
  hTimeSinceLast3->Draw(); c->SaveAs(plotname);
  hNHitsPerPMT20->Draw(); c->SaveAs(plotname);
  //hNHitsPerPMTM->Draw(); c->SaveAs(plotname); //this doesn't work properly. wrong logic
  hNHitsPerPMT3->Draw(); c->SaveAs(plotname);

  cout << hNHitsPerPMT3->GetEntries() << " entries in 3inch PMT nhit with mean " << hNHitsPerPMT3->GetMean() << endl;
  cout << hNHitsPerPMT20->GetEntries() << " entries in 20inch PMT nhit with mean " << hNHitsPerPMT20->GetMean() << endl;
  cout << hNHitsPerPMTM->GetEntries() << " entries in mPMT PMT nhit with mean " << hNHitsPerPMTM->GetMean() << endl;
  cout << hTimeSinceLast3->GetEntries() << " entries in 3inch PMT time since last with mean " << hTimeSinceLast3->GetMean() << endl;
  cout << hTimeSinceLast20->GetEntries() << " entries in 20inch PMT time since last with mean " << hTimeSinceLast20->GetMean() << endl;
  cout << hTimeSinceLastM->GetEntries() << " entries in mPMT PMT time since last with mean " << hTimeSinceLastM->GetMean() << endl;

  c->SaveAs(Form("%s]", plotname));

  cout << nev << " total events analysed" << endl;
}
