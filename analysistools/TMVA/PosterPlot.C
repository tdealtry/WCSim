#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TList.h"
#include "TCollection.h"
#include "TH1D.h"
#include "TLegend.h"

#include <iostream>

void PosterPlot(const char * pdfname = "tmvaplot.pdf", const char * filename = "TMVAPlots.root", const char * stackname = "hs_value_Bdt", const char * drawopt = "")
{
  TFile f(filename);
  THStack * hs = 0;
  f.GetObject(stackname, hs);

  TList * list = hs->GetHists();
  TIter next(list);
  TObject * o;
  TH1D * htemp;
  while(o = next.Next()) {
    htemp = (TH1D*)o;
    htemp->GetXaxis()->SetTitle("Separation Variable");
    //htemp->SetFillColor(kBlue);
    cout << htemp << endl;
  }


  TCanvas * c = new TCanvas();
  hs->SetTitle("bdt;Separation variable;Number of events");
  hs->SetMaximum(40000 * 0.97);
  hs->Draw(TString::Format("%s", drawopt));
  hs->GetYaxis()->SetNoExponent(false);
  hs->GetYaxis()->SetTitleOffset(1.1);
  hs->Draw(TString::Format("%s", drawopt));
  hs->Dump();

  /*
  // legend at the top
  TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), 1 - gPad->GetRightMargin(), 1.0, "");//, "F");
  l->SetNColumns(5);
  */

  /*
  // legend at the right
  c->SetRightMargin(0.25);
  TLegend * l = gPad->BuildLegend(1 - gPad->GetRightMargin(), 1 - gPad->GetTopMargin(), 1, 0, "");
  l->SetTextFont(132);
  cout << l->GetMargin() << endl;
  l->SetMargin(0.1);
  */

  //legend over text
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->SetBottomMargin(0.13);
  TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), 0.7, 0.62);
  //TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), gPad->GetLeftMargin() + 0.5, 1 - gPad->GetTopMargin() - 0.4, "");
  //TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), gPad->GetLeftMargin() + 0.3, 1 - gPad->GetTopMargin() - 0.5, "");
  //TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), gPad->GetLeftMargin() + 0.3, 1 - gPad->GetTopMargin() - 0.5, "");
  l->SetNColumns(2);
  l->SetTextFont(132);
  l->SetMargin(0.2);

  gPad->Update();
  c->SaveAs(pdfname);
}
