void PosterPlot(const char * pdfname = "tmvaplot.pdf", const char * filename = "TMVAPlots.root", const char * stackname = "hs_value_Bdt", const char * drawopt = "")
{
  TFile f(filename);
  THStack * hs = 0;
  f.GetObject(stackname, hs);

  TCanvas * c = new TCanvas();
  c->DrawFrame(-1,0,1,40000);
  hs->Draw(TString::Format("%s", drawopt));
  hs->GetYaxis()->SetNoExponent(false);
  hs->GetYaxis()->SetTitleOffset(1.1);

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
  c->SetRightMargin(0.01);
  c->SetTopMargin(0.01);
  c->SetBottomMargin(0.13);
  TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), gPad->GetLeftMargin() + 0.5, 1 - gPad->GetTopMargin() - 0.4, "");
  //TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), gPad->GetLeftMargin() + 0.3, 1 - gPad->GetTopMargin() - 0.5, "");
  //TLegend * l = gPad->BuildLegend(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin(), gPad->GetLeftMargin() + 0.3, 1 - gPad->GetTopMargin() - 0.5, "");
  l->SetNColumns(2);
  l->SetTextFont(132);
  l->SetMargin(0.2);
  hs->GetYaxis()->SetRangeUser(0,40000);

  gPad->Update();
  c->SaveAs(pdfname);
}
