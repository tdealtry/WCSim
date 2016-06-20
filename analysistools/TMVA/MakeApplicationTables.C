#include <iostream>
#include <iomanip>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"

void MakeApplicationTables(TString variables = "Bdt,BdtB,BdtG,BdtD,Cuts", int verbose = 1, const char * filename = "TMVAPlots.root")
{
  TFile f(filename);
  
  const int nenergies = 10;
  //book histograms
  TH1D * h_value [nenergies] = {0};
  TH1D * h_efficiency [nenergies] = {0};
  TH1D * h_purity [nenergies] = {0};
  TH1D * h_ep [nenergies] = {0};
  //loop variables
  int iv = 0;
  Ssiz_t from = 0;  //dummy
  TString variable;
  while(variables.Tokenize(variable, from, ",")) {
    cout << variable << endl;

    //table header
    ofstream str(TString::Format("TMVAPlots_%d.tex", iv));
    str << "\\begin{tabular}{| l | c c c  | c | c |}" << endl
	<< "\\hline" << endl
	<< "Energy & $\\eta\\rho$ & $\\eta$ & $\\rho$ & "
	<< variable << " & Noise rejection \\\\" << endl
	<< "\\hline" << endl;

    for(int ie = 0; ie < nenergies; ie++) {

      //get the histograms
      f.GetObject(TString::Format("h_value_%s_e%d",      variable.Data(), ie), h_value[ie]);
      f.GetObject(TString::Format("h_efficiency_%s_e%d", variable.Data(), ie), h_efficiency[ie]);
      f.GetObject(TString::Format("h_purity_%s_e%d",     variable.Data(), ie), h_purity[ie]);
      f.GetObject(TString::Format("h_ep_%s_e%d",         variable.Data(), ie), h_ep[ie]);
      if(!h_value[ie] || !h_efficiency[ie] || !h_purity[ie] || !h_ep[ie]) {
	cerr << "Not all histograms found for " << variable
	     << " energy bin " << ie << endl;
	continue;
      }

      //populate the table
      int maxbin = h_ep[ie]->GetMaximumBin();
      str << "$" << ie << " \\le E_{e_{\\textrm{tot}}} < " << ie+1 << "$ & "
	      << std::fixed << std::setprecision(1)
	      << 100 * h_ep        [ie]->GetBinContent(maxbin) << " & "
	      << 100 * h_efficiency[ie]->GetBinContent(maxbin) << " & "
	      << 100 * h_purity    [ie]->GetBinContent(maxbin) << " & "
	      << std::setprecision(2)
	      << h_value           [ie]->GetBinCenter(maxbin)  << " & "
	      << std::setprecision(1)
	      << 100 * (1 - h_efficiency[0]->GetBinContent(maxbin)) << " \\\\"
	      << endl;
    }//ie

    str << "\\hline" << endl
	<< "\\end{tabular}" << endl;
    str.close();

    //loop
    iv++;
  }//loop over variables

}
