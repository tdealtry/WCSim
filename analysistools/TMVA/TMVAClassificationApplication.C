/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TVector3.h"

#include "TMVAClassification.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void AnalyseEventWithMVA(TMVA::Reader * reader, TString method, TH1F * h, double & value, double & prob, double & rarity)
{
  value  = reader->EvaluateMVA(method);
  if(method.Contains("BDT")) {
    prob   = reader->GetProba(method);
    rarity = reader->GetRarity(method);
  }
  h->Fill(value);
}

void TMVAClassificationApplication( bool test = true,
				    TString myMethodList = "BDT",
				    TString tag = "test",
				    const char * filename = "$WCSIMDIR/runs/20160616_tmva_testevents/tmva_analysiswcsim_makekin_*e-*.0.root"
				    )
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 1;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 1; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 1;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	      std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   //create output file & tree
   TString outfilename = TString::Format("TMVApp_%s.root", tag.Data());
   TFile *target  = new TFile( outfilename,"CREATE" );
   if(!target || target->IsZombie()) {
     cerr << "File " << outfilename << " already exists. Exiting..." << endl;
     return;
   }
   TTree *tout = new TTree("tmva_results", "Result of trained TMVA applied to independent data");
   
   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   vector<pair<string, char> > vnames = GetVariableNames(test, 0);
   vector<pair<string, char> > vnames_plain = GetVariableNames(test, 2);
   Double_t var_d[vnames.size()];
   Float_t  var_f[vnames.size()];
   Int_t    var_i[vnames.size()];
   Bool_t   var_b[vnames.size()];
   for(size_t iv = 0; iv < vnames.size(); iv++) {
     reader->AddVariable(vnames[iv].first, &var_f[iv]);
   }//iv
   /*
   Float_t var1, var2;
   Float_t var3, var4;
   reader->AddVariable( "myvar1 := var1+var2", &var1 );
   reader->AddVariable( "myvar2 := var1-var2", &var2 );
   reader->AddVariable( "var3",                &var3 );
   reader->AddVariable( "var4",                &var4 );
   */

   // Spectator variables declared in the training have to be added to the reader, too
   vector<pair<string, char> > snames = GetSpectatorNames();
   Int_t    spec_b[snames.size()];
   Int_t    spec_i[snames.size()];
   Float_t  spec_f[snames.size()];
   Double_t spec_d[snames.size()];
   //TVector3 spec_v[snames.size()];
   for(size_t iv = 0; iv < snames.size(); iv++)
     reader->AddSpectator(snames[iv].first.c_str(), &spec_f[iv]);
   /*
   Float_t spec1,spec2;
   reader->AddSpectator( "spec1 := var1*2",   &spec1 );
   reader->AddSpectator( "spec2 := var1*3",   &spec2 );
   */

   Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }

   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = TString::Format("TMVAClassification_%s", tag.Data());

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
     if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }

   // Book output histograms
   UInt_t nbin = 100;
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtB(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);
   //and tree variables (return value)
   bool   bCuts(0), bCutsD(0), bCutsPCA(0), bCutsGA(0), bCutsSA(0);
   double bLk(0), bLkD(0), bLkPCA(0), bLkKDE(0), bLkMIX(0), bPD(0), bPDD(0);
   double bPDPCA(0), bPDEFoam(0), bKNN(0), bHm(0);
   double bFi(0), bFiG(0), bFiB(0), bLD(0), bNn(0),bNnbfgs(0),bNnbnn(0);
   double bNnC(0), bNnT(0), bBdt(0), bBdtG(0), bBdtB(0), bBdtD(0), bRf(0), bSVMG(0);
   double bSVMP(0), bSVML(0), bFDAMT(0), bFDAGA(0), bCat(0), bPBdt(0);
   //and tree variables (probablity)
   double bLk_prob(0), bLkD_prob(0), bLkPCA_prob(0), bLkKDE_prob(0), bLkMIX_prob(0), bPD_prob(0), bPDD_prob(0);
   double bPDPCA_prob(0), bPDEFoam_prob(0), bKNN_prob(0), bHm_prob(0);
   double bFi_prob(0), bFiG_prob(0), bFiB_prob(0), bLD_prob(0), bNn_prob(0),bNnbfgs_prob(0),bNnbnn_prob(0);
   double bNnC_prob(0), bNnT_prob(0), bBdt_prob(0), bBdtG_prob(0), bBdtB_prob(0), bBdtD_prob(0), bRf_prob(0), bSVMG_prob(0);
   double bSVMP_prob(0), bSVML_prob(0), bFDAMT_prob(0), bFDAGA_prob(0), bCat_prob(0), bPBdt_prob(0);
   //and tree variables (rarity)
   double bLk_rarity(0), bLkD_rarity(0), bLkPCA_rarity(0), bLkKDE_rarity(0), bLkMIX_rarity(0), bPD_rarity(0), bPDD_rarity(0);
   double bPDPCA_rarity(0), bPDEFoam_rarity(0), bKNN_rarity(0), bHm_rarity(0);
   double bFi_rarity(0), bFiG_rarity(0), bFiB_rarity(0), bLD_rarity(0), bNn_rarity(0),bNnbfgs_rarity(0),bNnbnn_rarity(0);
   double bNnC_rarity(0), bNnT_rarity(0), bBdt_rarity(0), bBdtG_rarity(0), bBdtB_rarity(0), bBdtD_rarity(0), bRf_rarity(0), bSVMG_rarity(0);
   double bSVMP_rarity(0), bSVML_rarity(0), bFDAMT_rarity(0), bFDAGA_rarity(0), bCat_rarity(0), bPBdt_rarity(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["BDTB"])          histBdtB    = new TH1F( "MVA_BDTB",          "MVA_BDTB",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_PBDT",          nbin, -0.8, 0.8 );

   if (Use["Cuts"])          tout->Branch("Cuts",    &bCuts    );
   if (Use["CutsD"])         tout->Branch("CutsD",   &bCutsD   );
   if (Use["CutsPCA"])       tout->Branch("CutsPCA", &bCutsPCA );
   if (Use["CutsGA"])        tout->Branch("CutsGA",  &bCutsGA  );
   if (Use["CutsSA"])        tout->Branch("CutsSA",  &bCutsSA  );
   if (Use["Likelihood"])    tout->Branch("Likelihood",      &bLk      );
   if (Use["LikelihoodD"])   tout->Branch("LikelihoodD",     &bLkD     );
   if (Use["LikelihoodPCA"]) tout->Branch("LikelihoodPCA",   &bLkPCA   );
   if (Use["LikelihoodKDE"]) tout->Branch("LikelihoodKDE",   &bLkKDE   );
   if (Use["LikelihoodMIX"]) tout->Branch("LikelihoodMIX",   &bLkMIX   );
   if (Use["PDERS"])         tout->Branch("PDERS",      &bPD      );
   if (Use["PDERSD"])        tout->Branch("PDERSD",     &bPDD     );
   if (Use["PDERSPCA"])      tout->Branch("PDERSPCA",   &bPDPCA   );
   if (Use["KNN"])           tout->Branch("KNN",     &bKNN     );
   if (Use["HMatrix"])       tout->Branch("HMatrix",      &bHm      );
   if (Use["Fisher"])        tout->Branch("Fisher",      &bFi      );
   if (Use["FisherG"])       tout->Branch("FisherG",     &bFiG     );
   if (Use["BoostedFisher"]) tout->Branch("BoostedFisher",     &bFiB     );
   if (Use["LD"])            tout->Branch("LD",      &bLD      );
   if (Use["MLP"])           tout->Branch("MLP",      &bNn      );
   if (Use["MLPBFGS"])       tout->Branch("MLPBFGS",  &bNnbfgs  );
   if (Use["MLPBNN"])        tout->Branch("MLPBNN",   &bNnbnn   );
   if (Use["CFMlpANN"])      tout->Branch("CFMlpANN",     &bNnC     );
   if (Use["TMlpANN"])       tout->Branch("TMlpANN",     &bNnT     );
   if (Use["BDT"])           tout->Branch("BDT",     &bBdt     );
   if (Use["BDTD"])          tout->Branch("BDTD",    &bBdtD    );
   if (Use["BDTG"])          tout->Branch("BDTG",    &bBdtG    );
   if (Use["BDTB"])          tout->Branch("BDTB",    &bBdtB    );
   if (Use["RuleFit"])       tout->Branch("RuleFit",      &bRf      );
   if (Use["SVM_Gauss"])     tout->Branch("SVM_Gauss",    &bSVMG    );
   if (Use["SVM_Poly"])      tout->Branch("SVM_Poly",    &bSVMP    );
   if (Use["SVM_Lin"])       tout->Branch("SVM_Lin",    &bSVML    );
   if (Use["FDA_MT"])        tout->Branch("FDA_MT",   &bFDAMT   );
   if (Use["FDA_GA"])        tout->Branch("FDA_GA",   &bFDAGA   );
   if (Use["Category"])      tout->Branch("Category",     &bCat     );
   if (Use["Plugin"])        tout->Branch("Plugin",    &bPBdt    );
   //prob
   if (Use["Likelihood"])    tout->Branch("Likelihood_prob",      &bLk_prob      );
   if (Use["LikelihoodD"])   tout->Branch("LikelihoodD_prob",     &bLkD_prob     );
   if (Use["LikelihoodPCA"]) tout->Branch("LikelihoodPCA_prob",   &bLkPCA_prob   );
   if (Use["LikelihoodKDE"]) tout->Branch("LikelihoodKDE_prob",   &bLkKDE_prob   );
   if (Use["LikelihoodMIX"]) tout->Branch("LikelihoodMIX_prob",   &bLkMIX_prob   );
   if (Use["PDERS"])         tout->Branch("PDERS_prob",      &bPD_prob      );
   if (Use["PDERSD"])        tout->Branch("PDERSD_prob",     &bPDD_prob     );
   if (Use["PDERSPCA"])      tout->Branch("PDERSPCA_prob",   &bPDPCA_prob   );
   if (Use["KNN"])           tout->Branch("KNN_prob",     &bKNN_prob     );
   if (Use["HMatrix"])       tout->Branch("HMatrix_prob",      &bHm_prob      );
   if (Use["Fisher"])        tout->Branch("Fisher_prob",      &bFi_prob      );
   if (Use["FisherG"])       tout->Branch("FisherG_prob",     &bFiG_prob     );
   if (Use["BoostedFisher"]) tout->Branch("BoostedFisher_prob",     &bFiB_prob     );
   if (Use["LD"])            tout->Branch("LD_prob",      &bLD_prob      );
   if (Use["MLP"])           tout->Branch("MLP_prob",      &bNn_prob      );
   if (Use["MLPBFGS"])       tout->Branch("MLPBFGS_prob",  &bNnbfgs_prob  );
   if (Use["MLPBNN"])        tout->Branch("MLPBNN_prob",   &bNnbnn_prob   );
   if (Use["CFMlpANN"])      tout->Branch("CFMlpANN_prob",     &bNnC_prob     );
   if (Use["TMlpANN"])       tout->Branch("TMlpANN_prob",     &bNnT_prob     );
   if (Use["BDT"])           tout->Branch("BDT_prob",     &bBdt_prob     );
   if (Use["BDTD"])          tout->Branch("BDTD_prob",    &bBdtD_prob    );
   if (Use["BDTG"])          tout->Branch("BDTG_prob",    &bBdtG_prob    );
   if (Use["BDTB"])          tout->Branch("BDTB_prob",    &bBdtB_prob    );
   if (Use["RuleFit"])       tout->Branch("RuleFit_prob",      &bRf_prob      );
   if (Use["SVM_Gauss"])     tout->Branch("SVM_Gauss_prob",    &bSVMG_prob    );
   if (Use["SVM_Poly"])      tout->Branch("SVM_Poly_prob",    &bSVMP_prob    );
   if (Use["SVM_Lin"])       tout->Branch("SVM_Lin_prob",    &bSVML_prob    );
   if (Use["FDA_MT"])        tout->Branch("FDA_MT_prob",   &bFDAMT_prob   );
   if (Use["FDA_GA"])        tout->Branch("FDA_GA_prob",   &bFDAGA_prob   );
   if (Use["Category"])      tout->Branch("Category_prob",     &bCat_prob     );
   if (Use["Plugin"])        tout->Branch("Plugin_prob",    &bPBdt_prob    );
   //rarity
   if (Use["Likelihood"])    tout->Branch("Likelihood_rarity",      &bLk_rarity      );
   if (Use["LikelihoodD"])   tout->Branch("LikelihoodD_rarity",     &bLkD_rarity     );
   if (Use["LikelihoodPCA"]) tout->Branch("LikelihoodPCA_rarity",   &bLkPCA_rarity   );
   if (Use["LikelihoodKDE"]) tout->Branch("LikelihoodKDE_rarity",   &bLkKDE_rarity   );
   if (Use["LikelihoodMIX"]) tout->Branch("LikelihoodMIX_rarity",   &bLkMIX_rarity   );
   if (Use["PDERS"])         tout->Branch("PDERS_rarity",      &bPD_rarity      );
   if (Use["PDERSD"])        tout->Branch("PDERSD_rarity",     &bPDD_rarity     );
   if (Use["PDERSPCA"])      tout->Branch("PDERSPCA_rarity",   &bPDPCA_rarity   );
   if (Use["KNN"])           tout->Branch("KNN_rarity",     &bKNN_rarity     );
   if (Use["HMatrix"])       tout->Branch("HMatrix_rarity",      &bHm_rarity      );
   if (Use["Fisher"])        tout->Branch("Fisher_rarity",      &bFi_rarity      );
   if (Use["FisherG"])       tout->Branch("FisherG_rarity",     &bFiG_rarity     );
   if (Use["BoostedFisher"]) tout->Branch("BoostedFisher_rarity",     &bFiB_rarity     );
   if (Use["LD"])            tout->Branch("LD_rarity",      &bLD_rarity      );
   if (Use["MLP"])           tout->Branch("MLP_rarity",      &bNn_rarity      );
   if (Use["MLPBFGS"])       tout->Branch("MLPBFGS_rarity",  &bNnbfgs_rarity  );
   if (Use["MLPBNN"])        tout->Branch("MLPBNN_rarity",   &bNnbnn_rarity   );
   if (Use["CFMlpANN"])      tout->Branch("CFMlpANN_rarity",     &bNnC_rarity     );
   if (Use["TMlpANN"])       tout->Branch("TMlpANN_rarity",     &bNnT_rarity     );
   if (Use["BDT"])           tout->Branch("BDT_rarity",     &bBdt_rarity     );
   if (Use["BDTD"])          tout->Branch("BDTD_rarity",    &bBdtD_rarity    );
   if (Use["BDTG"])          tout->Branch("BDTG_rarity",    &bBdtG_rarity    );
   if (Use["BDTB"])          tout->Branch("BDTB_rarity",    &bBdtB_rarity    );
   if (Use["RuleFit"])       tout->Branch("RuleFit_rarity",      &bRf_rarity      );
   if (Use["SVM_Gauss"])     tout->Branch("SVM_Gauss_rarity",    &bSVMG_rarity    );
   if (Use["SVM_Poly"])      tout->Branch("SVM_Poly_rarity",    &bSVMP_rarity    );
   if (Use["SVM_Lin"])       tout->Branch("SVM_Lin_rarity",    &bSVML_rarity    );
   if (Use["FDA_MT"])        tout->Branch("FDA_MT_rarity",   &bFDAMT_rarity   );
   if (Use["FDA_GA"])        tout->Branch("FDA_GA_rarity",   &bFDAGA_rarity   );
   if (Use["Category"])      tout->Branch("Category_rarity",     &bCat_rarity     );
   if (Use["Plugin"])        tout->Branch("Plugin_rarity",    &bPBdt_rarity    );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   

   TChain * theTree = new TChain("t_tmva");
   theTree->Add(filename);

   /*
   TFile *input(0);
   TString fname = "./tmva_example.root";   
   if (!gSystem->AccessPathName( fname )) 
      input = TFile::Open( fname ); // check if file in local directory exists
   else    
      input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
   
   // --- Event loop

   // Prepare the event tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("TreeS");
   */

   cout << "setting chain branch addresses" << endl;
   for(size_t iv = 0; iv < vnames_plain.size(); iv++) {
     var_b[iv] = false;
     if(vnames_plain[iv].second == 'F')
       theTree->SetBranchAddress(vnames_plain[iv].first.c_str(), &var_d[iv]);
     else if(vnames_plain[iv].second == 'I') {
       theTree->SetBranchAddress(vnames_plain[iv].first.c_str(), &var_i[iv]);
       var_b[iv] = true;
     }
   }

   //book some more branches
   Bool_t isSignal;
   tout->Branch("isSignal", &isSignal);

   for(size_t iv = 0; iv < snames.size(); iv++) {
     spec_b[iv] = 0;
     if(snames[iv].second == 'F') {
       theTree->SetBranchAddress(snames[iv].first.c_str(), &spec_d[iv]);
       tout->Branch(snames[iv].first.c_str(), &spec_d[iv]);
     }
     else if(snames[iv].second == 'I') {
       theTree->SetBranchAddress(snames[iv].first.c_str(), &spec_i[iv]);
       tout->Branch(snames[iv].first.c_str(), &spec_i[iv]);
       spec_b[iv] = 1;
     }
     /*
     else if(snames[iv].second == 'V') {
       theTree->SetBranchAddress(snames[iv].first.c_str(), &spec_v[iv]);
       tout->Branch(snames[iv].first.c_str(), &spec_v[iv]);
       spec_b[iv] = 2;
     }
     */
   }//iv
   cout << "chain branch addresses set" << endl;

   /*
   Float_t userVar1, userVar2;
   theTree->SetBranchAddress( "var1", &userVar1 );
   theTree->SetBranchAddress( "var2", &userVar2 );
   theTree->SetBranchAddress( "var3", &var3 );
   theTree->SetBranchAddress( "var4", &var4 );
   */

   // Efficiency calculator for cut method
   Int_t    nSelCuts(0), nSelCutsD(0), nSelCutsPCA(0), nSelCutsGA(0), nSelCutsSA(0);
   Double_t effS       = 0.99; //want a high efficiency for signal

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      //cout << "Getting event " << ievt << endl;
      theTree->GetEntry(ievt);
      //cout << "Got event" << endl;
      for(size_t iv = 0; iv < vnames.size(); iv++) {
	if(var_b[iv])
	  var_f[iv] = var_i[iv];
	else
	  var_f[iv] = var_d[iv];
	//cout << vnames[iv].first << ":" << var_f[iv] << " ";
      }//iv
      //cout << endl;
      for(size_t iv = 0; iv < snames.size(); iv++) {
	if(spec_b[iv] == 1)
	  spec_f[iv] = spec_i[iv];
	/*else if(spec_b[iv] == 2)
	  spec_f[iv] = spec_v[iv].Mag();*/
	else
	  spec_f[iv] = spec_d[iv];
      }//iv

      isSignal = (spec_f[2] > 0.1); //TRUE_e_energy

      /*
      var1 = userVar1 + userVar2;
      var2 = userVar1 - userVar2;
      */

      // --- Return the MVA outputs and fill into histograms

      if (Use["Cuts"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "Cuts method", effS );
         if (passed) nSelCuts++;
      }
      if (Use["CutsD"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsD method", effS );
         if (passed) nSelCutsD++;
      }
      if (Use["CutsPCA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsPCA method", effS );
         if (passed) nSelCutsPCA++;
      }
      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }
      if (Use["CutsSA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsSA method", effS );
         if (passed) nSelCutsSA++;
      }

      if (Use["Likelihood"   ]) AnalyseEventWithMVA(reader, "Likelihood method",    histLk,     bLk,     bLk_prob,     bLk_rarity     );
      if (Use["LikelihoodD"  ]) AnalyseEventWithMVA(reader, "LikelihoodD method",   histLkD,    bLkD,    bLkD_prob,    bLkD_rarity    );
      if (Use["LikelihoodPCA"]) AnalyseEventWithMVA(reader, "LikelihoodPCA method", histLkPCA,  bLkPCA,  bLkPCA_prob,  bLkPCA_rarity  );
      if (Use["LikelihoodKDE"]) AnalyseEventWithMVA(reader, "LikelihoodKDE method", histLkKDE,  bLkKDE,  bLkKDE_prob,  bLkKDE_rarity  );
      if (Use["LikelihoodMIX"]) AnalyseEventWithMVA(reader, "LikelihoodMIX method", histLkMIX,  bLkMIX,  bLkMIX_prob,  bLkMIX_rarity  );
      if (Use["PDERS"        ]) AnalyseEventWithMVA(reader, "PDERS method",         histPD,     bPD,     bPD_prob,     bPD_rarity     );
      if (Use["PDERSD"       ]) AnalyseEventWithMVA(reader, "PDERSD method",        histPDD,    bPDD,    bPDD_prob,    bPDD_rarity    );
      if (Use["PDERSPCA"     ]) AnalyseEventWithMVA(reader, "PDERSPCA method",      histPDPCA,  bPDPCA,  bPDPCA_prob,  bPDPCA_rarity  );
      if (Use["KNN"          ]) AnalyseEventWithMVA(reader, "KNN method",           histKNN,    bKNN,    bKNN_prob,    bKNN_rarity    );
      if (Use["HMatrix"      ]) AnalyseEventWithMVA(reader, "HMatrix method",       histHm,     bHm,     bHm_prob,     bHm_rarity     );
      if (Use["Fisher"       ]) AnalyseEventWithMVA(reader, "Fisher method",        histFi,     bFi,     bFi_prob,     bFi_rarity     );
      if (Use["FisherG"      ]) AnalyseEventWithMVA(reader, "FisherG method",       histFiG,    bFiG,    bFiG_prob,    bFiG_rarity    );
      if (Use["BoostedFisher"]) AnalyseEventWithMVA(reader, "BoostedFisher method", histFiB,    bFiB,    bFiB_prob,    bFiB_rarity    );
      if (Use["LD"           ]) AnalyseEventWithMVA(reader, "LD method",            histLD,     bLD,     bLD_prob,     bLD_rarity     );
      if (Use["MLP"          ]) AnalyseEventWithMVA(reader, "MLP method",           histNn,     bNn,     bNn_prob,     bNn_rarity     );
      if (Use["MLPBFGS"      ]) AnalyseEventWithMVA(reader, "MLPBFGS method",       histNnbfgs, bNnbfgs, bNnbfgs_prob, bNnbfgs_rarity );
      if (Use["MLPBNN"       ]) AnalyseEventWithMVA(reader, "MLPBNN method",        histNnbnn,  bNnbnn,  bNnbnn_prob,  bNnbnn_rarity  );
      if (Use["CFMlpANN"     ]) AnalyseEventWithMVA(reader, "CFMlpANN method",      histNnC,    bNnC,    bNnC_prob,    bNnC_rarity    );
      if (Use["TMlpANN"      ]) AnalyseEventWithMVA(reader, "TMlpANN method",       histNnT,    bNnT,    bNnT_prob,    bNnT_rarity    );
      if (Use["BDT"          ]) AnalyseEventWithMVA(reader, "BDT method",           histBdt,    bBdt,    bBdt_prob,    bBdt_rarity    );
      if (Use["BDTB"         ]) AnalyseEventWithMVA(reader, "BDTB method",          histBdtB,   bBdtB,   bBdtB_prob,   bBdtB_rarity   );
      if (Use["BDTD"         ]) AnalyseEventWithMVA(reader, "BDTD method",          histBdtD,   bBdtD,   bBdtD_prob,   bBdtD_rarity   );
      if (Use["BDTG"         ]) AnalyseEventWithMVA(reader, "BDTG method",          histBdtG,   bBdtG,   bBdtG_prob,   bBdtG_rarity   );
      if (Use["RuleFit"      ]) AnalyseEventWithMVA(reader, "RuleFit method",       histRf,     bRf,     bRf_prob,     bRf_rarity     );
      if (Use["SVM_Gauss"    ]) AnalyseEventWithMVA(reader, "SVM_Gauss method",     histSVMG,   bSVMG,   bSVMG_prob,   bSVMG_rarity   );
      if (Use["SVM_Poly"     ]) AnalyseEventWithMVA(reader, "SVM_Poly method",      histSVMP,   bSVMP,   bSVMP_prob,   bSVMP_rarity   );
      if (Use["SVM_Lin"      ]) AnalyseEventWithMVA(reader, "SVM_Lin method",       histSVML,   bSVML,   bSVML_prob,   bSVML_rarity   );
      if (Use["FDA_MT"       ]) AnalyseEventWithMVA(reader, "FDA_MT method",        histFDAMT,  bFDAMT,  bFDAMT_prob,  bFDAMT_rarity  );
      if (Use["FDA_GA"       ]) AnalyseEventWithMVA(reader, "FDA_GA method",        histFDAGA,  bFDAGA,  bFDAGA_prob,  bFDAGA_rarity  );
      if (Use["Category"     ]) AnalyseEventWithMVA(reader, "Category method",      histCat,    bCat,    bCat_prob,    bCat_rarity    );
      if (Use["Plugin"       ]) AnalyseEventWithMVA(reader, "P_BDT method",         histPBdt,   bPBdt,   bPBdt_prob,   bPBdt_rarity   );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
	probHistFi  ->Fill(bFi_prob  );
	rarityHistFi->Fill(bFi_rarity);
      }
      tout->Fill();
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   // --- Write histograms

   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
   if (Use["BDT"          ])   histBdt    ->Write();
   if (Use["BDTB"         ])   histBdtB   ->Write();
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write(); 
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();

   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
   tout->Write();
   target->Close();

   std::cout << "--- Created root file: \"" << outfilename << "\" containing the MVA output histograms" << std::endl;
  
   delete reader;
    
   std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 
