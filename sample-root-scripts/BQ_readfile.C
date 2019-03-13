
#include <iostream>
#include <vector>
#include <TH1F.h>
#include <TFile.h>
#include <TH1D.h>
#include <TBranch.h>
#include <TMath.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <TH2D.h>
#include <TTree.h>
#include <TSystem.h>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

using namespace std;

// Simple example of reading a generated Root file
const int nPMTtypes = 2;

//2D plots
TH2D * ChargeProfile2D[nPMTtypes];
TH2D * ChargeProfile2D_onlyFront[nPMTtypes];
TH2D * ChargeProfile2DTop[nPMTtypes];
TH2D * ChargeProfile2DCylindricCoordinate[nPMTtypes];
TH2D * HitProfile2DCylindricCoordinate[nPMTtypes];
TH2D * ChargeProfile2DRelativeAngle[nPMTtypes];
TH2D * TimeTOFProfileXTOF[nPMTtypes];

//1D plots of charge & hit
TH1D * ChargeProfile[nPMTtypes];
TH1D * HitProfile[nPMTtypes];
TH1D * HitProfilePMTId[nPMTtypes];

//1D time profiles
TH1D * TimeProfile[nPMTtypes];
TH1D * TimeHitProfile[nPMTtypes];
TH1D * TimeTOFProfile[nPMTtypes];
TH1D * HitTimeTOFProfile[nPMTtypes];

//Charge per PMT
TH1D * ChargePerPMT[nPMTtypes];

//Total charge
TH1D * TotalCharge[nPMTtypes];
TH1D * TotalHit[nPMTtypes];

TH2D * TotalChargeXdWall[nPMTtypes];

/*
void plotBQ(double PMTpos[3], double PMTdir[3],double phiAngle,double relativeAngle,double peForTube,int pmtType){
  if(TMath::Abs(PMTpos[2]) < 2500){//Temp
    ChargeProfile2D[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
    ChargeProfile2DTop[pmtType]->Fill(PMTpos[0],PMTpos[1],peForTube);
    ChargeProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],peForTube);
    HitProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],1);
    ChargeProfile2DRelativeAngle[pmtType]->Fill(relativeAngle,PMTpos[2],peForTube);
    
    ChargeProfile[pmtType]->Fill(relativeAngle,peForTube);
    HitProfilePMTId[pmtType]->Fill(tubeNumber,1);
    HitProfile[pmtType]->Fill(relativeAngle,1);
    
    TimeProfile[pmtType]->Fill(time,peForTube);
    
    ChargePerPMT[pmtType]->Fill(peForTube);
    
    //totalPe += peForTube;
    //totalHit ++;
  }
}
*/

void BQ_readfile(char *filename=NULL, bool verbose=true)
{
  bool hybrid = true;//false;//true;
  bool Gamma = false;//Is the mother particle a gamma or another particle?
  float cvacuum = 3e8 / 1e9;//speed of light, in meter per ns.
  float nindex = 1.373;//refraction index of water
  bool plotDigitized = true;//false;//true;//false;//true;//false;//true;//false;
  
  // Clear global scope
  //gROOT->Reset();
  /*
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetStatColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleSize(0.04);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPalette(1);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleX(.5);
  gStyle->SetTitleY(0.99);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetHatchesLineWidth(2);
  gStyle->SetLineWidth(1.5);
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetCanvasBorderMode(0);
  */
  // Load the library with class dictionary info
  // (create with "gmake shared")
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  }else{
    gSystem->Load("../libWCSimRoot.so");
  }

  bool HK=true;
  double TotalNHits = 1e4;
  double TankSize = 1100/2;//in cm, the maximal size of the plots for vertices etc...
  double TankRadius = 742/2;//in cm, the maximal size of the plots for vertices etc...
  double TankHalfHeight = 1042/2;//in cm, the maximal size of the plots for vertices etc...

  if(HK){
    TotalNHits = 3e4;
    TankSize = 7100/2;//in cm, the maximal size of the plots for vertices etc...
    TankRadius = 7080/2;//in cm, the maximal size of the plots for vertices etc...
    TankHalfHeight = 5480/2;//in cm, the maximal size of the plots for vertices etc...
    //nevt = std::min( 1e3 , ((double) WSTreeHits->GetEntries()) );
  }
    //nevt = std::min( 1e3 , ((double) WSTreeHits->GetEntries()) );

  
  TFile *file;
  // Open the file
  if (filename==NULL){
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hk_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hkmpmt_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkmpmt_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_trash.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_2muniformsphere.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_2muniformsphere.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_1000events.root","read");   
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_10hitstrigger.root","read");
    //wcsim_hkhybridmpmt_e10_center_nominal_fulltank_0hitstrigger_1000events.root
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_0hitstrigger_1000events.root","read");
    file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e10_center_nominal_fulltank_0hitstrigger_nodn.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374_e10_center_nominal_fulltank_0hitstrigger.root","read");

    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374_e10_center_nominal_fulltank_0hitstrigger_noqe.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt14374_e10_trash.root","read");
    
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_2muniformsphere_100events.root","read");
    //file = new TFile("/disk01/usr5/bquilain/trash.root","read");

    
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_onlybal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal_2.root","read");
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hkmpmt_mu500_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_100_0_0_0/11_100_0_0_0.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_900_0_6_0/11_900_0_6_0.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkmpmt100HzDN_e4_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkmpmt_e4_center_nominal_fulltank.root","read");
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hk_e3_center_nominal_fulltank.root","read");
    //file = new TFile("../test.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_900_0_3_0/11_900_0_3_0.root","read");
    //file = new TFile("/disk01/usr5/bquilain/fiTQun/HK_704cmx548cmID_mPMT_40perCent/timepdf/11_900_0_3_0/11_900_0_3_0.root","read");


    
  }else{
    file = new TFile(filename,"read");
  }
  if (!file->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return;
  }
  
  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");
  
  // Get the number of events
  int nevent = tree->GetEntries();
  if(verbose) printf("nevent %d\n",nevent);
  
  // Create a WCSimRootEvent to put stuff from the tree in

  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  TBranch *branch2;
  if(hybrid){
    branch2 = tree->GetBranch("wcsimrootevent2");
    branch2->SetAddress(&wcsimrootsuperevent2);
  // Force deletion to prevent memory leak 
    tree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
  }

  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = 0; 
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
      exit(9);
  }
  geotree->GetEntry(0);
  cout << "Number of PMTs of 1st type = " << geo->GetWCNumPMT() << endl;
  
  // Options tree - only need 1 "event"
  TTree *opttree = (TTree*)file->Get("wcsimRootOptionsT");
  WCSimRootOptions *opt = 0; 
  opttree->SetBranchAddress("wcsimrootoptions", &opt);
  if(verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
  if (opttree->GetEntries() == 0) {
    exit(9);
  }
  opttree->GetEntry(0);
  opt->Print();

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;
  WCSimRootTrigger* wcsimrootevent2;

  TH1F *h1 = new TH1F("PMT Hits", "PMT Hits", 8000, 0, 8000);
  TH1F *hvtx0 = new TH1F("Event VTX0", "Event VTX0", 200, -1500, 1500);
  TH1F *hvtx1 = new TH1F("Event VTX1", "Event VTX1", 200, -1500, 1500);
  TH1F *hvtx2 = new TH1F("Event VTX2", "Event VTX2", 200, -1500, 1500);
  
  int num_trig=0;

    //TFile * fOutput = new TFile(OutputFile,"recreate");

  for(int i=0;i<nPMTtypes;i++){
    ChargeProfile2D_onlyFront[i] = new TH2D(Form("ChargeProfile2D_onlyFront_pmtType%d",i),"",TankSize/2,-TankSize/2.,TankSize/2.,TankSize/2,-TankSize/2.,TankSize/2.);
    ChargeProfile2D[i] = new TH2D(Form("ChargeProfile2D_pmtType%d",i),"",TankSize/2,-TankSize/2.,TankSize/2.,TankSize/2,-TankSize/2.,TankSize/2.);
    ChargeProfile2DTop[i] = new TH2D(Form("ChargeProfile2DTop_pmtType%d",i),"",TankSize/2,-TankSize/2.,TankSize/2.,TankSize/2,-TankSize/2.,TankSize/2.);
    //ChargeProfile2DCylindricCoordinate[i] = new TH2D(Form("ChargeProfile2DCylindricCoordinate_pmtType%d",i),"",3600,0,360,TankSize,-TankSize,TankSize);
    ChargeProfile2DCylindricCoordinate[i] = new TH2D(Form("ChargeProfile2DCylindricCoordinate_pmtType%d",i),";Phi (deg);z (cm);Charge (p.e.) per event",720,0,360,TankSize/2,-TankSize/2.,TankSize/2.);
    HitProfile2DCylindricCoordinate[i] = new TH2D(Form("HitProfile2DCylindricCoordinate_pmtType%d",i),";Phi (deg);z (cm);NHits per event",720,0,360,TankSize/2,-TankSize/2.,TankSize/2.);
    ChargeProfile2DRelativeAngle[i] = new TH2D(Form("ChargeProfile2DRelativeAngle_pmtType%d",i),"",720,0,360,TankSize/2,-TankSize/2.,TankSize/2.);
    TimeTOFProfileXTOF[i] = new TH2D(Form("TimeTOFProfileXTOF_pmtType%d",i),"",1e4,-1e2,1e3,1e2,0,5e2);TimeTOFProfileXTOF[i]->Sumw2();
    
    ChargeProfile[i] = new TH1D(Form("ChargeProfile_pmtType%d",i),"",720,0,180);
    HitProfile[i] = new TH1D(Form("HitProfile_pmtType%d",i),"",720,0,180);
    HitProfilePMTId[i] = new TH1D(Form("HitProfilePMTId_pmtType%d",i),";PMT ID;NHits per event",10894,0,206986);

    TimeProfile[i] = new TH1D(Form("TimeProfile_pmtType%d",i),"",1e4,0,5e3);
    TimeHitProfile[i] = new TH1D(Form("TimeHitProfile_pmtType%d",i),"",1e4,0,5e3);

    TimeTOFProfile[i] = new TH1D(Form("TimeTOFProfile_pmtType%d",i),"",1e4,-1e2,1e3);TimeTOFProfile[i]->Sumw2();
    HitTimeTOFProfile[i] = new TH1D(Form("HitTimeTOFProfile_pmtType%d",i),"",1e4,-1e2,1e3);HitTimeTOFProfile[i]->Sumw2();

    ChargePerPMT[i] = new TH1D(Form("ChargePerPMT_pmtType%d",i),"",500,0,500);

    TotalCharge[i] = new TH1D(Form("TotalCharge_pmtType%d",i),"",1e4,0,3e4);
    TotalHit[i] = new TH1D(Form("TotalHit_pmtType%d",i),"",1e4,0,TotalNHits);

    TotalChargeXdWall[i] = new TH2D(Form("TotalChargeXdWall_pmtType%d",i),"",TankSize/2,0,TankSize/2.,1e4,0,3e4);

    ChargeProfile2D[i]->Sumw2();
    ChargeProfile2D_onlyFront[i]->Sumw2();
    ChargeProfile2DTop[i]->Sumw2();
    ChargeProfile2DCylindricCoordinate[i]->Sumw2();
    HitProfile2DCylindricCoordinate[i]->Sumw2();
    ChargeProfile2DRelativeAngle[i]->Sumw2();
    ChargeProfile[i]->Sumw2(); HitProfile[i]->Sumw2();TimeProfile[i]->Sumw2();TimeHitProfile[i]->Sumw2();ChargePerPMT[i]->Sumw2();
    TotalCharge[i]->Sumw2();TotalHit[i]->Sumw2();
    TotalChargeXdWall[i]->Sumw2();
  }
  /*
  
  TH2D * HitProfile2D = new TH2D("HitProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * TimeProfile2D = new TH2D("TimeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * ChargeXPositionXTime = new TH2D("ChargeXPositionXTime","",36,0,180,TankSize,0,TankSize);
  TH2D * ChargeProfileXdWall = new TH2D("ChargeProfileXdWall","",720,0,180,100,0,TankSize);
  TH2D * TotalChargeXdWall = new TH2D("TotalChargeXdWall","",100,0,TankSize,1e3,0,1e5);

  TH2D * TimeAngleProfile = new TH2D("TimeAngleHitProfile","",720,0,180,1e4,0,5e3);TimeAngleProfile->Sumw2();
  TH2D * HitTimeAngleProfile = new TH2D("HitTimeAngleHitProfile","",720,0,180,1e4,0,5e3);HitTimeAngleProfile->Sumw2();
  TH2D * TimeTOFAngleProfile = new TH2D("TimeTOFAngleProfile","",720,0,180,1e4,-1e2,5e3);TimeTOFAngleProfile->Sumw2();
  TH2D * HitTimeTOFAngleProfile = new TH2D("HitTimeTOFAngleProfile","",720,0,180,1e4,-1e2,5e3);HitTimeTOFAngleProfile->Sumw2();
  TH2D * TimeTOFTOF = new TH2D("TimeTOFTOF","",1e4,-1e2,5e3,1e2,0,5e2);TimeTOFTOF->Sumw2();
  TH2D * HitTimeTOFTOF = new TH2D("HitTimeTOFTOF","",1e4,-1e2,5e3,1e2,0,5e2);HitTimeTOFTOF->Sumw2();


  TH3D * VertexPosition = new TH3D("VertexPosition","",100,-TankSize,TankSize,100,-TankSize,TankSize,100,-TankSize,TankSize);
  TH3D * VertexDirection = new TH3D("VertexDirection","",100,-1.1,1.1,100,-1.1,1.1,100,-1.1,1.1);
  TH1D * VertexXdWall = new TH1D("VertexXdWall","",100,0,TankSize);
  TH1D * ParentFlyingDistance = new TH1D("ParentFlyingDistance","",10000,0,TankSize);

  TH2D * EGammaSeparationXdWall = new TH2D("EGammaSeparationxdWall","",100,0,TankSize,1000,0.,5.);
  TH1D * EGammaSeparation = new TH1D("EGammaSeparation","",1000,0.,5.);

  //Set the 2D histogram to 0 to increase lisibility of plots
  for(int ibinx=1;ibinx<=ChargeProfile2D->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=ChargeProfile2D->GetNbinsY();ibiny++){
      HitProfile2D->SetBinContent(ibinx,ibiny,0);
      ChargeProfile2D->SetBinContent(ibinx,ibiny,0);     
      ChargeProfile2DXZ->SetBinContent(ibinx,ibiny,0);     
    }
  }
  //
  double ConversionLength = 36;//cm
  TH2D * ChargeProfile2D_PerEvent[nevt];
  */
  
  //TH2D * ChargeProfile2D = new TH2D("ChargeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);

  // Now loop over events
  for (int ev=0; ev<nevent; ev++)
  {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);      
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    //wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    if(verbose){
      printf("********************************************************");
      printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	     wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n",
	     wcsimrootsuperevent->GetNumberOfSubEvents());
      
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    }
    hvtx0->Fill(wcsimrootevent->GetVtx(0));
    hvtx1->Fill(wcsimrootevent->GetVtx(1));
    hvtx2->Fill(wcsimrootevent->GetVtx(2));

    if(verbose){
      printf("Jmu %d\n", wcsimrootevent->GetJmu());
      printf("Npar %d\n", wcsimrootevent->GetNpar());
      printf("Ntrack %d\n", wcsimrootevent->GetNtrack());
      
    }

    std::vector<float> triggerInfo;
    triggerInfo.clear();
    triggerInfo = wcsimrootevent->GetTriggerInfo();

    std::vector<float> triggerInfo2;
    triggerInfo2.clear();
    if(hybrid) triggerInfo2 = wcsimrootevent2->GetTriggerInfo();

    if(verbose){
      for(size_t v=0;v<triggerInfo.size();v++){
	cout << "Trigger entry #" << v << ", info = " << triggerInfo[v] << endl;
      }
      if(hybrid){
	for(size_t v=0;v<triggerInfo2.size();v++){
	  cout << "Trigger2 entry #" << v << ", info = " << triggerInfo2[v] << endl;
	}
      }
    }

    double triggerShift[nPMTtypes];
    double triggerTime[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      triggerShift[pmtType]=0;
      triggerTime[pmtType]=0;
      if(triggerInfo.size()>=3){
	if(pmtType==0){
	  triggerShift[pmtType] = triggerInfo[1];
	  triggerTime[pmtType] = triggerInfo[2];
	}
      }
      if(triggerInfo2.size()>=3){
	if(pmtType==1 && hybrid){
	  triggerShift[pmtType] = triggerInfo2[1];
	  triggerTime[pmtType] = triggerInfo2[2];
	}
      }
    }    
    // Now read the tracks in the event
    
    // Get the number of tracks
    int ntrack = wcsimrootevent->GetNtrack();
    if(verbose) printf("ntracks=%d\n",ntrack);

    double particleStart[3];
    double particleStop[3];
    double particleDir[3];

    int i;
    // Loop through elements in the TClonesArray of WCSimTracks
    for (i=0; i<ntrack; i++)
    {
      TObject *element = (wcsimrootevent->GetTracks())->At(i);
      
      WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);

      if(i==(ntrack-1)){//Mother particle
	for(int j=0; j<3; j++){
	  particleStart[j] = wcsimroottrack->GetStart(j);
	  particleStop[j] = wcsimroottrack->GetStop(j);
	  particleDir[j] = wcsimroottrack->GetDir(j);
	}//j
	//if(particleDir[0] != 1 || particleDir[1] != 0 || particleDir[2] != 0) cout << endl << endl << endl << "###############################################" << endl << "Particle direction has an issue" << endl << "###############################################" << endl;
	//if(particleStart[0] != 0 || particleStart[1] != 0 || particleStart[2] != 0) cout << endl << endl << endl << "###############################################" << endl << "Particle position has an issue" << endl << "###############################################" << endl;
      }//i is mother mother

      if(verbose){
	printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
	printf("Track energy: %f\n", wcsimroottrack->GetE());
	printf("Track momentum: %f\n", wcsimroottrack->GetP());
	printf("Track mass: %f\n", wcsimroottrack->GetM());
      
	for (int j=0; j<3; j++) {
	  printf("Track start: %d %f\n",j, wcsimroottrack->GetStart(j));
	  printf("Track dir: %d %f\n",j, wcsimroottrack->GetDir(j));
	}//j
      }

      
    }  // End of loop over tracks
  
    double dwall, dwallBarrel, dwallCap;
    double StartRadius,StartHeight;
    if(Gamma){
      StartRadius = TMath::Sqrt(particleStop[0]*particleStop[0]+particleStop[1]*particleStop[1]);
      StartHeight = TMath::Abs(particleStop[2]);
    }
    else{
      StartRadius = TMath::Sqrt(particleStart[0]*particleStart[0]+particleStart[1]*particleStart[1]);
      StartHeight = TMath::Abs(particleStart[2]);
    }
    dwallBarrel = TankRadius - StartRadius;
    dwallCap = TankHalfHeight - StartHeight;
    dwall = dwallBarrel;//min(dwallBarrel,dwallCap);
    if(ev%10000 == 0){
      cout<<"Evt #"<<ev<<endl;
      //cout<<"**************************"<<endl<<"Vertex:"<<Start_x[0]<<", "<<Start_y[0]<<", "<<Start_z[0]<<endl;
      cout<<"Dwall = "<<dwall<<", radius = "<<StartRadius<<", Height = "<<StartHeight<<endl;
      //VertexXdWall->Fill(dwall);
    }
    // Now look at the Cherenkov hits
    
    // Get the number of Cherenkov hits.
    // Note... this is *NOT* the number of photons that hit tubes.
    // It is the number of tubes hit with Cherenkov photons.
    // The number of digitized tubes will be smaller because of the threshold.
    // Each hit "raw" tube has several photon hits.  The times are recorded.
    // See chapter 5 of ../doc/DetectorDocumentation.pdf
    // for more information on the structure of the root file.
    //  
    // The following code prints out the hit times for the first 10 tubes and also
    // adds up the total pe.
    // 
    // For digitized info (one time/charge tube after a trigger) use
    // the digitized information.
    //

    int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits(); 
    int ncherenkovhits2 = 0; if(hybrid) ncherenkovhits2 = wcsimrootevent2->GetNcherenkovhits();
    int ncherenkovdigihits2 = 0;if(hybrid) ncherenkovdigihits2 = wcsimrootevent2->GetNcherenkovdigihits(); 
    
    h1->Fill(ncherenkovdigihits);
    if(verbose){
      printf("node id: %i\n", ev);
      printf("Ncherenkovhits %d\n",     ncherenkovhits);
      printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      printf("Ncherenkovhits2 %d\n",     ncherenkovhits2);
      printf("Ncherenkovdigihits2 %d\n", ncherenkovdigihits2);
      cout << "RAW HITS:" << endl;
    }

  
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(verbose) cout << "PMT Type = " << pmtType << endl;
      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray;
      if(pmtType==0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
      else timeArray = wcsimrootevent2->GetCherenkovHitTimes();
      
      double particleRelativePMTpos[3];
      double totalPe = 0;
      int totalHit = 0;

      int nhits;
      if(pmtType == 0) nhits = ncherenkovhits;
      else nhits = ncherenkovhits2;
      
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for (i=0; i< nhits ; i++)
	{
	  TObject *Hit;
	  if(pmtType==0) Hit = (wcsimrootevent->GetCherenkovHits())->At(i);
	  else Hit = (wcsimrootevent2->GetCherenkovHits())->At(i);

	  WCSimRootCherenkovHit *wcsimrootcherenkovhit = 
	    dynamic_cast<WCSimRootCherenkovHit*>(Hit);
	  
	  int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	  int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	  double peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);

	  WCSimRootPMT pmt;
	  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
	  else pmt  = geo->GetPMT(tubeNumber-1,true); 
	  
	  double PMTpos[3];
	  double PMTdir[3];
	  for(int j=0;j<3;j++){
	    PMTpos[j] = pmt.GetPosition(j);
	    PMTdir[j] = pmt.GetOrientation(j);
	  }
	  
	  if ( verbose && i < 10 ) // Only print first XX=10 tubes
	    {
	      if(verbose) printf("Total pe: %d times( ",peForTube);
	      for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
		{
		  WCSimRootCherenkovHitTime * HitTime = 
		    dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(j));
		  
		  if(verbose) printf("%6.2f ", HitTime->GetTruetime() );	     
		}
	      if(verbose) cout << ")" << endl;
	      if(verbose) std::cout << "Position of PMT = " << PMTpos[0] << ", " << PMTpos[1] << ", " << PMTpos[2] << endl;
	    }
	  
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  double radius = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
	  double phiAngle;
	  if(PMTpos[1] >= 0) phiAngle = TMath::ACos(PMTpos[0]/radius);
	  else phiAngle = TMath::Pi() + (TMath::Pi() - TMath::ACos(PMTpos[0]/radius));
	  
	  if(Gamma){//Produce the profile at the gamma conversion point, since it takes ~36 cm to convert/be visible
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStop[j];
	  }
	  else{
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStart[j];
	  }
	  
	  double vDir[3];double vOrientation[3];
	  double vDir2D[3];double vOrientation2D[3];
	  for(int j=0;j<3;j++){
	    vDir[j] = particleRelativePMTpos[j];
	    vOrientation[j] = PMTdir[j];
	    if(j<2){
	      vDir2D[j] = particleRelativePMTpos[j];
	      vOrientation2D[j] = PMTdir[j];
	    }	  
	  }
	  double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
	  double tof = Norm*1e-2/(cvacuum/nindex);
	  //if(verbose) cout << "Time of flight = " << tof << endl;
	  double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
	  for(int j=0;j<3;j++){
	    vDir[j] /= Norm;
	    vOrientation[j] /= NormOrientation;
	  }
	  double Norm2D = TMath::Sqrt(vDir2D[0]*vDir2D[0]+vDir2D[1]*vDir2D[1]);
	  double NormOrientation2D = TMath::Sqrt(vOrientation2D[0]*vOrientation2D[0]+vOrientation2D[1]*vOrientation2D[1]);
	  for(int j=0;j<2;j++){
	    vDir2D[j] /= Norm2D;
	    vOrientation2D[j] /= NormOrientation2D;
	  }
	  double HorizontalScalarProduct = TMath::ACos(vDir2D[0]*vOrientation2D[0] + vDir2D[1]*vOrientation2D[1])*180./TMath::Pi();
	  double relativeAngle = TMath::ACos(vDir[0]*particleDir[0]+vDir[1]*particleDir[1]+vDir[2]*particleDir[2])*180./TMath::Pi();
	  WCSimRootCherenkovHitTime * HitTime = dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(timeArrayIndex));//Assumes that the earliest time the PMT is hit is in the first element -> Everything is ordered.
	  double time = HitTime->GetTruetime();
	  
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  
	  //std::cout << "Relative angle = " << relativeAngle << ", time = " << time << endl;
	  
	  if(!plotDigitized){
	    //plotBQ(PMTpos[3], PMTdir[3],phiAngle,relativeAngle,peForTube,pmtType);
	    //if(TMath::Abs(PMTpos[2]) < 4500){//Temp
	      if(phiAngle < TMath::Pi()/2 || phiAngle > 3*TMath::Pi()/2) ChargeProfile2D_onlyFront[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2D[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2DTop[pmtType]->Fill(PMTpos[0],PMTpos[1],peForTube);
	      ChargeProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],peForTube);
	      HitProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],1);
	      ChargeProfile2DRelativeAngle[pmtType]->Fill(relativeAngle,PMTpos[2],peForTube);
	      
	      ChargeProfile[pmtType]->Fill(relativeAngle,peForTube);
	      HitProfile[pmtType]->Fill(relativeAngle,1);
	      HitProfilePMTId[pmtType]->Fill(tubeNumber, 1);
	      
	      TimeProfile[pmtType]->Fill(time,peForTube);
	      TimeTOFProfile[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],peForTube);
	      TimeTOFProfileXTOF[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],tof,peForTube);
	      
	      ChargePerPMT[pmtType]->Fill(peForTube);

	      totalPe += peForTube;
	      totalHit ++;
	      //}
	   
 //if(i%50 == 0) cout << "Hit angle wrt particle direction = " << relativeAngle << ", Phi = " << phiAngle*180/TMath::Pi() << " Z position = " << PMTpos[2] << endl;
	  }
	} // End of loop over Cherenkov hits
      if(verbose) cout << "Total Pe : " << totalPe << endl;
      if(!plotDigitized){
	TotalCharge[pmtType]->Fill(totalPe);
	TotalHit[pmtType]->Fill(totalHit);
	TotalChargeXdWall[pmtType]->Fill(dwall,totalPe);
      }
    }
    // Get the number of digitized hits
    // Loop over sub events


    
    if(verbose) cout << "DIGITIZED HITS:" << endl;
    //for (int index = 0 ; index < wcsimrootsuperevent->GetNumberOfEvents(); index++) 
    //{
    //wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
    //if(verbose) cout << "Sub event number = " << index << "\n";
    
    //int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
    //if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);

    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(verbose) cout << "PMT Type = " << pmtType << endl;
      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray;
      if(pmtType==0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
      else timeArray = wcsimrootevent2->GetCherenkovHitTimes();
      
      double particleRelativePMTpos[3];
      double totalPe = 0;
      int totalHit = 0;

      int nhits;
      if(pmtType == 0) nhits = ncherenkovdigihits;
      else nhits = ncherenkovdigihits2;
      
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for (i=0; i< nhits ; i++)
	{
	  TObject *Hit;
	  if(pmtType==0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
	  else Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);

	  WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
	    dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
	  
	  int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId();
	  double peForTube      = wcsimrootcherenkovdigihit->GetQ();

	  //int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	  //int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	  //int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	  WCSimRootPMT pmt;
	  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
	  else pmt  = geo->GetPMT(tubeNumber-1,true); 
	  
	  double PMTpos[3];
	  double PMTdir[3];
	  for(int j=0;j<3;j++){
	    PMTpos[j] = pmt.GetPosition(j);
	    PMTdir[j] = pmt.GetOrientation(j);
	  }
	  
	  if ( i < 10 ) // Only print first XX=10 tubes
	    {
	      if(verbose) printf("Total pe: %1.1f times( ",peForTube);
	      /*	      for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
		{
		  WCSimRootCherenkovHitTime HitTime = 
		    dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(j));
		  
		  if(verbose) printf("%6.2f ", HitTime.GetTruetime() );	     
		  }*/
	      if(verbose) cout << ")" << endl;
	      if(verbose) std::cout << "Position of PMT = " << PMTpos[0] << ", " << PMTpos[1] << ", " << PMTpos[2] << endl;
	    }
	  
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  double radius = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
	  double phiAngle;
	  if(PMTpos[1] >= 0) phiAngle = TMath::ACos(PMTpos[0]/radius);
	  else phiAngle = TMath::Pi() + (TMath::Pi() - TMath::ACos(PMTpos[0]/radius));
	  
	  if(Gamma){//Produce the profile at the gamma conversion point, since it takes ~36 cm to convert/be visible
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStop[j];
	  }
	  else{
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStart[j];
	  }
	  
	  double vDir[3];double vOrientation[3];
	  double vDir2D[3];double vOrientation2D[3];
	  for(int j=0;j<3;j++){
	    vDir[j] = particleRelativePMTpos[j];
	    vOrientation[j] = PMTdir[j];
	    if(j<2){
	      vDir2D[j] = particleRelativePMTpos[j];
	      vOrientation2D[j] = PMTdir[j];
	    }	  
	  }
	  double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
	  double tof = Norm*1e-2/(cvacuum/nindex);
	  double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
	  for(int j=0;j<3;j++){
	    vDir[j] /= Norm;
	    vOrientation[j] /= NormOrientation;
	  }
	  double Norm2D = TMath::Sqrt(vDir2D[0]*vDir2D[0]+vDir2D[1]*vDir2D[1]);
	  double NormOrientation2D = TMath::Sqrt(vOrientation2D[0]*vOrientation2D[0]+vOrientation2D[1]*vOrientation2D[1]);
	  for(int j=0;j<2;j++){
	    vDir2D[j] /= Norm2D;
	    vOrientation2D[j] /= NormOrientation2D;
	  }
	  double HorizontalScalarProduct = TMath::ACos(vDir2D[0]*vOrientation2D[0] + vDir2D[1]*vOrientation2D[1])*180./TMath::Pi();
	  double relativeAngle = TMath::ACos(vDir[0]*particleDir[0]+vDir[1]*particleDir[1]+vDir[2]*particleDir[2])*180./TMath::Pi();
	  //WCSimRootCherenkovHitTime HitTime = dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(timeArrayIndex));//Assumes that the earliest time the PMT is hit is in the first element -> Everything is ordered.
	  //double time = HitTime.GetTruetime();
	  double time = wcsimrootcherenkovdigihit->GetT();
	 
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  
	  //std::cout << "Relative angle = " << relativeAngle << ", time = " << time << endl;
	  
	  if(plotDigitized){
	    //plotBQ(PMTpos[3], PMTdir[3],phiAngle,relativeAngle,peForTube,pmtType);
	    //if(TMath::Abs(PMTpos[2]) < 4500){//Temp
	      if(phiAngle < TMath::Pi()/2 || phiAngle > 3*TMath::Pi()/2) ChargeProfile2D_onlyFront[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2D[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2DTop[pmtType]->Fill(PMTpos[0],PMTpos[1],peForTube);
	      ChargeProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],peForTube);
	      HitProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],1);
	      ChargeProfile2DRelativeAngle[pmtType]->Fill(relativeAngle,PMTpos[2],peForTube);
	      
	      ChargeProfile[pmtType]->Fill(relativeAngle,peForTube);
	      HitProfile[pmtType]->Fill(relativeAngle,1);
	      HitProfilePMTId[pmtType]->Fill(tubeNumber, 1);
	      
	      TimeProfile[pmtType]->Fill(time,peForTube);
	      TimeTOFProfile[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],peForTube);
	      TimeTOFProfileXTOF[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],tof,peForTube);
	      
	      ChargePerPMT[pmtType]->Fill(peForTube);

	      totalPe += peForTube;
	      totalHit ++;
	      //}
	    //if(i%50 == 0) cout << "Hit angle wrt particle direction = " << relativeAngle << ", Phi = " << phiAngle*180/TMath::Pi() << " Z position = " << PMTpos[2] << endl;
	  }
	} // End of loop over Cherenkov hits
      if(verbose) cout << "Total Pe : " << totalPe << endl;
      if(plotDigitized){
	TotalCharge[pmtType]->Fill(totalPe);
	TotalHit[pmtType]->Fill(totalHit);
	TotalChargeXdWall[pmtType]->Fill(dwall,totalPe);
      }
    }
    // Get the number of digitized hits
    // Loop over sub events
    
    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();
    if(hybrid) wcsimrootsuperevent2->ReInitialize();
    
    } // End of loop over events
  //  TCanvas c1("c1"); 
  /*
  float win_scale = 0.75;
  int n_wide(2);
  int n_high(2);
  TCanvas* c1 = new TCanvas("c1", "First canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
  c1->Draw();
  c1->Divide(2,2);
  c1->cd(1); hvtx0->Draw();
  c1->cd(2); hvtx1->Draw();
  c1->cd(3); hvtx2->Draw();
  c1->cd(4); h1->Draw();
  */
  std::cout<<"num_trig "<<num_trig<<"\n";
  TFile * outfile = new TFile("out.root","RECREATE");

  const double scale = 1. / (double)nevent;
  TH2D * HitProfile2DCylindricCoordinate_ALL = new TH2D(*HitProfile2DCylindricCoordinate[0]);
  HitProfile2DCylindricCoordinate_ALL->SetName("HitProfile2DCylindricCoordinate_ALL");
  TH2D * ChargeProfile2DCylindricCoordinate_ALL = new TH2D(*ChargeProfile2DCylindricCoordinate[0]);
  ChargeProfile2DCylindricCoordinate_ALL->SetName("ChargeProfile2DCylindricCoordinate_ALL");
  for(int i=0;i<nPMTtypes;i++){
    
    HitProfile2DCylindricCoordinate[i]->Scale(scale);
    ChargeProfile2DCylindricCoordinate[i]->Scale(scale);
    HitProfilePMTId[i]->Scale(scale);

    ChargeProfile2D_onlyFront[i]->Write();
    ChargeProfile2D[i]->Write();
    ChargeProfile2DTop[i]->Write();
    ChargeProfile2DCylindricCoordinate[i]->Write();
    HitProfile2DCylindricCoordinate[i]->Write();
    ChargeProfile2DRelativeAngle[i]->Write();
    ChargeProfile[i]->Write();
    HitProfile[i]->Write();
    HitProfilePMTId[i]->Write();
    TimeProfile[i]->Write();
    TimeTOFProfile[i]->Write();
    TimeTOFProfileXTOF[i]->Write();
    ChargePerPMT[i]->Write();
    TotalCharge[i]->Write();
    TotalHit[i]->Write();
    TotalChargeXdWall[i]->Write();//Fill(dWall,totalPe);

    if(i > 0) {
      HitProfile2DCylindricCoordinate_ALL->Add(HitProfile2DCylindricCoordinate[i]);
      HitProfile2DCylindricCoordinate_ALL->Scale(1./2.);
      ChargeProfile2DCylindricCoordinate_ALL->Add(ChargeProfile2DCylindricCoordinate[i]);
      ChargeProfile2DCylindricCoordinate_ALL->Scale(1./2.);
    }
  }
  HitProfile2DCylindricCoordinate_ALL->Write();
  ChargeProfile2DCylindricCoordinate_ALL->Write();
  
  outfile->Close();
}
