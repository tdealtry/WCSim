#include <iostream>
#include <TH1F.h>
#include <stdio.h>     
#include <stdlib.h>    
// Simple example of reading a generated Root file
void Benjamin_readfile(char *filename=NULL, bool verbose=true)
{
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

  TFile *file;
  // Open the file
  if (filename==NULL){
    //file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt_e500_center_nominal.root","read");
    //file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hkmpmt_e500_center_nominal.root","read");
    file = new TFile("/disk01/usr5/bquilain/WCSimData/wcsim_hkmpmt_e500_center_nominal.root","read");
  }else{
    file = new TFile(filename,"read");
  }
  if (!file->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }
  
  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");
  
  // Get the number of events
  int nevent = tree->GetEntries();
  if(verbose){
    printf("nevent %d\n",nevent);
    tree->Print();
  }
  // Create a WCSimRootEvent to put stuff from the tree in

  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  //branch->
  //TBranch *branch2 = tree->GetBranch("wcsimrootevent2");
  //branch2->SetAddress(&wcsimrootsuperevent2);

  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
  //tree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);


  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = 0; 
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
      exit(9);
  }
  geotree->GetEntry(0);

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

  TH1F *h1 = new TH1F("PMT Hits", "PMT Hits", 8000, 0, 8000);
  TH1F *hvtx0 = new TH1F("Event VTX0", "Event VTX0", 200, -1500, 1500);
  TH1F *hvtx1 = new TH1F("Event VTX1", "Event VTX1", 200, -1500, 1500);
  TH1F *hvtx2 = new TH1F("Event VTX2", "Event VTX2", 200, -1500, 1500);
  TH1F *hnphotons = new TH1F("hnphotons", "hnphotons", 150, 0, 150000);
  TH1F *hnphotons2 = new TH1F("hnphotons2", "hnphotons2", 150, 0, 150000);
  TH1F *hnhits = new TH1F("hnhits", "hnhits", 150, 0, 150000);

  bool HK=true;;
  double TankSize = 500;//in cm, the maximal size of the plots for vertices etc...      
  double TotalNHits = 1e4;

  if(HK){
    TankSize = 5000;
    TotalNHits = 1e5;
    //nevt = std::min( 1e3 , ((double) WSTreeHits->GetEntries()) );
  }
  //cout<<"Number of events: "<<nevt<<endl;

  //TFile * fOutput = new TFile(OutputFile,"recreate");
  TH1D * ChargeProfile = new TH1D("ChargeProfile","",720,0,180);
  TH1D * ChargeProfileRadius = new TH1D("ChargeProfileRadius","",400,0,4000);
  TH1D * HitProfile = new TH1D("HitProfile","",720,0,180);
  TH1D * TimeProfile = new TH1D("TimeProfile","",TankSize,0,TankSize);
  TH1D * TimeHitProfile = new TH1D("TimeHitProfile","",TankSize,0,TankSize);
  TH1D * ChargePerPMT = new TH1D("ChargePerPMT","",500,0,500);
  ChargeProfile->Sumw2(); HitProfile->Sumw2();TimeProfile->Sumw2();TimeHitProfile->Sumw2();ChargePerPMT->Sumw2();
  //TH2D * ChargeProfileXTime = new TH2D("ChargeProfileXTime","",36,0,180,);

  TH2D * ChargeProfile2D = new TH2D("ChargeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * HitProfile2D = new TH2D("HitProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * TimeProfile2D = new TH2D("TimeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * ChargeXPositionXTime = new TH2D("ChargeXPositionXTime","",36,0,180,TankSize,0,TankSize);
  //for(int ibinx=1;ibinx<=)
  TH2D * ChargeProfile2DXDWall = new TH2D("ChargeProfile2Dbis","",720,0,180,100,0,TankSize);

  TH1D * TotalCharge = new TH1D("TotalCharge","",200,0,1.5e5);
  TH1D * TotalHit = new TH1D("TotalHit","",200,0,TotalNHits);
  TotalCharge->Sumw2();TotalHit->Sumw2();

  TH3D * VertexPosition = new TH3D("VertexPosition","",100,-TankSize,TankSize,100,-TankSize,TankSize,100,-TankSize,TankSize);
  TH3D * VertexDirection = new TH3D("VertexDirection","",100,-1.1,1.1,100,-1.1,1.1,100,-1.1,1.1);


  int num_trig=0;
  
  // Now loop over events
  for (int ev=0; ev<nevent; ev++)
  {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);      
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
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
    // Now read the tracks in the event
    
    // Get the number of tracks
    int ntrack = wcsimrootevent->GetNtrack();
    if(verbose) printf("ntracks=%d\n",ntrack);
    
    int i;
    // Loop through elements in the TClonesArray of WCSimTracks
    double particleDir[3];
    for (i=0; i<ntrack; i++)
    {
      TObject *element = (wcsimrootevent->GetTracks())->At(i);
      
      WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);

      if(verbose){
	printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
      
	for (int j=0; j<3; j++)
	  printf("Track dir: %d %f\n",j, wcsimroottrack->GetDir(j));
	printf("Track energy: %f\n", wcsimroottrack->GetE());
	printf("Track momentum: %f\n", wcsimroottrack->GetP());
	printf("Track mass: %f\n", wcsimroottrack->GetM());
	particleDir[j] = wcsimroottrack->GetDir(j);
      }

      
    }  // End of loop over tracks
    
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
    
    h1->Fill(ncherenkovdigihits);
    if(verbose){
      printf("node id: %i\n", ev);
      printf("Ncherenkovhits %d\n",     ncherenkovhits);
      printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      cout << "RAW HITS:" << endl;
    }



    // Grab the big arrays of times and parent IDs
    TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();
    hnphotons->Fill(wcsimrootevent->GetCherenkovHitTimes()->GetEntries());
    hnhits->Fill(ncherenkovhits);
    
    int totalPe = 0;
    int TotHit = 0;
    double TotCharge = 0;
    // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
    for (i=0; i< ncherenkovhits; i++)
    {
      TObject *Hit = (wcsimrootevent->GetCherenkovHits())->At(i);
      WCSimRootCherenkovHit *wcsimrootcherenkovhit = 
	dynamic_cast<WCSimRootCherenkovHit*>(Hit);

      int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
      int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
      int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
      WCSimRootPMT pmt   = geo->GetPMT(tubeNumber-1);
      totalPe += peForTube;

     
      if ( i < 10 ) // Only print first XX=10 tubes
      {
	if(verbose) printf("Total pe: %d times( ",peForTube);
	for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
	{
	  WCSimRootCherenkovHitTime HitTime = 
	    dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(j));
	  
	  if(verbose) printf("%6.2f ", HitTime.GetTruetime() );	     
	}
	if(verbose) cout << ")" << endl;
      }
      //What is needed: PMT position, PMT Orientation
      //What is needed: Particle Direction
      double PMTpos[3];
      double PMTdir[3];
      for(int i=0;i<3;i++){
	PMTpos[i] = pmt.GetPosition(i);
	PMTdir[i] = pmt.GetOrientation(i);
      }	
      
      double rPMT = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
      double cosPMT = (180./TMath::Pi())*TMath::ACos(PMTpos[1] / rPMT);
      
      double vDir[3] = {PMTpos[0],PMTpos[1],PMTpos[2]};
      double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
      double vOrientation[3] = {PMTdir[0],PMTdir[1],PMTdir[2]};
      double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
      for(int i=0;i<3;i++){
	vDir[i] /= Norm;
	vOrientation[i] /= NormOrientation;
      }
      double vDir2D[2] = {PMTpos[0],PMTpos[1]};
      double Norm2D = TMath::Sqrt(vDir2D[0]*vDir2D[0]+vDir2D[1]*vDir2D[1]);
      double vOrientation2D[2] = {PMTdir[0],PMTdir[1]};
      double NormOrientation2D = TMath::Sqrt(vOrientation2D[0]*vOrientation2D[0]+vOrientation2D[1]*vOrientation2D[1]);
      for(int i=0;i<2;i++){
	vDir2D[i] /= Norm2D;
	vOrientation2D[i] /= NormOrientation2D;
      }
      double HorizontalScalarProduct = TMath::ACos(vDir2D[0]*vOrientation2D[0] + vDir2D[1]*vOrientation2D[1])*180./TMath::Pi();
      //if(TMath::Abs(HorizontalScalarProduct - 179.422) > 0.1) continue;
      //if(TMath::Abs(PMTpos[2]) > 2700) continue;
      double relativeAngle = TMath::ACos(vDir[0]*particleDir[0]+vDir[1]*particleDir[0]+vDir[2]*particleDir[0])*180./TMath::Pi();
      if(i%1000 == 0){
	cout << "Tube number = " << tubeNumber <<endl;
	cout<<"PMT pos (x,y,z)= "<<PMTpos[0]<<","<<PMTpos[1]<<","<<PMTpos[2]<<endl;
	cout<<"PMT pos (r,theta,z)= "<<rPMT<<","<<cosPMT<<","<<PMTpos[2]<<endl;
	cout<<"PMT dir = "<<vDir[0]<<","<<vDir[1]<<","<<vDir[2]<<endl;
	cout<<"PMT orientation = "<<vOrientation[0]<<","<<vOrientation[1]<<","<<vOrientation[2]<<endl;
	cout<<"Product Scal.= "<<HorizontalScalarProduct<<endl;
	cout<<"Relative angle = "<<relativeAngle<<endl;
	cout<<"Time = "<<timeArrayIndex<<endl;
      }
      //double relativeDir[3] = {vDir[0]-Dirx[0],vDir[1]-Diry[0],vDir[2]-Dirz[0]};
      //cout<<"Relative dir = "<<relativeDir[0]<<", "<<relativeDir[1]<<", "<<relativeDir[1]<<endl;
      //cout<<"PMT radius = "<<rPMT<<", z = "<<PMTpos[2]<<", ang = "<<cosPMT<<endl;
      ChargeProfileRadius->Fill(rPMT,peForTube);
      ChargeProfile->Fill(relativeAngle,peForTube);
      HitProfile->Fill(relativeAngle);
      HitProfile2D->Fill(PMTpos[2],PMTpos[1]);
      ChargeProfile2D->Fill(PMTpos[2],PMTpos[1],peForTube);
      //ChargeProfile2D_PerEvent[ievt]->Fill(PMTpos[2],PMTpos[1],peForTube);
      //ChargeProfile2DXDWall->Fill(relativeAngle,dwall,peForTube);
      TimeProfile2D->Fill(PMTpos[2],PMTpos[1],timeArrayIndex);
      ChargeXPositionXTime->Fill(relativeAngle,timeArrayIndex,peForTube);
      TimeProfile->Fill(timeArrayIndex,peForTube);
      TimeHitProfile->Fill(timeArrayIndex);
      ChargePerPMT->Fill(peForTube);
      TotHit++;
      TotCharge += peForTube;


    } // End of loop over Cherenkov hits
    if(verbose) cout << "Total Pe : " << totalPe << endl;
    TotalCharge->Fill(TotCharge);
    TotalHit->Fill(TotHit);
    //cout<<"Total charge counted = "<<TotCharge<<" / expected = "<<QTot<<", and NHits = "<<TotHit<<" / expected = "<<NHits<<endl<<endl;
    hnphotons2->Fill(totalPe);
    
    // Look at digitized hit info

    // Get the number of digitized hits
    // Loop over sub events
   
    if(verbose) cout << "DIGITIZED HITS:" << endl;
    for (int index = 0 ; index < wcsimrootsuperevent->GetNumberOfEvents(); index++) 
    {
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
      if(verbose) cout << "Sub event number = " << index << "\n";
      
      int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
      if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
     
      if(ncherenkovdigihits>0)
	num_trig++;
      //for (i=0;i<(ncherenkovdigihits>4 ? 4 : ncherenkovdigihits);i++){
      for (i=0;i<ncherenkovdigihits;i++)
      {
    	// Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
	
    	TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
	
    	WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
    	  dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
	
	if(verbose){
	  if ( i < 10 ) // Only print first XX=10 tubes
	    printf("q, t, tubeid: %f %f %d \n",wcsimrootcherenkovdigihit->GetQ(),
		   wcsimrootcherenkovdigihit->GetT(),wcsimrootcherenkovdigihit->GetTubeId());
	}
      } // End of loop over Cherenkov digihits
    } // End of loop over trigger
    
    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();
    wcsimrootsuperevent2->ReInitialize();
    
  } // End of loop over events
  //  TCanvas c1("c1"); 
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
  TCanvas* c2 = new TCanvas("c2", "Second canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
  c2->Divide(2,2);
  c2->cd(1);
  hnphotons->Draw();
  c2->cd(2);
  hnphotons2->Draw();
  c2->cd(3);
  hnhits->Draw();

  cout<<"Writing output"<<endl;
  cout<<"Writing output"<<endl;
  ChargeProfileRadius->Write();
  ChargeProfile->Write();
  cout<<"Writing output"<<endl;
  HitProfile->Write();
  TimeProfile->Write();
  TimeHitProfile->Write();
  ChargeProfile2D->Write();
  ChargeProfile2DXDWall->Write();
  HitProfile2D->Write();
  TimeProfile2D->Write();
  ChargeXPositionXTime->Write();
  TotalCharge->Write();
  TotalHit->Write();
  VertexPosition->Write();
  VertexDirection->Write();
  ChargePerPMT->Write();
  fOutput->Close();
  
  std::cout<<"num_trig "<<num_trig<<"\n";
}
