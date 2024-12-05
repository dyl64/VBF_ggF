#include "HbbISRAnalysis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TRandom3.h"
#include "TTreeReaderValue.h"
#include "TMath.h"
using namespace std;

#include "SelectionTools.C"

void setParameters(){

  intLumi = 10.; //Equiv. data luminosity (fb-1)

  PtBin = "PtAll"; //PtAll >450; Pt02 450-650; Pt03 650-1000; Pt04 > 1000 GeV
  
  minMH   =    30.;
  maxMH   =   210.;

  minEta =   0;
  maxEta = 2.0;

  useMCWeight = true;

  iData = -1;
  hasData = false;

  minFJM  =  60;
  
  iQCD = -1;
  iVjets = -1;
  iTTbar = -1;
  iZbb = -1;
  iZqq = -1;
  iWqq = -1;

  iggH = -1;
  iVBF = -1;
  ittH = -1;
  iVH = -1;

}

void HbbISRRecoilAnalysis() {

  setParameters();
  
  std::vector<TFile*> vFiles;
  std::vector<TFile*> vWeiFiles;

  //Input file vector
  TString fileNames[] = {
    "/fs/ddn/sdf/group/ldmx/users/dongyi/Higgs/evttree-mc23_13p6TeV.601600.PhPy8_PDF4LHC21_ggH125_kt200_MiNLO_bb.merge.DAOD_PHYS.e8472_e8455_s3873_s3874_r14622_p6491.root",
    "/fs/ddn/sdf/group/ldmx/users/dongyi/Higgs/evttree-mc23_13p6TeV.601597.PhPy8EG_PDF4LHC21_VBFH125_bb.DAOD_PHYS.e8532_e8528_s4162_s4114_r14622_r14663_p6491.root",
  };

  //Weight file vector
  TString weiNames[] = {
    "/fs/ddn/sdf/group/ldmx/users/dongyi/Higgs/weitree-mc23_13p6TeV.601600.PhPy8_PDF4LHC21_ggH125_kt200_MiNLO_bb.merge.DAOD_PHYS.e8472_e8455_s3873_s3874_r14622_p6491.root",
    "/fs/ddn/sdf/group/ldmx/users/dongyi/Higgs/weitree-mc23_13p6TeV.601597.PhPy8EG_PDF4LHC21_VBFH125_bb.DAOD_PHYS.e8532_e8528_s4162_s4114_r14622_r14663_p6491.root",
  };

  //Sample Name vector
  std::string sampleNames[] = {
    "ggH",
    "VBF",
  };

  //----- please do not modify unless you are adding processes ------------
  
  int iProcMC[] = {
    601600, //ggH
    601597, //VBF
    604224, //ttH
    700599, //Zbb
    700598, //Zqq
    700441, //Wqq
    410471, //ttbar
    364703, //JZ3
    364704, //JZ4
    364705, //JZ5
    364706  //JZ6
  };

  Float_t xSec[] = {
    0.27386,     //ggH
    4.2207,     //VBF
    0.255,      //ttH    
    13.363,      //Zbb
    47.571,      //Zqq
    144.67,      //Wqq
    729.77,      //ttbar
    26450000.0,  //JZ3
      254610.0,  //JZ4
      4553.2,    //JZ5
      257.54,    //JZ6
  }; //pb
  
  Float_t filterEff[] = {
    1.0,
    1.0,
    1.0,
    1.0,          //Zbb
    1.0,          //Zqq
    1.0,          //Wqq
    4.5623E-01,   //ttbar
    1.165833E-02, //JZ3
    1.336553E-02, //JZ4
    1.452648E-02, //JZ5
    9.4734E-03,   //JZ6
  };

  Float_t kFactor[] = {
    1.0,
    1.0,
    1.0,
    1.0, //Zbb
    1.0, //Zqq
    1.0, //Wqq
    1.13974074379, //ttbar
    1.0, //JZ3
    1.0, //JZ4
    1.0, //JZ5
    1.0, //JZ6
  };

  //---------------------------------
  
  std::vector<TH1F*> vhCutFlow;

  std::vector<TH1F*> vhSRLMass;
  std::vector<TH1F*> vhSRSMass;
  std::vector<TH1F*> vdRL;
  std::vector<TH1F*> vdRS;

  std::vector<TH1F*> vhSRLPt;
  std::vector<TH1F*> vhSRSPt;

  std::vector<TH1F*> vhRecoilMass;
  std::vector<TH1F*> vhRecoilPt;    
  
  const int nFiles = sizeof(fileNames)/sizeof(fileNames[0]);
  const int nProc = sizeof(iProcMC)/sizeof(iProcMC[0]);

  Float_t sumOfWeightMC[nProc];
  
  float nev_FJSel[nFiles];
  float nev_RecSel[nFiles];
  float nev_CandFJ[nFiles];
  float nev_SRL0[nFiles];
  float nev_SRS0[nFiles];
  float nev_SRL[nFiles];
  float nev_SRS[nFiles];
  float nev_SRLM[nFiles];
  float nev_SRSM[nFiles];

  // for testing
  float num_event[nFiles];
  float num_event_0[nFiles];
  float num_event_1[nFiles];
  float num_event_2[nFiles];
  float num_event_3[nFiles];
  float num_leading_4[nFiles];
  float num_subleading_4[nFiles];
  float num_SRL0[nFiles];
  float num_SRS0[nFiles];
  float num_Recoil[nFiles];
  
  if(PtBin == "PtAll"){
    minFJPt =    450.; 
    maxFJPt =   1500.;
  }else  if(PtBin == "Pt01"){
    minFJPt =    300.; 
    maxFJPt =    450.;    
  }else  if(PtBin == "Pt02"){
    minFJPt =    450.; 
    maxFJPt =    650.;
  }else  if(PtBin == "Pt03"){
    minFJPt =    650.; 
    maxFJPt =   1000.;
  }else  if(PtBin == "Pt04"){
    minFJPt =   1000.; 
    maxFJPt =   2000.;
  }

  std::string mnpt = std::to_string((Int_t)minFJPt);
  std::string mxpt = std::to_string((Int_t)maxFJPt);
  std::string ptstr1 = " < p^{FJ}_{T} < ";
  std::string ptstr2 = " GeV ";
  std::string ptstr = mnpt + ptstr1 + mxpt + ptstr2;
  const char* headerPt = ptstr.c_str();

  const char * header;
    
  cout << "Processing " << nFiles << " files " << endl;

  cout << " " << endl;
  cout << "Recoil Analysis: " << endl;
  cout << "Int. Lum. = " << intLumi << " fb-1; " << minFJPt << " < FJ pT < " << maxFJPt << " GeV" << endl;

  Float_t pTGeVPerBin = 50.;
  Float_t massGeVPerBin = 5.;

  Int_t nBinsPt   = (Int_t)((maxFJPt-minFJPt)/pTGeVPerBin);
  Int_t nBinsMass = (Int_t)((maxMH-minMH)/massGeVPerBin);
  
  for (int i=0;i<nFiles;i++) { 

    //input files
    TFile *myFile = new TFile(fileNames[i]);
    vFiles.push_back(myFile);
    TFile *myWeiFile = new TFile(weiNames[i]);
    vWeiFiles.push_back(myWeiFile);

    TH1F * dRL = new TH1F(Form("dRL%d",i),"",50,0,5);
    dRL->Sumw2();
    vdRL.push_back(dRL);

    TH1F * dRS = new TH1F(Form("dRS%d",i),"",50,0,5);
    dRS->Sumw2();
    vdRS.push_back(dRL);

    TH1F * hSRLMass = new TH1F(Form("hSRLMass%d",i),"",nBinsMass,minMH,maxMH);
    hSRLMass->Sumw2();
    vhSRLMass.push_back(hSRLMass);

    TH1F * hSRLPt = new TH1F(Form("hSRLPt%d",i),"",nBinsPt,minFJPt,maxFJPt);
    hSRLPt->Sumw2();
    vhSRLPt.push_back(hSRLPt);
    
    TH1F * hSRSMass = new TH1F(Form("hSRSMass%d",i),"",nBinsMass,minMH,maxMH);
    hSRSMass->Sumw2();
    vhSRSMass.push_back(hSRSMass);

    TH1F * hSRSPt = new TH1F(Form("hSRSPt%d",i),"",nBinsPt,minFJPt,maxFJPt);
    hSRSPt->Sumw2();
    vhSRSPt.push_back(hSRSPt);

    TH1F * hRecoilMass = new TH1F(Form("hRecoilMass%d",i),"",nBinsMass,minMH,maxMH);
    hRecoilMass->Sumw2();
    vhRecoilMass.push_back(hRecoilMass);

    TH1F * hRecoilPt = new TH1F(Form("hRecoilPt%d",i),"",nBinsPt,minFJPt,maxFJPt);
    hRecoilPt->Sumw2();
    vhRecoilPt.push_back(hRecoilPt);    
    
    //add your histos here
    //...
  }

  for (Int_t i=0; i<nFiles; i++) {

    TFile * thisFile = vFiles[i];
    TFile * thisWeiFile = vWeiFiles[i];

    TTreeReader mySumWeightsReader("weitree", thisWeiFile);    
    TTreeReader myReader("evttree", thisFile);

    TTreeReaderValue< int > dsid(mySumWeightsReader, "DSID");
    TTreeReaderValue< float > totalEventsWeighted(mySumWeightsReader, "totalEventsWeighted");
    TTreeReaderValue< int > totalEvents(mySumWeightsReader, "totalEvents");
    
    TTreeReaderValue< int > isMC(myReader, "mcFlag");
    TTreeReaderValue< int > mcProcess(myReader, "mcProcess");
    TTreeReaderValue< float > eventWeightMC(myReader, "eventMCWeight");
    
    TTreeReaderValue< int > runNumber(myReader, "runNumber");

    //Triggers
    TTreeReaderValue< bool > HLT_j440_a10t_lcw_jes_L1J100(myReader, "HLT_j440_a10t_lcw_jes_L1J100");
    TTreeReaderValue< bool > HLT_j420_a10t_lcw_jes_35smcINF_L1SC111(myReader, "HLT_j420_a10t_lcw_jes_35smcINF_L1SC111");
    TTreeReaderValue< bool > HLT_j420_a10t_lcw_jes_35smcINF_L1J100(myReader, "HLT_j420_a10t_lcw_jes_35smcINF_L1J100");
    TTreeReaderValue< bool > HLT_j420_a10_lcw_L1J100(myReader, "HLT_j420_a10_lcw_L1J100");
    TTreeReaderValue< bool > HLT_j390_a10t_lcw_jes_30smcINF_L1J100(myReader, "HLT_j390_a10t_lcw_jes_30smcINF_L1J100");
    TTreeReaderValue< bool > HLT_j460_a10t_lcw_jes_L1J100(myReader, "HLT_j460_a10t_lcw_jes_L1J100");
    TTreeReaderValue< bool > HLT_j360_a10_lcw_sub_L1J100(myReader, "HLT_j360_a10_lcw_sub_L1J100");  

    //Missing ET
    TTreeReaderValue< float > MET(myReader, "MET");
    TTreeReaderValue< float > phiMET(myReader, "phiMET");

    //Boson Truth info
    TTreeReaderValue< vector<int> > bosonPdgId(myReader, "bosonPdgId");    
    TTreeReaderValue< vector<float> > bosonM(myReader, "bosonM");
    TTreeReaderValue< vector<float> > bosonPt(myReader, "bosonPt");
    TTreeReaderValue< vector<float> > bosonPx(myReader, "bosonPx");
    TTreeReaderValue< vector<float> > bosonPy(myReader, "bosonPy");
    TTreeReaderValue< vector<float> > bosonPz(myReader, "bosonPz");

    //Large-R Jet info
    TTreeReaderValue< vector<float> > fatJetPx(myReader, "fatJetPx");
    TTreeReaderValue< vector<float> > fatJetPy(myReader, "fatJetPy");
    TTreeReaderValue< vector<float> > fatJetPz(myReader, "fatJetPz");
    TTreeReaderValue< vector<float> > fatJetE(myReader, "fatJetE");    
    TTreeReaderValue< vector<float> > fatJetPt(myReader, "fatJetPt");		
    TTreeReaderValue< vector<float> > fatJetM(myReader, "fatJetM");

    //Large-R Jet truth level info    
    TTreeReaderValue< vector<float> > fatJetTruthPx(myReader, "fatJetTruthPx");
    TTreeReaderValue< vector<float> > fatJetTruthPy(myReader, "fatJetTruthPy");
    TTreeReaderValue< vector<float> > fatJetTruthPz(myReader, "fatJetTruthPz");
    TTreeReaderValue< vector<float> > fatJetTruthPt(myReader, "fatJetTruthPt");
    TTreeReaderValue< vector<float> > fatJetTruthE(myReader, "fatJetTruthE");
    TTreeReaderValue< vector<float> > fatJetTruthM(myReader, "fatJetTruthM");

    //Nb of b and c-hadrons in large-R Jet 
    TTreeReaderValue< vector<int> > fatJetNBHadrons(myReader, "fatJetNBHadrons");
    TTreeReaderValue< vector<int> > fatJetNCHadrons(myReader, "fatJetNCHadrons");

    //anti-kT4 Jet info
    TTreeReaderValue< vector<float> > antiKt4JetPx(myReader, "jetPx");
    TTreeReaderValue< vector<float> > antiKt4JetPy(myReader, "jetPy");
    TTreeReaderValue< vector<float> > antiKt4JetPz(myReader, "jetPz");
    TTreeReaderValue< vector<float> > antiKt4JetE(myReader, "jetE");    
    TTreeReaderValue< vector<float> > antiKt4JetPt(myReader, "jetPt");		
    TTreeReaderValue< vector<float> > antiKt4JetM(myReader, "jetM");
    
    Int_t nEntries = 0;
    Float_t sumOfEventWeights = 0.;
    Float_t sumOfFillWeights = 0.;
    
    nev_FJSel[i] = 0.;
    nev_RecSel[i] = 0.;
    nev_CandFJ[i] = 0.;
    nev_SRL0[i] = 0.;
    nev_SRS0[i] = 0.;    
    nev_SRL[i] = 0.;
    nev_SRS[i] = 0.;
    nev_SRLM[i] = 0.;
    nev_SRSM[i] = 0.;  

    // for testing
    num_event_0[i] = 0.;
    num_event_1[i] = 0.;
    num_event_2[i] = 0.;
    num_event_3[i] = 0.;
    num_leading_4[i] = 0.;
    num_subleading_4[i] = 0.;
    num_SRL0[i] = 0.;
    num_SRS0[i] = 0.;
    num_Recoil[i] = 0.;

    for(Int_t ip=0; ip < nProc; ip++){
      sumOfWeightMC[ip] = 0.;
    }
    
    cout << " " << endl;
    cout << "Reading from file " << fileNames[i] << endl;
    cout << " " << endl;

    while (mySumWeightsReader.Next()) {
      cout << *dsid << " " << *totalEventsWeighted << " " << *totalEvents << endl;
      for(Int_t ip=0; ip < nProc; ip++){
	if(iProcMC[ip] == *dsid){
	  sumOfWeightMC[ip] += *totalEventsWeighted;
	}
      }
    }

    // Loop over TTree entries    
    while (myReader.Next()) {

      nEntries++;    

      Float_t wei = 1;
      if(useMCWeight && *isMC){
	for(Int_t ip=0; ip < nProc; ip++){
	  if(iProcMC[ip] == *mcProcess){
	    wei = xSec[ip] * kFactor[ip] * filterEff[ip] * intLumi * 1000. * (*eventWeightMC)/(sumOfWeightMC[ip]);
	  }
	}
      }

      if(!*isMC){
	hasData = true;
	iData = i;
      }else if(*mcProcess == 601600){
	iggH = i;	
      }else if(*mcProcess == 601597){
	iVBF = i;
      }else if(*mcProcess == 604224){
	ittH = i;	
      }

      sumOfEventWeights += (*eventWeightMC);
      sumOfFillWeights += wei;

      //Higgs boson in event
      
      TVector3 boson(0.,0.,0.);
      if(*isMC){
	for(Int_t ib=0; ib<bosonPdgId->size(); ib++){
	  if((bosonPdgId->at(ib) == 23 && (*mcProcess == 700855 || *mcProcess == 700849)) || (abs(bosonPdgId->at(ib)) == 24 && (*mcProcess == 700843)) || (bosonPdgId->at(ib) == 25 && (*mcProcess == 601600 || *mcProcess == 601597 || *mcProcess == 604224))){
	    boson.SetXYZ(bosonPx->at(ib), bosonPy->at(ib), bosonPz->at(ib));
	    break;
	  }
	}
      }

      if (Nan4MomVec(*fatJetM)) continue;
      if (Nan4MomVec(*fatJetPt)) continue;
      
      Bool_t null = true;

      //PRD event selection
      
      std::vector<int> list_reorderd_fj = ListSortedFJ(*fatJetPt);
      
      if(list_reorderd_fj.size()<1) continue;
      //testing
      num_event_0[i] += 1;

      // >= 1 FJ pT>450 M>60
      
      int firstSelection = FirstEventSelection_ST(list_reorderd_fj, null, *fatJetPt, *fatJetM, minFJPt, minFJM);

      
      if(firstSelection==0) continue;
      num_event_1[i] += 1; // testing


      //>= 2 FJ pT>200 |Eta|<2.0
      if(list_reorderd_fj.size()<2) continue;
      num_event_2[i] += 1; // testing

      //recoil to leading
      
      TLorentzVector fatJet;
      fatJet.SetPxPyPzE(fatJetPx->at(list_reorderd_fj.at(1)), fatJetPy->at(list_reorderd_fj.at(1)), fatJetPz->at(list_reorderd_fj.at(1)), fatJetE->at(list_reorderd_fj.at(1)));

      if(fatJetPt->at(list_reorderd_fj.at(1)) < 350 || fabs(fatJet.Eta()) > maxEta)continue;
      num_event_3[i] += 1; // testing

      
      for (int iL = 0; iL < 2; iL++ ){

	Int_t iFJ = list_reorderd_fj.at(iL);

	TLorentzVector thisFatJet;
	
	thisFatJet.SetPxPyPzE(fatJetPx->at(iFJ), fatJetPy->at(iFJ), fatJetPz->at(iFJ), fatJetE->at(iFJ));

	if(fatJetPt->at(iFJ) > minFJPt && fatJetPt->at(iFJ) < maxFJPt && boson.Pt()>50.){
	  //testing
	  //cout << "######################boson Pt " << boson.Pt() << endl;
	  if(iL==0){ 
          num_leading_4[i] += 1;
          ((TH1F*)vdRL.at(i))->Fill(thisFatJet.Vect().DeltaR(boson));
          }
          if(iL==1){
          num_subleading_4[i] += 1;
          ((TH1F*)vdRS.at(i))->Fill(thisFatJet.Vect().DeltaR(boson));
          }
	  // signal large-R Jet bb matched to true boson
	  // recoil is anything not included in this jet (i.e. thisFatJet.Vect().DeltaR(XXX) > 1.0	  
	  if(fatJetNBHadrons->at(iFJ) == 2 && thisFatJet.Vect().DeltaR(boson) < 0.30){
	    
	    if(iL==0){
	      //leading	
              cout << "dR " << thisFatJet.Vect().DeltaR(boson) << endl; 
	      nev_SRL0[i] += wei;
              num_SRL0[i] += 1;
	      ((TH1F*)vhSRLMass.at(i))->Fill(thisFatJet.M(), wei);
	      ((TH1F*)vhSRLPt.at(i))->Fill(thisFatJet.Pt(), wei);
	      
	    }else if(iL==1){
	      //sub-leading	
	      nev_SRS0[i] += wei;
	      num_SRS0[i] += 1;
	      
	      ((TH1F*)vhSRSMass.at(i))->Fill(thisFatJet.M(), wei);
	      ((TH1F*)vhSRSPt.at(i))->Fill(thisFatJet.Pt(), wei);
	    }
	    
	  }else{
            num_Recoil[i] += 1;
	    
	  //recoil
	    ((TH1F*)vhRecoilMass.at(i))->Fill(thisFatJet.M(), wei);
	    ((TH1F*)vhRecoilPt.at(i))->Fill(thisFatJet.Pt(), wei);
	  }
	}
      }
    }
    
    cout<< "Sample: " << i << " " << sampleNames[i] << endl;
    cout << " " << endl;
    
    for(Int_t ip=0; ip < nProc; ip++){
      if(sumOfWeightMC[ip]>0)cout << "DSID " << iProcMC[ip] << ": xSec " << xSec[ip] << " pb; Filter Eff " << filterEff[ip] << "; Sum of MC Weights " << sumOfWeightMC[ip] << "; Sum of Event Weights " << sumOfEventWeights << "; <Weight> " << sumOfFillWeights/(float)nEntries << "; Nb of Entries " << nEntries << endl;
    }

    cout << " " << endl;
    cout<<"Lead " << nev_SRL0[i] << "; subLead " <<  nev_SRS0[i] << endl;
    cout << "fatjet size >= 1 cut: " << num_event_0[i] << "; firstselection cut: " << num_event_1[i] << "; fatjet size >= 2 cut: " << num_event_2[i] << "; subleading fatjet pT > 350 & eta < 2.0 cut: " << num_event_3[i] << "; 450 < fatjetpT < 1500 & bosonpT > 50 cut: leading jet # " << num_leading_4[i] << ", subleading jet # " << num_subleading_4[i] << "; # of SRL: " << num_SRL0[i] << ", # of SRS " << num_SRS0[i] << ", # of Recoil jet " << num_Recoil[i] << endl;
    
  }

   TCanvas * c0dRL = new TCanvas("c0dRL", "dR between boson and leading fatjet", 650, 600);

   TH1F * h0dRL = new TH1F("h0dRL","",50,0,5);
   h0dRL->SetMaximum(1.5*((TH1F*)vdRL.at(0))->GetMaximum());
   h0dRL->GetXaxis()->SetTitle("dR between boson and leading fatjet");	   
   h0dRL->Draw();
   
   for(Int_t i=0; i<vdRL.size();i++){
     ((TH1F*)vdRL.at(i))->SetMarkerStyle(20+i);
     ((TH1F*)vdRL.at(i))->Draw("sameE1");
   }
   
   TLegend * leg0dRL = new TLegend(0.175, 0.75, 0.675, 0.895);
   leg0dRL->SetBorderSize(0);
   leg0dRL->SetFillColor(0); 
   leg0dRL->SetTextSize(0.03);
   leg0dRL->SetTextFont(42);
   leg0dRL->SetHeader(headerPt);
   for(Int_t i=0; i<vdRL.size();i++){
     leg0dRL->AddEntry((TH1F*)vdRL.at(i),sampleNames[i].c_str(),"P");
   }
   leg0dRL->Draw();

   TCanvas * c0dRS = new TCanvas("c0dRS", "dR between boson and subleading fatjet", 650, 600);

   TH1F * h0dRS = new TH1F("h0dRS","",50,0,5);
   h0dRS->SetMaximum(1.5*((TH1F*)vdRS.at(0))->GetMaximum());
   h0dRS->GetXaxis()->SetTitle("dR between boson and fatjet");	   
   h0dRS->Draw();
   
   for(Int_t i=0; i<vdRS.size();i++){
     ((TH1F*)vdRS.at(i))->SetMarkerStyle(20+i);
     ((TH1F*)vdRS.at(i))->Draw("sameE1");
   }
   
   TLegend * leg0dRS = new TLegend(0.175, 0.75, 0.675, 0.895);
   leg0dRS->SetBorderSize(0);
   leg0dRS->SetFillColor(0); 
   leg0dRS->SetTextSize(0.03);
   leg0dRS->SetTextFont(42);
   leg0dRS->SetHeader(headerPt);
   for(Int_t i=0; i<vdRS.size();i++){
     leg0dRS->AddEntry((TH1F*)vdRS.at(i),sampleNames[i].c_str(),"P");
   }
   leg0dRS->Draw();


   TCanvas * c0M = new TCanvas("c0M", "Signal Mass", 650, 600);

   TH1F * h0M = new TH1F("h0M","",100,minMH,maxMH);
   h0M->SetMaximum(1.5*((TH1F*)vhSRLMass.at(0))->GetMaximum());
   h0M->GetXaxis()->SetTitle("Signal Leading Large-R Jet Mass [GeV]");	   
   h0M->Draw();
   
   for(Int_t i=0; i<vhSRLMass.size();i++){
     ((TH1F*)vhSRLMass.at(i))->SetMarkerStyle(20+i);
     ((TH1F*)vhSRLMass.at(i))->Draw("sameE1");
   }
   
   TLegend * leg0M = new TLegend(0.175, 0.75, 0.675, 0.895);
   leg0M->SetBorderSize(0);
   leg0M->SetFillColor(0); 
   leg0M->SetTextSize(0.03);
   leg0M->SetTextFont(42);
   leg0M->SetHeader(headerPt);
   for(Int_t i=0; i<vhSRLMass.size();i++){
     leg0M->AddEntry((TH1F*)vhSRLMass.at(i),sampleNames[i].c_str(),"P");
   }
   leg0M->Draw();

   TCanvas * c0RM = new TCanvas("c0RM", "Recoil Mass", 650, 600);

   TH1F * h0RM = new TH1F("h0RM","",100,minMH,maxMH);
   h0RM->SetMaximum(1.5*((TH1F*)vhSRLMass.at(0))->GetMaximum());
   h0RM->GetXaxis()->SetTitle("Recoil Large-R Jet Mass [GeV]");	   
   h0RM->Draw();
   
   for(Int_t i=0; i<vhRecoilMass.size();i++){
     ((TH1F*)vhRecoilMass.at(i))->SetMarkerStyle(20+i);
     ((TH1F*)vhRecoilMass.at(i))->Draw("sameE1");
   }
   
   TLegend * leg0RM = new TLegend(0.175, 0.75, 0.675, 0.895);
   leg0RM->SetBorderSize(0);
   leg0RM->SetFillColor(0); 
   leg0RM->SetTextSize(0.03);
   leg0RM->SetTextFont(42);
   leg0RM->SetHeader(headerPt);
   for(Int_t i=0; i<vhRecoilMass.size();i++){
     leg0RM->AddEntry((TH1F*)vhRecoilMass.at(i),sampleNames[i].c_str(),"P");
   }
   leg0RM->Draw();   
   
}

