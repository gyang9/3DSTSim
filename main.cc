#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TTree.h>
#include <TVector3.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TPolyMarker3D.h>
#include <string>
#include <iostream>
#include <fstream>

#include "ND280UpConst.hh"
#include "ND280UpRecoTrack.hh"


int main(int argc, char* argv[])
{

  TH2F * hMuon_CosTh_TrueVsReco = new TH2F("hMuon_CosTh_TrueVsReco","Muon CosTheta True Vs Reco",50,-1,+1,50,-1,+1);
  TH2F * hPion_CosTh_TrueVsReco = new TH2F("hPion_CosTh_TrueVsReco","Pion CosTheta True Vs Reco",50,-1,+1,50,-1,+1);
  TH2F * hProt_CosTh_TrueVsReco = new TH2F("hProt_CosTh_TrueVsReco","Prot CosTheta True Vs Reco",50,-1,+1,50,-1,+1);
  TH2F * hElec_CosTh_TrueVsReco = new TH2F("hElec_CosTh_TrueVsReco","Elec CosTheta True Vs Reco",50,-1,+1,50,-1,+1);

  TH1F * hMuon_CosTh_RecMinTr = new TH1F("hMuon_CosTh_RecMinTr","Muon CosTheta Reco-True",100,-1,+1);
  TH1F * hPion_CosTh_RecMinTr = new TH1F("hPion_CosTh_RecMinTr","Pion CosTheta Reco-True",100,-1,+1);
  TH1F * hProt_CosTh_RecMinTr = new TH1F("hProt_CosTh_RecMinTr","Prot CosTheta Reco-True",100,-1,+1);
  TH1F * hElec_CosTh_RecMinTr = new TH1F("hElec_CosTh_RecMinTr","Elec CosTheta Reco-True",100,-1,+1);

  TH2F * hMuon_Len_TrueVsReco = new TH2F("hMuon_Len_TrueVsReco","Muon Length True Vs Reco",50,0,+3000,50,0,+3000);
  TH2F * hPion_Len_TrueVsReco = new TH2F("hPion_Len_TrueVsReco","Pion Length True Vs Reco",50,0,+3000,50,0,+3000);
  TH2F * hProt_Len_TrueVsReco = new TH2F("hProt_Len_TrueVsReco","Prot Length True Vs Reco",50,0,+3000,50,0,+3000);
  TH2F * hElec_Len_TrueVsReco = new TH2F("hElec_Len_TrueVsReco","Elec Length True Vs Reco",50,0,+3000,50,0,+3000);

  TH1F * hMuon_Len_RecMinTr = new TH1F("hMuon_Len_RecMinTr","Muon Length Reco-True",100,-50,+50);
  TH1F * hPion_Len_RecMinTr = new TH1F("hPion_Len_RecMinTr","Pion Length Reco-True",100,-50,+50);
  TH1F * hProt_Len_RecMinTr = new TH1F("hProt_Len_RecMinTr","Prot Length Reco-True",100,-50,+50);
  TH1F * hElec_Len_RecMinTr = new TH1F("hElec_Len_RecMinTr","Elec Length Reco-True",100,-50,+50);

  TH2F * hMuon_Len_RecMinTr_Vs_TrLen = new TH2F("hMuon_Len_RecMinTr_Vs_TrLen","Muon Length (Reco-True) Vs True",100,-50,+50,50,0,+3000);
  TH2F * hPion_Len_RecMinTr_Vs_TrLen = new TH2F("hPion_Len_RecMinTr_Vs_TrLen","Pion Length (Reco-True) Vs True",100,-50,+50,50,0,+3000);
  TH2F * hProt_Len_RecMinTr_Vs_TrLen = new TH2F("hProt_Len_RecMinTr_Vs_TrLen","Prot Length (Reco-True) Vs True",100,-50,+50,50,0,+3000);
  TH2F * hElec_Len_RecMinTr_Vs_TrLen = new TH2F("hElec_Len_RecMinTr_Vs_TrLen","Elec Length (Reco-True) Vs True",100,-50,+50,50,0,+3000);

  TH2F * hMuon_TrMomVsRecoMom = new TH2F("hMuon_TrMomVsRecoMom","All Muons True Mom Vs reco. Mom",100,0,10000,100,0,10000);

  TH2F * hMuon_TrMomVsTrCosTh = new TH2F("hMuon_TrMomVsTrCosTh","All Muons True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hPion_TrMomVsTrCosTh = new TH2F("hPion_TrMomVsTrCosTh","All Pions True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hProt_TrMomVsTrCosTh = new TH2F("hProt_TrMomVsTrCosTh","All Protons True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hElec_TrMomVsTrCosTh = new TH2F("hElec_TrMomVsTrCosTh","All Electrons True Mom Vs CosTheta",100,0,10000,10,-1,+1);

  TH2F * hMuon_AllIso_TrMomVsTrCosTh = new TH2F("hMuon_AllIso_TrMomVsTrCosTh","All IsoTarget Muons True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hPion_AllIso_TrMomVsTrCosTh = new TH2F("hPion_AllIso_TrMomVsTrCosTh","All IsoTarget Pions True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hProt_AllIso_TrMomVsTrCosTh = new TH2F("hProt_AllIso_TrMomVsTrCosTh","All IsoTarget Protons True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hElec_AllIso_TrMomVsTrCosTh = new TH2F("hElec_AllIso_TrMomVsTrCosTh","All IsoTarget Electrons True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hMuon_EffIso_TrMomVsTrCosTh = new TH2F("hMuon_EffIso_TrMomVsTrCosTh","Efficiency IsoTarget Muons True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hPion_EffIso_TrMomVsTrCosTh = new TH2F("hPion_EffIso_TrMomVsTrCosTh","Efficiency IsoTarget Pions True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hProt_EffIso_TrMomVsTrCosTh = new TH2F("hProt_EffIso_TrMomVsTrCosTh","Efficiency IsoTarget Protons True Mom Vs CosTheta",100,0,10000,10,-1,+1);
  TH2F * hElec_EffIso_TrMomVsTrCosTh = new TH2F("hElec_EffIso_TrMomVsTrCosTh","Efficiency IsoTarget Electrons True Mom Vs CosTheta",100,0,10000,10,-1,+1);

  TH1F * hMuon_AllIso_TrMom = new TH1F("hMuon_AllIso_TrMom","All IsoTarget Muons True Mom",100,0,10000);
  TH1F * hPion_AllIso_TrMom = new TH1F("hPion_AllIso_TrMom","All IsoTarget Pions True Mom",100,0,10000);
  TH1F * hProt_AllIso_TrMom = new TH1F("hProt_AllIso_TrMom","All IsoTarget Protons True Mom",100,0,10000);
  TH1F * hElec_AllIso_TrMom = new TH1F("hElec_AllIso_TrMom","All IsoTarget Electrons True Mom",100,0,10000);
  TH1F * hMuon_EffIso_TrMom = new TH1F("hMuon_EffIso_TrMom","Efficiency IsoTarget Muons True Mom",100,0,10000);
  TH1F * hPion_EffIso_TrMom = new TH1F("hPion_EffIso_TrMom","Efficiency IsoTarget Pions True Mom",100,0,10000);
  TH1F * hProt_EffIso_TrMom = new TH1F("hProt_EffIso_TrMom","Efficiency IsoTarget Protons True Mom",100,0,10000);
  TH1F * hElec_EffIso_TrMom = new TH1F("hElec_EffIso_TrMom","Efficiency IsoTarget Electrons True Mom",100,0,10000);

  TH1F * hMuon_AllIso_TrCosTh = new TH1F("hMuon_AllIso_TrCosTh","All IsoTarget Muons True CosTheta",10,-1,+1);
  TH1F * hPion_AllIso_TrCosTh = new TH1F("hPion_AllIso_TrCosTh","All IsoTarget Pions True CosTheta",10,-1,+1);
  TH1F * hProt_AllIso_TrCosTh = new TH1F("hProt_AllIso_TrCosTh","All IsoTarget Protons True CosTheta",10,-1,+1);
  TH1F * hElec_AllIso_TrCosTh = new TH1F("hElec_AllIso_TrCosTh","All IsoTarget Electrons True CosTheta",10,-1,+1);
  TH1F * hMuon_EffIso_TrCosTh = new TH1F("hMuon_EffIso_TrCosTh","Efficiency IsoTarget Muons True CosTheta",10,-1,+1);
  TH1F * hPion_EffIso_TrCosTh = new TH1F("hPion_EffIso_TrCosTh","Efficiency IsoTarget Pions True CosTheta",10,-1,+1);
  TH1F * hProt_EffIso_TrCosTh = new TH1F("hProt_EffIso_TrCosTh","Efficiency IsoTarget Protons True CosTheta",10,-1,+1);
  TH1F * hElec_EffIso_TrCosTh = new TH1F("hElec_EffIso_TrCosTh","Efficiency IsoTarget Electrons True CosTheta",10,-1,+1);

  TH1F * hMuon_Len = new TH1F("hMuon_Len","Muon track length",100,0,2000); // mm 
  TH1F * hPion_Len = new TH1F("hPion_Len","Pion track length",100,0,2000);
  TH1F * hProt_Len = new TH1F("hProt_Len","Prot track length",100,0,2000);
  TH1F * hElec_Len = new TH1F("hElec_Len","Elec track length",100,0,2000);
  TH1F * hMuon_Stopped_Len = new TH1F("hMuon_Stopped_Len","Muon Stopped track length",100,0,2000); // mm 
  TH1F * hPion_Stopped_Len = new TH1F("hPion_Stopped_Len","Pion Stopped track length",100,0,2000);
  TH1F * hProt_Stopped_Len = new TH1F("hProt_Stopped_Len","Prot Stopped track length",100,0,2000);
  TH1F * hElec_Stopped_Len = new TH1F("hElec_Stopped_Len","Elec Stopped track length",100,0,2000);

  TH1F * hMuon_Edep = new TH1F("hMuon_Edep","Muon energy deposit",200,0,10000); // pe
  TH1F * hPion_Edep = new TH1F("hPion_Edep","Pion energy deposit",200,0,10000);
  TH1F * hProt_Edep = new TH1F("hProt_Edep","Proton energy deposit",200,0,10000);
  TH1F * hElec_Edep = new TH1F("hElec_Edep","Electron energy deposit",200,0,10000);
  TH1F * hMuon_Stopped_Edep = new TH1F("hMuon_Stopped_Edep","Muon Stopped energy deposit",200,0,10000); // pe
  TH1F * hPion_Stopped_Edep = new TH1F("hPion_Stopped_Edep","Pion Stopped energy deposit",200,0,10000);
  TH1F * hProt_Stopped_Edep = new TH1F("hProt_Stopped_Edep","Proton Stopped energy deposit",200,0,10000);
  TH1F * hElec_Stopped_Edep = new TH1F("hElec_Stopped_Edep","Electron Stopped energy deposit",200,0,10000);

  TH1F * hMuon_EdepOverLen = new TH1F("hMuon_EdepOverLen","Muon energy deposit over track length",200,0,200); // ~ pe / mm (cube)
  TH1F * hPion_EdepOverLen = new TH1F("hPion_EdepOverLen","Pion energy deposit over track length",200,0,200);
  TH1F * hProt_EdepOverLen = new TH1F("hProt_EdepOverLen","Prot energy deposit over track length",200,0,200);
  TH1F * hElec_EdepOverLen = new TH1F("hElec_EdepOverLen","Elec energy deposit over track length",200,0,200);
  TH1F * hMuon_Stopped_EdepOverLen = new TH1F("hMuon_Stopped_EdepOverLen","Stopping Muon energy deposit over track length",200,0,200);
  TH1F * hPion_Stopped_EdepOverLen = new TH1F("hPion_Stopped_EdepOverLen","Stopping Pion energy deposit over track length",200,0,200);
  TH1F * hProt_Stopped_EdepOverLen = new TH1F("hProt_Stopped_EdepOverLen","Stopping Prot energy deposit over track length",200,0,200);
  TH1F * hElec_Stopped_EdepOverLen = new TH1F("hElec_Stopped_EdepOverLen","Stopping Elec energy deposit over track length",200,0,200);

  TH2F * hMuon_EdepVsLen = new TH2F("hMuon_EdepVsLen","Muon energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)
  TH2F * hPion_EdepVsLen = new TH2F("hPion_EdepVsLen","Pion energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)
  TH2F * hProt_EdepVsLen = new TH2F("hProt_EdepVsLen","Prot energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)
  TH2F * hElec_EdepVsLen = new TH2F("hElec_EdepVsLen","Elec energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)
  TH2F * hMuon_Stopped_EdepVsLen = new TH2F("hMuon_Stopped_EdepVsLen","Muon Stopped energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)
  TH2F * hPion_Stopped_EdepVsLen = new TH2F("hPion_Stopped_EdepVsLen","Pion Stopped energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)
  TH2F * hProt_Stopped_EdepVsLen = new TH2F("hProt_Stopped_EdepVsLen","Prot Stopped energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)
  TH2F * hElec_Stopped_EdepVsLen = new TH2F("hElec_Stopped_EdepVsLen","Elec Stopped energy deposit Vs track length",200,0,20000,300,0,3000); // ~ pe Vs mm (cube)

  TH1F * hLikeRatio_MuProt_TrueMu = new TH1F("hLikeRatio_MuProt_TrueMu","hLikeRatio_MuProt_TrueMu",100,-10,10);
  TH1F * hLikeRatio_MuProt_TrueProt = new TH1F("hLikeRatio_MuProt_TrueProt","hLikeRatio_MuProt_TrueProt",100,-10,10);
  TH1F * hLikeRatio_PionProt_TruePion = new TH1F("hLikeRatio_PionProt_TruePion","hLikeRatio_PionProt_TruePion",100,-10,10);
  TH1F * hLikeRatio_PionProt_TrueProt = new TH1F("hLikeRatio_PionProt_TrueProt","hLikeRatio_PionProt_TrueProt",100,-10,10);
  TH1F * hLikeRatio_MuProt_TrueMu2 = new TH1F("hLikeRatio_MuProt_TrueMu","hLikeRatio_MuProt_TrueMu",100,-10,10);
  TH1F * hLikeRatio_MuProt_TrueProt2 = new TH1F("hLikeRatio_MuProt_TrueProt","hLikeRatio_MuProt_TrueProt",100,-10,10);
  TH1F * hLikeRatio_PionProt_TruePion2 = new TH1F("hLikeRatio_PionProt_TruePion","hLikeRatio_PionProt_TruePion",100,-10,10);
  TH1F * hLikeRatio_PionProt_TrueProt2 = new TH1F("hLikeRatio_PionProt_TrueProt","hLikeRatio_PionProt_TrueProt",100,-10,10);

int event;
double hitLocation[3],hitPE[3],hitT[3],adc[3],loadc[3],Q[3],hitLowQ[3];
int separateS=0,separateF=0;
string track1Name;
string track2Name;
double trueCos,trueC;
double trueLen,trueL;
double trueMom,trueM;
double trueCos_muon,trueC_muon;
double trueLen_muon,trueL_muon;
double trueMom_muon,trueM_muon;
double trueCos_pion,trueC_pion;
double trueLen_pion,trueL_pion;
double trueMom_pion,trueM_pion;
double trueCos_proton,trueC_proton;
double trueLen_proton,trueL_proton;
double trueMom_proton,trueM_proton;

TH1D* recoEff_ang = new TH1D("reco_eff_ang","reco_eff_ang", 10,-1,1);
TH1D* recoEff_mom = new TH1D("reco_eff_mom","reco_eff_mom", 21,0,2100);
TH1D* recoPass_ang = new TH1D("recoPass_ang","recoPass_ang", 10,-1,1);
TH1D* recoFail_ang = new TH1D("recoFail_ang","recoFail_ang", 10,-1,1);
TH1D* recoPass_mom = new TH1D("recoPass_mom","recoPass_mom", 21,0,2100);
TH1D* recoFail_mom = new TH1D("recoFail_mom","recoFail_mom", 21,0,2100);
TH1D* recoPass_momP = new TH1D("recoPass_momP","recoPass_momP", 21,0,2100);
TH1D* recoFail_momP = new TH1D("recoFail_momP","recoFail_momP", 21,0,2100);

int prim, PDG;
double ener;
int muonID = 0;

TFile file("testEvent.root");
TTree* c = (TTree*)file.Get("EDepSimTree");
c->SetBranchAddress("event",&event);
c->SetBranchAddress("hitLocation",&hitLocation);
c->SetBranchAddress("hitPE",&hitPE);
c->SetBranchAddress("hitT",&hitT);
c->SetBranchAddress("hitADC",&adc);
c->SetBranchAddress("hitLowADC",&loadc);
c->SetBranchAddress("hitQ",&Q);
c->SetBranchAddress("hitLowQ",&hitLowQ);
c->SetBranchAddress("hitPrim",&prim);
c->SetBranchAddress("hitPDG",&PDG);
c->SetBranchAddress("hitE",&ener);
c->SetBranchAddress("trueLen",&trueLen);
c->SetBranchAddress("trueMom",&trueMom);
c->SetBranchAddress("trueCos",&trueCos);

Int_t nevent = c->GetEntries();

Int_t nnevent =2000;
TH2F* hist2D_XY_Q[nnevent];
TH2F* hist2D_XZ_Q[nnevent];
TH2F* hist2D_YZ_Q[nnevent];
TH2F* hist2D_XY_PE[nnevent];
TH2F* hist2D_XZ_PE[nnevent];
TH2F* hist2D_YZ_PE[nnevent];
TH2F* hist2D_XY_ADC[nnevent];
TH2F* hist2D_XZ_ADC[nnevent];
TH2F* hist2D_YZ_ADC[nnevent];

for(Int_t i=0;i<nnevent;i++){
hist2D_XY_Q[i] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_XZ_Q[i] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_YZ_Q[i] = new TH2F("","",100,0,1000,100,0,1000);

hist2D_XY_PE[i] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_XZ_PE[i] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_YZ_PE[i] = new TH2F("","",100,0,1000,100,0,1000);

hist2D_XY_ADC[i] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_XZ_ADC[i] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_YZ_ADC[i] = new TH2F("","",100,0,1000,100,0,1000);
}

TH2F* hist2D_XY_E[10][nnevent];
TH2F* hist2D_XZ_E[10][nnevent];
TH2F* hist2D_YZ_E[10][nnevent];
for(Int_t i=0;i<5;i++){
for(Int_t j=0;j<nnevent;j++){
hist2D_XY_E[i][j] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_YZ_E[i][j] = new TH2F("","",100,0,1000,100,0,1000);
hist2D_XZ_E[i][j] = new TH2F("","",100,0,1000,100,0,1000);
}
}

int eventS=-1;
int eventE=-2;
int seeProton =0;
bool newEvent=false;
int list[100000]={};
int initA=1;
double hitInit[3]={};

for(Int_t ii=0;ii<nevent;ii++){

c->GetEntry(ii);

//if(event == 226 || event ==227) continue;

eventS = event;
//list[eventS] =1;
if(eventS != eventE){newEvent = true;}
else {newEvent = false;}
if(newEvent) {
cout<<"event "<<event<<endl;
cout<<hitLocation[0]<<" "<<hitLocation[1]<<" "<<hitLocation[2]<<endl;
}
//if(list[eventS]==0 && newEvent) cout<<"double new "<<endl;
if( !newEvent ){
//if( event ==1){
hist2D_XY_Q[event]->Fill(hitLocation[0],hitLocation[1],Q[0]+Q[1]);
hist2D_XZ_Q[event]->Fill(hitLocation[0],hitLocation[2],Q[0]+Q[2]);
hist2D_YZ_Q[event]->Fill(hitLocation[1],hitLocation[2],Q[1]+Q[2]);
hist2D_XY_PE[event]->Fill(hitLocation[0],hitLocation[1],hitPE[0]+hitPE[1]);
hist2D_XZ_PE[event]->Fill(hitLocation[0],hitLocation[2],hitPE[0]+hitPE[2]);
hist2D_YZ_PE[event]->Fill(hitLocation[1],hitLocation[2],hitPE[1]+hitPE[2]);
hist2D_XY_ADC[event]->Fill(hitLocation[0],hitLocation[1],adc[0]+adc[1]);
hist2D_XZ_ADC[event]->Fill(hitLocation[0],hitLocation[2],adc[0]+adc[2]);
hist2D_YZ_ADC[event]->Fill(hitLocation[1],hitLocation[2],adc[1]+adc[2]);
//cout<<"PDG and energy "<<PDG<<" "<<ener<<endl;

//proton2212  pion 211  muon 13
if(PDG == 13 || PDG == -13 ){
hist2D_XY_E[0][event]->Fill(hitLocation[0],hitLocation[1],ener);
hist2D_XZ_E[0][event]->Fill(hitLocation[0],hitLocation[2],ener);
hist2D_YZ_E[0][event]->Fill(hitLocation[1],hitLocation[2],ener);
trueCos_muon = trueCos;
trueLen_muon = trueLen;
trueMom_muon = trueMom;
if(initA==1) {
hitInit[0]=hitLocation[0];
hitInit[1]=hitLocation[1];
hitInit[2]=hitLocation[2];
}
initA=0;
}
if(PDG == 211 || PDG == -211 ){
hist2D_XY_E[1][event]->Fill(hitLocation[0],hitLocation[1],ener);
hist2D_XZ_E[1][event]->Fill(hitLocation[0],hitLocation[2],ener);
hist2D_YZ_E[1][event]->Fill(hitLocation[1],hitLocation[2],ener);
trueCos_pion = trueCos;
trueLen_pion = trueLen;
trueMom_pion = trueMom;
}
if( PDG == 2212 ){
hist2D_XY_E[2][event]->Fill(hitLocation[0],hitLocation[1],ener);
hist2D_XZ_E[2][event]->Fill(hitLocation[0],hitLocation[2],ener);
hist2D_YZ_E[2][event]->Fill(hitLocation[1],hitLocation[2],ener);
trueCos_proton = trueCos;
trueLen_proton = trueLen;
trueMom_proton = trueMom;
}
}
//list[ii]=1;
//eventE = event;
//}

if(newEvent && ii>0){

//      if(hist2D_XY_E[0][event-1]->Integral() <=0 || event>=1000 ) continue;
    std::cout << "Going to do tracking for event " <<event<< std::endl;

    //TFile input2D("2DExample.root");

      // Initialize the object for track reconstruction
      ND280UpRecoTrack nd280UpRecoTrack;
/*
      TH2F* pp1 = (TH2F*)input2D.Get("XY_E_muon");
      TH2F* pp2 = (TH2F*)input2D.Get("XZ_E_muon");
      TH2F* pp3 = (TH2F*)input2D.Get("YZ_E_muon");
      TH2F* pp4 = (TH2F*)input2D.Get("XY_E_pion");
      TH2F* pp5 = (TH2F*)input2D.Get("XZ_E_pion");
      TH2F* pp6 = (TH2F*)input2D.Get("YZ_E_pion");

      TH2F* f2d_xy = (TH2F*)input2D.Get("XY_E_muon");
      TH2F* f2d_xz = (TH2F*)input2D.Get("XZ_E_muon");
      TH2F* f2d_yz = (TH2F*)input2D.Get("YZ_E_muon");
      TH2F* f2d_xy_other = (TH2F*)input2D.Get("XY_E_pion");
      TH2F* f2d_xz_other = (TH2F*)input2D.Get("XZ_E_pion");
      TH2F* f2d_yz_other = (TH2F*)input2D.Get("YZ_E_pion");
*/
      nd280UpRecoTrack.SetMPPCXY(hist2D_XY_E[0][event-1]);
      nd280UpRecoTrack.SetMPPCXZ(hist2D_XZ_E[0][event-1]);
      nd280UpRecoTrack.SetMPPCYZ(hist2D_YZ_E[0][event-1]);

      for(int iiii=0;iiii<10;iiii++){
        for(int jjjj=0;jjjj<10;jjjj++){
               // cout<<f2d_xy->GetBinContent(iiii,jjjj)<<endl;
        }
		}

      nd280UpRecoTrack.SetMinPE(0.2); // it is 2. initially

        nd280UpRecoTrack.SetTrackSeparationMin(10); // (mm) min track distance is 1 cube

      nd280upconv::TargetType_t DetType = nd280upconv::kSuperFGD;
      nd280UpRecoTrack.DoTracking(nd280upconv::kSuperFGD); // run the tracking process



      // Check the track separation
      // if the reco track is separated from all the other tracks

      bool isseparated = true;


        // Check separation between current "itrk" and other tracks "itrkoth"

        nd280UpRecoTrack.SetMPPCXY_Other(hist2D_XY_E[1][event-1]);
        nd280UpRecoTrack.SetMPPCXZ_Other(hist2D_XZ_E[1][event-1]);
        nd280UpRecoTrack.SetMPPCYZ_Other(hist2D_YZ_E[1][event-1]);

        nd280UpRecoTrack.DoTrackSeparation();


         isseparated = nd280UpRecoTrack.IsSeparated();

        cout<<"well.. isseparted result is "<<isseparated<<endl;
        if (isseparated) separateS ++;
        if (!isseparated) separateF ++;       
      //
      // Get the Reco
      //

      // if the reco track is reconstructed
      bool isreco = nd280UpRecoTrack.IsReco();
      // if the reco track is OutFV (entering or exiting)
      bool isoutfv = nd280UpRecoTrack.IsOutFV();
      // reco track deposited energy
      double trkedep = nd280UpRecoTrack.GetEdep();
      // reco track length
      double trklen_reco = nd280UpRecoTrack.GetLength();
      double trklen_recoX = nd280UpRecoTrack.GetLengthX();
      double trklen_recoY = nd280UpRecoTrack.GetLengthY();
      double trklen_recoZ = nd280UpRecoTrack.GetLengthZ();
      // reco track costheta
      double costh_reco = nd280UpRecoTrack.GetRecoCosTheta();

  TFile inputPDF("PID_PDF.root");
  TH1F* muonPDF = (TH1F*)inputPDF.Get("Muon");
  TH1F* pionPDF = (TH1F*)inputPDF.Get("Pion");
  TH1F* protPDF = (TH1F*)inputPDF.Get("Prot");
  TH1F* elecPDF = (TH1F*)inputPDF.Get("Elec");

  // PDF is not binning-dependent! -Guang Oct.9 2017
  ND280UpPID nd280UpPID;
  nd280UpPID.SetPDF("Muon",muonPDF); // ok until TSpline3 (no TF1!)
  nd280UpPID.SetPDF("Pion",pionPDF);
  nd280UpPID.SetPDF("Prot",protPDF);
  nd280UpPID.SetPDF("Elec",elecPDF);

        cout<<"track length and energy deposit "<<trklen_reco<<" "<<trkedep<<endl;
        double logratio = nd280UpPID.CalcLogLikeRatio("Muon","Prot",trkedep/trklen_reco);
        double logratio2 = nd280UpPID.CalcLogLikeRatio("Muon","Pion",trkedep/trklen_reco);
        double logratio3 = nd280UpPID.CalcLogLikeRatio("Muon","Elec",trkedep/trklen_reco);
        //hLikeRatio_MuProt_TrueMu->Fill(logratio2);
        cout<<"dE/dx AND PID Muon/Proton and Muon/Pion "<<trkedep/trklen_reco<<" "<< logratio<<" "<<logratio2<<endl;

        hLikeRatio_MuProt_TrueMu->Fill(logratio);
        hLikeRatio_MuProt_TrueProt->Fill(logratio2);
        hLikeRatio_PionProt_TrueProt->Fill(logratio3);
        if(logratio>0 && logratio2>0 ) {track1Name = "muon"; }
        if(logratio<0 && logratio2>0 ) {track1Name = "proton"; }
        if(logratio2<0 ) {track1Name = "pion"; }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Proton

      if(trueMom_proton >0){
      ND280UpRecoTrack nd280UpRecoTrackP;

      nd280UpRecoTrackP.SetMPPCXY(hist2D_XY_E[2][event-1]);
      nd280UpRecoTrackP.SetMPPCXZ(hist2D_XZ_E[2][event-1]);
      nd280UpRecoTrackP.SetMPPCYZ(hist2D_YZ_E[2][event-1]);

      nd280UpRecoTrackP.SetMinPE(0.2); // it is 2. initially

      nd280UpRecoTrackP.SetTrackSeparationMin(10); // (mm) min track distance is 1 cube

      nd280UpRecoTrackP.DoTracking(nd280upconv::kSuperFGD); // run the tracking process


      // Check the track separation
      // if the reco track is separated from all the other tracks

      bool isseparatedP = true;


        // Check separation between current "itrk" and other tracks "itrkoth"

        nd280UpRecoTrackP.SetMPPCXY_Other(hist2D_XY_E[2][event-1]);
        nd280UpRecoTrackP.SetMPPCXZ_Other(hist2D_XZ_E[2][event-1]);
        nd280UpRecoTrackP.SetMPPCYZ_Other(hist2D_YZ_E[2][event-1]);

        nd280UpRecoTrackP.DoTrackSeparation();


        cout<<"well.. isseparted result for Proton is "<<isseparatedP<<endl;
      //
      // Get the Reco
      //

      // if the reco track is reconstructed
      bool isrecoP = nd280UpRecoTrackP.IsReco();
      // if the reco track is OutFV (entering or exiting)
      bool isoutfvP = nd280UpRecoTrackP.IsOutFV();
      // reco track deposited energy
      double trkedepP = nd280UpRecoTrackP.GetEdep();
      // reco track length
      double trklen_recoP = nd280UpRecoTrackP.GetLength();
      double trklen_recoXP = nd280UpRecoTrackP.GetLengthX();
      double trklen_recoYP = nd280UpRecoTrackP.GetLengthY();
      double trklen_recoZP = nd280UpRecoTrackP.GetLengthZ();
      // reco track costheta
      double costh_recoP = nd280UpRecoTrackP.GetRecoCosTheta();

  if(isrecoP) {recoPass_momP->Fill(trueM_proton);}
        else    {recoFail_momP->Fill(trueM_proton);}

       seeProton ++;
      }
//////////////////////////////////////////////////////////////////////////////////////////////////////////


        trueC_muon=trueCos_muon;
        trueL_muon=trueLen_muon;
        trueM_muon=trueMom_muon;
        trueC_pion=trueCos_pion;
        trueL_pion=trueLen_pion;
        trueM_pion=trueMom_pion;
        trueC_proton=trueCos_proton;
        trueL_proton=trueLen_proton;
        trueM_proton=trueMom_proton;

	if(track1Name == "muon"){ muonID++; }

if(hitInit[0]>100 && hitInit[0]<900 && hitInit[1]>100 && hitInit[1]<900 && hitInit[2]>100 && hitInit[2]<900){
        if(track1Name == "muon"){
 		hMuon_Len->Fill(trklen_reco);
		hMuon_Edep->Fill(trkedep);
		hMuon_AllIso_TrCosTh->Fill(costh_reco);

		hMuon_CosTh_TrueVsReco->Fill(costh_reco ,trueC_muon);
		hMuon_Len_TrueVsReco->Fill(trklen_reco ,trueL_muon);
		hMuon_TrMomVsRecoMom->Fill(trkedep ,trueM_muon);	
	}

  if(isreco) {recoPass_ang->Fill(trueC_muon); recoPass_mom->Fill(trueM_muon);}     
        else    {recoFail_ang->Fill(trueC_muon); recoFail_mom->Fill(trueM_muon);}

}
	trueCos_muon = -1;
        trueLen_muon = -1;
        trueMom_muon = -1;
        trueCos_pion = -1;
        trueLen_pion = -1;
        trueMom_pion = -1;
        trueCos_proton = -1;
        trueLen_proton = -1;
        trueMom_proton = -1;
	initA = 1;
/*
      ND280UpRecoTrack nd280UpRecoTrack2;

      nd280UpRecoTrack2.SetMPPCXY(hist2D_XY_E[1][event-1]);
      nd280UpRecoTrack2.SetMPPCXZ(hist2D_XZ_E[1][event-1]);
      nd280UpRecoTrack2.SetMPPCYZ(hist2D_YZ_E[1][event-1]);
      bool isreco2 = nd280UpRecoTrack2.IsReco();
      bool isoutfv2 = nd280UpRecoTrack2.IsOutFV();
      double trkedep2 = nd280UpRecoTrack2.GetEdep();
      double trklen_reco2 = nd280UpRecoTrack2.GetLength();
      double trklen_recoX2 = nd280UpRecoTrack2.GetLengthX();
      double trklen_recoY2 = nd280UpRecoTrack2.GetLengthY();
      double trklen_recoZ2 = nd280UpRecoTrack2.GetLengthZ();
      double costh_reco2 = nd280UpRecoTrack2.GetRecoCosTheta();

      ND280UpPID nd280UpPID2;
      nd280UpPID2.SetPDF("Muon",muonPDF); // ok until TSpline3 (no TF1!)
      nd280UpPID2.SetPDF("Pion",pionPDF);
      nd280UpPID2.SetPDF("Prot",protPDF);
      nd280UpPID2.SetPDF("Elec",elecPDF);
//if(isseparated){
        cout<<endl;
        cout<<"---------------------------------------------------------------------------"<<endl;
        cout<<"Now you are in the second track"<<endl;
        cout<<"track length and energy deposit "<<trklen_reco2<<" "<<trkedep2<<endl;
        double logratioo = nd280UpPID2.CalcLogLikeRatio("Pion","Muon",trkedep2/trklen_reco2);
        double logratioo2 = nd280UpPID2.CalcLogLikeRatio("Pion","Prot",trkedep2/trklen_reco2);
        double logratioo3 = nd280UpPID2.CalcLogLikeRatio("Pion","Elec",trkedep2/trklen_reco2);
        //hLikeRatio_MuProt_TrueMu->Fill(logratio2);
        cout<<"dE/dx AND PID Muon/Proton and Muon/Pion "<<trkedep2/trklen_reco2<<" "<< logratioo<<" "<<logratioo2<<endl;
        cout<<"---------------------------------------------------------------------------"<<endl;
        cout<<endl;

        if(logratioo>0 && logratioo2>0 ) {track1Name = "pion"; }
        if(logratioo2<0 && logratioo>0 ) {track1Name = "proton"; }
        if(logratioo<0 ) {track1Name = "muon"; }
//}
*/


  std::cout<<"summarize information for event: "<<event<<std::endl;
  std::cout<<"if there are multiple tracks: "<<isseparated<<std::endl;
  std::cout<<"if first track successfully reconstructed: "<<isreco<<std::endl;
  std::cout<<"PID for first track: "<<track1Name<<std::endl;
  std::cout<<"First track length, energy deposit and angle to z direction: "<<trklen_reco<<" "<<trkedep<<" "<<costh_reco<<std::endl;
//  std::cout<<"if second track successfully reconstructed: "<<isreco<<std::endl;
//  std::cout<<"PID for second track "<<track2Name<<std::endl;
//  std::cout<<"Second track length, energy deposit and angle to z direction: "<<trklen_reco2<<" "<<trkedep2<<" "<<costh_reco2<<std::endl;
  std::cout<<std::endl;
} 
//list[ii]=1;
eventE = event;
}

  std::cout<<"succesfully separated "<<separateS<<" while not "<<separateF<<std::endl;
  std::cout<<"Muon identified "<<muonID<<" out of "<<event<<std::endl;
  std::cout<<"See proton number: "<<seeProton<<std::endl;
  TFile *fileout = new TFile("outputFile_fid.root","RECREATE");
  fileout->cd();
        hLikeRatio_MuProt_TrueMu->Write("pid-muon-prot");
        hLikeRatio_MuProt_TrueProt->Write("pid-muon-pion");
        hLikeRatio_PionProt_TrueProt->Write("pid-muon-elec");
        hMuon_Len->Write("len");
        hMuon_Edep->Write("edep");
        hMuon_AllIso_TrCosTh->Write("cos");
        hMuon_CosTh_TrueVsReco->Write("cos2D");
        hMuon_Len_TrueVsReco->Write("len2D");
        hMuon_TrMomVsRecoMom->Write("mom2D");
	recoPass_ang->Write("recoPass_ang");
        recoPass_mom->Write("recoPass_mom");
        recoFail_ang->Write("recoFail_ang");
        recoFail_mom->Write("recoFail_mom");
        recoPass_momP->Write("recoPass_momP");
        recoFail_momP->Write("recoFail_momP");
}                                                       
