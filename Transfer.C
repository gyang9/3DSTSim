{
#include "/home/gyang/work/edep-sim/src/TG4Event.h"

Double_t E,startX,startY,startZ,stopX,stopY,stopZ,startT,stopT;
Int_t eventN;
string det;
Int_t iii=0;

  TH1F * hMuon_Stopped_EdepOverLen = new TH1F("hMuon_Stopped_EdepOverLen","Stopping Muon energy deposit over track length",200,0,40);
  TH1F * hPion_Stopped_EdepOverLen = new TH1F("hPion_Stopped_EdepOverLen","Stopping Pion energy deposit over track length",200,0,40);
  TH1F * hProt_Stopped_EdepOverLen = new TH1F("hProt_Stopped_EdepOverLen","Stopping Prot energy deposit over track length",200,0,40);
  TH1F * hElec_Stopped_EdepOverLen = new TH1F("hElec_Stopped_EdepOverLen","Stopping Elec energy deposit over track length",200,0,40);

//constants for energy calibration
const double CBIRKS = 0.00208; // mm/MeV
const double EdepToPhotConv_FGD = 70.8; // CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 
const double DistMPPCscint_FGD = 41; //*CLHEP::mm;
const double LongCompFrac_FGD = 0.816;
const double LongAtt_FGD = 11926.; //*CLHEP::mm;
const double ShortAtt_FGD = 312.; //*CLHEP::mm;
const double DecayLength_FGD = 0.0858; // CLHEP::mm;
const double Lbar_FGD = 1864.3; //* CLHEP::mm;
const double TransTimeInFiber = 1./28. *10.; //mm/ns
// SuperFGD constants
const double MPPCEff_SuperFGD = 0.38;

// Approximate collection factors from PDG2016 (Detectors and accelerators section)
const double CollFactor_SingleClad = 0.06;
//const double CollFactor_DoubleClad = 0.054; // from Licciardi's thesis  
const double CollFactor_DoubleClad = 0.10;

const double Pedestal = 0;//145;  // pedeltal of ADC counts
const double Gain = 10;  // Gain ADC counts of high gain channel
const double LowGain  = 1;  // Gain ADC counts of low gain channel
const double ElecNoise = 1.7;  // sigma of high gain electronics noise
const double LowElecNoise = 1.2;  // sigma of low gain electronics noise
const double PixelGainVari = 0.031;  // gain variation among pixels

  double a=0.;        // long attenuation component fraction
  double d=0.;        // distance MPPC-scint outside the bar
  double LongAtt=0.;  // long attenuation length
  double ShortAtt=0.; // short attenuation length
  double Ldecay=0.;   // decay length
  double Lbar=0.;     // bar length
  
double hitLocation[3]={},hitPE[3]={},hitT[3]={};
  double adc_tmp[3]={},loadc_tmp[3]={},Q[3]={},loQ[3]={},adc[3]={};
  double loadc[3]={};
Int_t prim,PDG;
double ener;
double trueMom,trueCos,trueLen;

TFile* outFile = TFile::Open(Form("testEvent.root"),"RECREATE");
c = new TTree("EDepSimTree","EDepSimTree");
c->Branch("event",&eventN,"event/I");
c->Branch("hitLocation",&hitLocation,"hitLocation[3]/D");
c->Branch("hitPE",&hitPE,"hitPE[3]/D");
c->Branch("hitT",&hitT,"hitT[3]/D");
c->Branch("hitADC",&adc,"hitADC[3]/D");
c->Branch("hitLowADC",&loadc,"hitLowADC[3]/D");
c->Branch("hitQ",&Q,"hitQ[3]/D");
c->Branch("hitLowQ",&loQ,"hitLowQ[3]/D");
c->Branch("hitPrim",&prim,"hitPrim/I");
c->Branch("hitPDG",&PDG,"hitPDG/I");
c->Branch("hitE",&ener,"hitE/D");
c->Branch("trueMom",&trueMom,"trueMom/D");
c->Branch("trueCos",&trueCos,"trueCos/D");
c->Branch("trueLen",&trueLen,"trueLen/D");

TFile g("../3DST_event.root");
TTree* events = (TTree*) g.Get("EDepSimEvents");

//TGeoManager* geo = new TGeoManager;
//geo->Import("../Juan.gdml");

TG4Event* event=NULL;
events->SetBranchAddress("Event",&event);

Int_t nevent = 2000;//events->GetEntries();

for(Int_t ii=0;ii<nevent;ii++){

events->GetEntry(ii);
cout<<"event number "<<ii<<"----------------- number of prim. particle "<<event->Primaries[0].Particles.size()<<endl;

eventN = ii;

for(auto sd : event->SegmentDetectors)
{
for(Int_t i=0;i<sd.second.size();i++){
//cout<<"evt size "<<sd.second.size()<<" det. number "<<sd.first<<" Energy deposit "<<sd.second[i].EnergyDeposit<<" "<<sd.second[i].Start.Z()<<" "<<sd.second[i].Stop.Z()<<endl;

Int_t detNumber = stoi(sd.first);
//double xlocation = (Int_t)(detNumber/1000000+0.00001)*10;
//double ylocation = (Int_t)(detNumber%1000000-detNumber%100+0.0000001)/10000.*10;
//double zlocation = (Int_t)(detNumber%100)*10;
//geo->FindNode( (sd.second[i].Stop.X()+sd.second[i].Start.X())/2., (sd.second[i].Stop.Y()+sd.second[i].Start.Y())/2.,  (sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2. );
//double Master[3]={},Local[3]={0,0,0};
//geo->LocalToMaster(Local,Master);
//double xlocation = Master[0]+500;
//double ylocation = Master[1]+500;
//double zlocation = Master[2]+500;

//<position name="volPrimary_pos" unit="cm" x="0.0" y="-5500.0" z="59500.0"/>

int aveX= ((sd.second[i].Stop.X()+sd.second[i].Start.X())/2. +500)/10 ;
int aveY= ((sd.second[i].Stop.Y()+sd.second[i].Start.Y())/2. +500)/10 ;
int aveZ= ((sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2. +500)/10 ;

double xlocation = aveX*10. + 5 + 0.;
double ylocation = aveY*10. + 5 + 55000.;
double zlocation = aveZ*10. + 5 - 595000.;

cout<<sd.second[i].Start.X()<<" "<<sd.second[i].Start.Y()<<" "<<sd.second[i].Start.Z()<<endl;
cout<<detNumber<<" "<< (sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2.<<" "<<xlocation<<" "<<ylocation<<" "<<zlocation<<endl;

prim = sd.second[i].PrimaryId;
PDG = event->Primaries[0].Particles[prim].PDGCode; 

trueMom = event->Primaries[0].Particles[prim].Momentum.Energy();
trueLen = 0;
trueCos = event->Primaries[0].Particles[prim].Momentum.CosTheta();

hitLocation[0]=xlocation;
hitLocation[1]=ylocation;
hitLocation[2]=zlocation;

ener = sd.second[i].EnergyDeposit;

Double_t dedx = sd.second[i].EnergyDeposit/sd.second[i].TrackLength;
Double_t edep= sd.second[i].EnergyDeposit/(1. + CBIRKS*dedx);

if(dedx){
if(PDG == 13)  hMuon_Stopped_EdepOverLen->Fill(dedx);
if(PDG == 211) hPion_Stopped_EdepOverLen->Fill(dedx);
if(PDG == 2212)hProt_Stopped_EdepOverLen->Fill(dedx);
if(PDG == 11)  hElec_Stopped_EdepOverLen->Fill(dedx);
}

    // Account for the 3 fibers in the same scintillator cube
    double collfact = CollFactor_DoubleClad;
    double fact_fib1 = collfact;
    double fact_fib2 = (1-fact_fib1)*collfact;
    double fact_fib3 = (1-fact_fib2)*collfact;
    double CollFactAve = (fact_fib1+fact_fib2+fact_fib3)/3.;
    double NormShadowLight = CollFactAve / collfact; // fraction 
    //cout << "NormShadowLight = " << NormShadowLight << endl;   
    double Nphot = edep * EdepToPhotConv_FGD * NormShadowLight;
    cout<<"before dividing Nphot: "<<Nphot<<endl;

  a = LongCompFrac_FGD;
  d = DistMPPCscint_FGD;
  LongAtt = LongAtt_FGD;
  ShortAtt = ShortAtt_FGD;
  Ldecay= DecayLength_FGD;
  Lbar = Lbar_FGD;

  double xx = 1000 - xlocation;
  double yy = 1000 - ylocation;
  double zz = 1000 - zlocation;  
  //cout<<"xx "<<xx<<endl;
  double xNphot = Nphot * ( a*exp((-xx-d)/LongAtt) + (1-a)*exp((-xx-d)/ShortAtt) );
  double yNphot = Nphot * ( a*exp((-yy-d)/LongAtt) + (1-a)*exp((-yy-d)/ShortAtt) );
  double zNphot = Nphot * ( a*exp((-zz-d)/LongAtt) + (1-a)*exp((-zz-d)/ShortAtt) );
  //cout<<"a xx d longAtt ShortAtt xNphot "<<a<<" "<<xx<<" "<<d<<" "<<LongAtt<<" "<<ShortAtt<<" "<<xNphot<<" "<<exp((-xx-d)/LongAtt)<<endl;
  //cout<<"short att "<<exp((-xx-d)/ShortAtt)<<" "<<(1-a)*exp((-xx-d)/ShortAtt)<<" "<<a*exp((-zz-d)/LongAtt)<<endl;
  double xTimeDelay =  sd.second[i].Start.T()+TransTimeInFiber * xx;
  double yTimeDelay =  sd.second[i].Start.T()+TransTimeInFiber * yy;
  double zTimeDelay =  sd.second[i].Start.T()+TransTimeInFiber * zz;

  double xpe = xNphot * MPPCEff_SuperFGD;
  double ype = yNphot * MPPCEff_SuperFGD;
  double zpe = zNphot * MPPCEff_SuperFGD;

  //cout<<xNphot<<" "<<xpe<<endl;

  hitT[0]=xTimeDelay;
  hitT[1]=yTimeDelay;
  hitT[2]=zTimeDelay;
  hitPE[0]=xpe;
  hitPE[1]=ype;
  hitPE[2]=zpe;

  for(Int_t dim =0;dim<3;dim++){
  //PE to ADC
  adc_tmp[dim] = Pedestal + (hitPE[dim])*Gain;
  loadc_tmp[dim] = Pedestal + (hitPE[dim])*LowGain*14.29/13.55;

  //Electronics noise
  adc_tmp[dim] = gRandom->Gaus(adc_tmp[dim],ElecNoise);
  loadc_tmp[dim] = gRandom->Gaus(loadc_tmp[dim],LowElecNoise);

  //ADC to Charge
  //Q=(adc_tmp+53)/217;
  //loQ=(loadc_tmp+82)/26;
  Q[dim]=(adc_tmp[dim])/135.5;
  loQ[dim]=(loadc_tmp[dim])/14.29;

  //Non linearlity of high gain ADC
  if(Q[dim]<0.65) adc[dim]=135.5*Q[dim];
  else if(Q[dim]<3.2)  adc[dim]=217*Q[dim]-53;
  else if(Q[dim]<4.2)  adc[dim]=158.6*Q[dim]+133.9;
  else if(Q[dim]<14)  adc[dim]=5.1*Q[dim]+778.6;
  else  adc[dim]=850;

  //Non linearlity of low gain ADC
  if(loQ[dim]<7)  loadc[dim]=14.29*loQ[dim];
  else if(loQ[dim]<27)  loadc[dim]=26*loQ[dim]-82;
  else if(loQ[dim]<35.5)  loadc[dim]=21.12*loQ[dim]+48.24;
  else if(loQ[dim]<178.4)  loadc[dim]=0.7*loQ[dim]+775.1;
  else  loadc[dim]=900;
  }
//  cout<<"location, PE, Q, ADC along x: "<<hitLocation[0]<<" "<<hitPE[0]<<" "<<Q[0]<<" "<<adc[0]<<" "<<endl; 
//  cout<<"location, PE, Q, ADC along y: "<<hitLocation[1]<<" "<<hitPE[1]<<" "<<Q[1]<<" "<<adc[1]<<" "<<endl;
//  cout<<"location, PE, Q, ADC along z: "<<hitLocation[2]<<" "<<hitPE[2]<<" "<<Q[2]<<" "<<adc[2]<<" "<<endl;
/*
det=sd.first;
startX=sd.second[i].Start.X();
startY=sd.second[i].Start.Y();
startZ=sd.second[i].Start.Z();
startT=sd.second[i].Start.T();
stopX=sd.second[i].Stop.X();
stopY=sd.second[i].Stop.Y();
stopZ=sd.second[i].Stop.Z();
stopT=sd.second[i].Stop.T();
*/
c->Fill();
//events->Fill();
}
}
}
c->AddFriend(events);
outFile->Write();

for(Int_t i=0;i<hMuon_Stopped_EdepOverLen->GetNbinsX();i++){
hMuon_Stopped_EdepOverLen->SetBinContent(i+1,hMuon_Stopped_EdepOverLen->GetBinContent(i+1)/hMuon_Stopped_EdepOverLen->Integral());
hPion_Stopped_EdepOverLen->SetBinContent(i+1,hPion_Stopped_EdepOverLen->GetBinContent(i+1)/hPion_Stopped_EdepOverLen->Integral());
hProt_Stopped_EdepOverLen->SetBinContent(i+1,hProt_Stopped_EdepOverLen->GetBinContent(i+1)/hProt_Stopped_EdepOverLen->Integral());
hElec_Stopped_EdepOverLen->SetBinContent(i+1,hElec_Stopped_EdepOverLen->GetBinContent(i+1)/hElec_Stopped_EdepOverLen->Integral());
}
TFile* outFile2 = TFile::Open(Form("PID_PDF.root"),"RECREATE");
hMuon_Stopped_EdepOverLen->Write("Muon");
hProt_Stopped_EdepOverLen->Write("Prot");
hPion_Stopped_EdepOverLen->Write("Pion");
hElec_Stopped_EdepOverLen->Write("Elec");
}
