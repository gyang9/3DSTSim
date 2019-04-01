#include "TG4Event.h"

using namespace std;

class Transfer {
public:
  Transfer();
  void SetUp();
  void Transferring(TString inputFileName, int startingEventStackNumber);
}	

Transfer::Transfer()
{
  this->SetUp(TString inputFileName);
}

void Transfer::SetUp(TString inputFileName)
{
  
  // Basic detector geometry variables, the Master volume location does not seem to work well..
  double halfXsize = 1200; // mm
  double halfYsize = 1200; // mm
  double halfZsize = 1000; // mm
  
  double detXlocation = 0; // mm
  double detYlocation = 55000; // mm, 55 m underground
  double detZlocation = 594000; // mm, 594 m from beam target

  //constants for energy calibration
  // things changed since first version of 3DSTSim: EdepToPhotConv_FGD; 
  const double CBIRKS = 0.00208; // mm/MeV
  const double EdepToPhotConv_FGD = 156.42; // CLHEP::MeV; // contains both collection in fiber and edep->gamma conversion 
  const double DistMPPCscint_FGD = 41; //*CLHEP::mm;
  const double LongCompFrac_FGD = 0.816;
  const double LongAtt_FGD = 11926.; //*CLHEP::mm;
  const double ShortAtt_FGD = 312.; //*CLHEP::mm;
  const double DecayLength_FGD = 0.0858; // CLHEP::mm;
  const double Lbar_FGD = 1864.3; //* CLHEP::mm;
  const double TransTimeInFiber = 1./280.; // speed in fiber: 280 mm/ns
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
  
  double hitLocation[3]={},hitPE[6]={},hitT[6]={};
  double adc_tmp[6]={},loadc_tmp[6]={},Q[6]={},loQ[6]={},adc[6]={};
  double loadc[6]={};
  Int_t prim,PDG;
  double ener;
  double trueMom,trueCos,trueLen;
  double true3Mom[3]={};
  Int_t eventN;
  string det;
  Int_t iii=0;

  TFile* outFile = TFile::Open(Form("output_%s.root", inputFileName),"RECREATE");
  c = new TTree("EDepSimTree","EDepSimTree");
  c->Branch("event",&eventN,"event/I");
  c->Branch("hitLocation",&hitLocation,"hitLocation[3]/D");
  c->Branch("hitPE",&hitPE,"hitPE[6]/D");
  c->Branch("hitT",&hitT,"hitT[6]/D");
  c->Branch("hitADC",&adc,"hitADC[6]/D");
  c->Branch("hitLowADC",&loadc,"hitLowADC[3]/D");
  c->Branch("hitQ",&Q,"hitQ[3]/D");
  c->Branch("hitLowQ",&loQ,"hitLowQ[3]/D");
  c->Branch("hitPrim",&prim,"hitPrim/I");
  c->Branch("hitPDG",&PDG,"hitPDG/I");
  c->Branch("hitE",&ener,"hitE/D");
  c->Branch("trueMom",&trueMom,"trueMom/D");
  c->Branch("trueCos",&trueCos,"trueCos/D");
  c->Branch("trueLen",&trueLen,"trueLen/D");
  c->Branch("true3Mom",&true3Mom,"true3Mom[3]/D");

}

void Transfer::Transferring(TString inputFileName, int startingEventStackNumber)
{
  // group 200 events for one stack
  TFile g(Form("%s.root", inputFileName));
  TTree* events = (TTree*) g.Get("EDepSimEvents");

  TG4Event* event=NULL;
  events->SetBranchAddress("Event",&event);

  Int_t nevent = 200;//events->GetEntries();

  for(Int_t ii=nevent*startingEventStackNumber;ii<nevent*startingEventStackNumber+1);ii++)
  {
    events->GetEntry(ii);
    cout<<"event number "<<ii<<"----------------- number of prim. particle "<<event->Primaries[0].Particles.size()<<endl;

    eventN = ii;

    for(auto sd : event->SegmentDetectors)
    {
	for(Int_t i=0;i<sd.second.size();i++)
	{

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

		if(sd.second[i].TrackLength>0)
		{

		int aveX= ((sd.second[i].Stop.X()+sd.second[i].Start.X())/2. +halfXsize)/10 ;
		int aveY= ((sd.second[i].Stop.Y()+sd.second[i].Start.Y())/2. +halfYsize)/10 ;
		int aveZ= ((sd.second[i].Stop.Z()+sd.second[i].Start.Z())/2. +halfZsize)/10 ;

		// change the global coordinate to local coordinate
		double xlocation = aveX*10. + 5 + detXlocation;
		double ylocation = aveY*10. + 5 + detYlocation;
		double zlocation = aveZ*10. + 5 - detZlocation;

		prim = sd.second[i].PrimaryId;
		PDG = event->Primaries[0].Particles[prim].PDGCode; 

		trueMom = event->Primaries[0].Particles[prim].Momentum.Energy();
		trueLen = 0;
		trueCos = event->Primaries[0].Particles[prim].Momentum.CosTheta();
		true3Mom[0]=event->Primaries[0].Particles[prim].Momentum.Px();
		true3Mom[1]=event->Primaries[0].Particles[prim].Momentum.Py();
		true3Mom[2]=event->Primaries[0].Particles[prim].Momentum.Pz();

		hitLocation[0]=xlocation;
		hitLocation[1]=ylocation;
		hitLocation[2]=zlocation;

		ener = sd.second[i].EnergyDeposit;

		Double_t dedx = sd.second[i].EnergyDeposit/sd.second[i].TrackLength;
		Double_t edep= sd.second[i].EnergyDeposit/(1. + CBIRKS*dedx);

		// Account for the 3 fibers in the same scintillator cube
		double collfact = CollFactor_DoubleClad;
		double fact_fib1 = collfact;
		double fact_fib2 = (1-fact_fib1)*collfact;
		double fact_fib3 = (1-fact_fib2)*collfact;
		double CollFactAve = (fact_fib1+fact_fib2+fact_fib3)/3.;
		double NormShadowLight = CollFactAve / collfact; // fraction 
		//cout << "NormShadowLight = " << NormShadowLight << endl;   
		double Nphot = edep * EdepToPhotConv_FGD * NormShadowLight;

		a = LongCompFrac_FGD;
		d = DistMPPCscint_FGD;
		LongAtt = LongAtt_FGD;
		ShortAtt = ShortAtt_FGD;
		Ldecay= DecayLength_FGD;
		Lbar = Lbar_FGD;

		double xx = 2 * halfXsize - xlocation;
		double yy = 2 * halfYsize - ylocation;
		double zz = 2 * halfZsize - zlocation;  

		double NphotXY = Nphot * ( a*exp((-zz-d)/LongAtt) + (1-a)*exp((-zz-d)/ShortAtt) );
		double NphotXZ = Nphot * ( a*exp((-yy-d)/LongAtt) + (1-a)*exp((-yy-d)/ShortAtt) );
		double NphotYZ = Nphot * ( a*exp((-xx-d)/LongAtt) + (1-a)*exp((-x-d)/ShortAtt) );

		double TimeDelayXY =  sd.second[i].Start.T()+TransTimeInFiber * zz;
		double TimeDelayXZ =  sd.second[i].Start.T()+TransTimeInFiber * yy;
		double TimeDelayYZ =  sd.second[i].Start.T()+TransTimeInFiber * xx;

		double peXY = NphotXY * MPPCEff_SuperFGD;
		double peXY = NphotXY * MPPCEff_SuperFGD;
		double peXZ = NphotXZ * MPPCEff_SuperFGD;

		hitT[0]=TimeDelayXY;
		hitT[1]=TimeDelayXZ;
		hitT[2]=TimeDelayYZ;
		hitPE[0]=peXY;
		hitPE[1]=peXZ;
		hitPE[2]=peYZ;

		for(Int_t dim =0;dim<3;dim++){
		  //PE to ADC
		  adc_tmp[dim] = Pedestal + (hitPE[dim])*Gain;
		  loadc_tmp[dim] = Pedestal + (hitPE[dim])*LowGain*14.29/13.55;

		  //Electronics noise
		  adc_tmp[dim] = gRandom->Gaus(adc_tmp[dim],ElecNoise);
		  loadc_tmp[dim] = gRandom->Gaus(loadc_tmp[dim],LowElecNoise);

		  //ADC to Charge
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
	c->Fill();

        }
      }
    }
  }
  c->AddFriend(events);
  outFile->Write();
}

void main(int argc, char* argv[])
{
  Transfer* module;
  module->Transferring("file",atoi(gApplication->Argv(1)));
}
