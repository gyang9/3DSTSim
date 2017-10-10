#include <stdio.h>
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TObject.h"
#include "TString.h"
#include "TMath.h"   
#include "TRandom.h"   
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TVector3.h"
#include "/home/gyang/work/edep-sim/src/TG4Event.h"

using namespace std;

class EventHeader { 
    public:
   	Int_t fEvtNum;
  	Int_t fRun;
   	Int_t fDate;
   	EventHeader() : fEvtNum(0), fRun(0), fDate(0) { };
   	//virtual ~EventHeader() { }
	~EventHeader() { };
   	void   Set(Int_t i, Int_t r, Int_t d) { fEvtNum = i; fRun = r; fDate = d; };
   	Int_t  GetEvtNum() const { return fEvtNum; };
   	Int_t  GetRun() const { return fRun; };
   	Int_t  GetDate() const { return fDate; };
   	void SetHeader(Int_t header1, Int_t header2, Int_t header3){fEvtNum = header1; fRun = header2; fDate = header3;} 
   	//ClassDef(EventHeader,1)  // Event Header
};

class Event : public TObject {
   public:
      char                 fType[20];
      Int_t                fNhit;
      Int_t                fNseg;
      Int_t                fip;
      std::string          fDet;
      UInt_t               fFlag;
      EventHeader          fEvtHdr;
      //TClonesArray         *fHits;            //->
      // ... list of methods
      //ClassDef(Event,1)  //Event structure
      //virtual ~Event(){};
};

class Hits : public TObject {
   public:
      Float_t   hitX;         //X component of the momentum
      Float_t   hitY;         //Y component of the momentum
      Float_t   hitZ;         //Z component of the momentum
      Float_t   hitT;     //A random track quantity
      Float_t   hitQ;      //The mass square of this particle
      TClonesArray         *fHits;
      // method definitions ...
      ClassDef(Hits,1);          //A track segment
      virtual ~Hits(){};
};

class HitAna {
        public:
        HitAna(std::string inputFileName, std::string outputFileName);
        int processFileName(std::string inputFileName, std::string &baseFileName);
        void doHitAnalysis(std::string inputFileName,std::string outputFileName);
        void parseFile();

	Int_t fNhit = 0;
        TClonesArray  *fHits;

        //Files
        ifstream infile;
        TFile *gOut;
	TTree *inTree;

        std::string inputFileName ;
	std::string outputFileName ;
        std::string fileName ;

        Hits *AddHit()
        {
        // Add a new hit to the list of hits in detector 
        TClonesArray &hits = *fHits;
        Hits *hit = new(hits[fNhit++]) Hits();
        return hit;
        }
};

int HitAna::processFileName(std::string inputFileName, std::string &baseFileName){
        //check if filename is empty
        if( inputFileName.size() == 0 ){
                std::cout << "processFileName : Invalid filename " << std::endl;
                return 0;
        }

        //remove path from name
        size_t pos = 0;
        std::string delimiter = "/";
        while ((pos = inputFileName.find(delimiter)) != std::string::npos)
                inputFileName.erase(0, pos + delimiter.length());

        if( inputFileName.size() == 0 ){
                std::cout << "processFileName : Invalid filename " << std::endl;
                return 0;
        }

        //replace / with _
        std::replace( inputFileName.begin(), inputFileName.end(), '/', '_'); // replace all 'x' to 'y'
        std::replace( inputFileName.begin(), inputFileName.end(), '-', '_'); // replace all 'x' to 'y'

        baseFileName = inputFileName;

        return 1;
}


int main(int argc, char* argv[])
{

  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  std::cout << "inputFileName " << inputFileName << std::endl;
  HitAna ana(inputFileName,outputFileName);
  ana.doHitAnalysis(inputFileName,outputFileName);

}

HitAna::HitAna(std::string inputFileName, std::string outputFileName){

        //get input file
        if( inputFileName.empty() ){
                std::cout << "Error invalid file name" << std::endl;
                exit(0);
        }

        infile.open(inputFileName);  //, std::ifstream::in | std::ifstream::binary);  
        if (infile.fail()) { 
                std::cout << "Error opening input file, exiting" << std::endl;
                exit(0);
        }


        if( processFileName( inputFileName, outputFileName ) )
                outputFileName = "output_processNtuple_" + outputFileName;
        else
                outputFileName = "output_processNtuple.root";

        gOut = new TFile(outputFileName.c_str() , "RECREATE");

}

void HitAna::doHitAnalysis(std::string inputFileName, std::string outputFileName){

   //Define the input tree
   TFile fint(inputFileName.c_str());
   inTree = (TTree*)fint.Get("EDepSimEvents"); 

   TFile fout(outputFileName.c_str(),"RECREATE");
   // create a ROOT Tree
   TTree t4("HitTree","A Tree with Events");
   // create a pointer to an Event object
   Event *event = new Event();
   //TG4Event* G4event=NULL;
   TG4Event* G4event;
   inTree->SetBranchAddress("Event",&G4event);
   Hits *hit;
   EventHeader *evtHeader;
   // create two branches, split one
   t4.Branch("event_branch", "Event", &event,16000,2);
   t4.Branch("event_not_split", "Event", &hit,16000,0);

   // a local variable for the event type
   char etype[20];
   Int_t nevent = inTree->GetEntries();

   // do the analysis and fill the tree
   for (Int_t ev = 0; ev <nevent; ev++) {

      	hit = AddHit();
      	inTree->GetEntry(ev);
	for(auto sd : G4event->SegmentDetectors){
		for(Int_t iii=0;iii<sd.second.size();iii++){
			cout<<"Event "<<ev<<" detector "<<sd.first<<" starting point "<<sd.second[iii].Start.Z()<<" Energy "<<sd.second[iii].EnergyDeposit<<endl;

      evtHeader->SetHeader(ev, 200, 960312);
      event->fip  = 1;
      event->fDet = sd.first; 
      hit->hitX = sd.second[iii].Start.X();
      hit->hitY = sd.second[iii].Start.Y();
      hit->hitZ = sd.second[iii].Start.Z();
      hit->hitT = 0;
      hit->hitQ = sd.second[iii].EnergyDeposit;
      //event->fHits->Add(hit);
      }
	      }

   t4.Fill();
   event->Clear();
   }
   fout.Write();
   t4.Print();
   
   //Add the true edep-sim results as friend tree
   t4.AddFriend(inTree,"friendTree");

}


