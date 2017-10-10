#include <stdio.h>
#include "TFile.h"
#include "TH2D.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "TMath.h"
#include "TObject.h"
#include "TString.h"
#include "TMath.h"   
#include "TRandom.h"   
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TVector3.h"
 
using namespace std;


class Event : public TObject {
   private:
      char                 fType[20];
      Int_t                fNhit;
      Int_t                fNseg;
      Int_t                fip;
      char                 fDet[50];
      UInt_t               fFlag;
      EventHeader          fEvtHdr;
      //TClonesArray         *fHits;            //->
      // ... list of methods
      ClassDef(Event,1)  //Event structure
};

class Hits : public TObject {
   private:
      Float_t   hitX;         //X component of the momentum
      Float_t   hitY;         //Y component of the momentum
      Float_t   hitZ;         //Z component of the momentum
      Float_t   hitT;     //A random track quantity
      Float_t   hitQ;      //The mass square of this particle
      TClonesArray         *fHits;
      // method definitions ...
      ClassDef(Track,1)          //A track segment
};


int Analyze::processFileName(std::string inputFileName, std::string &baseFileName){
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

   if (!TClassTable::GetDict("Event")) {
      gSystem->Load("$ROOTSYS/test/libEvent.so");
   }

  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];
  std::cout << "inputFileName " << inputFileName << std::endl;
  HitAna ana(inputFileName,outputFileName);


}

int Analyze::processFileName(std::string inputFileName, std::string &baseFileName){
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

HitAna::HitAna(std::string inputFileName, std::string outputFileName){

        //get input file
        if( inputFileName.empty() ){
                std::cout << "Error invalid file name" << std::endl;
                gSystem->Exit(0);
        }

        infile.open(inputFileName);  //, std::ifstream::in | std::ifstream::binary);  
        if (infile.fail()) { 
                std::cout << "Error opening input file, exiting" << std::endl;
                gSystem->Exit(0);
        }


        if( processFileName( inputFileName, outputFileName ) )
                outputFileName = "output_processNtuple_" + outputFileName;
        else
                outputFileName = "output_processNtuple.root";

        gOut = new TFile(outputFileName.c_str() , "RECREATE");

}

void HitAna::doHitAnalysis(){

   //Define the input tree
   TFile fint(inputFileName);
   TTree inTree = (TTree*)fint->Get("EDepSimEvents"); 


   TFile fout(outputFileName,"RECREATE");
   // create a ROOT Tree
   TTree t4("HitTree","A Tree with Events");
   // create a pointer to an Event object
   Event *event = new Event();
   // create two branches, split one
   t4.Branch("event_branch", "Event", &event,16000,2);
   t4.Branch("event_not_split", "Event", &event,16000,0);

   // a local variable for the event type
   char etype[20];
   Hit *hit;

   Int_t nevent = inTree->GetEntries()

   // do the analysis and fill the tree
   for (Int_t ev = 0; ev <nevent; ev++) {
      sprintf(etype,"type%d",ev%5);
      event->SetType(etype);
      hit->AddHit();



	events->GetEntry(ii);

	for(auto sd : event->SegmentDetectors){
		for(Int_t iii=0;iii<sd.second.size();iii++){
			cout<<"Event "<<ii<<" detector "<<sd.first<<" starting point "<<sd.second[iii].Start.Z()<<" Energy "<<sd.second[iii].EnergyDeposit<<endl;

      event->SetHeader(ev, 200, 960312, ev + 0.1);
      event->fip  = 1;
      event->fDet = sd.first; 
      hit->hitX = sd.second[iii].Start.X();
      hit->hitY = sd.second[iii].Start.Y();
      hit->hitZ = sd.second[iii].Start.Z();
      hit->hitT = 0;
      hit->hitQ = sd.second[iii].EnergyDeposit;
      event->fHits->Add(hit);
      }
	      }

   t4.Fill();
   event->Clear();
   }
   fout.Write();
   t4.Print();
   
   //Add the true edep-sim results as friend tree
   t4->AddFriend(inTree);

}








  Hits *AddHit()
  {
    // Add a new hit to the list of hits in detector 
    TClonesArray &hits = *fHits;
    Hit *hit = new(hits[fNhit++]) Hit();
    return hit;
  }








	



	




