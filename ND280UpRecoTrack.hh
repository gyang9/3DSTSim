#ifndef ND280UP_RECOTRACK
#define ND280UP_RECOTRACK

#include "ND280UpConst.hh"

#include <TH2F.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TGraphErrors.h>
#include <TF1.h>


using namespace std;

//
// Davide Sgalaberna 15/7/17
//
// Container of methods to reconstruct single tracks
//

class ND280UpRecoTrack {
public:
  
  ND280UpRecoTrack();
  ~ND280UpRecoTrack();

  // Event display of MPPC hits with single track
  TH2F *f2d_xy;
  TH2F *f2d_xz;
  TH2F *f2d_yz;
  TH2F *f2d_xy_other;
  TH2F *f2d_xz_other;
  TH2F *f2d_yz_other;
  //void SetMPPCXY(TH2F *h){f2d_xy = new TH2F(*h);}
  //void SetMPPCXZ(TH2F *h){f2d_xz = new TH2F(*h);} 
  //void SetMPPCYZ(TH2F *h){f2d_yz = new TH2F(*h);}
  //void SetMPPCXY(TH2F *h){h->Copy(*f2d_xy);}
  //void SetMPPCXZ(TH2F *h){h->Copy(*f2d_xz);} 
  //void SetMPPCYZ(TH2F *h){h->Copy(*f2d_yz);}  
  //void SetMPPCXY(TH2F *h){f2d_xy=(TH2F*)h->Clone("f2d_xy");} // TODO: Use Clone
  //void SetMPPCXZ(TH2F *h){f2d_xz=(TH2F*)h->Clone("f2d_xz");} // TODO: Use Clone
  //void SetMPPCYZ(TH2F *h){f2d_yz=(TH2F*)h->Clone("f2d_yz");} // TODO: Use Clone
  void SetMPPCXY(TH2F *h){f2d_xy=h;} // TODO: Use Clone
  void SetMPPCXZ(TH2F *h){f2d_xz=h;} // TODO: Use Clone
  void SetMPPCYZ(TH2F *h){f2d_yz=h;} // TODO: Use Clone
  
  void SetMPPCXY_Other(TH2F *h){f2d_xy_other=h;} // TODO: Use Clone
  void SetMPPCXZ_Other(TH2F *h){f2d_xz_other=h;} // TODO: Use Clone
  void SetMPPCYZ_Other(TH2F *h){f2d_yz_other=h;} // TODO: Use Clone

  //void SetMPPCXY(TH2F *h); // TODO: Use Clone
  //void SetMPPCXZ(TH2F *h); // TODO: Use Clone
  //void SetMPPCYZ(TH2F *h); // TODO: Use Clone

  //void SetMPPCXY_Other(TH2F *h); // TODO: Use Clone
  //void SetMPPCXZ_Other(TH2F *h); // TODO: Use Clone
  //void SetMPPCYZ_Other(TH2F *h); // TODO: Use Clone


  void DoTracking(nd280upconv::TargetType_t dettype);

  nd280upconv::TargetType_t GetDetectorType(){return fDetType;}

  double GetRecoCosTheta(){return fRecoCosTheta;}
  double GetRecoMom(){return fRecoMom;}
  
  double GetLength(){return fTrkLength;}
  double GetLengthX(){return fTrkLengthX;}
  double GetLengthY(){return fTrkLengthY;}
  double GetLengthZ(){return fTrkLengthZ;}

  double GetLengthStraight(){return fTrkLengthStraight;}
  double GetLengthStraightX(){return fTrkLengthStraightX;}
  double GetLengthStraightY(){return fTrkLengthStraightY;}
  double GetLengthStraightZ(){return fTrkLengthStraightZ;}

  double GetLengthFit(){return fTrkLengthFit;}
  double GetLengthFitX(){return fTrkLengthFitX;}
  double GetLengthFitY(){return fTrkLengthFitY;}
  double GetLengthFitZ(){return fTrkLengthFitZ;}
  
  double GetEdep(){return fTrkEdep;}

  void SetMinPE(double pe){fMinPE = pe;}
  double GetMinPE(){return fMinPE;}

  bool IsViewXY(){ return fUseView_xy;}
  bool IsViewXZ(){ return fUseView_xz;}
  bool IsViewYZ(){ return fUseView_yz;}
  int GetNumOfViews(){return fNumView;}

  int GetNHits(){return fNHits;}
  int GetNHitsXY(){return fNHitsXY;}
  int GetNHitsYZ(){return fNHitsYZ;}
  int GetNHitsXZ(){return fNHitsXZ;}

  bool IsOutFV(){return fIsOutFV;}
  bool IsReco(){return fIsReco;}

  TVector3 GetVertexReco(){return fVtxReco;}
   
  // Track separation

  void DoTrackSeparation(); // check if the loaded track is separated another track 
  
  void SetTrackSeparationMin(double len){fTrackSeparationMin = len;}
  double GetTrackSeparationMin(){return fTrackSeparationMin;}

  bool IsSeparated(){return fIsSeparated;}

  bool CalcOutFVTrue(TH2F *h2dmppchit,TGraph *gmchit); // based on the MC hits


  void Init();
  void Delete();

  // Methods for tracking
  void CleanMPPCHits();
  void CleanMPPCHits(TH2F *h2d);
  void CalcLength();
  void CalcLengthFit();
  void CalcLengthStraight(); // also gets the first and last MPPC hit
  double CalcLengthStraight1D(TH2F *h2D,bool axisX); // axisX=true --> X, axisX=false --> Y
  void CalcNHits();
  void CalcEdep();
  void CalcOutFV(); // based on the MPPC hits
  bool CalcOutFVView(TH2F *h2D);
  void DefineRecoTrack(); // cuts 
  void CalcRecoVar(); 

  void CalcVtxReco();

  //

  void SetIsReco(bool isreco){fIsReco = isreco;}

  void SetLength(double len){fTrkLength=len;}
  void SetLengthX(double len){fTrkLengthX=len;}
  void SetLengthY(double len){fTrkLengthY=len;}
  void SetLengthZ(double len){fTrkLengthZ=len;}

  double CalcPtDistance2D(double *pt1,double *pt2);
  double CalcPtDistance3D(double *pt1,double *pt2);
  
  // Track separation
  void FindTrackExtremes(TH2F *h2d,double *extreme,bool SelProjX);
  void FindTrackExtremes(TH2F *h2d,double *Pt1,double *Pt2,bool IsViewUsed);
  bool GetTrackSeparationInView(double *posdif_cur,double *posdif_oth);

  void GetPositionDiffView(double *Pt1_cur,double *Pt2_cur,double *Pt1_oth,double *Pt2_oth,double *posdif_cur,double *posdif_oth);

  double fTrkLengthStraight; // Track length (farthest hits)
  double fTrkLengthStraightX; // Track length (farthest hits)
  double fTrkLengthStraightY; // Track length (farthest hits)
  double fTrkLengthStraightZ; // Track length (farthest hits)

  double fTrkLengthFit; // from fit of hits
  double fTrkLengthFitX; // from fit of hits
  double fTrkLengthFitY; // from fit of hits
  double fTrkLengthFitZ; // from fit of hits

  double fTrkLength; // from fit of hits
  double fTrkLengthX; // from fit of hits
  double fTrkLengthY; // from fit of hits
  double fTrkLengthZ; // from fit of hits

  double fTrkEdep; // total energy deposited in single track MPPC hits (all views)

  bool fUseView_xy;
  bool fUseView_xz;
  bool fUseView_yz;
  int fNumView;
  bool fIsOutFV;
  bool fIsReco;

  int fNHitsXY;
  int fNHitsYZ;
  int fNHitsXZ;
  int fNHits;

  double fMinPE;

  double fRecoCosTheta;
  double fRecoMom;

  nd280upconv::TargetType_t fDetType;  

  // Track separation

  double fTrackSeparationMin;
  
  bool fIsSeparated;

  double fPt1_XY[2]; // first  track extreme point in XY
  double fPt2_XY[2]; // second track extreme point in XY
  double fPt1_XZ[2]; // first  track extreme point in XZ
  double fPt2_XZ[2]; // second track extreme point in XZ
  double fPt1_YZ[2]; // first  track extreme point in YZ
  double fPt2_YZ[2]; // second track extreme point in YZ

  // Vertex
  TVector3 fVtxReco;

  private:

};

class ND280UpPID {
public:
  ND280UpPID();
  ~ND280UpPID();

  void SetPDF(TString hyp,TH1F *hpdf);
  TH1F *GetPDF(TString hyp); 
  TGraph *GetGraph(TString hyp); 
  TSpline3 *GetSpline(TString hyp);
  
  double GetIntegral(TString hyp);

  void CalcIntegral(TString hyp,const int npoints=100000);

  double CalcLogLikeRatio(double prob1,double prob2);
  double CalcLogLikeRatio(TString hyp1,TString hyp2,double input);
  

  void Init();
  void Delete();

  double GetIntegral_Muon(){return fIntegral_Muon;}
  double GetIntegral_Pion(){return fIntegral_Pion;}
  double GetIntegral_Prot(){return fIntegral_Prot;}
  double GetIntegral_Elec(){return fIntegral_Elec;}

  void SetPDF_Muon(TH1F *hpdf);
  void SetPDF_Pion(TH1F *hpdf);
  void SetPDF_Prot(TH1F *hpdf);
  void SetPDF_Elec(TH1F *hpdf);

  TH1F *GetPDF_Muon(){return fPDF_Muon;}
  TH1F *GetPDF_Pion(){return fPDF_Pion;}
  TH1F *GetPDF_Prot(){return fPDF_Prot;}
  TH1F *GetPDF_Elec(){return fPDF_Elec;}

  TGraph *GetGraph_Muon(){return fGraph_Muon;}
  TGraph *GetGraph_Pion(){return fGraph_Pion;}
  TGraph *GetGraph_Prot(){return fGraph_Prot;}
  TGraph *GetGraph_Elec(){return fGraph_Elec;}

  TSpline3 *GetSpline_Muon(){return fSpline_Muon;}
  TSpline3 *GetSpline_Pion(){return fSpline_Pion;}
  TSpline3 *GetSpline_Prot(){return fSpline_Prot;}
  TSpline3 *GetSpline_Elec(){return fSpline_Elec;}

  TH1F *fPDF_Muon;
  TH1F *fPDF_Pion;
  TH1F *fPDF_Prot;
  TH1F *fPDF_Elec;

  TGraph *fGraph_Muon;
  TGraph *fGraph_Pion;
  TGraph *fGraph_Prot;
  TGraph *fGraph_Elec;
  
  TSpline3 *fSpline_Muon;
  TSpline3 *fSpline_Pion;
  TSpline3 *fSpline_Prot;
  TSpline3 *fSpline_Elec;

  double fIntegral_Muon;
  double fIntegral_Pion;
  double fIntegral_Prot;
  double fIntegral_Elec;

  private:
};

#endif
