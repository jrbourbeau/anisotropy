
#ifndef __SimpleDST_h__
#define __SimpleDST_h__

#include <Rtypes.h>

class TChain;
class TBranch;

class SimpleDST {

 public:

  SimpleDST() { }
  SimpleDST(TChain* chain, std::string config) { SetupChain(chain, config); }
  ~SimpleDST() { }

// Joint parameters
  Double_t ModJulDay;
  // IceCube parameters
  UShort_t NChannels;
  Bool_t IsGoodLineFit;
  Float_t LLHAzimuthDeg;
  Float_t LLHZenithDeg;
  // IceTop parameters
  Double_t Azimuth;
  Double_t Zenith;
  Int_t nStations;
  Bool_t isReco;
  Bool_t isSTA3;
  Bool_t isSTA3ii;
  Bool_t isSTA5;
  Bool_t isSTA8;
  // ShowerLLH parameters
  Double_t llhEnergy;
  Double_t pLLH;
  Double_t hLLH;
  Double_t oLLH;
  Double_t fLLH;
  Text_t comp;
  Double_t s125;
  Double_t ss125;


  // Joint parameters
  TBranch* b_ModJulDay;
  //IceCube parameters
  TBranch* b_nchan;
  TBranch* b_isGoodLineFit;
  TBranch* b_LLHAzimuth;
  TBranch* b_LLHZenith;
  //IceTop parameters
  TBranch* b_Azimuth;
  TBranch* b_Zenith;
  TBranch* b_nStations;
  TBranch* b_isReco;
  TBranch* b_isSTA3;
  TBranch* b_isSTA3ii;
  TBranch* b_isSTA5;
  TBranch* b_isSTA8;
  // ShowerLLH parameters
  TBranch* b_llhEnergy;
  TBranch* b_pLLH;
  TBranch* b_hLLH;
  TBranch* b_oLLH;
  TBranch* b_fLLH;
  TBranch* b_comp;
  TBranch* b_s125;
  TBranch* b_ss125;

  void SetupChain(TChain* chain, std::string config);

};

#endif // __SimpleDST_h__

