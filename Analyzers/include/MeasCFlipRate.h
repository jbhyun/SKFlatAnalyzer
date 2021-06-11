#ifndef MeasCFlipRate_h
#define MeasCFlipRate_h

#include "AnalyzerCore.h"

class MeasCFlipRate : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();

  bool SS2l, TriLep, TetraLep, OS2l, CFlip, CFMCClos, MCCFRate, MDists; 
  bool DblMu, DblEG, MuEG;
  bool FakeRun, ConvRun, FlipRun, SystRun;
  vector<TString> TrigList_DblMu, TrigList_DblEG, TrigList_MuEG;



  void GetMassDists(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                    vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                    vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void CheckMCCFClosure(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                        vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                        vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  float GetTmpCFRate(float PT, float fEta);
  void CheckChargeFlip(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                       vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                       vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void MeasMCCFRate(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                    vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                    vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);


  float CalcTestFakeWeight(vector<Muon>& MuColl, vector<Electron>& ElColl, TString MuIDT, TString MuIDL, TString ElIDT, TString ElIDL, int NBJet=0, int SystDir=0);
  float GetTestMuFR(Muon& Mu, TString Key, int SystDir=0);
  float GetTestElFR(Electron& El, TString Key, int SystDir=0);
  float GetCFRSF(Electron& El, TString Tag, TString Option="");
  float GetFlipCorrPT(Electron& El, TString Tag, TString Option="");


  MeasCFlipRate();
  ~MeasCFlipRate();

};



#endif
