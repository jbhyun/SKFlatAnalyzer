#ifndef FakeRateMeas_h
#define FakeRateMeas_h

#include "AnalyzerCore.h"

class FakeRateMeas : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool ElFR, MuFR, MeasFR, MeasPU, PrVal, METMTWCut, FkVal;
  bool SS2l;
  bool SystRun;
  vector<TString> TrigList_ElFR, TrigList_MuFR;

  void MeasFakeRate(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                    vector<Electron>& EleTColl, vector<Electron>& EleLColl, vector<Electron>& EleVColl,
                    vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label);
  void MeasPileUp(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                  vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                  vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label);
  void CheckPromptValid(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                  vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                  vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label);
  void OptimizeMETMTWCuts(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                          vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                          vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label);

  float GetResidualNormSF(TString Key);
  float NvtxReweight(TString Key);
  bool  IsNearBJet(Lepton& Lepton, vector<Jet>& BJetColl);

  FakeRateMeas();
  ~FakeRateMeas();

};


#endif
