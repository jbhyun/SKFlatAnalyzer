#ifndef IDOptimization_h
#define IDOptimization_h

#include "AnalyzerCore.h"

class IDOptimization : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MuID, ElID, ZData; 
  bool ElTrigCut, ElTrigEffect, MuTrigEffect, TrigEffCheck, SelEffCheck, FakeRateCheck;
  bool DblMu, DblEG, MuEG;
  bool MuMu, ElEl;
  bool SystRun;
  vector<TString> TrigList_DblMu, TrigList_DblMu1, TrigList_DblMu2, TrigList_DblEG;

  void CheckIDVar_MC(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                     vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void CheckIDVar_ZData(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                        vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void CheckTrigCut_MC(vector<Muon>& MuLColl, vector<Muon>& MuRawColl, vector<Electron>& ElLColl, vector<Electron>& ElRawColl,
                       vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void CheckTrigEffect(vector<Muon>& MuRawColl, vector<Electron>& ElRawColl,
                       vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void CheckTrigEfficiency(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                           vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void CheckSelectionEfficiency(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuRawColl,
                                vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElRawColl,
                                vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void CheckTrigEffonTightLooseRatio(vector<Electron>& ElRawColl, TString Str_ElLID, TString Str_ElTID,
                                     vector<Gen>& TruthColl, Event& ev, float weight, TString Label);

  void TEST(vector<Muon>& MuRawColl, vector<Electron>& ElRawColl,
            vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);





  IDOptimization();
  ~IDOptimization();

};



#endif

