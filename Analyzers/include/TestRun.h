#ifndef TestRun_h
#define TestRun_h

#include "AnalyzerCore.h"

class TestRun : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool SS2l, TriLep, TetraLep; 
  bool DblMu, DblEG, MuEG;
  bool SystRun;
  vector<TString> TrigList_DblMu, TrigList_DblEG, TrigList_MuEG;

  void TestThis(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);

  TestRun();
  ~TestRun();

};



#endif
