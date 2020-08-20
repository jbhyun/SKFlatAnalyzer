#ifndef TestRun_h
#define TestRun_h

#include "AnalyzerCore.h"

class TestRun : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MuMuMu, ElMuMu; 
  bool SystRun;
  vector<TString> TrigList_1e2mu, TrigList_3mu;

  void AnalyzeTriMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                     std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeElectronDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                             std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);



  TestRun();
  ~TestRun();

};



#endif

