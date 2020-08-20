#ifndef GenMatchingValid_h
#define GenMatchingValid_h

#include "AnalyzerCore.h"

class GenMatchingValid : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool PhotonMatch, LepTypeValid, HNTypeValid;
  bool DblMu, DblEG, MuEG;
  bool SystRun;
  vector<TString> TrigList;
  int NPrintMax, NCountPrint;

  void CheckLepTypeValidity(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                            vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, vector<Gen>& TruthColl, float weight, TString Label);
  void CheckPhotonMatching(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                           vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, vector<Gen>& TruthColl, float weight, TString Label);


  GenMatchingValid();
  ~GenMatchingValid();

};



#endif

