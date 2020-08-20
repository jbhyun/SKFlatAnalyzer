#ifndef TrigCheck_h
#define TrigCheck_h

#include "AnalyzerCore.h"

class TrigCheck : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MFilter, QFilter; 
  bool DblMu, DblEG, MuEG;
  bool SystRun;
  vector<TString> TrigList, TrigList_DblMu, TrigList_DblEG, TrigList_MuEG;

  void CheckMFilterBias(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                        std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event ev, float weight, TString Label);
  void CheckQFilterBias(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                        std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event ev, std::vector<Gen>& TruthColl, float weight, TString Label);


  TrigCheck();
  ~TrigCheck();

};



#endif

