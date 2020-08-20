#ifndef SyncYield_h
#define SyncYield_h

#include "AnalyzerCore.h"

class SyncYield : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MuMuMu, ElMuMu, HctoWASync; 
  bool DblMu, DblEG, MuEG;
  bool SystRun;
  vector<TString> TrigList_DblMu_BtoG, TrigList_DblMu_H, TrigList_MuEG_BtoG, TrigList_MuEG_H;

  void AnalyzeTriMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                     std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeElectronDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                             std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  int  TriMuChargeIndex_SyncYield(vector<Muon>& MuonColl, float MET, float METx, float METy, TString charge);


  SyncYield();
  ~SyncYield();

};



#endif

