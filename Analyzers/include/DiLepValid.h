#ifndef DiLepValid_h
#define DiLepValid_h

#include "AnalyzerCore.h"

class DiLepValid : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool ElEl, MuMu, ElMu; 
  bool TrigEffMeas, HEMCheck;
  bool SystRun;
  vector<TString> TrigList;
  TString SFKey_Trig;

  void AnalyzeDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                     std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeDiElectron(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                         std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeElectronMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                           std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);


  void MeasSiglEleTrigEff(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                          std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void CheckHEMIssue(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                     std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label, TString Option);



  DiLepValid();
  ~DiLepValid();

};



#endif

