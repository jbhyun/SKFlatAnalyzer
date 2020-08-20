#ifndef HNTopFeas_h
#define HNTopFeas_h

#include "AnalyzerCore.h"

class HNTopFeas : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MuMuMu, ElMuMu, SS2l, TetraLep; 
  bool DblMu, DblEG, MuEG;
  bool SystRun;
  vector<TString> TrigList_DblMu, TrigList_DblEG, TrigList_MuEG;

  void AnalyzeSSDiLepton(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                         std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeTriMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                     std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeElectronDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                             std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeTetraLepton(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                          std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);


  HNTopFeas();
  ~HNTopFeas();

};



#endif

