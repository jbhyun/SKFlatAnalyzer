#ifndef KinVarSearch_h
#define KinVarSearch_h

#include "AnalyzerCore.h"

class KinVarSearch : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MuMuMu, ElMuMu, SSMuMu, SSElEl, SSElMu, SS2l, TriLep, TetraLep; 
  bool DblMu, DblEG, MuEG;
  bool SystRun;
  vector<TString> TrigList_DblMu, TrigList_DblEG, TrigList_MuEG;

  void AnalyzeSSDiLepton(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                         std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeTriLepton(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                 vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);

  void AnalyzeTriMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                     std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeElectronDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                             std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeTetraLepton(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& EleTColl, std::vector<Electron>& EleLColl,
                          std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);


  KinVarSearch();
  ~KinVarSearch();

};



#endif

