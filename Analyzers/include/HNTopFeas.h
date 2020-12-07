#ifndef HNTopFeas_h
#define HNTopFeas_h

#include "AnalyzerCore.h"

class HNTopFeas : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MuMuMu, ElMuMu, SSMuMu, SSElEl, SSElMu, SS2l, TriLep, TetraLep; 
  bool DblMu, DblEG, MuEG;
  bool SystRun;
  vector<TString> TrigList_DblMu, TrigList_DblEG, TrigList_MuEG;

  void AnalyzeSSDiLepton(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                         vector<Electron>& EleTColl, vector<Electron>& EleLColl, vector<Electron>& EleVColl,
                         vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeTriLepton(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                        vector<Electron>& EleTColl, vector<Electron>& EleLColl, vector<Electron>& EleVColl,
                        vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  void AnalyzeTetraLepton(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                          vector<Electron>& EleTColl, vector<Electron>& EleLColl, vector<Electron>& EleVColl,
                          vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label);
  int  CheckDecMode(vector<Gen>& TruthColl);


  HNTopFeas();
  ~HNTopFeas();

};



#endif

