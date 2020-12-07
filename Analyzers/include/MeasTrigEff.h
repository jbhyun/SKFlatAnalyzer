#ifndef MeasTrigEff_h
#define MeasTrigEff_h

#include "AnalyzerCore.h"

class MeasTrigEff : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool ElEl, MuMu, ElMu; 
  bool DiMuTrig_MuLeg, DiMuTrig_DZ, DiMuMFilter, EMuTrig_ElLeg, EMuTrig_MuLeg, EMuTrig_DZ;
  bool SystRun;
  vector<TString> TrigList;
  TString SFKey_Trig;

  void MeasSiglEleTrigEff(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                          vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffDiMuTrig_MuLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                             vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffDiMuTrig_DZ(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                          vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffEMuTrig_ElLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                            vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffEMuTrig_MuLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                            vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffEMuTrig_DZ(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                         vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);


  MeasTrigEff();
  ~MeasTrigEff();

};



#endif

