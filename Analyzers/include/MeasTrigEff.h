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

  void MeasSiglEleTrigEff(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                          std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffDiMuTrig_MuLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                             vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffDiMuTrig_DZ(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                           std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffEMuTrig_ElLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                            vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);
  void MeasEffEMuTrig_MuLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                            vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label);



  MeasTrigEff();
  ~MeasTrigEff();

};



#endif

