#ifndef BTagEff_MC_h
#define BTagEff_MC_h

#include "AnalyzerCore.h"

class BTagEff_MC : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MCEffMeas; 
  bool SystRun;
  vector<TString> TrigList;
  TString SFKey_Trig;

  void MeasureMCBTagEfficiency(std::vector<Jet>& JetColl, JetTagging::Parameters jtp, float weight, TString Label);


  BTagEff_MC();
  ~BTagEff_MC();

};



#endif

