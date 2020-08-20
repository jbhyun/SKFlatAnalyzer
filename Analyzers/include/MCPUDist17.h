#ifndef MCPUDist17_h
#define MCPUDist17_h

#include "AnalyzerCore.h"

class MCPUDist17 : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool MeasMCPU; 
  bool SystRun;
  vector<TString> TrigList;
  TString SFKey_Trig;

  void MeasureMCBTagEfficiency(std::vector<Jet>& JetColl, JetTagging::Parameters jtp, float weight, TString Label);


  MCPUDist17();
  ~MCPUDist17();

};



#endif

