#ifndef SkimTree_TrigInfo_h
#define SkimTree_TrigInfo_h

#include "AnalyzerCore.h"

class SkimTree_TrigInfo : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_TrigInfo();
  ~SkimTree_TrigInfo();

  TTree *newtree;

  vector<TString> triggers;
  void WriteHist();

};



#endif

