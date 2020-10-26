#ifndef SkimRateCheck_h
#define SkimRateCheck_h

#include "AnalyzerCore.h"

class SkimRateCheck : public AnalyzerCore {

public:


  bool SS2lOR3lRate, TrigInfoRate; 
  void CheckSSOR3lSkim();
  void CheckTrigInfoSkim();

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimRateCheck();
  ~SkimRateCheck();

};



#endif

