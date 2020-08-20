#ifndef SkimRateCheck_h
#define SkimRateCheck_h

#include "AnalyzerCore.h"

class SkimRateCheck : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimRateCheck();
  ~SkimRateCheck();

};



#endif

