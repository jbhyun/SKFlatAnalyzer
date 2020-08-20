#ifndef NewAnalyzer_h
#define NewAnalyzer_h

#include "AnalyzerCore.h"

class NewAnalyzer : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();



  double weight_Prefire;


  NewAnalyzer();
  ~NewAnalyzer();

};



#endif

