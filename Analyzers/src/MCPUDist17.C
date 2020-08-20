#include "MCPUDist17.h"

void MCPUDist17::initializeAnalyzer(){


  //Set up the tagger map only for the param settings to be used.

}


void MCPUDist17::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  FillHist("NPileUp", nPileUp, weight, 99, 0., 99.); 


}



void MCPUDist17::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

MCPUDist17::MCPUDist17(){

}

MCPUDist17::~MCPUDist17(){

}


