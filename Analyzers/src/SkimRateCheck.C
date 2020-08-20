#include "SkimRateCheck.h"

void SkimRateCheck::initializeAnalyzer(){

}


void SkimRateCheck::executeEvent(){

  Event ev;

  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);

  int NEl  = electronPreColl.size();
  int NMu  = muonPreColl.size();
  int NLep = NEl+NMu;
  bool HasSS2lOR3l = false;
  if      ( NLep >= 3 ){ HasSS2lOR3l = true; }
  else if ( NLep == 2 ){
    int QTot = 0;
    for(unsigned int it_m=0; it_m<muonPreColl.size(); it_m++){ QTot+=muonPreColl.at(it_m).Charge(); }
    for(unsigned int it_e=0; it_e<electronPreColl.size(); it_e++){ QTot+=electronPreColl.at(it_e).Charge(); }
    if(abs(QTot)==2) HasSS2lOR3l = true;
  }

  if(!HasSS2lOR3l) return;


  bool CheckComposition=true;
  if(CheckComposition){
    float weight = 1.;

    std::vector<Muon>     muonPreColl8     = SelectMuons(muonPreColl, "NOCUT", 8., 2.4);
    std::vector<Electron> electronPreColl8 = SelectElectrons(electronPreColl, "NOCUT", 8., 2.5);

    int NEl8 = electronPreColl8.size();
    int NMu8 = muonPreColl8.size();
    int NLep8 = NEl8+NMu8;
    int QTot_E5M5 = 0, QTot_E8M5 = 0, QTot_E8M8 = 0;
    bool LeadMu17 = false, LeadEl23 = false, PassLeadLepPt = false;
    for(unsigned int it_m=0; it_m<muonPreColl.size(); it_m++){
      QTot_E5M5+=muonPreColl.at(it_m).Charge();
      QTot_E8M5+=muonPreColl.at(it_m).Charge();
      if(it_m==0){ LeadMu17=muonPreColl.at(it_m).Pt()>17; }
    }
    for(unsigned int it_m=0; it_m<muonPreColl8.size(); it_m++){
      QTot_E8M8+=muonPreColl8.at(it_m).Charge();
    }
    for(unsigned int it_e=0; it_e<electronPreColl.size(); it_e++){
      QTot_E5M5+=electronPreColl.at(it_e).Charge();
      if(it_e==0){ LeadEl23=electronPreColl.at(it_e).Pt()>23; }
    }
    for(unsigned int it_e=0; it_e<electronPreColl8.size(); it_e++){
      QTot_E8M5+=electronPreColl8.at(it_e).Charge();
      QTot_E8M8+=electronPreColl8.at(it_e).Charge();
    }
    PassLeadLepPt = LeadMu17 or LeadEl23;


    if(NLep==2 && abs(QTot_E5M5)==2){
      FillHist("Composition_E5M5", 1., weight, 10, 0., 10.);
      if(PassLeadLepPt) FillHist("Composition_E5M5", 3., weight, 10, 0., 10.);
    }
    else if(NLep>2){
      FillHist("Composition_E5M5", 2., weight, 10, 0., 10.);
      if(PassLeadLepPt) FillHist("Composition_E5M5", 4., weight, 10, 0., 10.);
    }
    if((NEl8+NMu)==2 && abs(QTot_E8M5)==2){
      FillHist("Composition_E8M5", 1., weight, 10, 0., 10.);
      if(PassLeadLepPt) FillHist("Composition_E8M5", 3., weight, 10, 0., 10.);
    }
    else if((NEl8+NMu)>2){
      FillHist("Composition_E8M5", 2., weight, 10, 0., 10.);
      if(PassLeadLepPt) FillHist("Composition_E8M5", 4., weight, 10, 0., 10.);
    }
    if(NLep8==2 && abs(QTot_E8M8)==2){
      FillHist("Composition_E8M8", 1., weight, 10, 0., 10.);
      if(PassLeadLepPt) FillHist("Composition_E8M8", 3., weight, 10, 0., 10.);
    }
    else if(NLep8>2){
      FillHist("Composition_E8M8", 2., weight, 10, 0., 10.);
      if(PassLeadLepPt) FillHist("Composition_E8M8", 4., weight, 10, 0., 10.);
    }
  }


}



void SkimRateCheck::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

SkimRateCheck::SkimRateCheck(){

}

SkimRateCheck::~SkimRateCheck(){

}


