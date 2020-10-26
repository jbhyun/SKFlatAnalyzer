#include "SkimRateCheck.h"

void SkimRateCheck::initializeAnalyzer(){

  SS2lOR3lRate=false, TrigInfoRate=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("SS2lOR3lRate")) SS2lOR3lRate=true; 
    if(Userflags.at(i).Contains("TrigInfoRate")) TrigInfoRate=true; 
  }

}


void SkimRateCheck::executeEvent(){


  if(SS2lOR3lRate){ CheckSSOR3lSkim(); }
  if(TrigInfoRate){ CheckTrigInfoSkim(); } 

}


void SkimRateCheck::CheckTrigInfoSkim(){

  Event ev = GetEvent();

  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  std::sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  std::sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);

  bool SiglElPD = DataStream.Contains("SingleElectron"), SiglMuPD = DataStream.Contains("SingleMuon");
  float MinPt = -1.;
  vector<TString> TagTrigList;

  if(SiglElPD){
    if(DataYear==2016){
      MinPt=27.; TagTrigList.push_back("HLT_Ele27_WPTight_Gsf_v"); TagTrigList.push_back("HLT_Ele27_eta2p1_WPTight_Gsf_v");
    }
    else if(DataYear==2017){
      MinPt=32.; TagTrigList.push_back("HLT_Ele32_WPTight_Gsf_v"); TagTrigList.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v");
    }
    else if(DataYear==2018){
      MinPt=32.; TagTrigList.push_back("HLT_Ele32_WPTight_Gsf_v");
    }
    else return;
  }
  else if(SiglMuPD){
    if(DataYear==2016){
      MinPt=24.; TagTrigList.push_back("HLT_IsoMu24_v"); TagTrigList.push_back("HLT_IsoTkMu24_v");
    }
    else if(DataYear==2017){
      MinPt=27.; TagTrigList.push_back("HLT_IsoMu27_v");
    }
    else if(DataYear==2018){
      MinPt=24.; TagTrigList.push_back("HLT_IsoMu24_v");
    }
    else return;
  }
  else return;


  if(SiglElPD){
    int it_cut=0;
    FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.); it_cut++;

    if(!ev.PassTrigger(TagTrigList)) return;
    FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.); it_cut++;

    if(!(electronPreColl.size()>0 && electronPreColl.at(0).Pt()>MinPt)) return;
    FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.); it_cut++;

    bool TagMatchHLT=false; int it_tag=-1;
    if( DataStream.Contains("TrigInfo") or MCSample.Contains("TrigInfo") ){
      for(unsigned int it_el=0; it_el<electronPreColl.size(); it_el++){
        if(electronPreColl.at(it_el).Pt()<MinPt) continue;
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(electronPreColl.at(it_el).Eta()-HLTObject_eta->at(it_obj),2)+pow(electronPreColl.at(it_el).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          for(unsigned int it_trig=0; it_trig<TagTrigList.size(); it_trig++){
            if(HLTObject_FiredPaths->at(it_obj).find(TagTrigList.at(it_trig))!=std::string::npos ){ TagMatchHLT=true; it_tag=it_el; break; }
          }
          if(TagMatchHLT) break;
        }
      }
    }
    if(!TagMatchHLT) return;
    FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.); it_cut++;

    if(electronPreColl.size()>1){
      FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.); it_cut++;

      bool HasZPair=false;
      for(unsigned int i=0; i<electronPreColl.size(); i++){
        for(unsigned int j=i+1; j<electronPreColl.size(); j++){
          float Melel = (electronPreColl.at(i)+electronPreColl.at(j)).M();
          HasZPair = Melel>50 && Melel<130;
          if(HasZPair){ break; break; }
        }
      }
      if(HasZPair){
        FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.);
      }
      it_cut++;
    }
    else{it_cut+=2;}

    if(muonPreColl.size()>0){
      FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.); it_cut++;
      
      int NCand=0;
      for(unsigned int it_m=0; it_m<muonPreColl.size(); it_m++){
        if(muonPreColl.at(it_m).DeltaR(electronPreColl.at(it_tag))<0.4) continue; 
        NCand++;
      }
      
     if(NCand>0){ FillHist("CutFlow_SiglEl", it_cut, 1, 20, 0., 20.); }
     it_cut++;
    }
  }

  //Single Muon PD
  if(SiglMuPD){
    int it_cut=0;
    FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.); it_cut++;

    if(!ev.PassTrigger(TagTrigList)) return;
    FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.); it_cut++;

    if(!(muonPreColl.size()>0 && muonPreColl.at(0).Pt()>MinPt)) return;
    FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.); it_cut++;

    bool TagMatchHLT=false; int it_tag=-1;
    if( DataStream.Contains("TrigInfo") or MCSample.Contains("TrigInfo") ){
      for(unsigned int it_mu=0; it_mu<muonPreColl.size(); it_mu++){
        if(muonPreColl.at(it_mu).Pt()<MinPt) continue;
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(muonPreColl.at(it_mu).Eta()-HLTObject_eta->at(it_obj),2)+pow(muonPreColl.at(it_mu).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          for(unsigned int it_trig=0; it_trig<TagTrigList.size(); it_trig++){
            if(HLTObject_FiredPaths->at(it_obj).find(TagTrigList.at(it_trig))!=std::string::npos ){ TagMatchHLT=true; it_tag=it_mu; break; }
          }
          if(TagMatchHLT) break;
        }
      }
    }
    if(!TagMatchHLT) return;
    FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.); it_cut++;

    if(muonPreColl.size()>1){
      FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.); it_cut++;

      bool HasZPair=false;
      for(unsigned int i=0; i<muonPreColl.size(); i++){
        for(unsigned int j=i+1; j<muonPreColl.size(); j++){
          float Mmumu = (muonPreColl.at(i)+muonPreColl.at(j)).M();
          HasZPair = Mmumu>50 && Mmumu<130;
          if(HasZPair){ break; break; }
        }
      }
      if(HasZPair){
        FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.);
      }
      it_cut++;
    }
    else{it_cut+=2;}

    if(electronPreColl.size()>0){
      FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.); it_cut++;
      
      int NCand=0;
      for(unsigned int it_e=0; it_e<electronPreColl.size(); it_e++){
        if(electronPreColl.at(it_e).DeltaR(muonPreColl.at(it_tag))<0.4) continue; 
        NCand++;
      }
     if(NCand>0){ FillHist("CutFlow_SiglMu", it_cut, 1, 20, 0., 20.); }
     it_cut++;
    }
  }


}


void SkimRateCheck::CheckSSOR3lSkim(){

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


