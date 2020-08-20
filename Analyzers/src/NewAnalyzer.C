#include "NewAnalyzer.h"

void NewAnalyzer::initializeAnalyzer(){


  //==== corresponding Muon ID SF Keys for mcCorr->MuonID_SF()


}

void NewAnalyzer::executeEvent(){

  if(IsDATA) return;

  if(!PassMETFilter()) return;

  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  if( !( muonPreColl.size()>1 && (electronPreColl.size()+muonPreColl.size())>2 ) ) return;


  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "POGIDTIsoM", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "POGMVAniWP90_Isop1", 10., 2.5);
  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "POGIDTIsoVVL", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "POGMVAniWP90_Isop1", 10., 2.5);


  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
  std::vector<Jet> jetNoVetoColl  = GetJets("tight", 25., 2.4);
  std::vector<Jet> bjetNoVetoColl = SelBJets(jetNoVetoColl, param_jets);
  std::vector<Jet> jetColl  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl, muonLooseColl, 0.4);
  std::vector<Jet> bjetColl = SelBJets(jetColl, param_jets);

  int NmuT=muonTightColl.size();
  int NelT=electronTightColl.size();
  //int NmuL=muonLooseColl.size();
  //int NelL=electronLooseColl.size();
  int NlepT=NmuT+NelT;

  if( !(NlepT>2 && NmuT>1) ) return;

  std::vector<Muon>   muonAColl;
  std::vector<Lepton> leptonWColl;
  std::vector<Lepton> leptonHcColl;
  std::vector<Gen>    truthColl = GetGens();
    for(unsigned int i=0; i<muonTightColl.size(); i++){
      int LeptonType = GetLeptonType_JH(muonTightColl.at(i), truthColl);
      if(LeptonType==2){
        muonAColl.push_back(muonTightColl.at(i));
      }
      else if(LeptonType==22){
        leptonHcColl.push_back(muonTightColl.at(i));
      }
      else if(LeptonType==1 or LeptonType==3){
        leptonWColl.push_back(muonTightColl.at(i));
      }
    }
    for(unsigned int i=0; i<electronTightColl.size(); i++){
      int LeptonType = GetLeptonType_JH(electronTightColl.at(i), truthColl);
      if(LeptonType==22){
        leptonHcColl.push_back(electronTightColl.at(i));
      }
      else if(LeptonType==1 or LeptonType==3){
        leptonWColl.push_back(electronTightColl.at(i));
      }
    }
  std::sort(leptonWColl.begin(), leptonWColl.end(), PtComparing);
  std::sort(leptonHcColl.begin(), leptonHcColl.end(), PtComparing);

  for(unsigned int i=0; i<muonAColl.size(); i++){
    FillHist("Pt_muA", muonAColl.at(i).Pt(), 1., 60, 0., 300.);
  }
  for(unsigned int i=0; i<leptonWColl.size(); i++){
    FillHist("Pt_lW", leptonWColl.at(i).Pt(), 1., 60, 0., 300.);
  }
  for(unsigned int i=0; i<leptonHcColl.size(); i++){
    FillHist("Pt_lHc", leptonHcColl.at(i).Pt(), 1., 60, 0., 300.);
  }


  FillHist("Ptmu1", muonTightColl.at(0).Pt(), 1., 60, 0., 300.);
  FillHist("Ptmu2", muonTightColl.at(1).Pt(), 1., 60, 0., 300.);
  if(NmuT>2) FillHist("Ptmu3", muonTightColl.at(2).Pt(), 1., 60, 0., 300.);
  if(NelT>0) FillHist("Ptel1", electronTightColl.at(0).Pt(), 1., 60., 0., 300);

  FillHist("Etamu1", muonTightColl.at(0).Eta(), 1., 20, -5., 5.);
  FillHist("Etamu2", muonTightColl.at(1).Eta(), 1., 20, -5., 5.);
  if(NmuT>2) FillHist("Etamu3", muonTightColl.at(2).Eta(), 1., 20, -5., 5.);
  if(NelT>0) FillHist("Etael1", electronTightColl.at(0).Eta(), 1., 20., -5., 5.);

  FillHist("Nj", jetColl.size(), 1., 10., 0., 10.);
  FillHist("Nb", bjetColl.size(), 1., 10., 0., 10.);

  bool Acc_1e2mu=false, Acc_3mu=false, Acc_4l=false;
  bool Acc_j_1e2mu=false, Acc_j_3mu=false, Acc_j_4l=false;
  bool Acc_2j=jetColl.size()>1;
  bool Acc_1b=bjetColl.size()>0;
  if(NelT==1 && NmuT==2){
    if(electronTightColl.at(0).Pt()>20 or muonTightColl.at(0).Pt()>20){
      Acc_1e2mu=true;
      if(Acc_2j && Acc_1b){ Acc_j_1e2mu=true; }
    }
  }
  else if(NelT==0 && NmuT==3){
    if(muonTightColl.at(0).Pt()>20){
      Acc_3mu=true;
      if(Acc_2j && Acc_1b){ Acc_j_3mu=true; }
    }
  }
  else if(NlepT==4 && NmuT>=2){
    if( (NelT>0 and electronTightColl.at(0).Pt()>20) or muonTightColl.at(0).Pt()>20){
      Acc_4l=true;
      if(Acc_2j && Acc_1b){ Acc_j_4l=true; }
    }
  }

  if(Acc_1e2mu){ FillHist("CountAcc", 0., 1., 10, 0., 10.); }
  if(Acc_3mu  ){ FillHist("CountAcc", 1., 1., 10, 0., 10.); }
  if(Acc_4l   ){ FillHist("CountAcc", 2., 1., 10, 0., 10.); }
  if(Acc_j_1e2mu){ FillHist("CountAcc", 3., 1., 10, 0., 10.); }
  if(Acc_j_3mu  ){ FillHist("CountAcc", 4., 1., 10, 0., 10.); }
  if(Acc_j_4l   ){ FillHist("CountAcc", 5., 1., 10, 0., 10.); }
  

//  AnalyzerParameter param;
//  executeEventFromParameter(param);

}

void NewAnalyzer::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

NewAnalyzer::NewAnalyzer(){

}

NewAnalyzer::~NewAnalyzer(){

}


