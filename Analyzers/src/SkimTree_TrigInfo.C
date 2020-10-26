#include "SkimTree_TrigInfo.h"

void SkimTree_TrigInfo::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_TrigInfo::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

}

void SkimTree_TrigInfo::executeEvent(){

  Event ev = GetEvent();

  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  std::sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  std::sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);

  bool SiglElPD = DataStream.Contains("SingleElectron"), SiglMuPD = DataStream.Contains("SingleMuon");
  float MinPt = -1.;
  vector<TString> TagTrigList;
  bool PassCriteria=false;

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
    if(!ev.PassTrigger(TagTrigList)) return;
    if(!(electronPreColl.size()>0 && electronPreColl.at(0).Pt()>MinPt)) return;

    bool TagMatchHLT=false;
    if( DataStream.Contains("TrigInfo") or MCSample.Contains("TrigInfo") ){
      for(unsigned int it_el=0; it_el<electronPreColl.size(); it_el++){
        if(electronPreColl.at(it_el).Pt()<MinPt) continue;
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(electronPreColl.at(it_el).Eta()-HLTObject_eta->at(it_obj),2)+pow(electronPreColl.at(it_el).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          for(unsigned int it_trig=0; it_trig<TagTrigList.size(); it_trig++){
            if(HLTObject_FiredPaths->at(it_obj).find(TagTrigList.at(it_trig))!=std::string::npos ){ TagMatchHLT=true; break; }
          }
          if(TagMatchHLT) break;
        }
      }
    }
    if(!TagMatchHLT) return;
    if(electronPreColl.size()>1){
      bool HasZPair=false;
      for(unsigned int i=0; i<electronPreColl.size(); i++){
        for(unsigned int j=i+1; j<electronPreColl.size(); j++){
          float Melel = (electronPreColl.at(i)+electronPreColl.at(j)).M();
          HasZPair = Melel>50 && Melel<130;
          if(HasZPair){ break; break; }
        }
      }
      if(HasZPair) PassCriteria=true;
    }
    if(muonPreColl.size()>0) PassCriteria=true;
  }

  //Single Muon PD
  if(SiglMuPD){
    if(!ev.PassTrigger(TagTrigList)) return;
    if(!(muonPreColl.size()>0 && muonPreColl.at(0).Pt()>MinPt)) return;

    bool TagMatchHLT=false;
    if( DataStream.Contains("TrigInfo") or MCSample.Contains("TrigInfo") ){
      for(unsigned int it_mu=0; it_mu<muonPreColl.size(); it_mu++){
        if(muonPreColl.at(it_mu).Pt()<MinPt) continue;
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(muonPreColl.at(it_mu).Eta()-HLTObject_eta->at(it_obj),2)+pow(muonPreColl.at(it_mu).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          for(unsigned int it_trig=0; it_trig<TagTrigList.size(); it_trig++){
            if(HLTObject_FiredPaths->at(it_obj).find(TagTrigList.at(it_trig))!=std::string::npos ){ TagMatchHLT=true; break; }
          }
          if(TagMatchHLT) break;
        }
      }
    }
    if(!TagMatchHLT) return;
    if(muonPreColl.size()>1){
      bool HasZPair=false;
      for(unsigned int i=0; i<muonPreColl.size(); i++){
        for(unsigned int j=i+1; j<muonPreColl.size(); j++){
          float Mmumu = (muonPreColl.at(i)+muonPreColl.at(j)).M();
          HasZPair = Mmumu>50 && Mmumu<130;
          if(HasZPair){ break; break; }
        }
      }
      if(HasZPair) PassCriteria=true;
    }
    if(electronPreColl.size()>0) PassCriteria=true;
  }

  if(PassCriteria) newtree->Fill();

}

void SkimTree_TrigInfo::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_TrigInfo::SkimTree_TrigInfo(){

  newtree = NULL;

}

SkimTree_TrigInfo::~SkimTree_TrigInfo(){

}

void SkimTree_TrigInfo::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}

