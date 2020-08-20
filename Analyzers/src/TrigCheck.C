#include "TrigCheck.h"

void TrigCheck::initializeAnalyzer(){

  MFilter=false, QFilter=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("MFilter")) MFilter=true;
    if(Userflags.at(i).Contains("QFilter")) QFilter=true;
    if(Userflags.at(i).Contains("SystRun")) SystRun=true; 
  }

  DblMu=false, DblEG=false, MuEG=false;
  if     (DataStream.Contains("DoubleMuon")) DblMu=true;
  else if(DataStream.Contains("MuonEG"))     MuEG =true;
  else if(DataStream.Contains("DoubleEG"))   DblEG=true;
  else if(DataYear==2018 and DataStream.Contains("EGamma")) DblEG=true;

  if(MFilter) TrigList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  if(DataYear==2016){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu.push_back("HLT_TripleMu_12_10_5_v");

    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  if(DataYear==2017){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    TrigList_DblMu.push_back("HLT_TripleMu_10_5_5_DZ_v");
    TrigList_DblMu.push_back("HLT_TripleMu_12_10_5_v");

    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  }
  if(DataYear==2018){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    TrigList_DblMu.push_back("HLT_TripleMu_10_5_5_DZ_v");
    TrigList_DblMu.push_back("HLT_TripleMu_12_10_5_v");

    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void TrigCheck::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  bool PassTrig=false;
  if(MFilter){ PassTrig = ev.PassTrigger(TrigList); } 
  if(QFilter) PassTrig=true;
  if(!PassTrig) return;
  float TmpW = IsDATA? 1:ev.MCweight()*weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
  FillHist("CutFlow", 0., weight*TmpW, 10, 0., 10.);
  if(!PassMETFilter()) return;
  FillHist("CutFlow", 1., weight*TmpW, 10, 0., 10.);

  bool PreCutPass=false;
  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  std::sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  std::sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);
  if((QFilter or MFilter) && muonPreColl.size()>1) PreCutPass=true; 
  if(!PreCutPass) return;
  FillHist("CutFlow", 2., weight*TmpW, 10, 0., 10.);


  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TESTT", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TESTT", 10., 2.5);
  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TESTL", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TESTL", 10., 2.5);


  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  std::vector<Jet> jetNoVetoColl  = GetJets("tight", 25., 2.4);
  std::sort(jetNoVetoColl.begin(), jetNoVetoColl.end(), PtComparing);
  std::vector<Jet> bjetNoVetoColl = SelBJets(jetNoVetoColl, param_jets);
  std::vector<Jet> jetColl  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl, muonLooseColl, 0.4);
  std::vector<Jet> bjetColl = SelBJets(jetColl, param_jets);


  Particle vMET = ev.GetMETVector();
  Particle vMET_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);


  std::vector<Gen> truthColl;


  bool EventCand = false;
  if(MFilter or QFilter){ EventCand = muonLooseColl.size()>1; }

  float w_gen = 1., w_filter = 1., w_topptrw = 1., w_lumi = 1., w_PU = 1., w_prefire = 1., sf_trig = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
  if((!IsDATA) and EventCand){
    w_gen     = ev.MCweight();
    w_lumi    = weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    //if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    //w_filter  = GetGenFilterEffCorr();
    //w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    w_PU      = GetPileUpWeight(nPileUp, 0);
    w_prefire = GetPrefireWeight(0);
    sf_muid   = GetMuonSF(muonTightColl, "POGTID_genTrk", "ID");
    sf_muiso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    sf_elreco = GetElectronSF(electronTightColl, "", "Reco");
    sf_elid   = GetElectronSF(electronTightColl, "POGMVAIsoWP90", "ID");
    sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    //sf_trig   = mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig, "");
    //cout<<"w_gen:"<<w_gen<<" w_lumi:"<<w_lumi<<" w_PU:"<<w_PU<<" w_prefire:"<<w_prefire<<" sf_trig:"<<sf_trig<<endl;
    //cout<<"sf_mutk"<<sf_mutk<<" sf_muid:"<<sf_muid<<" sf_muiso:"<<sf_muiso<<" sf_elreco:"<<sf_elreco<<" sf_elid:"<<sf_elid<<" sf_btag:"<<sf_btag<<endl;
  }
  weight *= w_gen * w_filter * w_topptrw * w_lumi * w_PU * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

 
  if(MFilter){
    CheckMFilterBias(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                     jetColl, bjetColl, vMET_xyCorr, ev, weight, "");
  }
  if(QFilter){
    CheckQFilterBias(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                     jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "");
  }



}


void TrigCheck::CheckMFilterBias(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                 std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event ev, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size();

  if( !(NMuT==2 && NElT==0) ) return;
  if( !(NMuL==2 && NElL==0) ) return;
  if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) return;

  float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();

  vector<TString> TrigListM3p8 = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
  vector<TString> TrigListM8p0 = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"};

  //LowMass
  FillHist("Mmumu_Low_TrIDIso"+Label, Mmumu, weight, 200, 0., 20.);
  if(ev.PassTrigger(TrigListM3p8)) FillHist("Mmumu_Low_TrIDIsoM3p8"+Label, Mmumu, weight, 200, 0., 20.);
  if(ev.PassTrigger(TrigListM8p0)) FillHist("Mmumu_Low_TrIDIsoM8p0"+Label, Mmumu, weight, 200, 0., 20.);

  //HighMass
  FillHist("Mmumu_High_TrIDIso"+Label, Mmumu, weight, 50, 0., 500.);
  if(ev.PassTrigger(TrigListM3p8)) FillHist("Mmumu_High_TrIDIsoM3p8"+Label, Mmumu, weight, 50, 0., 500.);
  if(ev.PassTrigger(TrigListM8p0)) FillHist("Mmumu_High_TrIDIsoM8p0"+Label, Mmumu, weight, 50, 0., 500.);


}


void TrigCheck::CheckQFilterBias(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                 std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event ev, std::vector<Gen>& TruthColl, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size();

  if( !(NMuT==2 && NElT==0) ) return;
  if( !(NMuL==2 && NElL==0) ) return;
  if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) return;

  TruthColl = GetGens();
  int LepType_Mu1 = GetLeptonType_JH(MuTColl.at(0), TruthColl);
  int LepType_Mu2 = GetLeptonType_JH(MuTColl.at(1), TruthColl);

  //Only prompt 2 leptons within trigger acceptance
  if( !(LepType_Mu1>0 && LepType_Mu1<4) ) return; 
  if( !(LepType_Mu2>0 && LepType_Mu2<4) ) return; 
  int NCountPr=0;
  for(unsigned int it_gen=2; it_gen<TruthColl.size(); it_gen++){
    int GenType = GetLeptonType_JH(it_gen, TruthColl);
    if( !(GenType>0 && GenType<4) ) continue;
    if(fabs(TruthColl.at(it_gen).Eta())>2.4) continue;
    if(TruthColl.at(it_gen).Pt()<8.) continue;
    NCountPr++;
  }
  if(NCountPr!=2) return;

  float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
  vector<TString> TrigListIsoDiMu = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"};
  vector<TString> TrigListSSDiMu  = {"HLT_Mu17_Mu8_SameSign_DZ_v"};
  vector<TString> TrigListDiMu    = {"HLT_Mu17_Mu8_DZ_v"};

  //TrigPerf
  int NPtBinEdges=10;
  double PtBinEdges[NPtBinEdges]   = {0., 10., 20., 30., 50., 80., 120., 200., 500., 1000.};
  FillHist("Mmumu_ID"+Label, Mmumu, weight, NPtBinEdges-1, PtBinEdges);
  if(ev.PassTrigger(TrigListIsoDiMu)) FillHist("Mmumu_IDTrigIso"+Label, Mmumu, weight, NPtBinEdges-1, PtBinEdges);
  if(Mmumu<10){
    FillHist("PtMu1_Mlt10_Incl", MuTColl.at(0).Pt(), weight, 20, 0., 200.);
    FillHist("PtMu2_Mlt10_Incl", MuTColl.at(1).Pt(), weight, 20, 0., 200.);
    FillHist("EtaMu1_Mlt10_Incl", MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("EtaMu2_Mlt10_Incl", MuTColl.at(1).Eta(), weight, 20, -5., 5.);
    FillHist("dRM1M2_Mlt10_Incl", MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 50, 0., 5.);
  }
  else{
    FillHist("PtMu1_Mgt10_Incl", MuTColl.at(0).Pt(), weight, 20, 0., 200.);
    FillHist("PtMu2_Mgt10_Incl", MuTColl.at(1).Pt(), weight, 20, 0., 200.);
    FillHist("EtaMu1_Mgt10_Incl", MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("EtaMu2_Mgt10_Incl", MuTColl.at(1).Eta(), weight, 20, -5., 5.);
    FillHist("dRM1M2_Mgt10_Incl", MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 50, 0., 5.);
  }

  if(MuTColl.at(0).Charge()==MuTColl.at(1).Charge()){
    FillHist("MmumuSS_ID"+Label, Mmumu, weight, NPtBinEdges-1, PtBinEdges);
    if(ev.PassTrigger(TrigListSSDiMu))  FillHist("MmumuSS_IDTrigSS" +Label, Mmumu, weight, NPtBinEdges-1, PtBinEdges);
    if(ev.PassTrigger(TrigListDiMu)){
      FillHist("MmumuSS_IDTrigDiMu"+Label, Mmumu, weight, NPtBinEdges-1, PtBinEdges);
      if(ev.PassTrigger(TrigListSSDiMu)) FillHist("MmumuSS_IDTrigDiMuSS" +Label, Mmumu, weight, NPtBinEdges-1, PtBinEdges);
    }

    if(Mmumu<10){
      FillHist("PtMu1_Mlt10_SS", MuTColl.at(0).Pt(), weight, 20, 0., 200.);
      FillHist("PtMu2_Mlt10_SS", MuTColl.at(1).Pt(), weight, 20, 0., 200.);
      FillHist("EtaMu1_Mlt10_SS", MuTColl.at(0).Eta(), weight, 20, -5., 5.);
      FillHist("EtaMu2_Mlt10_SS", MuTColl.at(1).Eta(), weight, 20, -5., 5.);
      FillHist("dRM1M2_Mlt10_SS", MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 50, 0., 5.);
    }
    else{
      FillHist("PtMu1_Mgt10_SS", MuTColl.at(0).Pt(), weight, 20, 0., 200.);
      FillHist("PtMu2_Mgt10_SS", MuTColl.at(1).Pt(), weight, 20, 0., 200.);
      FillHist("EtaMu1_Mgt10_SS", MuTColl.at(0).Eta(), weight, 20, -5., 5.);
      FillHist("EtaMu2_Mgt10_SS", MuTColl.at(1).Eta(), weight, 20, -5., 5.);
      FillHist("dRM1M2_Mgt10_SS", MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 50, 0., 5.);
    }
  }

}




void TrigCheck::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

TrigCheck::TrigCheck(){

}

TrigCheck::~TrigCheck(){

}


