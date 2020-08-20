#include "TestRun.h"

void TestRun::initializeAnalyzer(){

  ElMuMu=false, MuMuMu=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("ElMuMu")) ElMuMu=true; 
    if(Userflags.at(i).Contains("MuMuMu")) MuMuMu=true; 
    if(Userflags.at(i).Contains("SystRun")) SystRun=true; 
  }

  if(MuMuMu){
    if(DataYear==2016){
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
      TrigList_3mu.push_back("HLT_TripleMu_12_10_5_v");
    }
    if(DataYear==2017){
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
      TrigList_3mu.push_back("HLT_TripleMu_10_5_5_DZ_v");
      TrigList_3mu.push_back("HLT_TripleMu_12_10_5_v");
    }
    if(DataYear==2018){
      TrigList_3mu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
      TrigList_3mu.push_back("HLT_TripleMu_10_5_5_DZ_v");
      TrigList_3mu.push_back("HLT_TripleMu_12_10_5_v");
    }
  }
  if(ElMuMu){
    if(DataYear==2016){
      TrigList_1e2mu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
      TrigList_1e2mu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TrigList_1e2mu.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
      TrigList_1e2mu.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TrigList_1e2mu.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");
    }
    if(DataYear==2017){
      TrigList_1e2mu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TrigList_1e2mu.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
      TrigList_1e2mu.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v");
    }
    if(DataYear==2018){
      TrigList_1e2mu.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TrigList_1e2mu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TrigList_1e2mu.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v");
    }
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void TestRun::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  if(MuMuMu && !ev.PassTrigger(TrigList_3mu) ) return; 
  if(ElMuMu && !ev.PassTrigger(TrigList_1e2mu) ) return;
  float TmpW = IsDATA? 1:ev.MCweight()*weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
  FillHist("CutFlow", 0., weight*TmpW, 10, 0., 10.);
  if(!PassMETFilter()) return;
  FillHist("CutFlow", 1., weight*TmpW, 10, 0., 10.);

  bool PreCutPass=false;
  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  if(MuMuMu and muonPreColl.size()>2    ) PreCutPass=true;
  if(ElMuMu and electronPreColl.size()>0 && muonPreColl.size()>1) PreCutPass=true;
  if(!PreCutPass) return;
  FillHist("CutFlow", 2., weight*TmpW, 10, 0., 10.);


  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TESTT", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TESTT", 10., 2.5);
  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TESTL", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TESTL", 10., 2.5);

  std::vector<Muon>     muonTightColl2     = SelectMuons(muonPreColl, "TEST2T", 10., 2.4);
  std::vector<Electron> electronTightColl2 = SelectElectrons(electronPreColl, "TEST2T", 10., 2.5);
  std::vector<Muon>     muonLooseColl2     = SelectMuons(muonPreColl, "TEST2L", 10., 2.4);
  std::vector<Electron> electronLooseColl2 = SelectElectrons(electronPreColl, "TEST2L", 10., 2.5);

  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  std::vector<Jet> jetNoVetoColl  = GetJets("tight", 25., 2.4);
  std::vector<Jet> bjetNoVetoColl = SelBJets(jetNoVetoColl, param_jets);
  std::vector<Jet> jetColl  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl, muonLooseColl, 0.4);
  std::vector<Jet> bjetColl = SelBJets(jetColl, param_jets);

  std::vector<Jet> jetColl2  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl2, muonLooseColl2, 0.4);
  std::vector<Jet> bjetColl2 = SelBJets(jetColl2, param_jets);



  Particle vMET = ev.GetMETVector();
  Particle vMET_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);


  std::vector<Gen> truthColl;


  bool EventCand = false;
  if(MuMuMu){ EventCand = muonLooseColl.size()>2; }
  if(ElMuMu){ EventCand = electronLooseColl.size()>0 && muonLooseColl.size()>1; }

  float w_gen = 1., w_filter = 1., w_topptrw = 1., w_lumi = 1., w_PU = 1., w_prefire = 1., sf_trig = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
  if((!IsDATA) and EventCand){
    if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    w_gen     = ev.MCweight();
//    w_filter  = GetGenFilterEffCorr();
    w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    w_lumi    = weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
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

 
//  if(ElMu){
//      CheckHEMIssue(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET_xyCorr,
//                    weight, "", "");
//  }

  if(MuMuMu){
    AnalyzeTriMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                   jetColl, bjetColl, vMET_xyCorr, weight, "_Iso04");
    AnalyzeTriMuon(muonTightColl2, muonLooseColl2, electronTightColl2, electronLooseColl2,
                   jetColl2, bjetColl2, vMET_xyCorr, weight, "_MiniIso");
  }
  if(ElMuMu){
    AnalyzeElectronDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                          jetColl, bjetColl, vMET_xyCorr, weight, "_Iso04");
    AnalyzeElectronDiMuon(muonTightColl2, muonLooseColl2, electronTightColl2, electronLooseColl2,
                          jetColl2, bjetColl2, vMET_xyCorr, weight, "_MiniIso");
  }


}


void TestRun::AnalyzeTriMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  if( !(MuTColl.size()==3 && MuLColl.size()==3) ) return;
  if( ElLColl.size()!=0 ) return;
  if( !(MuTColl.at(0).Pt()>20. && MuTColl.at(2).Pt()>10) ) return;
  if( abs(SumCharge(MuTColl,ElTColl)) != 1 ) return;
  FillHist("CutFlow"+Label, 3., weight, 10, 0., 10.);
  
  
  int IdxOS = TriMuChargeIndex(MuTColl, "OS");
  int IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
  int IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
  float Mmumu1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
  float Mmumu2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
  if( !(Mmumu1>12 && Mmumu2>12) ) return;
  FillHist("CutFlow"+Label, 4., weight, 10, 0., 10.);
  if( !(fabs(Mmumu1-91.2)>10 && fabs(Mmumu2-91.2)>10) ) return;
  FillHist("CutFlow"+Label, 5., weight, 10, 0., 10.);

  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, 6., weight, 10, 0., 10.);

  FillHist("Mmumu1_3l1b1j_OnShell"+Label, Mmumu1, weight, 680, 12., 80.);
  FillHist("Mmumu1_3l1b1j_OffShell"+Label, Mmumu1, weight, 400, 100., 500.);
  FillHist("Mmumu2_3l1b1j_OnShell"+Label, Mmumu2, weight, 680, 12., 80.);
  FillHist("Mmumu2_3l1b1j_OffShell"+Label, Mmumu2, weight, 400, 100., 500.);
  if(JetColl.size()>1){
    FillHist("Mmumu1_3l1b2j_OnShell"+Label, Mmumu1, weight, 680, 12., 80.);
    FillHist("Mmumu1_3l1b2j_OffShell"+Label, Mmumu1, weight, 400, 100., 500.);
    FillHist("Mmumu2_3l1b2j_OnShell"+Label, Mmumu2, weight, 680, 12., 80.);
    FillHist("Mmumu2_3l1b2j_OffShell"+Label, Mmumu2, weight, 400, 100., 500.);
  }
  FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  FillHist("NPV"   +Label, nPV, weight, 80, 0., 80.);
  FillHist("Rho"   +Label, Rho, weight, 70, 0., 70.);

} 


void TestRun::AnalyzeElectronDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                    std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  if( !(ElTColl.size()==1 && ElLColl.size()==1) ) return;
  if( !(MuTColl.size()==2 && MuLColl.size()==2) ) return;
  if( !( (ElTColl.at(0).Pt()>25 && MuTColl.at(1).Pt()>10)
      or (ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>25 && MuTColl.at(1).Pt()>10) ) ) return;
  if( MuTColl.at(0).Charge() == MuTColl.at(1).Charge() ) return;
  if( !(fabs(ElTColl.at(0).Eta())<2.5) ) return;
  FillHist("CutFlow"+Label, 3., weight, 10, 0., 10.);
  float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
  if(Mmumu<12) return;
  FillHist("CutFlow"+Label, 4., weight, 10, 0., 10.);
  if(fabs(Mmumu-91.2)<10) return;
  FillHist("CutFlow"+Label, 5., weight, 10, 0., 10.);

  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, 6., weight, 10, 0., 10.);

  FillHist("Mmumu_3l1b1j_OnShell"+Label, Mmumu, weight, 680, 12., 80.);
  FillHist("Mmumu_3l1b1j_OffShell"+Label, Mmumu, weight, 400, 100., 500.);
  if(JetColl.size()>1){
    FillHist("Mmumu_3l1b2j_OnShell"+Label, Mmumu, weight, 680, 12., 80.);
    FillHist("Mmumu_3l1b2j_OffShell"+Label, Mmumu, weight, 400, 100., 500.);
  }
  FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  FillHist("NPV"   +Label, nPV, weight, 80, 0., 80.);
  FillHist("Rho"   +Label, Rho, weight, 70, 0., 70.);

}



void TestRun::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

TestRun::TestRun(){

}

TestRun::~TestRun(){

}


