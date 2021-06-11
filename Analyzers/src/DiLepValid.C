#include "DiLepValid.h"

void DiLepValid::initializeAnalyzer(){

  ElEl=false, MuMu=false, ElMu=false, SystRun=false; 
  TrigEffMeas=false, HEMCheck=false;
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("ElEl")) ElEl=true; 
    if(Userflags.at(i).Contains("MuMu")) MuMu=true; 
    if(Userflags.at(i).Contains("ElMu")) ElMu=true; 
    if(Userflags.at(i).Contains("TrigEffMeas")) TrigEffMeas=true;
    if(Userflags.at(i).Contains("HEMCheck")) HEMCheck=true;
    if(Userflags.at(i).Contains("SystRun")) SystRun=true; 
  }

  if(MuMu){
    //if(DataYear==2016){ TrigList.push_back("HLT_IsoMu24_v"); TrigList.push_back("HLT_IsoTkMu24_v"); SFKey_Trig="IsoORTkIsoMu24_POGTight";}
    //if(DataYear==2017){ TrigList.push_back("HLT_IsoMu27_v"); SFKey_Trig="IsoMu27_POGTight";}
    //if(DataYear==2018){ TrigList.push_back("HLT_IsoMu24_v"); SFKey_Trig="IsoMu24_POGTight";}
    if(DataYear==2017){ TrigList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
                        TrigList.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
                        TrigList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"); SFKey_Trig="DiMuIso_HNTopID";}
  }
  if(ElEl){
//    if(DataYear==2016){ TrigList.push_back("HLT_Ele27_WPTight_Gsf_v"); SFKey_Trig="Ele27WPTight_POGMVAIsoWP90";}
//    if(DataYear==2017){ TrigList.push_back("HLT_Ele32_WPTight_Gsf_v"); TrigList.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"); SFKey_Trig="Ele32WPTight1OR2_POGMVAIsoWP90";}
//    if(DataYear==2018){ TrigList.push_back("HLT_Ele32_WPTight_Gsf_v"); SFKey_Trig="Ele32WPTight_POGMVAIsoWP90";}
    if(DataYear>=2017){ TrigList.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"); SFKey_Trig="DiElIso_HNTopIDSS"; }
  }
  if(ElMu){
//    if(DataYear==2016){ TrigList.push_back("HLT_IsoMu24_v"); TrigList.push_back("HLT_IsoTkMu24_v"); SFKey_Trig="IsoORTkIsoMu24_POGTight";}
//    if(DataYear==2017){ TrigList.push_back("HLT_IsoMu27_v"); SFKey_Trig="IsoMu27_POGTight";}
//    if(DataYear==2018){ TrigList.push_back("HLT_IsoMu24_v"); SFKey_Trig="IsoMu24_POGTight";}
    if(DataYear>=2017){
      TrigList.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TrigList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      SFKey_Trig="EMuIso_HNTopIDSS"; }
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void DiLepValid::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  if(!ev.PassTrigger(TrigList)) return;
  if(!PassMETFilter()) return;


  bool PreCutPass=false;
  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  std::sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  std::sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);
  if(ElEl and electronPreColl.size()>1) PreCutPass=true;
  if(MuMu and muonPreColl.size()>1    ) PreCutPass=true;
  if(ElMu and electronPreColl.size()>0 && muonPreColl.size()>0) PreCutPass=true;
  if(!PreCutPass) return;


  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TopHN17T", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TopHN17SST", 10., 2.5);
  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TopHN17SSL", 10., 2.5);


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
  if(MuMu){ EventCand = muonLooseColl.size()>1; }
  if(ElEl){ EventCand = electronLooseColl.size()>1; }
  if(ElMu){ EventCand = electronLooseColl.size()>0 && muonLooseColl.size()>0; }

  float w_gen = 1., w_filter = 1., w_topptrw = 1., w_lumi = 1., w_PU = 1., w_prefire = 1., sf_trig = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
  float pogsf=1., pogsf2=1.;
  if((!IsDATA) and EventCand){
    if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    w_gen     = ev.MCweight();
    //w_filter  = GetGenFilterEffCorr();
    w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    w_lumi    = weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    w_PU      = GetPileUpWeight(nPileUp, 0);
    w_prefire = GetPrefireWeight(0);
    sf_muid   = GetMuonSF(muonTightColl, "TopHNID_TkMu", "ID");
    //pogsf     = GetMuonSF(muonTightColl, "POGMID_genTrk", "ID")*GetMuonSF(muonTightColl, "POGTIso_POGMID", "Iso");
    //pogsf2    = GetMuonSF(muonTightColl, "POGTID_genTrk", "ID")*GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    sf_elreco = GetElectronSF(electronTightColl, "", "Reco");
    sf_elid   = GetElectronSF(electronTightColl, "TopHNIDSS", "ID");
    sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    sf_trig   = mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig, "");
    //cout<<"w_gen:"<<w_gen<<" w_lumi:"<<w_lumi<<" w_PU:"<<w_PU<<" w_prefire:"<<w_prefire<<" sf_trig:"<<sf_trig<<endl;
    //cout<<"sf_mutk"<<sf_mutk<<" sf_muid:"<<sf_muid<<" sf_muiso:"<<sf_muiso<<" sf_elreco:"<<sf_elreco<<" sf_elid:"<<sf_elid<<" sf_btag:"<<sf_btag<<endl;
  }
  weight *= w_gen * w_filter * w_topptrw * w_lumi * w_PU * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

 
  if(MuMu){
    if(HEMCheck){
      CheckHEMIssue(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET_xyCorr, weight, "", "");
    }
    else{
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, weight, "");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi, "_NoW");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire, "_Pref");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU, "_PrefPU");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*sf_muid*sf_muiso, "_PrefPUID");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*sf_muid*sf_muiso*sf_trig, "_PrefPUIDTrig");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*pogsf*sf_trig, "_PrefPUIDTrig_POG");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*pogsf2*sf_trig, "_PrefPUIDTrig_POG2");
      AnalyzeDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*sf_trig, "_PrefPUTrig");
    }
  }
  if(ElEl){
    if(HEMCheck){
      CheckHEMIssue(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET_xyCorr, weight, "", "");
    }
    else{
      AnalyzeDiElectron(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, weight, "");
      AnalyzeDiElectron(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi, "_NoW");
      AnalyzeDiElectron(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire, "_Pref");
      AnalyzeDiElectron(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU, "_PrefPU");
      AnalyzeDiElectron(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*sf_elreco, "_PrefPUReco");
      AnalyzeDiElectron(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*sf_elreco*sf_elid, "_PrefPURecoID");
      AnalyzeDiElectron(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, w_gen*w_lumi*w_prefire*w_PU*sf_elreco*sf_elid*sf_trig, "_PrefPURecoIDTrig");
    }
  }
  if(ElMu){
    if(TrigEffMeas){
      MeasSiglEleTrigEff(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, w_gen*w_lumi*w_prefire*w_PU, "");
    }
    else if(HEMCheck){
      CheckHEMIssue(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET_xyCorr,
                    weight, "", "");
    }
    else{
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET_xyCorr,
                          weight, "");
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET,
                          w_gen*w_lumi, "_NoW");
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET,
                          w_gen*w_lumi*w_prefire, "_Pref");
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET,
                          w_gen*w_lumi*w_prefire*w_PU, "_PrefPU");
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET,
                          w_gen*w_lumi*w_prefire*w_PU*sf_muid*sf_muiso*sf_trig, "_PrefPUIDTrigMu");
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET,
                          w_gen*w_lumi*w_prefire*w_PU*sf_muid*sf_muiso*sf_elreco*sf_elid*sf_trig, "_PrefPUIDTrigMuEl");
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET_xyCorr,
                          w_gen*w_lumi*w_prefire*w_PU*sf_muid*sf_muiso*sf_elreco*sf_elid*sf_trig, "_PrefPUIDTrigMuElMETPhi");
      AnalyzeElectronMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET_xyCorr,
                          w_gen*w_lumi*w_prefire*w_PU*sf_muid*sf_muiso*sf_elreco*sf_elid*sf_trig*sf_btag, "_PrefPUIDTrigMuElMETPhiBTag");
    }
  }

}


void DiLepValid::AnalyzeDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  if( !(MuTColl.size()==2 && MuLColl.size()==2) ) return;
  if( ElLColl.size()!=0 ) return;
  //if( DataYear==2016 or DataYear==2018 ){ if(MuTColl.at(0).Pt()<26) return; }
  //else{                                   if(MuTColl.at(0).Pt()<29) return; }
  if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) return;
  if( MuTColl.at(0).Charge() == MuTColl.at(1).Charge() ) return;
  
  float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
  if(Mmumu<50) return;
       
  FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  FillHist("NPV"   +Label, nPV, weight, 80, 0., 80.);
  FillHist("Rho"   +Label, Rho, weight, 70, 0., 70.);
  FillHist("Mmumu" +Label, Mmumu, weight, 60, 60., 120.);
  FillHist("PtMu1" +Label, MuTColl.at(0).Pt(), weight, 40, 0., 200.);
  FillHist("PtMu2" +Label, MuTColl.at(1).Pt(), weight, 40, 0., 200.);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
  FillHist("EtaMu2"+Label, MuTColl.at(1).Eta(), weight, 20, -5., 5.);
  FillHist("NJet"  +Label, JetColl.size(), weight, 10, 0., 10.);
  FillHist("NBJet" +Label, BJetColl.size(), weight, 10, 0., 10.);
  FillHist("MET"   +Label, vMET.Pt(), weight, 40, 0., 200.);

  float HT=0; for(unsigned int itj=0; itj<JetColl.size(); itj++){ HT+=JetColl.at(itj).Pt(); }
  float ST=HT+vMET.Pt()+MuTColl.at(0).Pt()+MuTColl.at(1).Pt();
  float MET2ST=pow(vMET.Pt(),2.)/ST;
  float MET2HT=HT==0? 0.:pow(vMET.Pt(),2.)/HT;
  FillHist("ST"+Label, ST, weight, 60, 0., 300.);
  FillHist("HT"+Label, HT, weight, 40, 0., 200.);
  FillHist("MET2ST"+Label, MET2ST, weight, 40, 0., 200.);
  FillHist("MET2HT"+Label, MET2HT, weight, 40, 0., 200.);

} 


void DiLepValid::AnalyzeDiElectron(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  if( !(ElTColl.size()==2 && ElLColl.size()==2) ) return;
  if( MuLColl.size()!=0 ) return;
  if( !(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15) ) return; 
  if( ElTColl.at(0).Charge() == ElTColl.at(1).Charge() ) return;
  //if( DataYear==2016 ){ if(ElTColl.at(0).Pt()<30) return; }
  //else                { if(ElTColl.at(0).Pt()<35) return; }

  float Melel = (ElTColl.at(0)+ElTColl.at(1)).M();
  if(Melel<50) return;

  FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  FillHist("NPV"   +Label, nPV, weight, 80, 0., 80.);
  FillHist("Rho"   +Label, Rho, weight, 70, 0., 70.);
  FillHist("Melel" +Label, Melel, weight, 60, 60., 120.);
  FillHist("PtEl1" +Label, ElTColl.at(0).Pt(), weight, 40, 0., 200.);
  FillHist("PtEl2" +Label, ElTColl.at(1).Pt(), weight, 40, 0., 200.);
  FillHist("EtaEl1"+Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
  FillHist("EtaEl2"+Label, ElTColl.at(1).Eta(), weight, 20, -5., 5.);
  FillHist("NJet"  +Label, JetColl.size(), weight, 10, 0., 10.);
  FillHist("NBJet" +Label, BJetColl.size(), weight, 10, 0., 10.);
  FillHist("MET"   +Label, vMET.Pt(), weight, 40, 0., 200.);

}


void DiLepValid::AnalyzeElectronMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                    std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)

{
  if( !(ElTColl.size()==1 && ElLColl.size()==1) ) return;
  if( !(MuTColl.size()==1 && MuLColl.size()==1) ) return;
  if( !(JetColl.size()>1 && BJetColl.size()>0) ) return;
  if( MuTColl.at(0).Charge() == ElTColl.at(0).Charge() ) return;
  if( MuTColl.at(0).DeltaR(ElTColl.at(0))<0.4 ) return;
  //if( DataYear==2016 or DataYear==2018 ){ if(MuTColl.at(0).Pt()<26) return; }
  //else{                                   if(MuTColl.at(0).Pt()<29) return; }
  if(!(MuTColl.at(0).Pt()>10 && ElTColl.at(0).Pt()>15)) return; 
  if(!(MuTColl.at(0).Pt()>25 || ElTColl.at(0).Pt()>25)) return; 
  //if(!(MuTColl.at(0).Pt()>25 && ElTColl.at(0).Pt()>15)) return; 
  //if(!(MuTColl.at(0).Pt()>10 && ElTColl.at(0).Pt()>25)) return; 


  FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
  FillHist("Rho"+Label, Rho, weight, 70, 0., 70.);
  FillHist("PtMu1" +Label, MuTColl.at(0).Pt(), weight, 30, 0., 300.);
  FillHist("PtEl1" +Label, ElTColl.at(0).Pt(), weight, 30, 0., 300.);
  FillHist("PtJet1"+Label, JetColl.at(0).Pt(), weight, 50, 0., 500.);
  FillHist("PtJet2"+Label, JetColl.at(1).Pt(), weight, 40, 0., 400.);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
  FillHist("EtaEl1"+Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
  FillHist("EtaJet1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
  FillHist("EtaJet2"+Label, JetColl.at(1).Eta(), weight, 20, -5., 5.);
  FillHist("MET"   +Label, vMET.Pt(), weight, 40, 0., 400.);
  FillHist("METPhi"+Label, vMET.Phi(), weight, 10, -3.15, 3.15);
  FillHist("NJet"  +Label, JetColl.size(), weight, 10, 0., 10.);
  FillHist("NBJet" +Label, BJetColl.size(), weight, 10, 0., 10.);
  FillHist("dRElMu"+Label, MuTColl.at(0).DeltaR(ElTColl.at(0)), weight, 50, 0., 5.);

  float HT=0; for(unsigned int itj=0; itj<JetColl.size(); itj++){ HT+=JetColl.at(itj).Pt(); }
  float ST=HT+vMET.Pt()+MuTColl.at(0).Pt()+ElTColl.at(0).Pt();
  float MET2ST=pow(vMET.Pt(),2.)/ST;
  float MET2HT=pow(vMET.Pt(),2.)/HT;
  FillHist("ST"+Label, ST, weight, 48, 0., 1200.);
  FillHist("HT"+Label, HT, weight, 32, 0., 800.);
  FillHist("MET2ST"+Label, MET2ST, weight, 12, 0., 120.);
  FillHist("MET2HT"+Label, MET2HT, weight, 30, 0., 300.);

}


void DiLepValid::CheckHEMIssue(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                               std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label, TString Option)
{

  if(DataYear!=2018) return;

  if(MuMu){
    if( !(MuTColl.size()==2 && MuLColl.size()==2) ) return;
    if( ElLColl.size()!=0 ) return;
    if( DataYear==2016 or DataYear==2018 ){ if(MuTColl.at(0).Pt()<26) return; }
    else{                                   if(MuTColl.at(0).Pt()<29) return; }
    if( MuTColl.at(0).Charge() == MuTColl.at(1).Charge() ) return;
    float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<50) return;

    const int NEtaEdges=5, NPhiEdges=9;
    double EtaEdges_p[NEtaEdges]={0., 0.6, 1.3, 1.9, 2.5};
    double PhiEdges[NPhiEdges]={-3.15, -2.4, -1.6, -0.8, 0., 0.8, 1.6, 2.4, 3.15};
    float mu2_feta=fabs(MuTColl.at(1).Eta()), mu2_phi = MuTColl.at(1).Phi();

    if(mu2_feta>0){ FillHist("Mu2_fEtaPhi_p"+Label, mu2_feta, mu2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }
    else{           FillHist("Mu2_fEtaPhi_m"+Label, mu2_feta, mu2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }

    if(IsDATA && IsHEMIssueRun()){
      if(mu2_feta>0){ FillHist("Mu2_fEtaPhi_p_HEMRun"+Label, mu2_feta, mu2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }
      else{           FillHist("Mu2_fEtaPhi_m_HEMRun"+Label, mu2_feta, mu2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }
    }
  }
  if(ElEl){
    if( !(ElTColl.size()==2 && ElLColl.size()==2) ) return;
    if( MuLColl.size()!=0 ) return;
    if( ElTColl.at(0).Charge() == ElTColl.at(1).Charge() ) return;
    if( DataYear==2016 ){ if(ElTColl.at(0).Pt()<30) return; }
    else                { if(ElTColl.at(0).Pt()<35) return; }
    if( !(fabs(ElTColl.at(0).Eta())<2.5 && fabs(ElTColl.at(1).Eta())<2.5) ) return;
    float Melel = (ElTColl.at(0)+ElTColl.at(1)).M();
    if(Melel<50) return;

    const int NEtaEdges=5, NPhiEdges=9;
    double EtaEdges_p[NEtaEdges]={0., 0.6, 1.3, 1.9, 2.5};
    double PhiEdges[NPhiEdges]={-3.15, -2.4, -1.6, -0.8, 0., 0.8, 1.6, 2.4, 3.15};
    float el2_feta=fabs(ElTColl.at(1).Eta()), el2_phi = ElTColl.at(1).Phi();
    if(el2_feta>0){ FillHist("El2_fEtaPhi_p"+Label, el2_feta, el2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }
    else{         FillHist("El2_fEtaPhi_m"+Label, el2_feta, el2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }
    if(IsDATA && IsHEMIssueRun()){
      if(el2_feta>0){ FillHist("El2_fEtaPhi_p_HEMRun"+Label, el2_feta, el2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }
      else{           FillHist("El2_fEtaPhi_m_HEMRun"+Label, el2_feta, el2_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges); }
    }
  }
  if(ElMu){
    if( !(ElTColl.size()==1 && ElLColl.size()==1) ) return;
    if( !(MuTColl.size()==1 && MuLColl.size()==1) ) return;
    if( !(JetColl.size()>1 && BJetColl.size()>0) ) return;
    if( MuTColl.at(0).Charge() == ElTColl.at(0).Charge() ) return;
    if( fabs(ElTColl.at(0).Eta())>2.5 ) return;
    if( MuTColl.at(0).DeltaR(ElTColl.at(0))<0.4 ) return;
    if( DataYear==2016 or DataYear==2018 ){ if(MuTColl.at(0).Pt()<26) return; }
    else{                                   if(MuTColl.at(0).Pt()<29) return; }

    const int NEtaEdges=5, NPhiEdges=9;
    double EtaEdges_p[NEtaEdges]={0., 0.6, 1.3, 1.9, 2.5};
    double PhiEdges[NPhiEdges]={-3.15, -2.4, -1.6, -0.8, 0., 0.8, 1.6, 2.4, 3.15};
    for(unsigned int it_j=0; it_j<JetColl.size(); it_j++){
      float j_pt = JetColl.at(it_j).Pt(), j_eta = JetColl.at(it_j).Eta();
      float j_phi = JetColl.at(it_j).Phi(), j_feta = fabs(j_eta);
      if( (IsDATA && IsHEMIssueRun()) or (!IsDATA) ){
        if(IsHEMIssueReg(JetColl.at(it_j))){
          FillHist("CHEF_HEMReg" +Label, JetColl.at(it_j).GetIDVariable("CHEF"), weight, 20, 0., 1.);
          FillHist("NHEF_HEMReg" +Label, JetColl.at(it_j).GetIDVariable("NHEF"), weight, 20, 0., 1.);
          FillHist("CEMEF_HEMReg"+Label, JetColl.at(it_j).GetIDVariable("CEMEF"), weight, 20, 0., 1.);
          FillHist("NEMEF_HEMReg"+Label, JetColl.at(it_j).GetIDVariable("NEMEF"), weight, 20, 0., 1.);
          FillHist("MuEF_HEMReg" +Label, JetColl.at(it_j).GetIDVariable("MuEF"), weight, 20, 0., 1.);
          FillHist("CMult_HEMReg"+Label, JetColl.at(it_j).GetIDVariable("CMult"), weight, 25, 0., 50.);
          FillHist("NMult_HEMReg"+Label, JetColl.at(it_j).GetIDVariable("NMult"), weight, 25, 0., 50.);
        }
        else if(IsHEMCRReg(JetColl.at(it_j), "SameEta") or IsHEMCRReg(JetColl.at(it_j), "SamePhi")){
          FillHist("CHEF_NonHEMReg" +Label, JetColl.at(it_j).GetIDVariable("CHEF"), weight, 20, 0., 1.);
          FillHist("NHEF_NonHEMReg" +Label, JetColl.at(it_j).GetIDVariable("NHEF"), weight, 20, 0., 1.);
          FillHist("CEMEF_NonHEMReg"+Label, JetColl.at(it_j).GetIDVariable("CEMEF"), weight, 20, 0., 1.);
          FillHist("NEMEF_NonHEMReg"+Label, JetColl.at(it_j).GetIDVariable("NEMEF"), weight, 20, 0., 1.);
          FillHist("MuEF_NonHEMReg" +Label, JetColl.at(it_j).GetIDVariable("MuEF"), weight, 20, 0., 1.);
          FillHist("CMult_NonHEMReg"+Label, JetColl.at(it_j).GetIDVariable("CMult"), weight, 25, 0., 50.);
          FillHist("NMult_NonHEMReg"+Label, JetColl.at(it_j).GetIDVariable("NMult"), weight, 25, 0., 50.);
        }
      }


      if(j_eta>0){ FillHist("Nj_fEtaPhi_p" +Label, j_feta, j_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                   FillHist("Ptj_fEtaPhi_p"+Label, j_feta, j_phi, weight*j_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
      }
      else{        FillHist("Nj_fEtaPhi_m" +Label, j_feta, j_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                   FillHist("Ptj_fEtaPhi_m"+Label, j_feta, j_phi, weight*j_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
      }
      if(IsDATA && IsHEMIssueRun()){
        if(j_eta>0){ FillHist("Nj_fEtaPhi_p_HEMRun"+Label, j_feta, j_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptj_fEtaPhi_p_HEMRun"+Label, j_feta, j_phi, weight*j_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
        else{        FillHist("Nj_fEtaPhi_m_HEMRun"+Label, j_feta, j_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptj_fEtaPhi_m_HEMRun"+Label, j_feta, j_phi, weight*j_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
      }
      else if(IsDATA && !IsHEMIssueRun()){
        if(j_eta>0){ FillHist("Nj_fEtaPhi_p_NonHEMRun"+Label, j_feta, j_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptj_fEtaPhi_p_NonHEMRun"+Label, j_feta, j_phi, weight*j_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
        else{        FillHist("Nj_fEtaPhi_m_NonHEMRun"+Label, j_feta, j_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptj_fEtaPhi_m_NonHEMRun"+Label, j_feta, j_phi, weight*j_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
      }
    }//End of jet loop
    for(unsigned int it_b=0; it_b<BJetColl.size(); it_b++){
      float b_pt=BJetColl.at(it_b).Pt(), b_eta=BJetColl.at(it_b).Eta();
      float b_feta=fabs(b_eta), b_phi = BJetColl.at(it_b).Phi();
      if(b_eta>0){ FillHist("Nb_fEtaPhi_p"+Label, b_feta, b_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                   FillHist("Ptb_fEtaPhi_p"+Label, b_feta, b_phi, weight*b_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
      }
      else{        FillHist("Nb_fEtaPhi_m"+Label, b_feta, b_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                   FillHist("Ptb_fEtaPhi_m"+Label, b_feta, b_phi, weight*b_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
      }
      if(IsDATA && IsHEMIssueRun()){
        if(b_eta>0){ FillHist("Nb_fEtaPhi_p_HEMRun"+Label, b_feta, b_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptb_fEtaPhi_p_HEMRun"+Label, b_feta, b_phi, weight*b_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
        else{        FillHist("Nb_fEtaPhi_m_HEMRun"+Label, b_feta, b_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptb_fEtaPhi_m_HEMRun"+Label, b_feta, b_phi, weight*b_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
      }
      else if(IsDATA && !IsHEMIssueRun()){
        if(b_eta>0){ FillHist("Nb_fEtaPhi_p_NonHEMRun"+Label, b_feta, b_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptb_fEtaPhi_p_NonHEMRun"+Label, b_feta, b_phi, weight*b_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
        else{        FillHist("Nb_fEtaPhi_m_NonHEMRun"+Label, b_feta, b_phi, weight, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
                     FillHist("Ptb_fEtaPhi_m_NonHEMRun"+Label, b_feta, b_phi, weight*b_pt, NEtaEdges-1, EtaEdges_p, NPhiEdges-1, PhiEdges);
        }
      }
    }//End of b-jet loop
  }//End of ElMu

}



void DiLepValid::MeasSiglEleTrigEff(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                    std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)

{
  if( !(ElTColl.size()==1 && ElLColl.size()==1) ) return;
  if( !(MuTColl.size()==1 && MuLColl.size()==1) ) return;
  if( MuTColl.at(0).Charge() == ElTColl.at(0).Charge() ) return;
  if( DataYear==2016 or DataYear==2018 ){ 
    if(MuTColl.at(0).Pt()<26) return;
  }
  else{
    if(MuTColl.at(0).Pt()<29) return;
  }
  if( ElTColl.at(0).DeltaR(MuTColl.at(0))<0.4 ) return;

  //Sanity check
  FillHist("PtMu1" +Label, MuTColl.at(0).Pt(), weight, 30, 0., 300.);
  FillHist("PtEl1" +Label, ElTColl.at(0).Pt(), weight, 30, 0., 300.);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
  FillHist("EtaEl1"+Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
  FillHist("MET"   +Label, vMET.Pt(), weight, 40, 0., 400.);


  vector<TString> TrigListToMeas;
  if(DataYear==2016){ TrigListToMeas.push_back("HLT_Ele27_WPTight_Gsf_v"); }
  if(DataYear==2017){ TrigListToMeas.push_back("HLT_Ele32_WPTight_Gsf_v"); TrigListToMeas.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"); }
  if(DataYear==2018){ TrigListToMeas.push_back("HLT_Ele32_WPTight_Gsf_v"); }

  const int NPtBinEdges=11, NfEtaBinEdges=4;
  double PtBinEdges[NPtBinEdges];
  double PtBinEdges_16[NPtBinEdges]   = {0., 15., 24., 27., 30., 35., 45., 70., 100., 200., 500};
  double PtBinEdges_1718[NPtBinEdges] = {0., 20., 29., 32., 35., 40., 50., 70., 100., 200., 500};
  if(DataYear==2016){ std::copy(PtBinEdges_16, PtBinEdges_16+NPtBinEdges, PtBinEdges); }
  else              { std::copy(PtBinEdges_1718, PtBinEdges_1718+NPtBinEdges, PtBinEdges); }

  double fEtaBinEdges[NfEtaBinEdges]={0., 0.8, 1.479, 2.5};
  double PTEle   = ElTColl.at(0).Pt();
  double fEtaEle = fabs(ElTColl.at(0).Eta());

  FillHist("NEle_ALL_Pt_1D", PTEle, weight, NPtBinEdges-1, PtBinEdges);
  FillHist("NEle_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
  if(ev.PassTrigger(TrigListToMeas)){
    FillHist("NEleTrig_ALL_Pt_1D", PTEle, weight, NPtBinEdges-1, PtBinEdges);
    FillHist("NEleTrig_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
  }

  //Syst:QCD contamination 
  if(MuTColl.at(0).Pt()>35 && MuTColl.at(0).RelIso()<0.1){
    FillHist("NEle_AltTag_ALL_Pt_1D", PTEle, weight, NPtBinEdges-1, PtBinEdges);
    FillHist("NEle_AltTag_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
    if(ev.PassTrigger(TrigListToMeas)){
      FillHist("NEleTrig_AltTag_ALL_Pt_1D", PTEle, weight, NPtBinEdges-1, PtBinEdges);
      FillHist("NEleTrig_AltTag_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
    }
  }
  
}





void DiLepValid::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

DiLepValid::DiLepValid(){

}

DiLepValid::~DiLepValid(){

}


