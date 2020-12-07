#include "MeasTrigEff.h"

void MeasTrigEff::initializeAnalyzer(){

  ElEl=false, MuMu=false, ElMu=false, SystRun=false; 
  DiMuTrig_MuLeg=false, DiMuTrig_DZ=false, EMuTrig_ElLeg=false, EMuTrig_MuLeg=false, EMuTrig_DZ=false;
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("ElEl")) ElEl=true; 
    if(Userflags.at(i).Contains("MuMu")) MuMu=true; 
    if(Userflags.at(i).Contains("ElMu")) ElMu=true; 
    if(Userflags.at(i).Contains("DiMuTrig_MuLeg")) DiMuTrig_MuLeg=true; 
    if(Userflags.at(i).Contains("DiMuTrig_DZ"))    DiMuTrig_DZ=true; 
    if(Userflags.at(i).Contains("EMuTrig_ElLeg"))  EMuTrig_ElLeg=true; 
    if(Userflags.at(i).Contains("EMuTrig_MuLeg"))  EMuTrig_MuLeg=true; 
    if(Userflags.at(i).Contains("EMuTrig_DZ"))     EMuTrig_DZ=true; 
    if(Userflags.at(i).Contains("SystRun"))        SystRun=true; 
  }

  if(MuMu){
    if(DiMuTrig_MuLeg){ TrigList.push_back(DataYear==2017? "HLT_IsoMu27_v":"HLT_IsoMu24_v"); }
    //if(DiMuTrig_DZ){ TrigList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");}
    if(DiMuTrig_DZ){ TrigList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"); TrigList.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");}
  }
  if(ElEl){
    if(DataYear==2016){ TrigList.push_back("HLT_Ele27_WPTight_Gsf_v"); SFKey_Trig="Ele27WPTight_POGMVAIsoWP90";}
    if(DataYear==2017){ TrigList.push_back("HLT_Ele32_WPTight_Gsf_v"); TrigList.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"); SFKey_Trig="Ele32WPTight1OR2_POGMVAIsoWP90";}
    if(DataYear==2018){ TrigList.push_back("HLT_Ele32_WPTight_Gsf_v"); SFKey_Trig="Ele32WPTight_POGMVAIsoWP90";}
  }
  if(ElMu){
    if(EMuTrig_ElLeg){ TrigList.push_back(DataYear==2017? "HLT_IsoMu27_v":"HLT_IsoMu24_v"); }
    if(EMuTrig_MuLeg){
      TrigList.push_back(DataYear>2016? "HLT_Ele32_WPTight_Gsf_v":"HLT_Ele27_WPTight_Gsf_v"); 
      if(DataYear==2017){ TrigList.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"); }
    }
    if(EMuTrig_DZ){ 
      if(DataYear==2017){TrigList.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"); TrigList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");}
    }
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void MeasTrigEff::executeEvent(){


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
  if((!IsDATA) and EventCand){
    if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    w_gen     = ev.MCweight();
    //w_filter  = GetGenFilterEffCorr();
    //w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    w_lumi    = weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    w_PU      = GetPileUpWeight(nPileUp, 0);
    w_prefire = GetPrefireWeight(0);
    //sf_muid   = GetMuonSF(muonTightColl, "TopHNID_TkMu", "ID");
    //sf_muiso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    //sf_elreco = GetElectronSF(electronTightColl, "", "Reco");
    //sf_elid   = GetElectronSF(electronTightColl, "POGMVAIsoWP90", "ID");
    //sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    //sf_trig   = mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig, "");
    //cout<<"w_gen:"<<w_gen<<" w_lumi:"<<w_lumi<<" w_PU:"<<w_PU<<" w_prefire:"<<w_prefire<<" sf_trig:"<<sf_trig<<endl;
    //cout<<"sf_mutk"<<sf_mutk<<" sf_muid:"<<sf_muid<<" sf_muiso:"<<sf_muiso<<" sf_elreco:"<<sf_elreco<<" sf_elid:"<<sf_elid<<" sf_btag:"<<sf_btag<<endl;
  }
  weight *= w_gen * w_filter * w_topptrw * w_lumi * w_PU * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

 
  if(MuMu){
    if(DiMuTrig_MuLeg) MeasEffDiMuTrig_MuLeg(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, weight, "");
    if(DiMuTrig_DZ) MeasEffDiMuTrig_DZ(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, weight, "");
  }
  if(ElMu){
    if(EMuTrig_ElLeg){
      MeasEffEMuTrig_ElLeg(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, weight, "");
    }
    if(EMuTrig_MuLeg){
      MeasEffEMuTrig_MuLeg(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, weight, "");
    }
    if(EMuTrig_DZ){
      MeasEffEMuTrig_DZ(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, weight, "");
    }
  }


}



void MeasTrigEff::MeasEffDiMuTrig_MuLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                        vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label){

  if( !(MuTColl.size()==2 && ElTColl.size()==0) ) return;
  if( !(MuLColl.size()==2 && ElLColl.size()==0) ) return;
  if( MuTColl.at(0).Pt()<(DataYear!=2017? 26:29) ) return;
  if( MuTColl.at(0).Charge()==MuTColl.at(1).Charge() ) return;
  if( MuTColl.at(0).DeltaR(MuTColl.at(1))<0.4 ) return;

  const int NPtBinEdges1=10, NPtBinEdges2=10, NfEtaBinEdges=9;
  double PtBinEdges1[NPtBinEdges1], PtBinEdges2[NPtBinEdges2];
  double PtBinEdges1_17[NPtBinEdges1] = {0.,15.,17.,20.,25.,30.,40.,50.,100.,200.};
  double PtBinEdges2_17[NPtBinEdges2] = {0.,5.,8.,10.,15.,20.,30.,50.,100.,200.};
  std::copy(PtBinEdges1_17, PtBinEdges1_17+NPtBinEdges1, PtBinEdges1); 
  std::copy(PtBinEdges2_17, PtBinEdges2_17+NPtBinEdges2, PtBinEdges2); 

  double fEtaBinEdges[NfEtaBinEdges] = {-2.4,-2.1,-1.2,-0.9,0.,0.9,1.2,2.1,2.4};
  double PTMu   = MuTColl.at(0).Pt();
  //double fEtaMu = fabs(MuTColl.at(0).Eta());
  double fEtaMu = MuTColl.at(0).Eta();

  vector<TString> TrigListToMeas1={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"};
  vector<TString> TrigListToMeas2={"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"};
  bool PassLeg1=false, PassLeg2=false;
  if(MCSample.Contains("Trig") or DataStream.Contains("Trig")){
    for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
      float dR = sqrt(pow(MuTColl.at(0).Eta()-HLTObject_eta->at(it_obj),2)+pow(MuTColl.at(0).Phi()-HLTObject_phi->at(it_obj),2));
      if(dR>0.2) continue;
      if(HLTObject_FiredPaths->at(it_obj).find(TrigListToMeas1.at(0))!=std::string::npos ){ PassLeg1=true; }
      if(HLTObject_FiredPaths->at(it_obj).find(TrigListToMeas2.at(0))!=std::string::npos ){ PassLeg2=true; }
      if(PassLeg1 && PassLeg2) break;
    }
  }
  else{
    PassLeg1=ev.PassTrigger(TrigListToMeas1); PassLeg2=ev.PassTrigger(TrigListToMeas2);
  }


  FillHist("NMu1_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
  FillHist("NMu1_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
  FillHist("NMu2_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
  FillHist("NMu2_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
  if(PassLeg1){
    FillHist("NMu1Trig_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
    FillHist("NMu1Trig_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
  }
  if(PassLeg2){
    FillHist("NMu2Trig_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
    FillHist("NMu2Trig_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
  }

  //Syst:QCD contamination 
  if(ElTColl.at(0).Pt()>40){
    FillHist("NMu1_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
    FillHist("NMu1_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
    FillHist("NMu2_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
    FillHist("NMu2_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
    if(PassLeg1){
      FillHist("NMu1Trig_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
      FillHist("NMu1Trig_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
    }
    if(PassLeg2){
      FillHist("NMu2Trig_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
      FillHist("NMu2Trig_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
    }
  }

}



void MeasTrigEff::MeasEffEMuTrig_MuLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                       vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label){

  if( !(ElTColl.size()==1 && MuTColl.size()==1) ) return;
  if( !(ElLColl.size()==1 && MuLColl.size()==1) ) return;
  if( ElTColl.at(0).Pt()<(DataYear<2017? 29:32) ) return;
  if( MuTColl.at(0).Charge()==ElTColl.at(0).Charge() ) return;
  if( MuTColl.at(0).DeltaR(ElTColl.at(0))<0.4 ) return;

  const int NPtBinEdges1=10, NPtBinEdges2=10, NfEtaBinEdges=4;
  double PtBinEdges1[NPtBinEdges1], PtBinEdges2[NPtBinEdges2];
  double PtBinEdges1_17[NPtBinEdges1] = {0.,10.,20.,23.,25.,30.,40.,50.,100.,200.};
  double PtBinEdges2_17[NPtBinEdges2] = {0.,5.,8.,10.,15.,20.,30.,50.,100.,200.};
  std::copy(PtBinEdges1_17, PtBinEdges1_17+NPtBinEdges1, PtBinEdges1); 
  std::copy(PtBinEdges2_17, PtBinEdges2_17+NPtBinEdges2, PtBinEdges2); 

  double fEtaBinEdges[NfEtaBinEdges] = {0.,0.8,1.6,2.4};
  double PTMu   = MuTColl.at(0).Pt();
  double fEtaMu = fabs(MuTColl.at(0).Eta());

  vector<TString> TrigListToMeas1={"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
  vector<TString> TrigListToMeas2={"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"};
  bool PassLeg1=false, PassLeg2=false;
  if(MCSample.Contains("Trig") or DataStream.Contains("Trig")){
    for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
      float dR = sqrt(pow(MuTColl.at(0).Eta()-HLTObject_eta->at(it_obj),2)+pow(MuTColl.at(0).Phi()-HLTObject_phi->at(it_obj),2));
      if(dR>0.2) continue;
      if(HLTObject_FiredPaths->at(it_obj).find(TrigListToMeas1.at(0))!=std::string::npos ){ PassLeg1=true; }
      if(HLTObject_FiredPaths->at(it_obj).find(TrigListToMeas2.at(0))!=std::string::npos ){ PassLeg2=true; }
      if(PassLeg1 && PassLeg2) break;
    }
  }
  else{
    PassLeg1=ev.PassTrigger(TrigListToMeas1); PassLeg2=ev.PassTrigger(TrigListToMeas2);
  }


  FillHist("NMu1_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
  FillHist("NMu1_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
  FillHist("NMu2_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
  FillHist("NMu2_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
  if(PassLeg1){
    FillHist("NMu1Trig_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
    FillHist("NMu1Trig_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
  }
  if(PassLeg2){
    FillHist("NMu2Trig_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
    FillHist("NMu2Trig_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
  }

  //Syst:QCD contamination 
  if(ElTColl.at(0).Pt()>40){
    FillHist("NMu1_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
    FillHist("NMu1_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
    FillHist("NMu2_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
    FillHist("NMu2_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
    if(PassLeg1){
      FillHist("NMu1Trig_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges1-1, PtBinEdges1);
      FillHist("NMu1Trig_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
    }
    if(PassLeg2){
      FillHist("NMu2Trig_AltTag_ALL_Pt_1D", PTMu, weight, NPtBinEdges2-1, PtBinEdges2);
      FillHist("NMu2Trig_AltTag_ALL_PtEta_2D", PTMu, fEtaMu, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
    }
  }

}



void MeasTrigEff::MeasEffEMuTrig_ElLeg(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                        vector<Jet>& JetColl,  vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label){
  if( !(MuTColl.size()==1 && ElTColl.size()==1) ) return;
  if( !(MuLColl.size()==1 && ElLColl.size()==1) ) return;
  if( MuTColl.at(0).Pt()<(DataYear!=2017? 26:29) ) return;
  if( MuTColl.at(0).Charge()==ElTColl.at(0).Charge() ) return;
  if( MuTColl.at(0).DeltaR(ElTColl.at(0))<0.4 ) return;


  const int NPtBinEdges1=10, NPtBinEdges2=9, NfEtaBinEdges=4;
  double PtBinEdges1[NPtBinEdges1], PtBinEdges2[NPtBinEdges2];
  double PtBinEdges1_17[NPtBinEdges1] = {0.,10.,20.,23.,25.,30.,40.,50.,100.,200.};
  double PtBinEdges2_17[NPtBinEdges2] = {0.,10.,12.,15.,20.,30.,50.,100.,200.};
  std::copy(PtBinEdges1_17, PtBinEdges1_17+NPtBinEdges1, PtBinEdges1); 
  std::copy(PtBinEdges2_17, PtBinEdges2_17+NPtBinEdges2, PtBinEdges2); 

  double fEtaBinEdges[NfEtaBinEdges]={0., 0.8, 1.479, 2.5};
  double PTEle   = ElTColl.at(0).Pt();
  double fEtaEle = fabs(ElTColl.at(0).Eta());

  vector<TString> TrigListToMeas1={"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"};
  vector<TString> TrigListToMeas2={"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
  bool PassLeg1=false, PassLeg2=false;
  if(MCSample.Contains("Trig") or DataStream.Contains("Trig")){
    for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
      float dR = sqrt(pow(ElTColl.at(0).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElTColl.at(0).Phi()-HLTObject_phi->at(it_obj),2));
      if(dR>0.2) continue;
      if(HLTObject_FiredPaths->at(it_obj).find(TrigListToMeas1.at(0))!=std::string::npos ){ PassLeg1=true; }
      if(HLTObject_FiredPaths->at(it_obj).find(TrigListToMeas2.at(0))!=std::string::npos ){ PassLeg2=true; }
      if(PassLeg1 && PassLeg2) break;
    }
  }
  else{
    PassLeg1=ev.PassTrigger(TrigListToMeas1); PassLeg2=ev.PassTrigger(TrigListToMeas2);
  }


  FillHist("NEle1_ALL_Pt_1D", PTEle, weight, NPtBinEdges1-1, PtBinEdges1);
  FillHist("NEle1_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
  FillHist("NEle2_ALL_Pt_1D", PTEle, weight, NPtBinEdges2-1, PtBinEdges2);
  FillHist("NEle2_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
  if(PassLeg1){
    FillHist("NEle1Trig_ALL_Pt_1D", PTEle, weight, NPtBinEdges1-1, PtBinEdges1);
    FillHist("NEle1Trig_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
  }
  if(PassLeg2){
    FillHist("NEle2Trig_ALL_Pt_1D", PTEle, weight, NPtBinEdges2-1, PtBinEdges2);
    FillHist("NEle2Trig_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
  }

  //Syst:QCD contamination 
  if(MuTColl.at(0).Pt()>35 && MuTColl.at(0).RelIso()<0.15){
    FillHist("NEle1_AltTag_ALL_Pt_1D", PTEle, weight, NPtBinEdges1-1, PtBinEdges1);
    FillHist("NEle1_AltTag_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
    FillHist("NEle2_AltTag_ALL_Pt_1D", PTEle, weight, NPtBinEdges2-1, PtBinEdges2);
    FillHist("NEle2_AltTag_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
    if(PassLeg1){
      FillHist("NEle1Trig_AltTag_ALL_Pt_1D", PTEle, weight, NPtBinEdges1-1, PtBinEdges1);
      FillHist("NEle1Trig_AltTag_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges1-1, PtBinEdges1, NfEtaBinEdges-1, fEtaBinEdges);
    }
    if(PassLeg2){
      FillHist("NEle2Trig_AltTag_ALL_Pt_1D", PTEle, weight, NPtBinEdges2-1, PtBinEdges2);
      FillHist("NEle2Trig_AltTag_ALL_PtEta_2D", PTEle, fEtaEle, weight, NPtBinEdges2-1, PtBinEdges2, NfEtaBinEdges-1, fEtaBinEdges);
    }
  }

}


void MeasTrigEff::MeasEffEMuTrig_DZ(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                    vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)
{
  if( !(MuTColl.size()==1 && MuLColl.size()==1) ) return;
  if( !(ElTColl.size()==1 && ElLColl.size()==1) ) return;
  if( !(ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>10) ) return;
  if( !(ElTColl.at(0).Pt()>25 || MuTColl.at(0).Pt()>25) ) return;
  if(  ElTColl.at(0).Charge() == MuTColl.at(0).Charge() ) return;
  if(  ElTColl.at(0).DeltaR(MuTColl.at(0))<0.4  ) return;
  
  vector<TString> TrigListDen1    = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"};
  vector<TString> TrigListDen2    = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
  vector<TString> TrigListToMeas1 = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"};
  vector<TString> TrigListToMeas2 = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};

  bool IsCase1 = ElTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && ev.PassTrigger(TrigListDen1);
  bool IsCase2 = ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>25 && ev.PassTrigger(TrigListDen2);

  int NMatchIsoFilt=0, NMatchDZFilt=0;
  if(IsCase1){ NMatchIsoFilt=2; NMatchDZFilt=ev.PassTrigger(TrigListToMeas1)? 2:0;}
  if(IsCase2){ NMatchIsoFilt=2; NMatchDZFilt=ev.PassTrigger(TrigListToMeas2)? 2:0;}

  if(NMatchIsoFilt!=2) return;

  if(IsCase1){
    TString Label1=Label+"_HLT1";
    FillHist("NEvt"     +Label1, 0., weight, 1, 0., 1.);
    FillHist("NEvt_DZ1" +Label1, fabs(MuTColl.at(0).dZ()), weight, 20, 0., 0.1);
    FillHist("NEvt_DZ2" +Label1, fabs(ElTColl.at(0).dZ()), weight, 20, 0., 0.1);
    FillHist("NEvt_DZ12"+Label1, fabs(MuTColl.at(0).dZ()-ElTColl.at(0).dZ()), weight, 40, 0., 0.2);
    if(NMatchDZFilt==2){
      FillHist("NEvtTrig"     +Label1, 0., weight, 1, 0., 1.);
      FillHist("NEvtTrig_DZ1" +Label1, fabs(MuTColl.at(0).dZ()), weight, 20, 0., 0.1);
      FillHist("NEvtTrig_DZ2" +Label1, fabs(ElTColl.at(0).dZ()), weight, 20, 0., 0.1);
      FillHist("NEvtTrig_DZ12"+Label1, fabs(MuTColl.at(0).dZ()-ElTColl.at(0).dZ()), weight, 40, 0., 0.2);
    }
  }
  if(IsCase2){
    TString Label2=Label+"_HLT2";
    FillHist("NEvt"     +Label2, 0., weight, 1, 0., 1.);
    FillHist("NEvt_DZ1" +Label2, fabs(MuTColl.at(0).dZ()), weight, 20, 0., 0.1);
    FillHist("NEvt_DZ2" +Label2, fabs(ElTColl.at(0).dZ()), weight, 20, 0., 0.1);
    FillHist("NEvt_DZ12"+Label2, fabs(MuTColl.at(0).dZ()-ElTColl.at(0).dZ()), weight, 40, 0., 0.2);
    if(NMatchDZFilt==2){
      FillHist("NEvtTrig"     +Label2, 0., weight, 1, 0., 1.);
      FillHist("NEvtTrig_DZ1" +Label2, fabs(MuTColl.at(0).dZ()), weight, 20, 0., 0.1);
      FillHist("NEvtTrig_DZ2" +Label2, fabs(ElTColl.at(0).dZ()), weight, 20, 0., 0.1);
      FillHist("NEvtTrig_DZ12"+Label2, fabs(MuTColl.at(0).dZ()-ElTColl.at(0).dZ()), weight, 40, 0., 0.2);
    }
  }

} 




void MeasTrigEff::MeasEffDiMuTrig_DZ(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                      std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)
{
  if( !(MuTColl.size()==2 && MuLColl.size()==2) ) return;
  if( ElLColl.size()!=0 ) return;
  if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) return;
  if( MuTColl.at(0).Charge() == MuTColl.at(1).Charge() ) return;
  
  float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
  if(Mmumu<10) return;
  //if(Mmumu<4) return;

  int NMatchIsoFilt=0, NMatchDZFilt=0;
  if(MCSample.Contains("Trig") or DataStream.Contains("Trig")){
    for(unsigned int it_m=0; it_m<MuTColl.size(); it_m++){
      for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
        float dR = sqrt(pow(MuTColl.at(it_m).Eta()-HLTObject_eta->at(it_obj),2)+pow(MuTColl.at(it_m).Phi()-HLTObject_phi->at(it_obj),2));
        if(dR>0.2) continue;
        if(HLTObject_FiredPaths->at(it_obj).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")!=std::string::npos ){ NMatchIsoFilt++; break; }
      }
    }
    for(unsigned int it_m=0; it_m<MuTColl.size(); it_m++){
      for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
        float dR = sqrt(pow(MuTColl.at(it_m).Eta()-HLTObject_eta->at(it_obj),2)+pow(MuTColl.at(it_m).Phi()-HLTObject_phi->at(it_obj),2));
        if(dR>0.2) continue;
        if(HLTObject_FiredPaths->at(it_obj).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")!=std::string::npos ){ NMatchDZFilt++; break; }
      }
    }
  }
  else{NMatchIsoFilt=2; NMatchDZFilt= ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v")? 2:0;}
  //else{NMatchIsoFilt=2; NMatchDZFilt= ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v")? 2:0;}
  //else{NMatchIsoFilt=2; NMatchDZFilt= ev.PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")? 2:0;}

  if(NMatchIsoFilt!=2) return;
  FillHist("NEvt"     +Label, 0., weight, 1, 0., 1.);
  FillHist("NEvt_DZ1" +Label, fabs(MuTColl.at(0).dZ()), weight, 20, 0., 0.1);
  FillHist("NEvt_DZ2" +Label, fabs(MuTColl.at(1).dZ()), weight, 20, 0., 0.1);
  FillHist("NEvt_DZ12"+Label, fabs(MuTColl.at(0).dZ()-MuTColl.at(1).dZ()), weight, 40, 0., 0.2);
  if(NMatchDZFilt==2){
    FillHist("NEvtTrig"     +Label, 0., weight, 1, 0., 1.);
    FillHist("NEvtTrig_DZ1" +Label, fabs(MuTColl.at(0).dZ()), weight, 20, 0., 0.1);
    FillHist("NEvtTrig_DZ2" +Label, fabs(MuTColl.at(1).dZ()), weight, 20, 0., 0.1);
    FillHist("NEvtTrig_DZ12"+Label, fabs(MuTColl.at(0).dZ()-MuTColl.at(1).dZ()), weight, 40, 0., 0.2);
  }

} 



void MeasTrigEff::MeasSiglEleTrigEff(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
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





void MeasTrigEff::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

MeasTrigEff::MeasTrigEff(){

}

MeasTrigEff::~MeasTrigEff(){

}


