#include "HNTopFeas.h"

void HNTopFeas::initializeAnalyzer(){

  ElMuMu=false, MuMuMu=false, SSMuMu=false, SSElEl=false, SSElMu=false;
  TriLep=false, TetraLep=false, SS2l=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("ElMuMu")) ElMuMu=true; 
    if(Userflags.at(i).Contains("MuMuMu")) MuMuMu=true; 
    if(Userflags.at(i).Contains("SSMuMu")) SSMuMu=true;
    if(Userflags.at(i).Contains("SSElEl")) SSElEl=true;
    if(Userflags.at(i).Contains("SSElMu")) SSElMu=true;
    if(Userflags.at(i).Contains("TriLep")) TriLep=true;
    if(Userflags.at(i).Contains("SS2l"))   SS2l=true;
    if(Userflags.at(i).Contains("TetraLep")) TetraLep=true; 
    if(Userflags.at(i).Contains("SystRun")) SystRun=true; 
  }

  DblMu=false, DblEG=false, MuEG=false;
  if     (DataStream.Contains("DoubleMuon")) DblMu=true;
  else if(DataStream.Contains("MuonEG"))     MuEG =true;
  else if(DataStream.Contains("DoubleEG"))   DblEG=true;
  else if(DataYear==2018 and DataStream.Contains("EGamma")) DblEG=true;

  if(DataYear==2016){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  if(DataYear==2017){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");

    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  }
  if(DataYear==2018){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");

    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void HNTopFeas::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;
  if(!IsDATA){
    weight*=ev.MCweight()*weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    weight *= GetGenFilterEffCorr();
    weight*=GetPileUpWeight(nPileUp, 0);
  }
  FillHist("CutFlow", 0., weight, 20, 0., 20.);

  bool PassTrig=false;
  if     (MuMuMu){ PassTrig = ev.PassTrigger(TrigList_DblMu); }
  else if(ElMuMu){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu);
    else{
      if     (MuEG)  PassTrig = ev.PassTrigger(TrigList_MuEG);
      else if(DblMu) PassTrig = (!ev.PassTrigger(TrigList_MuEG)) and ev.PassTrigger(TrigList_DblMu);
    }
  }
  else if(TetraLep or TriLep){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu) or ev.PassTrigger(TrigList_DblEG);
    else{
      if     (MuEG)  PassTrig = ev.PassTrigger(TrigList_MuEG);
      else if(DblMu) PassTrig = (!ev.PassTrigger(TrigList_MuEG)) and ev.PassTrigger(TrigList_DblMu);
      else if(DblEG) PassTrig = (!(ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu))) and ev.PassTrigger(DblEG);
    }
  }
  else if(SS2l){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_DblMu) or ev.PassTrigger(TrigList_DblEG) or ev.PassTrigger(TrigList_MuEG);
    else{
      if     (SSMuMu) PassTrig = ev.PassTrigger(TrigList_DblMu);
      else if(SSElEl) PassTrig = ev.PassTrigger(TrigList_DblEG);
      else if(SSElMu) PassTrig = ev.PassTrigger(TrigList_MuEG);
    }
  }
  if(!PassTrig) return;
  FillHist("CutFlow", 1., weight, 20, 0., 20.);
  if(!PassMETFilter()) return;
  FillHist("CutFlow", 2., weight, 20, 0., 20.);

  bool PreCutPass=false;
  std::vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  std::sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  std::sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);
  if(MuMuMu and muonPreColl.size()>2    ) PreCutPass=true;
  if(ElMuMu and electronPreColl.size()>0 && muonPreColl.size()>1) PreCutPass=true;
  if(TriLep and (muonPreColl.size()+electronPreColl.size())>2 ) PreCutPass=true;
  if(TetraLep and (muonPreColl.size()+electronPreColl.size())>3 ) PreCutPass=true;
  if(SS2l){
    int NPreMu=muonPreColl.size(), NPreEl=electronPreColl.size();
    if( (NPreMu+NPreEl)>2 ) PreCutPass=true;
    else if(NPreMu==2 and SumCharge(muonPreColl)!=0) PreCutPass=true;
    else if(NPreEl==2 and SumCharge(electronPreColl)!=0) PreCutPass=true;
  }
  if(!PreCutPass) return;


  TString IDSSLabel = SS2l? "SS":"";
  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TopHN17T", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"T", 10., 2.5);
  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"L", 10., 2.5);


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
  if(SS2l){ EventCand = muonTightColl.size()==2 or electronTightColl.size()==2 or (muonTightColl.size()==1 && electronTightColl.size()==1); }
  if(TriLep){ EventCand = (muonTightColl.size()+electronTightColl.size())==3; }
  if(MuMuMu){ EventCand = muonLooseColl.size()>2; }
  if(ElMuMu){ EventCand = electronLooseColl.size()>0 && muonLooseColl.size()>1; }
  if(TetraLep){ EventCand = (muonLooseColl.size()+electronLooseColl.size())>3; }

  float w_topptrw = 1., w_prefire = 1., sf_trig = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
  if((!IsDATA) and EventCand){
    //if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    //w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    sf_muid   = GetMuonSF(muonTightColl, "TopHNID_TkMu", "ID");
    //sf_muiso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    sf_elreco = GetElectronSF(electronTightColl, "", "Reco");
    sf_elid   = GetElectronSF(electronTightColl, "TopHNID"+IDSSLabel, "ID");
    sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    TString SFKey_Trig_All = muonTightColl.size()==2? "DiMuIso_HNTopID":electronTightColl.size()==2? "DiElIso_HNTopIDSS":"EMuIso_HNTopIDSS"; 
    sf_trig   = SS2l? mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig_All, ""):1.;
    w_prefire = GetPrefireWeight(0);
    //cout<<"w_gen:"<<w_gen<<" w_lumi:"<<w_lumi<<" w_PU:"<<w_PU<<" w_prefire:"<<w_prefire<<" sf_trig:"<<sf_trig<<endl;
    //cout<<"sf_mutk"<<sf_mutk<<" sf_muid:"<<sf_muid<<" sf_muiso:"<<sf_muiso<<" sf_elreco:"<<sf_elreco<<" sf_elid:"<<sf_elid<<" sf_btag:"<<sf_btag<<endl;
  }
  weight *= w_topptrw * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

 
  if(SS2l){
    AnalyzeSSDiLepton(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                      jetColl, bjetColl, vMET_xyCorr, weight, "");
  }
  if(MuMuMu){
    AnalyzeTriMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                   jetColl, bjetColl, vMET_xyCorr, weight, "");
  }
  if(ElMuMu){
    AnalyzeElectronDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                          jetColl, bjetColl, vMET_xyCorr, weight, "");
  }
  if(TriLep){
    AnalyzeTriLepton(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                     jetColl, bjetColl, vMET_xyCorr, weight, "");
  }
  if(TetraLep){
    AnalyzeTetraLepton(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                       jetColl, bjetColl, vMET_xyCorr, weight, "");
  }

}


void HNTopFeas::AnalyzeSSDiLepton(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                  std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size();
  int it_cut=3;
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if(NMuT==2){
    Label+="_2M";
    if(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    float MSSSF = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(DataYear>2016 && MSSSF<4) return; 
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    FillHist("NBJet"+Label, BJetColl.size(), weight, 10, 0., 10.);
    if(BJetColl.size()<1) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  
    FillHist("NJet"+Label, JetColl.size(), weight, 10, 0., 10.); 
    if(JetColl.size()<3) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  
    FillHist("MSSSF"+Label, MSSSF, weight, 30, 0., 300.);
  
    float minSumdRql = -1., IdxJ_minSumdRql=-1;
    for(unsigned int it_j=0; it_j<JetColl.size(); it_j++){
      float SumdRql = MuTColl.at(0).DeltaR(JetColl.at(it_j))+MuTColl.at(1).DeltaR(JetColl.at(it_j));
      if(minSumdRql<0 or minSumdRql>SumdRql){ minSumdRql = SumdRql; IdxJ_minSumdRql=it_j; }
    }
    float dRll = MuTColl.at(0).DeltaR(MuTColl.at(1)), dEtall = abs(MuTColl.at(0).Eta()-MuTColl.at(1).Eta());
    FillHist("dRll"+Label, dRll, weight, 50, 0., 5.);
    FillHist("dEtall"+Label, dEtall, weight, 50, 0., 5.);
    FillHist("dPhill"+Label, MuTColl.at(0).DeltaPhi(MuTColl.at(1)), weight, 10, -3.15, 3.15);
    FillHist("minSumdRql"+Label, minSumdRql, weight, 50, 0., 10.);
    FillHist("Ml1q_minSumdRql"+Label, (MuTColl.at(0)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);
    FillHist("Ml2q_minSumdRql"+Label, (MuTColl.at(1)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);
    FillHist("Mllq_minSumdRql"+Label, (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);

    float Pzv = GetvPz(MuTColl.at(0), vMET);
    Particle vReco; vReco.SetXYZM(vMET.Px(),vMET.Py(),Pzv,0.);
    FillHist("MET"+Label, vMET.Pt(), weight, 30, 0., 300.);
    FillHist("MT_vl1"+Label, MT(MuTColl.at(0),vMET), weight, 30, 0., 300.); 
    FillHist("Mlvb"+Label, (MuTColl.at(0)+vReco+BJetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllb"+Label, (MuTColl.at(0)+MuTColl.at(1)+BJetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllbq1"+Label, (MuTColl.at(0)+MuTColl.at(1)+BJetColl.at(0)+JetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllbq2"+Label, (MuTColl.at(0)+MuTColl.at(1)+BJetColl.at(0)+JetColl.at(1)).M(), weight, 50, 0., 500.);
    FillHist("Mllbq3"+Label, (MuTColl.at(0)+MuTColl.at(1)+BJetColl.at(0)+JetColl.at(2)).M(), weight, 50, 0., 500.);
    FillHist("Mllq1"+Label, (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllq2"+Label, (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(1)).M(), weight, 50, 0., 500.);
    FillHist("Mllq3"+Label, (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(2)).M(), weight, 50, 0., 500.);


    FillHist("NPresel", 0., weight, 12, 0., 12.);
    FillHist("NPresel_Tot", 0., weight, 3, 0., 3.);

    if(vMET.Pt()>100) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    if(MSSSF>70.) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
    FillHist("dRll_atMETM"+Label, dRll, weight, 50, 0., 5.);
    FillHist("dEtall_atMETM"+Label, dEtall, weight, 50, 0., 5.);

    if(dRll>3.2) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
    FillHist("dEtall_atMETMdR"+Label, dEtall, weight, 50, 0., 5.);

    if(dEtall>2.) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
    FillHist("Mllq_minSumdRql_atMETMdREta"+Label, (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);

    FillHist("NBJet_atMETMdREta"+Label, BJetColl.size(), weight, 10, 0., 10.);
    FillHist("NJet_atMETMdREta"+Label, JetColl.size(), weight, 10, 0., 10.);
    FillHist("MSSSF_2Mu_atMETMdREta"+Label, MSSSF, weight, 30, 0., 300.);
  }
  else if(NElT==2){
    Label+="_2E";
    if(ElTColl.at(0).Charge()!=ElTColl.at(1).Charge()) return;
    if(!(ElTColl.at(0).Pt()>25 and ElTColl.at(1).Pt()>15)) return;
    float MSSSF = (ElTColl.at(0)+ElTColl.at(1)).M();
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    FillHist("NBJet"+Label, BJetColl.size(), weight, 10, 0., 10.);
    if(BJetColl.size()<1) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  
    FillHist("NJet"+Label, JetColl.size(), weight, 10, 0., 10.); 
    if(JetColl.size()<3) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  
    FillHist("MSSSF"+Label, MSSSF, weight, 30, 0., 300.);
  
    float minSumdRql = -1., IdxJ_minSumdRql=-1;
    for(unsigned int it_j=0; it_j<JetColl.size(); it_j++){
      float SumdRql = ElTColl.at(0).DeltaR(JetColl.at(it_j))+ElTColl.at(1).DeltaR(JetColl.at(it_j));
      if(minSumdRql<0 or minSumdRql>SumdRql){ minSumdRql = SumdRql; IdxJ_minSumdRql=it_j; }
    }
    float dRll = ElTColl.at(0).DeltaR(ElTColl.at(1)), dEtall = abs(ElTColl.at(0).Eta()-ElTColl.at(1).Eta());
    FillHist("dRll"+Label, dRll, weight, 50, 0., 5.);
    FillHist("dEtall"+Label, dEtall, weight, 50, 0., 5.);
    FillHist("dPhill"+Label, ElTColl.at(0).DeltaPhi(ElTColl.at(1)), weight, 10, -3.15, 3.15);
    FillHist("minSumdRql"+Label, minSumdRql, weight, 50, 0., 10.);
    FillHist("Ml1q_minSumdRql"+Label, (ElTColl.at(0)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);
    FillHist("Ml2q_minSumdRql"+Label, (ElTColl.at(1)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);
    FillHist("Mllq_minSumdRql"+Label, (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);

    float Pzv = GetvPz(ElTColl.at(0), vMET);
    Particle vReco; vReco.SetXYZM(vMET.Px(),vMET.Py(),Pzv,0.);
    FillHist("MET"+Label, vMET.Pt(), weight, 30, 0., 300.);
    FillHist("MT_vl1"+Label, MT(ElTColl.at(0),vMET), weight, 30, 0., 300.); 
    FillHist("Mlvb"+Label, (ElTColl.at(0)+vReco+BJetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllb"+Label, (ElTColl.at(0)+ElTColl.at(1)+BJetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllbq1"+Label, (ElTColl.at(0)+ElTColl.at(1)+BJetColl.at(0)+JetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllbq2"+Label, (ElTColl.at(0)+ElTColl.at(1)+BJetColl.at(0)+JetColl.at(1)).M(), weight, 50, 0., 500.);
    FillHist("Mllbq3"+Label, (ElTColl.at(0)+ElTColl.at(1)+BJetColl.at(0)+JetColl.at(2)).M(), weight, 50, 0., 500.);
    FillHist("Mllq1"+Label, (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(0)).M(), weight, 50, 0., 500.);
    FillHist("Mllq2"+Label, (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(1)).M(), weight, 50, 0., 500.);
    FillHist("Mllq3"+Label, (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(2)).M(), weight, 50, 0., 500.);


    FillHist("NPresel", 1., weight, 12, 0., 12.);

    if(vMET.Pt()>100) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    if(MSSSF>70.) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
    FillHist("dRll_atMETM"+Label, dRll, weight, 50, 0., 5.);
    FillHist("dEtall_atMETM"+Label, dEtall, weight, 50, 0., 5.);

    if(dRll>3.2) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
    FillHist("dEtall_atMETMdR"+Label, dEtall, weight, 50, 0., 5.);

    if(dEtall>2.) return;
    FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
    FillHist("Mllq_minSumdRql_atMETMdREta"+Label, (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(IdxJ_minSumdRql)).M(), weight, 30, 0., 300.);

    FillHist("NBJet_atMETMdREta"+Label, BJetColl.size(), weight, 10, 0., 10.);
    FillHist("NJet_atMETMdREta"+Label, JetColl.size(), weight, 10, 0., 10.);
    FillHist("MSSSF_2El_atMETMdREta"+Label, MSSSF, weight, 30, 0., 300.);
  }

}



void HNTopFeas::AnalyzeTriMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  if( !(MuTColl.size()==3 && MuLColl.size()==3) ) return;
  if( ElLColl.size()!=0 ) return;
  if( !(MuTColl.at(0).Pt()>20. && MuTColl.at(2).Pt()>10) ) return;
  if( abs(SumCharge(MuTColl,ElTColl)) != 1 ) return;
  FillHist("CutFlow"+Label, 3., weight, 10, 0., 10.);
  
  
  int IdxOS  = TriMuChargeIndex(MuTColl, "OS");
  int IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
  int IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
  float Mmumu1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
  float Mmumu2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
  float MSSSF = (MuTColl.at(IdxSS1)+MuTColl.at(IdxSS2)).M();
  if( !(Mmumu1>12 && Mmumu2>12) ) return;
  FillHist("CutFlow"+Label, 4., weight, 10, 0., 10.);
  if( !(fabs(Mmumu1-91.2)>10 && fabs(Mmumu2-91.2)>10) ) return;
  FillHist("CutFlow"+Label, 5., weight, 10, 0., 10.);

  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, 6., weight, 10, 0., 10.);

  if(JetColl.size()>1){
    float M3l = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M();
    FillHist("CutFlow"+Label, 7., weight, 10, 0., 10.);
    FillHist("Composition_3l1b2jNoZ2l"+Label, 2., weight, 10, 0., 10.);
    FillHist("M3l"+Label, M3l, weight, 30, 0., 300.);
    FillHist("MOSSF1"+Label, Mmumu1, weight, 30., 0., 300.);
    FillHist("MOSSF2"+Label, Mmumu2, weight, 30., 0., 300.);
    FillHist("MSSSF"+Label, MSSSF, weight, 30., 0., 300.);
    if(M3l<85){
      FillHist("SigBins_1b"+Label, 4., weight, 20, 0., 20.);
    }
    else{
      float mindR=-1, mindRM=-1;
      for(unsigned int i=0; i<MuTColl.size(); i++){
        for(unsigned int j=i+1; j<MuTColl.size(); j++){
           float TmpdR = MuTColl.at(i).DeltaR(MuTColl.at(j));
           float TmpM  = (MuTColl.at(i)+MuTColl.at(j)).M();
           if(mindR<0 or TmpdR<mindR){ mindR=TmpdR; mindRM=TmpM; }
        }
      }
      if(mindRM>0 and mindRM<85){
        FillHist("SigBins_1b"+Label, 5., weight, 20, 0., 20.);
      }
    }
    if(BJetColl.size()>1){
      FillHist("CutFlow"+Label, 8., weight, 10, 0., 10.);
      FillHist("Composition_3l2b2jNoZ2l"+Label, 2., weight, 10, 0., 10.);
      if(M3l<85){
        FillHist("SigBins_2b"+Label, 4., weight, 20, 0., 20.);
      }
      else{
        float mindR=-1, mindRM=-1;
        for(unsigned int i=0; i<MuTColl.size(); i++){
          for(unsigned int j=i+1; j<MuTColl.size(); j++){
             float TmpdR = MuTColl.at(i).DeltaR(MuTColl.at(j));
             float TmpM  = (MuTColl.at(i)+MuTColl.at(j)).M();
             if(mindR<0 or TmpdR<mindR){ mindR=TmpdR; mindRM=TmpM; }
          }
        }
        if(mindRM>0 and mindRM<85){
          FillHist("SigBins_2b"+Label, 5., weight, 20, 0., 20.);
        }
      }
    }
  }

} 


void HNTopFeas::AnalyzeTriLepton(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                 vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NElT=ElTColl.size(), NMuT=MuTColl.size(), NElL=ElLColl.size(), NMuL=MuLColl.size();
  int NLepT=NElT+NMuT, NLepL=NElL+NMuL;
  int it_cut=3;
  if( !(NLepT==3 && NLepL==3) ) return;
  if( !(NMuT==NMuL && NElT==NElL) ) return;



  bool PassTrigAccept=false;
  if( NMuT>=2 && MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 ) PassTrigAccept=true;
  else if( NElT>=2 && ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15 ) PassTrigAccept=true;
  else if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 ) PassTrigAccept=true;
  else if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>25 ) PassTrigAccept=true;
  if(!PassTrigAccept) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

  int SumQ_M = SumCharge(MuTColl);
  int SumQ_E = SumCharge(ElTColl);
  int SumQ_L = SumQ_M+SumQ_E;
  if(abs(SumQ_L)!=1) return;
  bool IsMu3El0 = NMuT==3 && NElT==0, IsMu2El1 = NMuT==2 && NElT==1, IsMu1El2 = NMuT==1 && NElT==2, IsMu0El3 = NMuT==0 && NElT==3;
  float MOSSF1=-1, MOSSF2=-1, MSSSF=-1, M3l=-1;
  int IdxOS=-1, IdxSS1=-1, IdxSS2=-1;
  if(IsMu3El0){
    IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    MOSSF1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    MOSSF2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    MSSSF  = (MuTColl.at(IdxSS1)+MuTColl.at(IdxSS2)).M();
    M3l    = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M(); 
  }
  else if(IsMu0El3){
    IdxOS  = TriElChargeIndex(ElTColl, "OS");
    IdxSS1 = TriElChargeIndex(ElTColl, "SS1");
    IdxSS2 = TriElChargeIndex(ElTColl, "SS2");
    MOSSF1 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS1)).M();
    MOSSF2 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS2)).M();
    MSSSF  = (ElTColl.at(IdxSS1)+ElTColl.at(IdxSS2)).M();
    M3l    = (ElTColl.at(0)+ElTColl.at(1)+ElTColl.at(2)).M(); 
  }
  else if(IsMu2El1){
    if(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) MOSSF1=(MuTColl.at(0)+MuTColl.at(1)).M();
    else MSSSF = (MuTColl.at(0)+MuTColl.at(1)).M();
  }
  else if(IsMu1El2){
    if(ElTColl.at(0).Charge()!=ElTColl.at(1).Charge()) MOSSF1=(ElTColl.at(0)+ElTColl.at(1)).M();
    else MSSSF = (ElTColl.at(0)+ElTColl.at(1)).M();
  }

  float IsQCDLike = (MOSSF1>0 && MOSSF1<12) or (MOSSF2>0 && MOSSF2<12) or (NMuT>1 && MSSSF>0 && MSSSF<4);
  if(IsQCDLike) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

  float IsZLike = (MOSSF1>0 && fabs(MOSSF1-91.2)<10) or (MOSSF2>0 && fabs(MOSSF2-91.2)<10);
  if(MOSSF1>0) FillHist("MOSSF1", MOSSF1, weight, 30, 0., 300.);
  if(MOSSF2>0) FillHist("MOSSF2", MOSSF1, weight, 30, 0., 300.);
  if(IsZLike) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  FillHist("NBJet"+Label, BJetColl.size(), weight, 10, 0., 10.);

  if(BJetColl.size()<1) return; 
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  FillHist("NJet"+Label, JetColl.size(), weight, 10, 0., 10.);

  if(JetColl.size()<2) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  FillHist("M3l_PreSel"+Label, M3l, weight, 30, 0., 300.);
  if(MOSSF1>0) FillHist("MOSSF1_PreSel", MOSSF1, weight, 30, 0., 300.);
  if(MOSSF2>0) FillHist("MOSSF2_PreSel", MOSSF1, weight, 30, 0., 300.);


  if     (IsMu3El0) FillHist("NPresel", 2., weight, 12, 0., 12.);
  else if(IsMu2El1) FillHist("NPresel", 3., weight, 12, 0., 12.);
  else if(IsMu1El2) FillHist("NPresel", 4., weight, 12, 0., 12.);
  else if(IsMu0El3) FillHist("NPresel", 5., weight, 12, 0., 12.);

  if(IsMu3El0 or IsMu2El1) FillHist("NPresel_Tot", 1., weight, 3, 0., 3.);

}



void HNTopFeas::AnalyzeElectronDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                      std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  if( !(ElTColl.size()==1 && ElLColl.size()==1) ) return;
  if( !(MuTColl.size()==2 && MuLColl.size()==2) ) return;
  if( !(ElTColl.at(0).Pt()>25 or MuTColl.at(0).Pt()>20)  ) return;
  FillHist("CutFlow"+Label, 3., weight, 10, 0., 10.);

  bool SSmm = MuTColl.at(0).Charge()==MuTColl.at(1).Charge(), OSmm = !SSmm; 
  float MOSSF = OSmm? (MuTColl.at(0)+MuTColl.at(1)).M():-1.;
  if(OSmm && MOSSF<12) return;
  FillHist("CutFlow"+Label, 4., weight, 10, 0., 10.);

  if(OSmm && fabs(MOSSF-91.2)<10) return;
  FillHist("CutFlow"+Label, 5., weight, 10, 0., 10.);

  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, 6., weight, 10, 0., 10.);

  if(JetColl.size()>1){
    FillHist("CutFlow"+Label, 7., weight, 10, 0., 10.);
    float M3l     = (ElTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    float Mem_OS1 = -1., Mem_OS2 = -1., Mem_SS = -1.;
    float MSSSF   = SSmm? (MuTColl.at(0)+MuTColl.at(1)).M():-1.;
    if(OSmm){
      Mem_OS1 = MuTColl.at(0).Charge()==ElTColl.at(0).Charge()? (ElTColl.at(0)+MuTColl.at(1)).M():(ElTColl.at(0)+MuTColl.at(0)).M();
      Mem_SS  = MuTColl.at(0).Charge()==ElTColl.at(0).Charge()? (ElTColl.at(0)+MuTColl.at(0)).M():(ElTColl.at(0)+MuTColl.at(1)).M();
      FillHist("Composition_3l1b2jNoZ2l"+Label, 1., weight, 10, 0., 10.);
      FillHist("M3l_OSmm"+Label, M3l, weight, 30, 0., 300.);
      FillHist("MemOS1_OSmm"+Label, Mem_OS1, weight, 30, 0., 300.);
      FillHist("MemSS_OSmm"+Label, Mem_SS, weight, 30, 0., 300.);
      FillHist("MOSSF_OSmm"+Label, MOSSF, weight, 30., 0., 300.);
      if     (M3l<85  ) FillHist("SigBins_1b"+Label, 2., weight, 20, 0., 20.);
      else if(MOSSF<85) FillHist("SigBins_1b"+Label, 3., weight, 20, 0., 20.);
    }
    if(SSmm){
      Mem_OS1 = (ElTColl.at(0)+MuTColl.at(0)).M();
      Mem_OS2 = (ElTColl.at(0)+MuTColl.at(1)).M();
      FillHist("Composition_3l1b2jNoZ2l"+Label, 3., weight, 10, 0., 10.);
      FillHist("M3l_SSmm"+Label, M3l, weight, 30, 0., 300.);
      FillHist("MemOS1_SSmm"+Label, Mem_OS1, weight, 30, 0., 300.);
      FillHist("MemOS2_SSmm"+Label, Mem_OS2, weight, 30, 0., 300.);
      FillHist("MSSSF_SSmm"+Label, MSSSF, weight, 30., 0., 300.);
      if     (M3l<85  ) FillHist("SigBins_1b"+Label, 0., weight, 20, 0., 20.);
      else if(MSSSF<85) FillHist("SigBins_1b"+Label, 1., weight, 20, 0., 20.);
    }
    if(BJetColl.size()>1){
      FillHist("CutFlow"+Label, 8., weight, 10, 0., 10.);
      if(OSmm){
        FillHist("Composition_3l2b2jNoZ2l"+Label, 1., weight, 10, 0., 10.);
        if     (M3l<85  ) FillHist("SigBins_2b"+Label, 2., weight, 20, 0., 20.);
        else if(MOSSF<85) FillHist("SigBins_2b"+Label, 3., weight, 20, 0., 20.);
      }
      if(SSmm){
        FillHist("Composition_3l2b2jNoZ2l"+Label, 3., weight, 10, 0., 10.);
        if     (M3l<85  ) FillHist("SigBins_2b"+Label, 0., weight, 20, 0., 20.);
        else if(MSSSF<85) FillHist("SigBins_2b"+Label, 1., weight, 20, 0., 20.);
      }
    }
  }

}


void HNTopFeas::AnalyzeTetraLepton(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                   vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NElT = ElTColl.size(), NMuT=MuTColl.size(), NElL = ElLColl.size(), NMuL=MuLColl.size();
  int it_cut=3;
  if( (NElT+NMuT)!=4 ) return;
  if( !(NElT==NElL && NMuT==NMuL) ) return;

  bool PassTrigAccept=false;
  if     ( NMuT>=2 && MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 ) PassTrigAccept=true;
  else if( NElT>=2 && ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15 ) PassTrigAccept=true;
  else if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 ) PassTrigAccept=true;
  else if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>25 ) PassTrigAccept=true;
  if(!PassTrigAccept) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

  int SumQ_M = SumCharge(MuTColl);
  int SumQ_E = SumCharge(ElTColl);
  int SumQ_L = SumQ_M+SumQ_E;
  if(SumQ_L!=0) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  
  vector<int> IdxMupColl, IdxMumColl, IdxElpColl, IdxElmColl;
  for(unsigned int i=0; i<MuTColl.size(); i++){
    if(MuTColl.at(i).Charge()>0){IdxMupColl.push_back(i);} 
    else {IdxMumColl.push_back(i);}
  }
  for(unsigned int i=0; i<ElTColl.size(); i++){
    if(ElTColl.at(i).Charge()>0){IdxElpColl.push_back(i);} 
    else {IdxElmColl.push_back(i);}
  }

  bool IsZlike=false, IsQCDlike=false, IsZTo4l=false, FailMFilter=false, OSSF=false, SSSF=false, OddF=false;
  for(unsigned int i=0; i<MuTColl.size(); i++){
    for(unsigned int j=i+1; j<MuTColl.size(); j++){
      float TmpM = (MuTColl.at(i)+MuTColl.at(j)).M();
      if(TmpM<4) FailMFilter=true;
    }
  }
  for(unsigned int it_p=0; it_p<IdxMupColl.size(); it_p++){
    for(unsigned int it_m=0; it_m<IdxMumColl.size(); it_m++){
      float TmpM = (MuTColl.at(IdxMupColl.at(it_p))+MuTColl.at(IdxMumColl.at(it_m))).M();
      if(TmpM<12) IsQCDlike=true;
      else if(fabs(TmpM-91.2)<10) IsZlike=true;
    }
  }
  for(unsigned int it_p=0; it_p<IdxElpColl.size(); it_p++){
    for(unsigned int it_m=0; it_m<IdxElmColl.size(); it_m++){
      float TmpM = (ElTColl.at(IdxElpColl.at(it_p))+ElTColl.at(IdxElmColl.at(it_m))).M();
      if(TmpM<12) IsQCDlike=true;
      else if(fabs(TmpM-91.2)<10) IsZlike=true;
    }
  }
  if(FailMFilter) FillHist("MOSSFType"+Label, 0., weight, 10, 0., 10.);
  if(IsQCDlike)   FillHist("MOSSFType"+Label, 1., weight, 10, 0., 10.);
  if(IsZlike)     FillHist("MOSSFType"+Label, 2., weight, 10, 0., 10.);
  if(IsQCDlike) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

  float px4l=0, py4l=0, pz4l=0, E4l=0;
  for(unsigned int it_e=0; it_e<ElTColl.size(); it_e++){
    px4l+=ElTColl.at(it_e).Px(); py4l+=ElTColl.at(it_e).Py(); pz4l+=ElTColl.at(it_e).Pz(); E4l+=ElTColl.at(it_e).E();
  }
  for(unsigned int it_m=0; it_m<MuTColl.size(); it_m++){
    px4l+=MuTColl.at(it_m).Px(); py4l+=MuTColl.at(it_m).Py(); pz4l+=MuTColl.at(it_m).Pz(); E4l+=MuTColl.at(it_m).E();
  }
  TLorentzVector V4l(px4l,py4l,pz4l,E4l);
  float M4l = V4l.M();
  cout<<M4l<<endl;

  SSSF = (IdxMupColl.size()==2 and IdxElmColl.size()==2) or (IdxMumColl.size()==2 and IdxElpColl.size()==2);
  OSSF = (IdxMupColl.size()==IdxMumColl.size()) and (IdxElpColl.size()==IdxElmColl.size());
  OddF = NElT==1 or NMuT==1;
  IsZTo4l = OSSF and fabs(M4l-91.2)<10;
  if(IsZTo4l)   FillHist("MOSSFType"+Label, 3., weight, 10, 0., 10.);
  
  if     (OSSF and NElT==4) FillHist("Composition_4l"+Label, 0., weight, 10, 0., 10.);
  else if(OSSF and NElT==2) FillHist("Composition_4l"+Label, 1., weight, 10, 0., 10.);
  else if(OSSF and NElT==0) FillHist("Composition_4l"+Label, 2., weight, 10, 0., 10.);
  else if(SSSF and NElT==2) FillHist("Composition_4l"+Label, 3., weight, 10, 0., 10.);
  else if(OddF and NElT==1) FillHist("Composition_4l"+Label, 4., weight, 10, 0., 10.);
  else if(OddF and NElT==3) FillHist("Composition_4l"+Label, 5., weight, 10, 0., 10.);
  if(!IsZlike){
    if     (OSSF and NElT==4) FillHist("Composition_4lNoZ2l"+Label, 0., weight, 10, 0., 10.);
    else if(OSSF and NElT==2) FillHist("Composition_4lNoZ2l"+Label, 1., weight, 10, 0., 10.);
    else if(OSSF and NElT==0) FillHist("Composition_4lNoZ2l"+Label, 2., weight, 10, 0., 10.);
    else if(SSSF and NElT==2) FillHist("Composition_4lNoZ2l"+Label, 3., weight, 10, 0., 10.);
    else if(OddF and NElT==1) FillHist("Composition_4lNoZ2l"+Label, 4., weight, 10, 0., 10.);
    else if(OddF and NElT==3) FillHist("Composition_4lNoZ2l"+Label, 5., weight, 10, 0., 10.);
    if(!IsZTo4l){
      if     (OSSF and NElT==4) FillHist("Composition_4lNoZ2l4l"+Label, 0., weight, 10, 0., 10.);
      else if(OSSF and NElT==2) FillHist("Composition_4lNoZ2l4l"+Label, 1., weight, 10, 0., 10.);
      else if(OSSF and NElT==0) FillHist("Composition_4lNoZ2l4l"+Label, 2., weight, 10, 0., 10.);
      else if(SSSF and NElT==2) FillHist("Composition_4lNoZ2l4l"+Label, 3., weight, 10, 0., 10.);
      else if(OddF and NElT==1) FillHist("Composition_4lNoZ2l4l"+Label, 4., weight, 10, 0., 10.);
      else if(OddF and NElT==3) FillHist("Composition_4lNoZ2l4l"+Label, 5., weight, 10, 0., 10.);
    }
  }
  

  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  if     (OSSF and NElT==4) FillHist("Composition_4l1b"+Label, 0., weight, 10, 0., 10.);
  else if(OSSF and NElT==2) FillHist("Composition_4l1b"+Label, 1., weight, 10, 0., 10.);
  else if(OSSF and NElT==0) FillHist("Composition_4l1b"+Label, 2., weight, 10, 0., 10.);
  else if(SSSF and NElT==2) FillHist("Composition_4l1b"+Label, 3., weight, 10, 0., 10.);
  else if(OddF and NElT==1) FillHist("Composition_4l1b"+Label, 4., weight, 10, 0., 10.);
  else if(OddF and NElT==3) FillHist("Composition_4l1b"+Label, 5., weight, 10, 0., 10.);
  if(!IsZlike){
    if     (OSSF and NElT==4) FillHist("Composition_4l1bNoZ2l"+Label, 0., weight, 10, 0., 10.);
    else if(OSSF and NElT==2) FillHist("Composition_4l1bNoZ2l"+Label, 1., weight, 10, 0., 10.);
    else if(OSSF and NElT==0) FillHist("Composition_4l1bNoZ2l"+Label, 2., weight, 10, 0., 10.);
    else if(SSSF and NElT==2) FillHist("Composition_4l1bNoZ2l"+Label, 3., weight, 10, 0., 10.);
    else if(OddF and NElT==1) FillHist("Composition_4l1bNoZ2l"+Label, 4., weight, 10, 0., 10.);
    else if(OddF and NElT==3) FillHist("Composition_4l1bNoZ2l"+Label, 5., weight, 10, 0., 10.);
    if(!IsZTo4l){
      if     (OSSF and NElT==4) FillHist("Composition_4l1bNoZ2l4l"+Label, 0., weight, 10, 0., 10.);
      else if(OSSF and NElT==2) FillHist("Composition_4l1bNoZ2l4l"+Label, 1., weight, 10, 0., 10.);
      else if(OSSF and NElT==0) FillHist("Composition_4l1bNoZ2l4l"+Label, 2., weight, 10, 0., 10.);
      else if(SSSF and NElT==2) FillHist("Composition_4l1bNoZ2l4l"+Label, 3., weight, 10, 0., 10.);
      else if(OddF and NElT==1) FillHist("Composition_4l1bNoZ2l4l"+Label, 4., weight, 10, 0., 10.);
      else if(OddF and NElT==3) FillHist("Composition_4l1bNoZ2l4l"+Label, 5., weight, 10, 0., 10.);
    }
  }


  if(JetColl.size()<2) return;
  FillHist("CutFlow"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
  if     (OSSF and NElT==4) FillHist("Composition_4l1b2j"+Label, 0., weight, 10, 0., 10.);
  else if(OSSF and NElT==2) FillHist("Composition_4l1b2j"+Label, 1., weight, 10, 0., 10.);
  else if(OSSF and NElT==0) FillHist("Composition_4l1b2j"+Label, 2., weight, 10, 0., 10.);
  else if(SSSF and NElT==2) FillHist("Composition_4l1b2j"+Label, 3., weight, 10, 0., 10.);
  else if(OddF and NElT==1) FillHist("Composition_4l1b2j"+Label, 4., weight, 10, 0., 10.);
  else if(OddF and NElT==3) FillHist("Composition_4l1b2j"+Label, 5., weight, 10, 0., 10.);
  if(!IsZlike){
    if     (OSSF and NElT==4) FillHist("NPresel"+Label, 6., weight, 12, 0., 12.);
    else if(OSSF and NElT==2) FillHist("NPresel"+Label, 7., weight, 12, 0., 12.);
    else if(OSSF and NElT==0) FillHist("NPresel"+Label, 8., weight, 12, 0., 12.);
    else if(SSSF and NElT==2) FillHist("NPresel"+Label, 9., weight, 12, 0., 12.);
    else if(OddF and NElT==1) FillHist("NPresel"+Label, 10., weight, 12, 0., 12.);
    else if(OddF and NElT==3) FillHist("NPresel"+Label, 11., weight, 12, 0., 12.);
    if(NElT<3) FillHist("NPresel_Tot", 2., weight, 3, 0., 3.);

    if     (OSSF and NElT==4) FillHist("Composition_4l1b2jNoZ2l"+Label, 0., weight, 10, 0., 10.);
    else if(OSSF and NElT==2) FillHist("Composition_4l1b2jNoZ2l"+Label, 1., weight, 10, 0., 10.);
    else if(OSSF and NElT==0) FillHist("Composition_4l1b2jNoZ2l"+Label, 2., weight, 10, 0., 10.);
    else if(SSSF and NElT==2) FillHist("Composition_4l1b2jNoZ2l"+Label, 3., weight, 10, 0., 10.);
    else if(OddF and NElT==1) FillHist("Composition_4l1b2jNoZ2l"+Label, 4., weight, 10, 0., 10.);
    else if(OddF and NElT==3) FillHist("Composition_4l1b2jNoZ2l"+Label, 5., weight, 10, 0., 10.);
    if(!IsZTo4l){
      if     (OSSF and NElT==4) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 0., weight, 10, 0., 10.);
      else if(OSSF and NElT==2) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 1., weight, 10, 0., 10.);
      else if(OSSF and NElT==0) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 2., weight, 10, 0., 10.);
      else if(SSSF and NElT==2) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 3., weight, 10, 0., 10.);
      else if(OddF and NElT==1) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 4., weight, 10, 0., 10.);
      else if(OddF and NElT==3) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 5., weight, 10, 0., 10.);
    }
  }


}




void HNTopFeas::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

HNTopFeas::HNTopFeas(){

}

HNTopFeas::~HNTopFeas(){

}


