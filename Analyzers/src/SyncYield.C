#include "SyncYield.h"

void SyncYield::initializeAnalyzer(){

  ElMuMu=false, MuMuMu=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("ElMuMu"))     ElMuMu=true; 
    if(Userflags.at(i).Contains("MuMuMu"))     MuMuMu=true; 
    if(Userflags.at(i).Contains("SystRun"))    SystRun=true; 
  }

  DblMu=false, MuEG=false;
  if     (DataStream.Contains("DoubleMuon")) DblMu=true;
  else if(DataStream.Contains("MuonEG"))     MuEG =true;

  if(DataYear==2016){
    TrigList_DblMu_BtoG.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
    TrigList_DblMu_BtoG.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
    TrigList_DblMu_H.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu_H.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

    TrigList_MuEG_BtoG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG_H.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::CSVv2, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void SyncYield::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  bool PassTrig=false;
  bool IsPeriodH = IsDATA and run>280385;
  //bool IsPeriodH = IsDATA and run>=280919;
  if     (MuMuMu){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_DblMu_BtoG) or ev.PassTrigger(TrigList_DblMu_H);
    else{
      if(DblMu) PassTrig = !IsPeriodH? ev.PassTrigger(TrigList_DblMu_BtoG):ev.PassTrigger(TrigList_DblMu_H);
    }
  }
  else if(ElMuMu && MuEG){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_MuEG_BtoG) or ev.PassTrigger(TrigList_MuEG_H);
    else{
      if(MuEG)  PassTrig = !IsPeriodH? ev.PassTrigger(TrigList_MuEG_BtoG):ev.PassTrigger(TrigList_MuEG_H);
    }
  }
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
  if(MuMuMu and muonPreColl.size()>2    ) PreCutPass=true;
  if(ElMuMu and electronPreColl.size()>0 && muonPreColl.size()>1) PreCutPass=true;
  if(!PreCutPass) return;
  FillHist("CutFlow", 2., weight*TmpW, 10, 0., 10.);


  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "HctoWA16T", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "HctoWA16T", 25., 2.5);
  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "HctoWA16L", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "HctoWA16L", 10., 2.5);

  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::CSVv2, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  std::vector<Jet> jetNoVetoColl  = GetJets("tight", 25., 2.4);
  std::sort(jetNoVetoColl.begin(), jetNoVetoColl.end(), PtComparing);
  std::vector<Jet> bjetNoVetoColl = SelBJets(jetNoVetoColl, param_jets);
  std::vector<Jet> jetColl  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl, muonLooseColl, 0.4);
  std::vector<Jet> bjetColl = SelBJets(jetColl, param_jets);

  Particle vMET = ev.GetMETVector();
  Particle vMET_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);


  std::vector<Gen> truthColl;


  bool EventCand = false;
  if(MuMuMu){ EventCand = muonLooseColl.size()>2; }
  if(ElMuMu){ EventCand = electronLooseColl.size()>0 && muonLooseColl.size()>1; }

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
    //sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    //sf_trig   = mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig, "");
    //cout<<"w_gen:"<<w_gen<<" w_lumi:"<<w_lumi<<" w_PU:"<<w_PU<<" w_prefire:"<<w_prefire<<" sf_trig:"<<sf_trig<<endl;
    //cout<<"sf_mutk"<<sf_mutk<<" sf_muid:"<<sf_muid<<" sf_muiso:"<<sf_muiso<<" sf_elreco:"<<sf_elreco<<" sf_elid:"<<sf_elid<<" sf_btag:"<<sf_btag<<endl;
  }
  weight *= w_gen * w_filter * w_topptrw * w_lumi * w_PU * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

 
  if(MuMuMu){
    AnalyzeTriMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                   jetColl, bjetColl, vMET_xyCorr, weight, "");
    if(muonTightColl.size()!=3){
      AnalyzeTriMuon(muonLooseColl, muonLooseColl, electronLooseColl, electronLooseColl,
                     jetColl, bjetColl, vMET_xyCorr, weight, "_AppReg");
    }
  }
  if(ElMuMu){
    AnalyzeElectronDiMuon(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                          jetColl, bjetColl, vMET_xyCorr, weight, "");
    if(!(muonTightColl.size()==2 and electronTightColl.size()==1)){
      AnalyzeElectronDiMuon(muonLooseColl, muonLooseColl, electronLooseColl, electronLooseColl,
                            jetColl, bjetColl, vMET_xyCorr, weight, "_AppReg");
    }
  }


}



void SyncYield::AnalyzeTriMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
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
  float MmumuZ = fabs(Mmumu1-91.2)<fabs(Mmumu2-91.2)? Mmumu1:Mmumu2;
  float M3l    = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M();
  bool IsZLike = fabs(Mmumu1-91.2)<10 or fabs(Mmumu2-91.2)<10;
  bool IsQCDLike = Mmumu1<12 or Mmumu2<12;
  if( IsQCDLike ) return;
  FillHist("CutFlow"+Label, 4., weight, 10, 0., 10.);

  if(IsZLike){
    FillHist("NJet_WZReg"+Label, JetColl.size(), weight, 10, 0., 10.);
    FillHist("NB_WZReg"+Label, BJetColl.size(), weight, 10, 0., 10.);
    FillHist("PTMu1_WZReg"+Label, MuTColl.at(0).Pt(), weight, 20, 0., 200.);
    FillHist("PTMu2_WZReg"+Label, MuTColl.at(1).Pt(), weight, 20, 0., 200.);
    FillHist("PTMu3_WZReg"+Label, MuTColl.at(2).Pt(), weight, 20, 0., 200.);
    FillHist("Mmumu_WZReg"+Label, MmumuZ, weight, 30, 60., 120.);
  }
  else if( (Mmumu1<81.2 or Mmumu2<81.2) and fabs(M3l-91.2)<10 and vMET.Pt()<50 ){
    FillHist("NJet_ZGReg"+Label, JetColl.size(), weight, 10, 0., 10.);
    FillHist("NB_ZGReg"+Label, BJetColl.size(), weight, 10, 0., 10.);
    FillHist("M3l_ZGReg"+Label, M3l, weight, 30, 60., 120.);
  }

  if(IsZLike) return;
  FillHist("CutFlow"+Label, 5., weight, 10, 0., 10.);

  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, 6., weight, 10, 0., 10.);

  if(JetColl.size()<2) return;
  FillHist("CutFlow"+Label, 7., weight, 10, 0., 10.);

  int IdxSSA = TriMuChargeIndex_SyncYield(MuTColl, vMET.Pt(), vMET.Px(), vMET.Py(), "SSA");
  float MmumuA = (MuTColl.at(IdxSSA)+MuTColl.at(IdxOS)).M();

  if(MmumuA>80) return;
  FillHist("CutFlow"+Label, 8., weight, 10, 0., 10.);
  FillHist("MmumuA_SR"+Label, MmumuA, weight, 80., 0., 80.);

} 


void SyncYield::AnalyzeElectronDiMuon(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                      std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  if( !(ElTColl.size()==1 && ElLColl.size()==1) ) return;
  if( !(MuTColl.size()==2 && MuLColl.size()==2) ) return;
  if( !(ElTColl.at(0).Pt()>25)  ) return;
  if( !(SumCharge(MuTColl)==0) ) return;
  FillHist("CutFlow"+Label, 3., weight, 10, 0., 10.);

  float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
  float M3l   = (ElTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
  bool  IsQCDLike = Mmumu<12;
  bool  IsZLike   = fabs(Mmumu-91.2)<10;
  if(IsQCDLike) return;
  FillHist("CutFlow"+Label, 4., weight, 10, 0., 10.);

  if(IsZLike){
    FillHist("NJet_WZReg"+Label, JetColl.size(), weight, 10, 0., 10.);
    FillHist("NB_WZReg"+Label, BJetColl.size(), weight, 10, 0., 10.);
    FillHist("PTMu1_WZReg"+Label, MuTColl.at(0).Pt(), weight, 20, 0., 200.);
    FillHist("PTMu2_WZReg"+Label, MuTColl.at(1).Pt(), weight, 20, 0., 200.);
    FillHist("PTEl1_WZReg"+Label, ElTColl.at(0).Pt(), weight, 20, 0., 200.);
    FillHist("Mmumu_WZReg"+Label, Mmumu, weight, 30, 60., 120.);
  }
  else if( Mmumu<81.2 and fabs(M3l-91.2)<10 and vMET.Pt()<50 ){
    FillHist("NJet_ZGReg"+Label, JetColl.size(), weight, 10, 0., 10.);
    FillHist("NB_ZGReg"+Label, BJetColl.size(), weight, 10, 0., 10.);
    FillHist("M3l_ZGReg"+Label, M3l, weight, 30, 60., 120.);
  }

  if(IsZLike) return;
  FillHist("CutFlow"+Label, 5., weight, 10, 0., 10.);

  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, 6., weight, 10, 0., 10.);

  if(JetColl.size()<2) return;
  FillHist("CutFlow"+Label, 7., weight, 10, 0., 10.);
  
  if(Mmumu>80) return;
  FillHist("CutFlow"+Label, 8., weight, 10, 0., 10.);
  FillHist("MmumuA_SR"+Label, Mmumu, weight, 80., 0., 80.);

}



int SyncYield::TriMuChargeIndex_SyncYield(vector<Muon>& MuonColl, float MET, float METx, float METy, TString charge){
    //First Choose 2SS, 1OS muons(++- or --+) SS means 2 of them having same sign, OS means 1 of them having different sign from others.
    // charge="OS" will return the index of muon that having different charge from other 2,
    // charge="SS1" will return the index of first muon that having same sign
    // charge="SS2" will return the index of second muon that having same sign

    if(MuonColl.size()!=3) return -1;
    if(fabs(SumCharge(MuonColl))>1) return -1;
    
    int IdxOS=-1, IdxSS1=-1, IdxSS2=-1;
    if     (MuonColl.at(0).Charge()==MuonColl.at(1).Charge()){ IdxSS1=0, IdxSS2=1, IdxOS=2; }
    else if(MuonColl.at(0).Charge()==MuonColl.at(2).Charge()){ IdxSS1=0, IdxSS2=2, IdxOS=1; }
    else if(MuonColl.at(1).Charge()==MuonColl.at(2).Charge()){ IdxSS1=1, IdxSS2=2, IdxOS=0; }

    int ReturnIdx=-1;
    if     (charge.Contains("OS") ) ReturnIdx=IdxOS;
    else if(charge.Contains("SS1")) ReturnIdx=IdxSS1;
    else if(charge.Contains("SS2")) ReturnIdx=IdxSS2;

    int IdxSSW=-1, IdxSSA=-1;
    if(charge.Contains("SSW") || charge.Contains("SSA")){
      if(IdxSS1<0 || IdxSS2<0) return ReturnIdx;
      float dPt   = fabs(MuonColl.at(IdxSS1).Pt()-MuonColl.at(IdxSS2).Pt());
      float MTW1  = sqrt(2)*sqrt(MET*MuonColl.at(IdxSS1).Pt()-METx*MuonColl.at(IdxSS1).Px()-METy*MuonColl.at(IdxSS1).Py());
      float MTW2  = sqrt(2)*sqrt(MET*MuonColl.at(IdxSS2).Pt()-METx*MuonColl.at(IdxSS2).Px()-METy*MuonColl.at(IdxSS2).Py());
      bool  SS1_W = MTW1>50. && MTW1<120.;
      bool  SS2_W = MTW2>50. && MTW2<120.;
      if( dPt<25. && SS2_W && (!SS1_W) ){ IdxSSW=IdxSS2; IdxSSA=IdxSS1; }
      else{ IdxSSW=IdxSS1; IdxSSA=IdxSS2; }
    }
    if     (charge.Contains("SSW")) ReturnIdx=IdxSSW;
    else if(charge.Contains("SSA")) ReturnIdx=IdxSSA;
    
    return ReturnIdx;

}


void SyncYield::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

SyncYield::SyncYield(){

}

SyncYield::~SyncYield(){

}


