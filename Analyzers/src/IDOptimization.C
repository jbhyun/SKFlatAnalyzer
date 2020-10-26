#include "IDOptimization.h"

void IDOptimization::initializeAnalyzer(){

  MuMu=false, ElEl=false;
  MuID=false, ElID=false, ZData=false, SystRun=false; 
  ElTrigCut=false, ElTrigEffect=false, MuTrigEffect=false, TrigEffCheck=false, SelEffCheck=false;
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("MuMu"))      MuMu    = true; 
    if(Userflags.at(i).Contains("ElEl"))      ElEl    = true; 
    if(Userflags.at(i).Contains("ElID"))      ElID    = true; 
    if(Userflags.at(i).Contains("MuID"))      MuID    = true; 
    if(Userflags.at(i).Contains("ZData"))     ZData   = true; 
    if(Userflags.at(i).Contains("ElTrigCut"))    ElTrigCut    = true; 
    if(Userflags.at(i).Contains("ElTrigEffect")) ElTrigEffect = true; 
    if(Userflags.at(i).Contains("MuTrigEffect")) MuTrigEffect = true; 
    if(Userflags.at(i).Contains("TrigEffCheck")) TrigEffCheck = true; 
    if(Userflags.at(i).Contains("SelEffCheck"))  SelEffCheck  = true;
    if(Userflags.at(i).Contains("SystRun"))   SystRun   = true; 
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

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  if(DataYear==2017){
    TrigList_DblMu1.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu1.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu2.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  }
  if(DataYear==2018){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void IDOptimization::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  bool PassTrig=false;
  bool IsPeriod17B = DataYear==2017 && GetDataPeriod()=="B", IsPeriod17CtoF= DataYear==2017 && !IsPeriod17B;
  if(ZData){
    if     (ElID) PassTrig = ev.PassTrigger(TrigList_DblEG);
    else if(MuID){
      if(!IsDATA){
        if(DataYear==2017) PassTrig = ev.PassTrigger(TrigList_DblMu1) or ev.PassTrigger(TrigList_DblMu2);
        else               PassTrig = ev.PassTrigger(TrigList_DblMu);
      }
      else{
        if(DataYear==2017) PassTrig = (IsPeriod17B && ev.PassTrigger(TrigList_DblMu1)) or
                                      (IsPeriod17CtoF && ev.PassTrigger(TrigList_DblMu2));
        else               PassTrig = ev.PassTrigger(TrigList_DblMu);
      }
    }
  }
  if(!ZData && (MuID or ElID)) PassTrig=true;
  if(ElTrigCut or ElTrigEffect or MuTrigEffect or TrigEffCheck or SelEffCheck) PassTrig=true;
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
  if( ZData && MuID && muonPreColl.size()>1    ) PreCutPass=true;
  if( ZData && ElID && electronPreColl.size()>1) PreCutPass=true;
  if(!ZData && MuID && muonPreColl.size()>0    ) PreCutPass=true;
  if(!ZData && ElID && electronPreColl.size()>0) PreCutPass=true; 
  if( ElTrigCut && electronPreColl.size()>1    ) PreCutPass=true;
  if( MuTrigEffect && muonPreColl.size()>0     ) PreCutPass=true;
  if( ElTrigEffect && electronPreColl.size()>0     ) PreCutPass=true;
  if( TrigEffCheck && MuMu && muonPreColl.size()>1 ) PreCutPass=true;
  if( TrigEffCheck && ElEl && electronPreColl.size()>1 ) PreCutPass=true;
  if( SelEffCheck ) PreCutPass=true; 
  if(!PreCutPass) return;
  FillHist("CutFlow", 2., weight*TmpW, 10, 0., 10.);


  vector<Muon>     muonRawColl       = GetMuons("NOCUT", 0., 2.4);
  vector<Electron> electronRawColl   = GetElectrons("NOCUT", 0., 2.5);
  std::sort(muonRawColl.begin(), muonRawColl.end(), PtComparing);
  std::sort(electronRawColl.begin(), electronRawColl.end(), PtComparing);

  vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TopHN17T", 10., 2.4);
  vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TESTT", 10., 2.5);
  vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TESTL", 10., 2.5);


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

  float w_gen = 1., w_filter = 1., w_topptrw = 1., w_lumi = 1., w_PU = 1., w_prefire = 1., sf_trig = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
  if((!IsDATA)){
    w_gen     = ev.MCweight();
    w_lumi    = weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    //if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    //w_filter  = GetGenFilterEffCorr();
    //w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    w_PU      = GetPileUpWeight(nPileUp, 0);
    w_prefire = GetPrefireWeight(0);
    if(EventCand){
    sf_muid   = GetMuonSF(muonTightColl, "POGTID_genTrk", "ID");
    sf_muiso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    sf_elreco = GetElectronSF(electronTightColl, "", "Reco");
    sf_elid   = GetElectronSF(electronTightColl, "POGMVAIsoWP90", "ID");
    sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    //sf_trig   = mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig, "");
    }
  }
  weight *= w_gen * w_filter * w_topptrw * w_lumi * w_PU * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

  if(MuID){ 
    vector<Muon>  muonPOGMColl = SelectMuons(muonPreColl, "POGMedium", 10., 2.4);
    CheckIDVar_MC(muonPOGMColl, muonRawColl, electronLooseColl, electronRawColl, jetColl, bjetColl, vMET, weight, "");
  }
  if(ElID){
    vector<Electron> ElPOGMHLTCvColl = SelectElectrons(electronPreColl, "POGMVAni90HLT172lCv", 10., 2.5);
    CheckIDVar_MC(muonLooseColl, muonRawColl, ElPOGMHLTCvColl, electronRawColl, jetColl, bjetColl, vMET, weight, "");
  }
  if(ElTrigCut){
    vector<Electron> ElPOGMVAMColl = SelectElectrons(electronPreColl, "passMVAID_noIso_WP90", 10., 2.5);
    CheckTrigCut_MC(muonLooseColl, muonRawColl, electronLooseColl, electronRawColl, jetColl, bjetColl, vMET, ev, weight, "");
  }
  if(MuTrigEffect or ElTrigEffect){
    CheckTrigEffect(muonRawColl, electronRawColl, jetColl, bjetColl, vMET, ev, weight, "");
  }
  if(TrigEffCheck){
    if(MuMu){
      vector<Muon> MuMiniTColl = SelectMuons(muonPreColl, "TopHN17T", 10., 2.4);
      vector<Muon> MuMiniLColl = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
      vector<Muon> MuMiniTkTColl = SelectMuons(muonPreColl, "TopHN17T_TkIso", 10., 2.4);
      vector<Muon> MuMiniTkLColl = SelectMuons(muonPreColl, "TopHN17L_TkIso", 10., 2.4);
      CheckTrigEfficiency(MuMiniTColl, MuMiniLColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, weight, "");
      CheckTrigEfficiency(MuMiniTkTColl, MuMiniTkLColl, electronTightColl, electronLooseColl, jetColl, bjetColl, vMET, ev, weight, "_Tkltp4");
    }
    if(ElEl){
      vector<Electron> ElTopHNTColl = SelectElectrons(electronPreColl, "TopHN17T", 10., 2.5);
      vector<Electron> ElTopHNLColl = SelectElectrons(electronPreColl, "TopHN17L", 10., 2.5);
      vector<Electron> ElTopHNTColl_NoHLT = SelectElectrons(electronPreColl, "TopHN17T_NoHLT", 10., 2.5);
      vector<Electron> ElTopHNLColl_NoHLT = SelectElectrons(electronPreColl, "TopHN17L_NoHLT", 10., 2.5);
      CheckTrigEfficiency(muonTightColl, muonLooseColl, ElTopHNTColl, ElTopHNLColl, jetColl, bjetColl, vMET, ev, weight, "");
      CheckTrigEfficiency(muonTightColl, muonLooseColl, ElTopHNTColl_NoHLT, ElTopHNLColl_NoHLT, jetColl, bjetColl, vMET, ev, weight, "_NoHLT");
    }
  }
  if(SelEffCheck){
    vector<Muon> MuMiniTColl = SelectMuons(muonPreColl, "TopHN17T", 10., 2.4);
    vector<Muon> MuMiniLColl = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
    vector<Muon> MuMiniTkTColl = SelectMuons(muonPreColl, "TopHN17T_TkIso", 10., 2.4);
    vector<Muon> MuMiniTkLColl = SelectMuons(muonPreColl, "TopHN17L_TkIso", 10., 2.4);
    vector<Muon> MuPOGMIsoTColl = SelectMuons(muonPreColl, "POGMIPIsoT", 10., 2.4);
    vector<Muon> MuPOGMIsoVVLColl = SelectMuons(muonPreColl, "POGMIPIsoVVL", 10., 2.4);


    vector<Electron> ElTopHNTColl = SelectElectrons(electronPreColl, "TopHN17T", 10., 2.5);
    vector<Electron> ElTopHNLColl = SelectElectrons(electronPreColl, "TopHN17L", 10., 2.5);
    vector<Electron> ElTopHNTColl_NoHLT = SelectElectrons(electronPreColl, "TopHN17T_NoHLT", 10., 2.5);
    vector<Electron> ElTopHNLColl_NoHLT = SelectElectrons(electronPreColl, "TopHN17L_NoHLT", 10., 2.5);

    if(ElEl){
      CheckSelectionEfficiency(MuMiniTColl, MuMiniLColl, muonRawColl, ElTopHNTColl, ElTopHNLColl, electronRawColl, jetNoVetoColl, bjetNoVetoColl, vMET, ev, weight, "_HLT");
      CheckSelectionEfficiency(MuMiniTColl, MuMiniLColl, muonRawColl, ElTopHNTColl_NoHLT, ElTopHNLColl_NoHLT, electronRawColl, jetNoVetoColl, bjetNoVetoColl, vMET, ev, weight, "_NoHLT");
    }
    if(MuMu){
      CheckSelectionEfficiency(MuMiniTColl, MuMiniLColl, muonRawColl, ElTopHNTColl, ElTopHNLColl, electronRawColl, jetNoVetoColl, bjetNoVetoColl, vMET, ev, weight, "_NoHLT");
      CheckSelectionEfficiency(MuMiniTkTColl, MuMiniTkLColl, muonRawColl, ElTopHNTColl, ElTopHNLColl, electronRawColl, jetNoVetoColl, bjetNoVetoColl, vMET, ev, weight, "_HLT");
      CheckSelectionEfficiency(MuPOGMIsoTColl, MuPOGMIsoVVLColl, muonRawColl, ElTopHNTColl, ElTopHNLColl, electronRawColl, jetNoVetoColl, bjetNoVetoColl, vMET, ev, weight, "_POG");
    }
  }


}

void IDOptimization::TEST(vector<Muon>& MuRawColl, vector<Electron>& ElRawColl,
vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)
{

 return;
}


void IDOptimization::CheckSelectionEfficiency(
vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuRawColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElRawColl,
vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)
{
  //Check if LVeto in VVV search is useful as they claimed. It seems not.
  int it_cut=0; 
  vector<Gen> TruthColl = GetGens();

  if(false){
  //if(MuMu){
    FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    if( !(MuTColl.size()==2 && ElTColl.size()==0) ) return;
    FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    if( !(MuLColl.size()==2 && ElLColl.size()==0) ) return;
    FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    if( MuTColl.at(0).Charge()!=MuTColl.at(1).Charge() ) return;
    FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    if( !(MuTColl.at(0).Pt()>20) ) return;
    FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    FillHist("NBJets"+Label, BJetColl.size(), weight, 10, 0., 10.);
    if(BJetColl.size()==0) return;
    FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    float Mu1IsoIncLR03 = MuTColl.at(0).MiniRelIso();
    float Mu2IsoIncLR03 = MuTColl.at(1).MiniRelIso();
    float Mu1IsoIncLMini = MuTColl.at(0).MiniRelIso();
    float Mu2IsoIncLMini = MuTColl.at(1).MiniRelIso();

    for(unsigned int it_mr=0; it_mr<MuRawColl.size(); it_mr++){
      float dR1 = MuRawColl.at(it_mr).DeltaR(MuTColl.at(0));
      float dR2 = MuRawColl.at(it_mr).DeltaR(MuTColl.at(1));
      float PTMu1 = MuTColl.at(0).Pt(), PTMu2 = MuTColl.at(1).Pt();
      float dRmini1 = PTMu1<50.? 0.2: PTMu1<200.? 0.25-1E-3*PTMu1:0.05;
      float dRmini2 = PTMu2<50.? 0.2: PTMu2<200.? 0.25-1E-3*PTMu2:0.05;
      if(dR1<1E-3 or dR2<1E-3) continue;
      if(dR1<dRmini1) Mu1IsoIncLMini+= MuRawColl.at(it_mr).Pt()/PTMu1;
      if(dR2<dRmini2) Mu2IsoIncLMini+= MuRawColl.at(it_mr).Pt()/PTMu2;
      if(dR1<0.3){
        Mu1IsoIncLR03+= MuRawColl.at(it_mr).Pt()/PTMu1;
        int MuType = GetLeptonType_JH(MuTColl.at(0), TruthColl);
        int NearMuIdx = GenMatchedIdx(MuRawColl.at(it_mr),TruthColl);
        int NearMuType = GetLeptonType_JH(MuRawColl.at(it_mr),TruthColl);
        if(MuType<0 && MuType>-5){
          FillHist("NearMuType"+Label, NearMuType, weight, 20, -10., 10.);
          if(NearMuIdx!=-1){
            int MIdx = FirstNonSelfMotherIdx(NearMuIdx, TruthColl);
            int MPID = MIdx!=-1? abs(TruthColl.at(MIdx).PID()):0;
            FillHist("NearMuMPID"+Label, MPID, weight, 6E3, 0., 6E3);
            printf("NearMuMPID: %d\n", MPID);
          }
        }
      }
      if(dR2<0.3){
        Mu2IsoIncLR03+= MuRawColl.at(it_mr).Pt()/PTMu2;
        int MuType = GetLeptonType_JH(MuTColl.at(1), TruthColl);
        int NearMuIdx = GenMatchedIdx(MuRawColl.at(it_mr),TruthColl);
        int NearMuType = GetLeptonType_JH(MuRawColl.at(it_mr),TruthColl);
        if(MuType<0 && MuType>-5){
          FillHist("NearMuType"+Label, NearMuType, weight, 20, -10., 10.);
          if(NearMuIdx!=-1){
            int MIdx = FirstNonSelfMotherIdx(NearMuIdx, TruthColl);
            int MPID = MIdx!=-1? abs(TruthColl.at(MIdx).PID()):0;
            FillHist("NearMuMPID"+Label, MPID, weight, 6E3, 0., 6E3);
            printf("NearMuMPID: %d\n", MPID);
          }
        }
      }
    }
    for(unsigned int it_mr=0; it_mr<ElRawColl.size(); it_mr++){
      float dR1 = ElRawColl.at(it_mr).DeltaR(MuTColl.at(0));
      float dR2 = ElRawColl.at(it_mr).DeltaR(MuTColl.at(1));
      float PTMu1 = MuTColl.at(0).Pt(), PTMu2 = MuTColl.at(1).Pt();
      float dRmini1 = PTMu1<50.? 0.2: PTMu1<200.? 0.25-1E-3*PTMu1:0.05;
      float dRmini2 = PTMu2<50.? 0.2: PTMu2<200.? 0.25-1E-3*PTMu2:0.05;
      if(dR1<1E-3 or dR2<1E-3) continue;
      if(dR1<dRmini1) Mu1IsoIncLMini+= ElRawColl.at(it_mr).Pt()/PTMu1;
      if(dR2<dRmini2) Mu2IsoIncLMini+= ElRawColl.at(it_mr).Pt()/PTMu2;
      if(dR1<0.3){
        Mu1IsoIncLR03+= ElRawColl.at(it_mr).Pt()/PTMu1;
        int MuType = GetLeptonType_JH(MuTColl.at(0), TruthColl);
        //int NearElIdx = GenMatchedIdx(ElRawColl.at(it_mr),TruthColl);
        int NearElType = GetLeptonType_JH(ElRawColl.at(it_mr),TruthColl);
        if(MuType<0 && MuType>-5){
          FillHist("NearElType"+Label, NearElType, weight, 20, -10., 10.);
        }
      }
      if(dR2<0.3){
        Mu2IsoIncLR03+= ElRawColl.at(it_mr).Pt()/PTMu2;
        int MuType = GetLeptonType_JH(MuTColl.at(1), TruthColl);
        //int NearElIdx = GenMatchedIdx(ElRawColl.at(it_mr),TruthColl);
        int NearElType = GetLeptonType_JH(ElRawColl.at(it_mr),TruthColl);
        if(MuType<0 && MuType>-5){
          FillHist("NearElType"+Label, NearElType, weight, 20, -10., 10.);
        }
      }
    }
    FillHist("Mu1IsoInclLMini"+Label, Mu1IsoIncLMini, weight, 100, 0., 1.);
    FillHist("Mu2IsoInclLMini"+Label, Mu2IsoIncLMini, weight, 100, 0., 1.);
    FillHist("Mu1IsoInclLR03"+Label, Mu1IsoIncLR03, weight, 100, 0., 1.);
    FillHist("Mu2IsoInclLR03"+Label, Mu2IsoIncLR03, weight, 100, 0., 1.);
    if(Mu1IsoIncLMini<0.1 && Mu2IsoIncLMini<0.1){ FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++; }
    if(Mu1IsoIncLR03<0.1 && Mu2IsoIncLR03<0.1){ FillHist("CutFlow_Sel"+Label, it_cut, weight, 20, 0., 20.); it_cut++; }
  }
  if(MuMu){
    vector<TString> ListDiMuM = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
    bool PassTrig = ev.PassTrigger(ListDiMuM);
 
    if(!PassTrig) return;
    if(!(MuTColl.size()==2 && MuLColl.size()==2)) return;
    if(ElLColl.size()!=0) return;
    if(MuTColl.at(0).Charge()!=MuLColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    int it_cut=0;
    FillHist("CutFlow"+Label, it_cut, weight, 10, 0., 10.); it_cut++;

    float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<4) return;
    FillHist("CutFlow"+Label, it_cut, weight, 10, 0., 10.); it_cut++;

    std::vector<Jet> JetVetoColl  = JetsVetoLeptonInside(JetColl, ElLColl, MuLColl, 0.4);
    std::vector<Jet> BJetVetoColl = JetsVetoLeptonInside(BJetColl, ElLColl, MuLColl, 0.4);

    if(BJetVetoColl.size()==0) return;
    FillHist("CutFlow"+Label, it_cut, weight, 10, 0., 10.); it_cut++;

    if(JetVetoColl.size()<2) return;
    FillHist("CutFlow"+Label, it_cut, weight, 10, 0., 10.); it_cut++;
  }
  if(ElEl){
    vector<TString> ListDiEl = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    bool PassTrig = ev.PassTrigger(ListDiEl);
 
    if(!PassTrig) return;
    if(!(ElTColl.size()==2 && ElLColl.size()==2)) return;
    if(MuLColl.size()!=0) return;
    if(ElTColl.at(0).Charge()!=ElLColl.at(1).Charge()) return;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    int it_cut=0;
    FillHist("CutFlow"+Label, it_cut, weight, 10, 0., 10.); it_cut++;


    std::vector<Jet> JetVetoColl  = JetsVetoLeptonInside(JetColl, ElLColl, MuLColl, 0.4);
    std::vector<Jet> BJetVetoColl = JetsVetoLeptonInside(BJetColl, ElLColl, MuLColl, 0.4);


    if(BJetVetoColl.size()==0) return;
    FillHist("CutFlow"+Label, it_cut, weight, 10, 0., 10.); it_cut++;

    if(JetVetoColl.size()<2) return;
    FillHist("CutFlow"+Label, it_cut, weight, 10, 0., 10.); it_cut++;
  }
}



void IDOptimization::CheckTrigEfficiency(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                         vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)
{

  if(MuMu){
    int NMuT=MuTColl.size(), NMuL=MuLColl.size();
    if(NMuT!=2) return;
    if(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
    vector<TString> ListDiMu  = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"};
    vector<TString> ListDiMuM = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"};
    bool PassDiMuTrig = ev.PassTrigger(ListDiMu), PassDiMuMTrig = ev.PassTrigger(ListDiMuM);

    FillHist("TrigCount"+Label, 0., weight, 2, 0., 2.);
    if(PassDiMuTrig) FillHist("TrigCount"+Label, 1., weight, 2, 0., 2.);
    if(Mmumu>4){
      FillHist("TrigCountM"+Label, 0., weight, 2, 0., 2.);
      if(PassDiMuMTrig) FillHist("TrigCountM"+Label, 1., weight, 2, 0., 2.);
    }
    if(NMuL==2){
      FillHist("TrigCount_VetoL"+Label, 0., weight, 2, 0., 2.);
      if(PassDiMuTrig) FillHist("TrigCount_VetoL"+Label, 1., weight, 2, 0., 2.);
      if(Mmumu>4){
        FillHist("TrigCountM_VetoL"+Label, 0., weight, 2, 0., 2.);
        if(PassDiMuMTrig) FillHist("TrigCountM_VetoL"+Label, 1., weight, 2, 0., 2.);
      }
    }
  }
  if(ElEl){
    int NElT=ElTColl.size(), NElL=ElLColl.size();
    if(NElT!=2) return;
    if(!(ElTColl.at(0).Charge()!=ElTColl.at(1).Charge())) return;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    vector<TString> ListDiEl = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"};
    bool PassDiElTrig   = ev.PassTrigger(ListDiEl);
    bool MatchTagSiglEl = false, MatchLeg1=false, MatchLeg2=false;
    TString TagTrigName = "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v";
    TString FiltName_Leg1 = "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter";
    TString FiltName_Leg2 = "hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter";
    int IdxProbe=-1;
    if(DataStream.Contains("TrigInfo") or MCSample.Contains("TrigInfo")){
      for(unsigned int it_e=0; it_e<ElTColl.size(); it_e++){
        if(ElTColl.at(it_e).Pt()<35) continue;
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(ElTColl.at(it_e).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElTColl.at(it_e).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          if(HLTObject_FiredPaths->at(it_obj).find(TagTrigName)!=std::string::npos ){ MatchTagSiglEl=true; IdxProbe=1-it_e; break; }
        }
      }
      if(MatchTagSiglEl){
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(ElTColl.at(IdxProbe).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElTColl.at(IdxProbe).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          if(HLTObject_FiredFilters->at(it_obj).find(FiltName_Leg1)!=std::string::npos ){ MatchLeg1=true; }
          if(HLTObject_FiredFilters->at(it_obj).find(FiltName_Leg2)!=std::string::npos ){ MatchLeg2=true; }
          if(MatchLeg1 && MatchLeg2) break;
        }
      }
      for(unsigned int it_e=0; it_e<ElTColl.size(); it_e++){
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(ElTColl.at(it_e).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElTColl.at(it_e).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          if(HLTObject_FiredFilters->at(it_obj).find(FiltName_Leg1)!=std::string::npos ){ MatchLeg1=true; break; }
        }
      }
    }


    const int NPtEdges1=7, NPtEdges2=7, NEtaEdges=9;
    double PtEdges_El1[NPtEdges1]={25., 30., 40., 50., 100., 200., 500.};
    double PtEdges_El2[NPtEdges2]={15., 25., 35., 50., 100., 200., 500.};
    double EtaEdges[NEtaEdges]   ={-2.5,-2.0,-1.5,-0.8,0.,0.8,1.5,2.0,2.5};

    FillHist("TrigCount"+Label, 0., weight, 2, 0., 2.);
    FillHist("NEl_PT2D"+Label, ElTColl.at(0).Pt(), ElTColl.at(1).Pt(), weight, NPtEdges1-1, PtEdges_El1, NPtEdges2-1, PtEdges_El2);
    if(PassDiElTrig){
      FillHist("TrigCount"+Label, 1., weight, 2, 0., 2.);
      FillHist("NElTrig_PT2D"+Label, ElTColl.at(0).Pt(), ElTColl.at(1).Pt(), weight, NPtEdges1-1, PtEdges_El1, NPtEdges2-1, PtEdges_El2);
    }
    if(MatchTagSiglEl){
      vector<Gen> TruthColl = GetGens();
      int ProbeType = GetLeptonType_JH(ElTColl.at(IdxProbe), TruthColl);
      bool IsPrompt = ProbeType>0, IsFake=ProbeType<0 && ProbeType>-5;

      if(IsPrompt){
        FillHist("NElPr_PTEta"+Label, ElTColl.at(IdxProbe).Pt(), ElTColl.at(IdxProbe).Eta(), weight, NPtEdges2-1, PtEdges_El2, NEtaEdges-1, EtaEdges);
        if(MatchLeg1) FillHist("NElPrLeg1_PTEta"+Label, ElTColl.at(IdxProbe).Pt(), ElTColl.at(IdxProbe).Eta(), weight, NPtEdges1-1, PtEdges_El1, NEtaEdges-1, EtaEdges);
        if(MatchLeg2) FillHist("NElPrLeg2_PTEta"+Label, ElTColl.at(IdxProbe).Pt(), ElTColl.at(IdxProbe).Eta(), weight, NPtEdges2-1, PtEdges_El2, NEtaEdges-1, EtaEdges);
      }
      if(IsFake){
        FillHist("NElFk_PTEta"+Label, ElTColl.at(IdxProbe).Pt(), ElTColl.at(IdxProbe).Eta(), weight, NPtEdges2-1, PtEdges_El2, NEtaEdges-1, EtaEdges);
        if(MatchLeg1) FillHist("NElFkLeg1_PTEta"+Label, ElTColl.at(IdxProbe).Pt(), ElTColl.at(IdxProbe).Eta(), weight, NPtEdges1-1, PtEdges_El1, NEtaEdges-1, EtaEdges);
        if(MatchLeg2) FillHist("NElFkLeg2_PTEta"+Label, ElTColl.at(IdxProbe).Pt(), ElTColl.at(IdxProbe).Eta(), weight, NPtEdges2-1, PtEdges_El2, NEtaEdges-1, EtaEdges);
      }
    }


    if(NElL==2){//Observed to have negligible impact, no additional lepton impact for 2l final state.
      FillHist("TrigCount_VetoL"+Label, 0., weight, 2, 0., 2.);
      FillHist("NEl_PT2D_VetoL"+Label, ElTColl.at(0).Pt(), ElTColl.at(1).Pt(), weight, NPtEdges1-1, PtEdges_El1, NPtEdges2-1, PtEdges_El2);
      if(PassDiElTrig){
        FillHist("TrigCount_VetoL"+Label, 1., weight, 2, 0., 2.);
        FillHist("NElTrig_PT2D_VetoL"+Label, ElTColl.at(0).Pt(), ElTColl.at(1).Pt(), weight, NPtEdges1-1, PtEdges_El1, NPtEdges2-1, PtEdges_El2);
      }
    }
  }

}


void IDOptimization::CheckTrigEffect(vector<Muon>& MuRawColl, vector<Electron>& ElRawColl,
                                     vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)
{


  if(ElTrigEffect){

    vector<Gen> TruthColl = GetGens();
    vector<Electron> ElTColl, ElLColl, ElLNoIsoMVAColl;
    for(unsigned int it_e=0; it_e<ElRawColl.size(); it_e++){
      int LeptonType = GetLeptonType_JH(ElRawColl.at(it_e), TruthColl);
      bool IsPrompt = LeptonType>0 && LeptonType<3;
      float SIP2D = ElRawColl.at(it_e).dXYerr()!=0?   fabs(ElRawColl.at(it_e).dXY()/ElRawColl.at(it_e).dXYerr()):0.;
      float SIP3D = ElRawColl.at(it_e).IP3Derr()!=0.? fabs(ElRawColl.at(it_e).IP3D()/ElRawColl.at(it_e).IP3Derr()):0.;
      bool PassHLTDZ = fabs(ElRawColl.at(it_e).dZ())<0.1, PassHLTIDIso = ElRawColl.at(it_e).Pass_CaloIdL_TrackIdL_IsoVL17();
      bool PassMVA = ElRawColl.at(it_e).passMVAID_noIso_WP90(), PassConv = ElRawColl.at(it_e).PassConversionVeto();
      bool PassNMissHit = ElRawColl.at(it_e).NMissingHits()<2;
      bool PassMiniIsoT  = ElRawColl.at(it_e).MiniRelIso()<0.1;
      if(MCSample.Contains("_MN")) IsPrompt = LeptonType>0 && LeptonType!=3;
      if(ElRawColl.at(it_e).Pt()<10) continue;

      if(SIP2D<3 && SIP3D<5 && PassHLTDZ && PassHLTIDIso && PassConv && PassNMissHit){
        ElLNoIsoMVAColl.push_back(ElRawColl.at(it_e));
      }

      if(!IsPrompt) continue;
      int it_cut=0;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!PassMVA) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!(SIP2D<3)) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!(SIP3D<5)) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!PassHLTDZ) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!PassConv) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!PassNMissHit) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!PassHLTIDIso) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!PassMiniIsoT) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

    }

    if(!(ElLNoIsoMVAColl.size()>1 && ElLNoIsoMVAColl.at(0).Pt()>35)) return;

    bool FireRefTrig=false; int IdxRefEl=-1;
    for(unsigned int it_e=0; it_e<ElLNoIsoMVAColl.size(); it_e++){
      if(ElLNoIsoMVAColl.at(it_e).Pt()<35) continue;
      if(MCSample.Contains("Trig")){
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(ElLNoIsoMVAColl.at(it_e).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElLNoIsoMVAColl.at(it_e).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          if(HLTObject_FiredPaths->at(it_obj).find("HLT_Ele32_WPTight_Gsf_v")!=std::string::npos ){ FireRefTrig=true; IdxRefEl=it_e; break; }
        }
      }
    }

    for(unsigned int it_e=0; it_e<ElLNoIsoMVAColl.size(); it_e++){
      if(!FireRefTrig) break;
      if(ElLNoIsoMVAColl.at(it_e).Pt()<15) continue;
      if((int) it_e==IdxRefEl) continue;
      if(ElLNoIsoMVAColl.at(it_e).DeltaR(ElLNoIsoMVAColl.at(IdxRefEl))<0.4) continue;

      float MiniIso = ElLNoIsoMVAColl.at(it_e).MiniRelIso();
      float MVA     = ElLNoIsoMVAColl.at(it_e).MVANoIso();
      float IsEB1 = fabs(ElLNoIsoMVAColl.at(it_e).scEta())<0.8;
      float IsEB2 = !IsEB1 && fabs(ElLNoIsoMVAColl.at(it_e).scEta())<1.5;
      float IsEE  = !IsEB1 && !IsEB2 && fabs(ElLNoIsoMVAColl.at(it_e).scEta())<2.5;
      if(IsEB1) FillHist("NElID_Iso_EB1"+Label, MiniIso, weight, 20, 0., 1.);
      if(IsEB2) FillHist("NElID_Iso_EB2"+Label, MiniIso, weight, 20, 0., 1.);
      if(IsEE ) FillHist("NElID_Iso_EE"+Label, MiniIso, weight, 20, 0., 1.);
      if(IsEB1) FillHist("NElID_MVA_EB1"+Label, MVA, weight, 40, -1., 1.);
      if(IsEB2) FillHist("NElID_MVA_EB2"+Label, MVA, weight, 40, -1., 1.);
      if(IsEE ) FillHist("NElID_MVA_EE"+Label, MVA, weight, 40, -1., 1.);

      bool FireTrig=false;
      if(MCSample.Contains("Trig")){
        for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
          float dR = sqrt(pow(ElLNoIsoMVAColl.at(it_e).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElLNoIsoMVAColl.at(it_e).Phi()-HLTObject_phi->at(it_obj),2));
          if(dR>0.2) continue;
          if(HLTObject_FiredPaths->at(it_obj).find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")!=std::string::npos ){FireTrig=true; break;}
        }
      }
      if(FireTrig){
        if(IsEB1) FillHist("NElIDTrig_Iso_EB1"+Label, MiniIso, weight, 20, 0., 1.);
        if(IsEB2) FillHist("NElIDTrig_Iso_EB2"+Label, MiniIso, weight, 20, 0., 1.);
        if(IsEE ) FillHist("NElIDTrig_Iso_EE"+Label, MiniIso, weight, 20, 0., 1.);
        if(IsEB1) FillHist("NElIDTrig_MVA_EB1"+Label, MVA, weight, 40, -1., 1.);
        if(IsEB2) FillHist("NElIDTrig_MVA_EB2"+Label, MVA, weight, 40, -1., 1.);
        if(IsEE ) FillHist("NElIDTrig_MVA_EE"+Label, MVA, weight, 40, -1., 1.);
      } 
    }
  }
  if(MuTrigEffect){

    vector<Gen> TruthColl = GetGens();
    vector<Muon> MuTColl, MuTColl2, MuLColl, MuLNoIsoColl;
    for(unsigned int it_m=0; it_m<MuRawColl.size(); it_m++){
      int LeptonType = GetLeptonType_JH(MuRawColl.at(it_m), TruthColl);
      bool IsPrompt = LeptonType>0 && LeptonType<3;
      float SIP2D = MuRawColl.at(it_m).dXYerr()!=0?   fabs(MuRawColl.at(it_m).dXY()/MuRawColl.at(it_m).dXYerr()):0.;
      float SIP3D = MuRawColl.at(it_m).IP3Derr()!=0.? fabs(MuRawColl.at(it_m).IP3D()/MuRawColl.at(it_m).IP3Derr()):0.;
      if(MCSample.Contains("_MN")) IsPrompt = LeptonType>0 && LeptonType!=3;
      if(!IsPrompt) continue;
      if(MuRawColl.at(it_m).Pt()<10) continue;

      int it_cut=0;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!MuRawColl.at(it_m).isPOGMedium()) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!(SIP2D<3)) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!(SIP3D<4)) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;

      if(!(fabs(MuRawColl.at(it_m).dZ())<0.1)) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
      FillHist("dXY_IDIP"+Label, fabs(MuRawColl.at(it_m).dXY()), weight, 100, 0., 0.1);

      if(MuRawColl.at(it_m).MiniRelIso()<0.4){ MuLColl.push_back(MuRawColl.at(it_m)); }

      if(!(MuRawColl.at(it_m).MiniRelIso()<0.1)) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
      MuTColl.push_back(MuRawColl.at(it_m));

      if(!(MuRawColl.at(it_m).TrkIso()/MuRawColl.at(it_m).Pt()<0.4)) continue;
      FillHist("CutFlow_ID"+Label, it_cut, weight, 20, 0., 20.); it_cut++;
      MuTColl2.push_back(MuRawColl.at(it_m));
    }

    if(MuTColl.size()==2 && MuLColl.size()==2 && MuTColl.at(0).Pt()>20){
      //tight-dimuon trigger efficiency w/o trk iso cut
      int it_cut=0;
      FillHist("CutFlow_Trig"+Label, it_cut, weight, 10, 0., 10.); it_cut++;

      bool PassDiMuTrig = DataYear==2017? ev.PassTrigger(TrigList_DblMu1):ev.PassTrigger(TrigList_DblMu);
      if(PassDiMuTrig){ FillHist("CutFlow_Trig"+Label, it_cut, weight, 10, 0., 10.); it_cut++; }
    }

    if(MuTColl2.size()==2 && MuLColl.size()==2 && MuTColl2.at(0).Pt()>20){
      //tight-dimuon trigger efficiency w/ trk iso cut
      int it_cut=0;
      FillHist("CutFlow_Trig2"+Label, it_cut, weight, 10, 0., 10.); it_cut++;

      bool PassDiMuTrig = DataYear==2017? ev.PassTrigger(TrigList_DblMu1):ev.PassTrigger(TrigList_DblMu);
      if(PassDiMuTrig){ FillHist("CutFlow_Trig2"+Label, it_cut, weight, 10, 0., 10.); it_cut++; }
    }

    for(unsigned int it_m=0; it_m<MuRawColl.size(); it_m++){
      float SIP2D = MuRawColl.at(it_m).dXYerr()!=0?   fabs(MuRawColl.at(it_m).dXY()/MuRawColl.at(it_m).dXYerr()):0.;
      float SIP3D = MuRawColl.at(it_m).IP3Derr()!=0.? fabs(MuRawColl.at(it_m).IP3D()/MuRawColl.at(it_m).IP3Derr()):0.;
      if(MuRawColl.at(it_m).Pt()<10) continue;
      if(!MuRawColl.at(it_m).isPOGMedium()) continue;
      if(!(SIP2D<3)) continue;
      if(!(SIP3D<4)) continue;
      if(!(fabs(MuRawColl.at(it_m).dZ())<0.1)) continue;
      MuLNoIsoColl.push_back(MuRawColl.at(it_m));
    }
    if(MuLNoIsoColl.size()==1){
      float MiniIso = MuLNoIsoColl.at(0).MiniRelIso();
      FillHist("NMuID_Iso"+Label, MiniIso, weight, 20, 0., 1.);

      if(ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v")) FillHist("NMuIDTrig_Iso"+Label, MiniIso, weight, 20, 0., 1.);
    }

    for(unsigned int it_m=0; it_m<MuLNoIsoColl.size(); it_m++){
      if(!MCSample.Contains("TrigInfo")) continue;
      float MiniIso = MuLNoIsoColl.at(it_m).MiniRelIso();
      float RelIso  = MuLNoIsoColl.at(it_m).RelIso();
      float TkIso   = MuLNoIsoColl.at(it_m).TrkIso()/MuLNoIsoColl.at(it_m).Pt();
      FillHist("NMuID_MiniIso_Filter"+Label, MiniIso, weight, 20, 0., 1.);
      FillHist("NMuID_RelIso_Filter"+Label, RelIso, weight, 20, 0., 1.);
      FillHist("NMuID_TkIso_Filter"+Label, TkIso, weight, 20, 0., 1.);
      if(TkIso<0.4) FillHist("NMuID_MiniIso_Tkltp4_Filter"+Label, MiniIso, weight, 20, 0., 1.);

      bool FireLeg1=false, FireLeg2=false;
      for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
        float dR = sqrt(pow(MuLNoIsoColl.at(it_m).Eta()-HLTObject_eta->at(it_obj),2)+pow(MuLNoIsoColl.at(it_m).Phi()-HLTObject_phi->at(it_obj),2));
        if(dR>0.2) continue;
        //if(HLTObject_FiredPaths->at(it_obj).find("HLT_Mu8_TrkIsoVVL_v")!=std::string::npos ) {FireLeg1=true;}
        if(HLTObject_FiredFilters->at(it_obj).find("hltL3fL1sMu5L1f0L2f5L3Filtered8TkIsoFiltered0p4")!=std::string::npos ){
          FireLeg1=true; FillHist("dR_muHLT"+Label, dR, weight, 50, 0., 0.5);
        }
        //if(HLTObject_FiredPaths->at(it_obj).find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v")!=std::string::npos ) {FireLeg2=true;}
        if(HLTObject_FiredFilters->at(it_obj).find("hltDiMuon178Mass3p8Filtered")!=std::string::npos ){ FireLeg2=true; }
        if(FireLeg1 && FireLeg2) break;
      }
      if(FireLeg1){
        FillHist("NMuIDTrig1_MiniIso_Filter"+Label, MiniIso, weight, 20, 0., 1.);
        FillHist("NMuIDTrig1_RelIso_Filter"+Label, RelIso, weight, 20, 0., 1.);
        FillHist("NMuIDTrig1_TkIso_Filter"+Label, TkIso, weight, 20, 0., 1.);
        if(TkIso<0.4) FillHist("NMuIDTrig1_MiniIso_Tkltp4_Filter"+Label, MiniIso, weight, 20, 0., 1.);
      }
      if(FireLeg2){
        FillHist("NMuIDTrig2_MiniIso_Filter"+Label, MiniIso, weight, 20, 0., 1.);
        FillHist("NMuIDTrig2_RelIso_Filter"+Label, RelIso, weight, 20, 0., 1.);
        FillHist("NMuIDTrig2_TkIso_Filter"+Label, TkIso, weight, 20, 0., 1.);
        if(TkIso<0.4) FillHist("NMuIDTrig2_MiniIso_Tkltp4_Filter"+Label, MiniIso, weight, 20, 0., 1.);
      }
    }
  }

}



void IDOptimization::CheckTrigCut_MC(vector<Muon>& MuLColl, vector<Muon>& MuRawColl, vector<Electron>& ElLColl, vector<Electron>& ElRawColl,
                                     vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, float weight, TString Label)
{

  vector<Gen> TruthColl = GetGens();
  if(ElTrigCut){
    TString SiglTrigName = DataYear>2016? "HLT_Ele32_WPTight_Gsf_v":"HLT_Ele27_WPTight_Gsf_v";
    TString DblTrigName  = DataYear>2016? "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v":"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

    if(ElLColl.size()!=2) return;
    if(!(ElLColl.at(0).Pt()>35 && ElLColl.at(1).Pt()>15)) return;
    if(ElLColl.at(0).DeltaR(ElLColl.at(1))<0.4) return;
    if(!(fabs(ElLColl.at(0).dZ())<0.1 && fabs(ElLColl.at(1).dZ())<0.1)) return;

    bool PassSiglTrig=false, PassElLeg=false;
    if(MCSample.Contains("TrigInfo")){
      for(unsigned int it_obj=0; it_obj<HLTObject_eta->size(); it_obj++){
        float dR1 = sqrt(pow(ElLColl.at(0).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElLColl.at(0).Phi()-HLTObject_phi->at(it_obj),2));
        float dR2 = sqrt(pow(ElLColl.at(1).Eta()-HLTObject_eta->at(it_obj),2)+pow(ElLColl.at(1).Phi()-HLTObject_phi->at(it_obj),2));
        
        if( dR1<0.2 && HLTObject_FiredPaths->at(it_obj).find(SiglTrigName)!=std::string::npos ){
          PassSiglTrig=true; FillHist("dR_el1SiglHLT"+Label, dR1, weight, 40, 0., 0.4);
        }
        if( dR2<0.2 && HLTObject_FiredPaths->at(it_obj).find(DblTrigName)!=std::string::npos ){
          PassElLeg=true; FillHist("dR_el2LegHLT"+Label, dR2, weight, 50, 0., 0.5);
        }
        if(PassSiglTrig && PassElLeg) break;
      }
    }
    if(!PassSiglTrig) return;

    float ElPt = ElLColl.at(1).Pt();
    //float RhoEA = ElLColl.at(1).Rho()*ElLColl.at(1).EA();
    float IsEB = fabs(ElLColl.at(1).scEta())<1.479;

    if(IsEB){
      FillHist("ElEB_SigIetaIeta"+Label, ElLColl.at(1).Full5x5_sigmaIetaIeta(), weight, 50, 0., 0.05);
      FillHist("ElEB_HoverE"+Label, ElLColl.at(1).HoverE(), weight, 50, 0., 1.);
      FillHist("ElEB_dEtaInSeed"+Label, fabs(ElLColl.at(1).dEtaSeed()), weight, 50, 0., 0.05);
      FillHist("ElEB_dPhiIn"+Label, fabs(ElLColl.at(1).dPhiIn()), weight, 40, 0., 0.2);
      FillHist("ElEB_ooEooP"+Label, fabs(ElLColl.at(1).InvEminusInvP()), weight, 50, 0., 1.);
      FillHist("ElEB_PFEcalIso"+Label, ElLColl.at(1).ecalPFClusterIso()/ElPt, weight, 100, 0., 1.);
      FillHist("ElEB_PFHcalIso"+Label, ElLColl.at(1).hcalPFClusterIso()/ElPt, weight, 100, 0., 1.);
      FillHist("ElEB_TrkIso"+Label, ElLColl.at(1).dr03TkSumPt()/ElPt, weight, 50, 0., 1.);
    }
    else{
      FillHist("ElEE_SigIetaIeta"+Label, ElLColl.at(1).Full5x5_sigmaIetaIeta(), weight, 50, 0., 0.05);
      FillHist("ElEE_HoverE"+Label, ElLColl.at(1).HoverE(), weight, 50, 0., 1.);
      FillHist("ElEE_dEtaInSeed"+Label, fabs(ElLColl.at(1).dEtaSeed()), weight, 50, 0., 0.05);
      FillHist("ElEE_dPhiIn"+Label, fabs(ElLColl.at(1).dPhiIn()), weight, 40, 0., 0.2);
      FillHist("ElEE_ooEooP"+Label, fabs(ElLColl.at(1).InvEminusInvP()), weight, 50, 0., 1.);
      FillHist("ElEE_PFEcalIso"+Label, ElLColl.at(1).ecalPFClusterIso()/ElPt, weight, 100, 0., 1.);
      FillHist("ElEE_PFHcalIso"+Label, ElLColl.at(1).hcalPFClusterIso()/ElPt, weight, 100, 0., 1.);
      FillHist("ElEE_TrkIso"+Label, ElLColl.at(1).dr03TkSumPt()/ElPt, weight, 50, 0., 1.);
    }

    if(PassElLeg){
      if(IsEB){
        FillHist("ElEBTrig_SigIetaIeta"+Label, ElLColl.at(1).Full5x5_sigmaIetaIeta(), weight, 50, 0., 0.05);
        FillHist("ElEBTrig_HoverE"+Label, ElLColl.at(1).HoverE(), weight, 50, 0., 1.);
        FillHist("ElEBTrig_dEtaInSeed"+Label, fabs(ElLColl.at(1).dEtaSeed()), weight, 50, 0., 0.05);
        FillHist("ElEBTrig_dPhiIn"+Label, fabs(ElLColl.at(1).dPhiIn()), weight, 40, 0., 0.2);
        FillHist("ElEBTrig_ooEooP"+Label, fabs(ElLColl.at(1).InvEminusInvP()), weight, 50, 0., 1.);
        FillHist("ElEBTrig_PFEcalIso"+Label, ElLColl.at(1).ecalPFClusterIso()/ElPt, weight, 100, 0., 1.);
        FillHist("ElEBTrig_PFHcalIso"+Label, ElLColl.at(1).hcalPFClusterIso()/ElPt, weight, 100, 0., 1.);
        FillHist("ElEBTrig_TrkIso"+Label, ElLColl.at(1).dr03TkSumPt()/ElPt, weight, 50, 0., 1.);
      }
      else{
        FillHist("ElEETrig_SigIetaIeta"+Label, ElLColl.at(1).Full5x5_sigmaIetaIeta(), weight, 50, 0., 0.05);
        FillHist("ElEETrig_HoverE"+Label, ElLColl.at(1).HoverE(), weight, 50, 0., 1.);
        FillHist("ElEETrig_dEtaInSeed"+Label, fabs(ElLColl.at(1).dEtaSeed()), weight, 50, 0., 0.05);
        FillHist("ElEETrig_dPhiIn"+Label, fabs(ElLColl.at(1).dPhiIn()), weight, 40, 0., 0.2);
        FillHist("ElEETrig_ooEooP"+Label, fabs(ElLColl.at(1).InvEminusInvP()), weight, 50, 0., 1.);
        FillHist("ElEETrig_PFEcalIso"+Label, ElLColl.at(1).ecalPFClusterIso()/ElPt, weight, 100, 0., 1.);
        FillHist("ElEETrig_PFHcalIso"+Label, ElLColl.at(1).hcalPFClusterIso()/ElPt, weight, 100, 0., 1.);
        FillHist("ElEETrig_TrkIso"+Label, ElLColl.at(1).dr03TkSumPt()/ElPt, weight, 50, 0., 1.);
      }
    }
  }//End of ElTrig

}



void IDOptimization::CheckIDVar_MC(vector<Muon>& MuLColl, vector<Muon>& MuRawColl, vector<Electron>& ElLColl, vector<Electron>& ElRawColl,
                                   vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  vector<Gen> TruthColl = GetGens();

  if(MuID){
    for(unsigned int it_m=0; it_m<MuLColl.size(); it_m++){
      int LeptonType = GetLeptonType_JH(MuLColl.at(it_m), TruthColl);
      bool IsPrompt = LeptonType>0 && LeptonType<3, IsHadronicFake = LeptonType<0 && LeptonType>-5;
      if(MCSample.Contains("_MN")) IsPrompt = LeptonType>0 && LeptonType!=3;
      float SIP2D = MuLColl.at(it_m).dXYerr()!=0?   fabs(MuLColl.at(it_m).dXY()/MuLColl.at(it_m).dXYerr()):0.;
      float SIP3D = MuLColl.at(it_m).IP3Derr()!=0.? fabs(MuLColl.at(it_m).IP3D()/MuLColl.at(it_m).IP3Derr()):0.;
      float RelIso04 = MuLColl.at(it_m).RelIso();
      float LepIso04 = 0.;
      for(unsigned int it_mr=0; it_mr<MuRawColl.size(); it_mr++){
        if(MuLColl.at(it_m).DeltaR(MuRawColl.at(it_mr))<1E-3) continue;
        if(MuLColl.at(it_m).DeltaR(MuRawColl.at(it_mr))<0.4) LepIso04+=MuRawColl.at(it_mr).Pt();
      }
      for(unsigned int it_er=0; it_er<ElRawColl.size(); it_er++){
        if(MuLColl.at(it_m).DeltaR(ElRawColl.at(it_er))<1E-3) continue;
        if(MuLColl.at(it_m).DeltaR(ElRawColl.at(it_er))<0.4) LepIso04+=ElRawColl.at(it_er).Pt();
      }
      LepIso04 *= 1./MuLColl.at(it_m).Pt();
      LepIso04 += RelIso04;

      if(IsPrompt){
        FillHist("MiniRelIsoPF_Type12" +Label, MuLColl.at(it_m).MiniRelIso(), weight, 20, 0., 0.4);
        FillHist("RelIsoPF_Type12" +Label, RelIso04, weight, 20, 0., 0.4);
        FillHist("RelIsoPFL_Type12" +Label, LepIso04, weight, 20, 0., 0.4);
        FillHist("dxy_Type12" +Label, fabs(MuLColl.at(it_m).dXY()) , weight, 50, 0., 0.05);
        FillHist("dz_Type12"  +Label, fabs(MuLColl.at(it_m).dZ())  , weight, 50, 0., 0.1);
        FillHist("IP3D_Type12"+Label, fabs(MuLColl.at(it_m).IP3D()), weight, 100, 0., 0.2);
        FillHist("SIP2D_Type12"+Label, SIP2D, weight, 20, 0., 10.);
        FillHist("SIP3D_Type12"+Label, SIP3D, weight, 20, 0., 10.);
        FillHist("h2_SIP2D_SIP3D_Type12"+Label, SIP2D, SIP3D, weight, 20, 0., 10., 20, 0., 10.); 
      }
      if(IsHadronicFake){
        FillHist("MiniRelIsoPF_Typem1234" +Label, MuLColl.at(it_m).MiniRelIso(), weight, 20, 0., 0.4);
        FillHist("RelIsoPF_Typem1234" +Label, RelIso04, weight, 20, 0., 0.4);
        FillHist("RelIsoPFL_Typem1234" +Label, LepIso04, weight, 20, 0., 0.4);
        FillHist("dxy_Typem1234" +Label, fabs(MuLColl.at(it_m).dXY()) , weight, 50, 0., 0.05);
        FillHist("dz_Typem1234"  +Label, fabs(MuLColl.at(it_m).dZ())  , weight, 50, 0., 0.1);
        FillHist("IP3D_Typem1234"+Label, fabs(MuLColl.at(it_m).IP3D()), weight, 100, 0., 0.2);
        FillHist("SIP2D_Typem1234"+Label, SIP2D, weight, 20, 0., 10.);
        FillHist("SIP3D_Typem1234"+Label, SIP3D, weight, 20, 0., 10.);
        FillHist("h2_SIP2D_SIP3D_Typem1234"+Label, SIP2D, SIP3D, weight, 20, 0., 10., 20, 0., 10.); 
      }
    }
  }
  if(ElID){
    for(unsigned int it_e=0; it_e<ElLColl.size(); it_e++){
      int LeptonType = GetLeptonType_JH(ElLColl.at(it_e), TruthColl);
      bool IsPrompt = LeptonType>0 && LeptonType<3, IsHadronicFake = LeptonType<0 && LeptonType>-5;
      if(MCSample.Contains("_MN")) IsPrompt = LeptonType>0 && LeptonType!=3;
      bool IsExConv = LeptonType<-4 && LeptonType>-6;
      float SIP2D = ElLColl.at(it_e).dXYerr()!=0?   fabs(ElLColl.at(it_e).dXY()/ElLColl.at(it_e).dXYerr()):0.;
      float SIP3D = ElLColl.at(it_e).IP3Derr()!=0.? fabs(ElLColl.at(it_e).IP3D()/ElLColl.at(it_e).IP3Derr()):0.;
      float RelIso03 = ElLColl.at(it_e).RelIso();
      int TightQPass = ElLColl.at(it_e).IsGsfCtfScPixChargeConsistent()? 1:0;      

      if(IsPrompt){
        FillHist("RelIsoPF_Type12" +Label, RelIso03, weight, 20, 0., 0.4);
        FillHist("MiniRelIsoPF_Type12" +Label, ElLColl.at(it_e).MiniRelIso(), weight, 20, 0., 0.4);
        FillHist("dxy_Type12" +Label, fabs(ElLColl.at(it_e).dXY()) , weight, 50, 0., 0.05);
        FillHist("dz_Type12"  +Label, fabs(ElLColl.at(it_e).dZ())  , weight, 50, 0., 0.1);
        FillHist("IP3D_Type12"+Label, fabs(ElLColl.at(it_e).IP3D()), weight, 100, 0., 0.2);
        FillHist("SIP2D_Type12"+Label, SIP2D, weight, 20, 0., 10.);
        FillHist("SIP3D_Type12"+Label, SIP3D, weight, 20, 0., 10.);
        FillHist("h2_SIP2D_SIP3D_Type12"+Label, SIP2D, SIP3D, weight, 20, 0., 10., 20, 0., 10.); 
        if(SIP2D<3 && SIP3D<5){
          FillHist("NMissHits_atSIP2D3D_Type12"+Label, ElLColl.at(it_e).NMissingHits(), weight, 10, 0., 10.);    
          FillHist("TightQ_atSIP2D3D_Type12"+Label, TightQPass, weight, 2, 0., 2.);
        }
      }
      if(IsHadronicFake){
        FillHist("RelIsoPF_Typem1234" +Label, RelIso03, weight, 20, 0., 0.4);
        FillHist("MiniRelIsoPF_Typem1234" +Label, ElLColl.at(it_e).MiniRelIso(), weight, 20, 0., 0.4);
        FillHist("dxy_Typem1234" +Label, fabs(ElLColl.at(it_e).dXY()) , weight, 50, 0., 0.05);
        FillHist("dz_Typem1234"  +Label, fabs(ElLColl.at(it_e).dZ())  , weight, 50, 0., 0.1);
        FillHist("IP3D_Typem1234"+Label, fabs(ElLColl.at(it_e).IP3D()), weight, 100, 0., 0.2);
        FillHist("SIP2D_Typem1234"+Label, SIP2D, weight, 20, 0., 10.);
        FillHist("SIP3D_Typem1234"+Label, SIP3D, weight, 20, 0., 10.);
        FillHist("h2_SIP2D_SIP3D_Typem1234"+Label, SIP2D, SIP3D, weight, 20, 0., 10., 20, 0., 10.); 
        if(SIP2D<3 && SIP3D<5){
          FillHist("NMissHits_atSIP2D3D_Typem1234"+Label, ElLColl.at(it_e).NMissingHits(), weight, 10, 0., 10.);
          FillHist("TightQ_atSIP2D3D_Typem1234"+Label, TightQPass, weight, 2, 0., 2.);
        }
      }
      if(IsExConv){
        FillHist("RelIsoPF_Typem56" +Label, RelIso03, weight, 20, 0., 0.4);
        FillHist("MiniRelIsoPF_Typem56" +Label, ElLColl.at(it_e).MiniRelIso(), weight, 20, 0., 0.4);
        FillHist("dxy_Typem56" +Label, fabs(ElLColl.at(it_e).dXY()) , weight, 50, 0., 0.05);
        FillHist("dz_Typem56"  +Label, fabs(ElLColl.at(it_e).dZ())  , weight, 50, 0., 0.1);
        FillHist("IP3D_Typem56"+Label, fabs(ElLColl.at(it_e).IP3D()), weight, 100, 0., 0.2);
        FillHist("SIP2D_Typem56"+Label, SIP2D, weight, 20, 0., 10.);
        FillHist("SIP3D_Typem56"+Label, SIP3D, weight, 20, 0., 10.);
        FillHist("h2_SIP2D_SIP3D_Typem56"+Label, SIP2D, SIP3D, weight, 20, 0., 10., 20, 0., 10.); 
        if(SIP2D<3 && SIP3D<5){
          FillHist("NMissHits_atSIP2D3D_Typem56"+Label, ElLColl.at(it_e).NMissingHits(), weight, 10, 0., 10.);
          FillHist("TightQ_atSIP2D3D_Typem56"+Label, TightQPass, weight, 2, 0., 2.);
        }
      }
    }
  }
}


void IDOptimization::CheckIDVar_ZData(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
                                      vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  if(MuID){
    if(!(MuTColl.size()==2 && MuLColl.size()==2)) return;
    if(!(ElTColl.size()==0 && ElLColl.size()==0)) return;
    if(MuTColl.at(0).Charge()==MuTColl.at(1).Charge()) return;
    float M2l = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(fabs(M2l-91.2)>10) return;

  }
  if(ElID){
    if(!(ElTColl.size()==2 && ElLColl.size()==2)) return;
    if(!(MuTColl.size()==0 && MuLColl.size()==0)) return;
    if(ElTColl.at(0).Charge()==ElTColl.at(1).Charge()) return;
    float M2l = (ElTColl.at(0)+ElTColl.at(1)).M();
    if(fabs(M2l-91.2)>10) return;
  }
}




void IDOptimization::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

IDOptimization::IDOptimization(){

}

IDOptimization::~IDOptimization(){

}


