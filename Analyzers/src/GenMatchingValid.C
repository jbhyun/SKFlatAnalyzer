#include "GenMatchingValid.h"

void GenMatchingValid::initializeAnalyzer(){

  PhotonMatch=false, LepTypeValid=false, HNTypeValid=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("PhotonMatch"))  PhotonMatch  = true;
    if(Userflags.at(i).Contains("LepTypeValid")) LepTypeValid = true;
    if(Userflags.at(i).Contains("HNTypeValid"))  HNTypeValid  = true;
    if(Userflags.at(i).Contains("SystRun"))      SystRun      = true; 
  }

  DblMu=false, DblEG=false, MuEG=false;
  if     (DataStream.Contains("DoubleMuon")) DblMu=true;
  else if(DataStream.Contains("MuonEG"))     MuEG =true;
  else if(DataStream.Contains("DoubleEG"))   DblEG=true;
  else if(DataYear==2018 and DataStream.Contains("EGamma")) DblEG=true;

  NPrintMax=20, NCountPrint=0;

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void GenMatchingValid::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;

  bool PassTrig=false;
  if(PhotonMatch or LepTypeValid) PassTrig=true;
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
  if(PhotonMatch  &&  electronPreColl.size()!=0) PreCutPass=true;
  if(LepTypeValid && (electronPreColl.size()!=0 or muonPreColl.size()!=0)) PreCutPass=true;
  if(!PreCutPass) return;
  FillHist("CutFlow", 2., weight*TmpW, 10, 0., 10.);


  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TESTT", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TESTT", 10., 2.5);
  std::vector<Muon>     muonCBPOGMColl     = SelectMuons(muonPreColl, "POGIDMIsoM", 10., 2.4);
  std::vector<Electron> electronCBPOGMColl = SelectElectrons(electronPreColl, "passMediumID", 10., 2.5);
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

 
  if(PhotonMatch){
    truthColl = GetGens();
    CheckPhotonMatching(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                        jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_TESTID");
//    CheckPhotonMatching(muonTightColl, muonLooseColl, electronPreColl, electronPreColl,
//                        jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_TESTIDNoConv");
//    CheckPhotonMatching(muonTightColl, muonLooseColl, electronPreColl, electronPreColl,
//                        jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_TESTIDNoIP");
//    CheckPhotonMatching(muonTightColl, muonLooseColl, electronPreColl, electronPreColl,
//                        jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_CBPOGT");
//    CheckPhotonMatching(muonTightColl, muonLooseColl, electronPreColl, electronPreColl,
//                        jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_CBPOGM");
//    CheckPhotonMatching(muonTightColl, muonLooseColl, electronPreColl, electronPreColl,
//                        jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_NoID");
  }
  if(LepTypeValid){
    truthColl = GetGens();
    CheckLepTypeValidity(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                         jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_TESTID");
    CheckLepTypeValidity(muonCBPOGMColl, muonLooseColl, electronCBPOGMColl, electronLooseColl,
                         jetColl, bjetColl, vMET_xyCorr, ev, truthColl, weight, "_CBPOGM");
  }


}


void GenMatchingValid::CheckLepTypeValidity(
vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, vector<Gen>& TruthColl, float weight, TString Label)
{

  vector<int> Idx_FakeMu, Idx_ExConvMu, Idx_InConvMu, Idx_PrMu;
  vector<int> Idx_FakeEl, Idx_ExConvEl, Idx_InConvEl, Idx_PrEl;
  int Idx_v=-1, Idx_vx=-1;
  for(unsigned int it_m=0; it_m<MuTColl.size(); it_m++){
    int LepType = GetLeptonType_JH(MuTColl.at(it_m),TruthColl);
    FillHist("LepType_Mu"+Label, LepType, 1, 20, -10., 10.);
    if     (LepType>0  && LepType<4){ Idx_PrMu.push_back(it_m); }
    else if(LepType>3              ){ Idx_InConvMu.push_back(it_m); }
    else if(LepType>-5 && LepType<0){ Idx_FakeMu.push_back(it_m); }
    else if(LepType<-4             ){ Idx_ExConvMu.push_back(it_m); }
  }

  vector<int> Idx_Stm1, Idx_Stm234;
  for(unsigned int it_e=0; it_e<ElTColl.size(); it_e++){
    int LepType = GetLeptonType_JH(ElTColl.at(it_e),TruthColl);
    FillHist("LepType_El"+Label, LepType, 1, 20, -10., 10.);
    if     (LepType>0  && LepType<4){ Idx_PrEl.push_back(it_e); }
    else if(LepType>3              ){ Idx_InConvEl.push_back(it_e); }
    else if(LepType>-5 && LepType<0){ Idx_FakeEl.push_back(it_e); }
    else if(LepType<-4             ){ Idx_ExConvEl.push_back(it_e); }

    if(LepType==-1) Idx_Stm1.push_back(it_e);
    if(LepType<=-2 && LepType>=-4) Idx_Stm234.push_back(it_e);
  }

  if(MCSample.Contains("_MN")){
    FillHist("NPr_Mu"+Label, Idx_PrMu.size(), 1, 10, 0., 10.);
    FillHist("NInConv_Mu"+Label, Idx_InConvMu.size(), 1, 10, 0., 10.);
    FillHist("NFake_Mu"+Label, Idx_InConvMu.size(), 1, 10, 0., 10.);
    if(Idx_PrMu.size()==2){
      FillHist("M2P_Mu"+Label, (MuTColl.at(Idx_PrMu.at(0))+MuTColl.at(Idx_PrMu.at(1))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrMu.size()==1 && Idx_InConvMu.size()==1){
      FillHist("M1P1inC_Mu"+Label, (MuTColl.at(Idx_PrMu.at(0))+MuTColl.at(Idx_InConvMu.at(0))).M(), 1, 200, 0., 200.);
    }
    if(Idx_InConvMu.size()==2){
      FillHist("M2inC_Mu"+Label, (MuTColl.at(Idx_InConvMu.at(0))+MuTColl.at(Idx_InConvMu.at(1))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrMu.size()==1 && Idx_FakeMu.size()>0){
      FillHist("M1P1F_Mu"+Label, (MuTColl.at(Idx_PrMu.at(0))+MuTColl.at(Idx_FakeMu.at(0))).M(), 1, 200, 0., 200.);
    }
  }

  if(MCSample.Contains("ZGTo2LG") or MCSample.Contains("ZGToLLG")){
    FillHist("NPr_Mu"+Label, Idx_PrMu.size(), 1, 10, 0., 10.);
    FillHist("NInConv_Mu"+Label, Idx_InConvMu.size(), 1, 10, 0., 10.);
    FillHist("NFake_Mu"+Label, Idx_InConvMu.size(), 1, 10, 0., 10.);
    if(Idx_PrMu.size()==2 && Idx_InConvMu.size()==1){
      FillHist("M2P1C_Mu"+Label, (MuTColl.at(Idx_PrMu.at(0))+MuTColl.at(Idx_PrMu.at(1))+MuTColl.at(Idx_InConvMu.at(0))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrMu.size()==2 && Idx_InConvMu.size()==2){
      FillHist("M2P2C_Mu"+Label, (MuTColl.at(Idx_PrMu.at(0))+MuTColl.at(Idx_PrMu.at(1))+MuTColl.at(Idx_InConvMu.at(0))+MuTColl.at(Idx_InConvMu.at(1))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrMu.size()==2 && Idx_FakeMu.size()==1){
      FillHist("M2P1F_Mu"+Label, (MuTColl.at(Idx_PrMu.at(0))+MuTColl.at(Idx_PrMu.at(1))+MuTColl.at(Idx_FakeMu.at(0))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrMu.size()==2 && Idx_FakeMu.size()==2){
      FillHist("M2P2F_Mu"+Label, (MuTColl.at(Idx_PrMu.at(0))+MuTColl.at(Idx_PrMu.at(1))+MuTColl.at(Idx_FakeMu.at(0))+MuTColl.at(Idx_FakeMu.at(1))).M(), 1, 200, 0., 200.);
    }


    FillHist("NPr_El"+Label, Idx_PrEl.size(), 1, 10, 0., 10.);
    FillHist("NInConv_El"+Label, Idx_InConvEl.size(), 1, 10, 0., 10.);
    FillHist("NFake_El"+Label, Idx_InConvEl.size(), 1, 10, 0., 10.);
    if(Idx_PrEl.size()==2 && Idx_InConvEl.size()==1){
      FillHist("M2P1inC_El"+Label, (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_InConvEl.at(0))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrEl.size()==2 && Idx_InConvEl.size()==2){
      FillHist("M2P2inC_El"+Label, (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_InConvEl.at(0))+ElTColl.at(Idx_InConvEl.at(1))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrEl.size()==2 && Idx_ExConvEl.size()==1){
      FillHist("M2P1exC_El"+Label, (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_ExConvEl.at(0))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrEl.size()==2 && Idx_ExConvEl.size()==2){
      FillHist("M2P2exC_El"+Label, (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_ExConvEl.at(0))+ElTColl.at(Idx_ExConvEl.at(1))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrEl.size()==2 && Idx_FakeEl.size()==1){
      FillHist("M2P1F_El"+Label, (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_FakeEl.at(0))).M(), 1, 200, 0., 200.);
      if(Idx_Stm1.size()>0)   FillHist("M2P1Fm1_El"+Label,   (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_Stm1.at(0))).M(), 1, 200, 0., 200.);
      if(Idx_Stm234.size()>0) FillHist("M2P1Fm234_El"+Label, (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_Stm234.at(0))).M(), 1, 200, 0., 200.);
    }
    if(Idx_PrEl.size()==2 && Idx_FakeEl.size()==2){
      FillHist("M2P2F_El"+Label, (ElTColl.at(Idx_PrEl.at(0))+ElTColl.at(Idx_PrEl.at(1))+ElTColl.at(Idx_FakeEl.at(0))+ElTColl.at(Idx_FakeEl.at(1))).M(), 1, 200, 0., 200.);
    }
  }

  if(MCSample.Contains("TT") and MCSample.Contains("powheg")){
    for(unsigned int it_gen=2; it_gen<TruthColl.size(); it_gen++){
      int pid = TruthColl.at(it_gen).PID(), apid = abs(pid);
      if( !(apid==12 or apid==14 or apid==16) ) continue;
      int MIdx_last = FirstNonSelfMotherIdx(it_gen,TruthColl);
      int MPID = MIdx_last>0? TruthColl.at(MIdx_last).PID():-1, aMPID=abs(MPID);
      if(aMPID!=24) continue;
      if(pid>0){ Idx_v=it_gen; } else { Idx_vx=it_gen; }
    }

    for(unsigned int it_m=0; it_m<Idx_PrMu.size(); it_m++){
      if(MuTColl.at(Idx_PrMu.at(it_m)).Charge()>0 && Idx_v>0){
        FillHist("Mlv_PrMu"+Label, (MuTColl.at(Idx_PrMu.at(it_m))+TruthColl.at(Idx_v)).M(), 1, 200., 0., 200.);
      }
      else if(MuTColl.at(Idx_PrMu.at(it_m)).Charge()<0 && Idx_vx>0){
        FillHist("Mlv_PrMu"+Label, (MuTColl.at(Idx_PrMu.at(it_m))+TruthColl.at(Idx_vx)).M(), 1, 200., 0., 200.);
      }
      else{ FillHist("ElseCount_PrMu"+Label, 0., 1, 1, 0., 1.); }//charge-flip rate

      FillHist("Iso_PrMu"+Label, MuTColl.at(Idx_PrMu.at(it_m)).RelIso(), 1, 50, 0., 10.);
    }

    for(unsigned int it_m=0; it_m<Idx_FakeMu.size(); it_m++){
      if(MuTColl.at(Idx_FakeMu.at(it_m)).Charge()>0 && Idx_v>0){
        FillHist("Mlv_FakeMu"+Label, (MuTColl.at(Idx_FakeMu.at(it_m))+TruthColl.at(Idx_v)).M(), 1, 200., 0., 200.);
      }
      else if(MuTColl.at(Idx_FakeMu.at(it_m)).Charge()<0 && Idx_vx>0){
        FillHist("Mlv_FakeMu"+Label, (MuTColl.at(Idx_FakeMu.at(it_m))+TruthColl.at(Idx_vx)).M(), 1, 200., 0., 200.);
      }
      FillHist("Iso_FakeMu"+Label, MuTColl.at(Idx_FakeMu.at(it_m)).RelIso(), 1, 50, 0., 10.);
    }

    for(unsigned int it_m=0; it_m<Idx_InConvMu.size(); it_m++){
      if(MuTColl.at(Idx_InConvMu.at(it_m)).Charge()>0 && Idx_v>0){
        FillHist("Mlv_InConvMu"+Label, (MuTColl.at(Idx_InConvMu.at(it_m))+TruthColl.at(Idx_v)).M(), 1, 200., 0., 200.);
      }
      else if(MuTColl.at(Idx_InConvMu.at(it_m)).Charge()<0 && Idx_vx>0){
        FillHist("Mlv_InConvMu"+Label, (MuTColl.at(Idx_InConvMu.at(it_m))+TruthColl.at(Idx_vx)).M(), 1, 200., 0., 200.);
      }
      FillHist("Iso_InConvMu"+Label, MuTColl.at(Idx_InConvMu.at(it_m)).RelIso(), 1, 50, 0., 10.);
    }


    for(unsigned int it_e=0; it_e<Idx_PrEl.size(); it_e++){
      if(ElTColl.at(Idx_PrEl.at(it_e)).Charge()>0 && Idx_v>0){
        FillHist("Mlv_PrEl"+Label, (ElTColl.at(Idx_PrEl.at(it_e))+TruthColl.at(Idx_v)).M(), 1, 200., 0., 200.);
      }
      else if(ElTColl.at(Idx_PrEl.at(it_e)).Charge()<0 && Idx_vx>0){
        FillHist("Mlv_PrEl"+Label, (ElTColl.at(Idx_PrEl.at(it_e))+TruthColl.at(Idx_vx)).M(), 1, 200., 0., 200.);
      }
      else{ FillHist("ElseCount_PrEl"+Label, 0., 1, 1, 0., 1.); }//charge-flip rate
      FillHist("Iso_PrEl"+Label, ElTColl.at(Idx_PrEl.at(it_e)).RelIso(), 1, 50, 0., 10.);
    }

    for(unsigned int it_e=0; it_e<Idx_FakeEl.size(); it_e++){
      if(ElTColl.at(Idx_FakeEl.at(it_e)).Charge()>0 && Idx_v>0){
        FillHist("Mlv_FakeEl"+Label, (ElTColl.at(Idx_FakeEl.at(it_e))+TruthColl.at(Idx_v)).M(), 1, 200., 0., 200.);
      }
      else if(ElTColl.at(Idx_FakeEl.at(it_e)).Charge()<0 && Idx_vx>0){
        FillHist("Mlv_FakeEl"+Label, (ElTColl.at(Idx_FakeEl.at(it_e))+TruthColl.at(Idx_vx)).M(), 1, 200., 0., 200.);
      }
      FillHist("Iso_FakeEl"+Label, ElTColl.at(Idx_FakeEl.at(it_e)).RelIso(), 1, 50, 0., 10.);
    }

    for(unsigned int it_e=0; it_e<Idx_InConvEl.size(); it_e++){
      if(ElTColl.at(Idx_InConvEl.at(it_e)).Charge()>0 && Idx_v>0){
        FillHist("Mlv_InConvEl"+Label, (ElTColl.at(Idx_InConvEl.at(it_e))+TruthColl.at(Idx_v)).M(), 1, 200., 0., 200.);

        float Mlv = (ElTColl.at(Idx_InConvEl.at(it_e))+TruthColl.at(Idx_v)).M();
        if(fabs(Mlv-81)<1.){
          int GenIdx_El = GenMatchedIdx(ElTColl.at(Idx_InConvEl.at(it_e)), TruthColl);
          PrintGen(TruthColl);
          printf("IdxEl:%d, Idxv:%d, ElPt:%.2f, ElEta:%.2f, ElPhi:%.2f\n", GenIdx_El, Idx_v,
                  ElTColl.at(Idx_InConvEl.at(it_e)).Pt(), ElTColl.at(Idx_InConvEl.at(it_e)).Eta(), ElTColl.at(Idx_InConvEl.at(it_e)).Phi());
        }
      }
      else if(ElTColl.at(Idx_InConvEl.at(it_e)).Charge()<0 && Idx_vx>0){
        FillHist("Mlv_InConvEl"+Label, (ElTColl.at(Idx_InConvEl.at(it_e))+TruthColl.at(Idx_vx)).M(), 1, 200., 0., 200.);

        float Mlv = (ElTColl.at(Idx_InConvEl.at(it_e))+TruthColl.at(Idx_vx)).M();
        if(fabs(Mlv-81)<1.){
          int GenIdx_El = GenMatchedIdx(ElTColl.at(Idx_InConvEl.at(it_e)), TruthColl);
          PrintGen(TruthColl);
          printf("IdxEl:%d, Idxv:%d, ElPt:%.2f, ElEta:%.2f, ElPhi:%.2f\n", GenIdx_El, Idx_vx,
                  ElTColl.at(Idx_InConvEl.at(it_e)).Pt(), ElTColl.at(Idx_InConvEl.at(it_e)).Eta(), ElTColl.at(Idx_InConvEl.at(it_e)).Phi());
        }
      }
      FillHist("Iso_InConvEl"+Label, ElTColl.at(Idx_InConvEl.at(it_e)).RelIso(), 1, 50, 0., 10.);
    }

    for(unsigned int it_e=0; it_e<Idx_ExConvEl.size(); it_e++){
      if(ElTColl.at(Idx_ExConvEl.at(it_e)).Charge()>0 && Idx_v>0){
        FillHist("Mlv_ExConvEl"+Label, (ElTColl.at(Idx_ExConvEl.at(it_e))+TruthColl.at(Idx_v)).M(), 1, 200., 0., 200.);
      }
      else if(ElTColl.at(Idx_ExConvEl.at(it_e)).Charge()<0 && Idx_vx>0){
        FillHist("Mlv_ExConvEl"+Label, (ElTColl.at(Idx_ExConvEl.at(it_e))+TruthColl.at(Idx_vx)).M(), 1, 200., 0., 200.);
      }
      FillHist("Iso_ExConvEl"+Label, ElTColl.at(Idx_ExConvEl.at(it_e)).RelIso(), 1, 50, 0., 10.);
    }
  }

}



void GenMatchingValid::CheckPhotonMatching(
vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Electron>& ElTColl, vector<Electron>& ElLColl,
vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& ev, vector<Gen>& TruthColl, float weight, TString Label)
{
  bool CheckIsFinalPhotonSt23=true, CheckdRGmEl=false, CheckPhotonStatus=false;

  if(CheckPhotonStatus){

    //bool IsPhotonFinalSt23=IsFinalPhotonSt23(TruthColl);
    bool HasSt20sGm=false, HasVHtGm=false, HasLepMom=false, HasHadMom=false;
    int PhotonIdx=-1;
    for(unsigned int it_gen=2; it_gen<TruthColl.size(); it_gen++){
      if(TruthColl.at(it_gen).PID()!=22) continue;
      if(!(TruthColl.at(it_gen).Status()==1 or TruthColl.at(it_gen).Status()==23)) continue;
      int LastSelfIdx = LastSelfMotherIdx(it_gen,TruthColl);
      int MotherIdx   = FirstNonSelfMotherIdx(it_gen,TruthColl);
      int aMPID       = abs(TruthColl.at(MotherIdx).PID());
      int Status_Orig  = TruthColl.at(LastSelfIdx).Status();
      //int Status_MLast = TruthColl.at(MotherIdx).Status();
 
      PhotonIdx = it_gen;
      if     (Status_Orig>20 && Status_Orig<30) HasSt20sGm=true;
      else if((aMPID>=23 && aMPID<=25) or aMPID==6) HasVHtGm=true;
      else if(aMPID>10 && aMPID<20) HasLepMom=true;
      else if(aMPID>50) HasHadMom=true;

    }
    if     (HasSt20sGm){ FillHist("PhotonCategory"+Label, 0., 1, 10, 0., 10.); }
    else if(HasVHtGm)  { FillHist("PhotonCategory"+Label, 1., 1, 10, 0., 10.); }
    else if(HasLepMom) { FillHist("PhotonCategory"+Label, 2., 1, 10, 0., 10.); }
    else if(HasHadMom) { FillHist("PhotonCategory"+Label, 3., 1, 10, 0., 10.); }
    else if(PhotonIdx==-1){ FillHist("PhotonCategory"+Label, 4., 1, 10, 0., 10.); }
    else{
      FillHist("PhotonCategory"+Label, 5., 1, 10, 0., 10.);
      PrintGen(TruthColl);
      printf("PhotonIdx: %d\n", PhotonIdx);
    }
  }

  if(CheckdRGmEl){
    for(unsigned int it_e=0; it_e<ElTColl.size(); it_e++){
      if( !(ElTColl.at(it_e).Pt()>10. && fabs(ElTColl.at(it_e).Eta())<2.5) ) continue;
      if     (Label=="_CBPOGM" && !ElTColl.at(it_e).passMediumID()) continue;
      else if(Label=="_CBPOGT" && !ElTColl.at(it_e).passTightID()) continue;
      else if(Label=="_TESTIDNoConv"){
        if(! ElTColl.at(it_e).passMVAID_noIso_WP90() ) continue;
        if(! (ElTColl.at(it_e).RelIso()<0.1)         ) continue;
        if(! (fabs(ElTColl.at(it_e).dXY())<0.025 && fabs(ElTColl.at(it_e).dZ())<0.1) ) continue;
        if(! (ElTColl.at(it_e).dXYerr()>0. && fabs(ElTColl.at(it_e).dXY()/ElTColl.at(it_e).dXYerr())<4.) ) continue;
        if(! (ElTColl.at(it_e).IP3Derr()>0. && fabs(ElTColl.at(it_e).IP3D()/ElTColl.at(it_e).IP3Derr())<4.) ) continue;
        if(! ElTColl.at(it_e).Pass_CaloIdL_TrackIdL_IsoVL16() ) continue;
      }
      else if(Label=="_TESTIDNoIP"){
        if(! ElTColl.at(it_e).passMVAID_noIso_WP90() ) continue;
        if(! (ElTColl.at(it_e).RelIso()<0.1)         ) continue;
        if(! ElTColl.at(it_e).PassConversionVeto() ) continue;
        if(! ElTColl.at(it_e).Pass_CaloIdL_TrackIdL_IsoVL16() ) continue;
      }

      int MatchedTruthIdx = GenMatchedIdx(ElTColl.at(it_e),TruthColl);
      if(MatchedTruthIdx>0) continue;
  
      float PTthreshold=8., dRmax=0.2, dRmin=999.; //float dPtRel=0.5;
      int NearPhotonIdx=-1;
      for(unsigned int it_gen=2; it_gen<TruthColl.size(); it_gen++){
        if( TruthColl.at(it_gen).MotherIndex()<0   ) continue;
        if( !(TruthColl.at(it_gen).PID()==22 && (TruthColl.at(it_gen).Status()==1 || TruthColl.at(it_gen).Status()==23)) ) continue;
        if( TruthColl.at(it_gen).Status()==23 && !IsFinalPhotonSt23(TruthColl) ) continue;
        if( TruthColl.at(it_gen).Pt()<PTthreshold  ) continue;
        if( ElTColl.at(it_e).DeltaR(TruthColl.at(it_gen))>dRmax ) continue;
        //if( !(ElTColl.at(it_e).Pt()/TruthColl.at(it_gen).Pt()>(1.-dPtRel) && ElTColl.at(it_e).Pt()/TruthColl.at(it_gen).Pt()<(1.+dPtRel)) ) continue;
        if( ElTColl.at(it_e).DeltaR(TruthColl.at(it_gen))<dRmin ){ dRmin=ElTColl.at(it_e).DeltaR(TruthColl.at(it_gen)); NearPhotonIdx=it_gen; }
      }

      int NearPhotonType = GetPhotonType_JH(NearPhotonIdx, TruthColl);
      if(NearPhotonIdx!=-1){
        FillHist("PtGm"+Label, TruthColl.at(NearPhotonIdx).Pt(), 1, 20, 0., 100.);
        FillHist("EtaGm"+Label, TruthColl.at(NearPhotonIdx).Eta(), 1, 20, -5., 5.);
        FillHist("mindRElGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
        FillHist("PtRelElGm"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        float ElPt=ElTColl.at(it_e).Pt();
        if(ElPt>10 && ElPt<25){
          FillHist("mindRElGm_10To25"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_10To25"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else if(ElPt>=25 && ElPt<50){
          FillHist("mindRElGm_25To50"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_25To50"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else if(ElPt>=50 && ElPt<100){
          FillHist("mindRElGm_50To100"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_50To100"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else if(ElPt>=100 && ElPt<200){
          FillHist("mindRElGm_100To200"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_100To200"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else{
          FillHist("mindRElGm_200ToInf"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_200ToInf"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
      }

      if(NearPhotonIdx!=-1 && NearPhotonType>0){
        FillHist("mindRElGm_PrGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
        FillHist("PtRelElGm_PrGm"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        float ElPt=ElTColl.at(it_e).Pt();
        if(ElPt>10 && ElPt<25){
          FillHist("mindRElGm_10To25_PrGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_10To25_PrGm"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else if(ElPt>=25 && ElPt<50){
          FillHist("mindRElGm_25To50_PrGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_25To50_PrGm"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else if(ElPt>=50 && ElPt<100){
          FillHist("mindRElGm_50To100_PrGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_50To100_PrGm"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else if(ElPt>=100 && ElPt<200){
          FillHist("mindRElGm_100To200_PrGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_100To200_PrGm"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
        else{
          FillHist("mindRElGm_200ToInf_PrGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(NearPhotonIdx)), 1, 20, 0., 0.2);
          FillHist("PtRelElGm_200ToInf_PrGm"+Label, ElTColl.at(it_e).Pt()/TruthColl.at(NearPhotonIdx).Pt(), 1, 200, 0., 2.);
        }
      }

    }
  }

  if(CheckIsFinalPhotonSt23){
    if(!IsFinalPhotonSt23(TruthColl)){
      FillHist("Bool_IsFinalPhotonSt23"+Label, 0., 1., 2, 0., 2.); 
    }
    else{
      FillHist("Bool_IsFinalPhotonSt23"+Label, 1., 1., 2, 0., 2.);
      int St23GmIdx=-1;
      for(unsigned int it_gen=2; it_gen<TruthColl.size(); it_gen++){
        if( TruthColl.at(it_gen).Status()!=23 ) continue;
        if( TruthColl.at(it_gen).PID()!=22 ) continue;
        St23GmIdx=it_gen;
      } 
      int NEl_R04=0, NEl_R02=0, NEl_R01=0;
      for(unsigned int it_e=0; it_e<ElTColl.size(); it_e++){
        if(ElTColl.at(it_e).DeltaR(TruthColl.at(St23GmIdx))<0.1){ NEl_R01++; }
        if(ElTColl.at(it_e).DeltaR(TruthColl.at(St23GmIdx))<0.2){ NEl_R02++; }
        if(ElTColl.at(it_e).DeltaR(TruthColl.at(St23GmIdx))<0.4){ NEl_R04++; }
        FillHist("dRElGm"+Label, ElTColl.at(it_e).DeltaR(TruthColl.at(St23GmIdx)), 1, 50, 0., 5.);
      }
      if(TruthColl.at(St23GmIdx).Pt()>20. && fabs(TruthColl.at(St23GmIdx).Eta())<2.5){
        FillHist("NEl_R01"+Label, NEl_R01, 1, 10, 0., 10.);
        FillHist("NEl_R02"+Label, NEl_R02, 1, 10, 0., 10.);
        FillHist("NEl_R04"+Label, NEl_R04, 1, 10, 0., 10.);
      }
        
      bool PrintLog=true;
      if(PrintLog){ PrintGen(TruthColl); printf("St23PhotonIdx: %d\n", St23GmIdx);}
    }
  }
}



void GenMatchingValid::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

GenMatchingValid::GenMatchingValid(){

}

GenMatchingValid::~GenMatchingValid(){

}


