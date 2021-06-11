#include "ControlPlots.h"

void ControlPlots::initializeAnalyzer(){

  TriLep=false, TetraLep=false, SS2l=false, TopCR_OS2l=false;
  TopBSrc=false, SB_SS2L=false, CFlip=false, ConvCR=false, FkCR3l=false, ConvVar=false;
  FakeRun=false, ConvRun=false, FlipRun=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("SS2l"))       SS2l       = true;
    if(Userflags.at(i).Contains("TriLep"))     TriLep     = true;
    if(Userflags.at(i).Contains("TetraLep"))   TetraLep   = true; 
    if(Userflags.at(i).Contains("TopCR_OS2l")) TopCR_OS2l = true;
    if(Userflags.at(i).Contains("SB_SS2L"))    SB_SS2L    = true;
    if(Userflags.at(i).Contains("TopBSrc"))    TopBSrc    = true;
    if(Userflags.at(i).Contains("CFlip"))      CFlip      = true;
    if(Userflags.at(i).Contains("FkCR3l"))     FkCR3l     = true;
    if(Userflags.at(i).Contains("ConvCR"))     ConvCR     = true;
    if(Userflags.at(i).Contains("ConvVar"))    ConvVar    = true;
    if(Userflags.at(i).Contains("FakeRun"))    FakeRun    = true; 
    if(Userflags.at(i).Contains("ConvRun"))    ConvRun    = true; 
    if(Userflags.at(i).Contains("FlipRun"))    FlipRun    = true; 
    if(Userflags.at(i).Contains("SystRun"))    SystRun    = true; 
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
  vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

  InitializeTreeVars();
  InitializeReader(MVAreader_Mu, "TMVAClassification_MN20_Mu_BDTG.weights.xml");
  InitializeReader(MVAreader_El, "TMVAClassification_MN20_El_BDTG.weights.xml");
  outfile->cd();

}


void ControlPlots::executeEvent(){


  Event ev = GetEvent();
  float weight = 1., w_GenNorm=1., w_BR=1., w_PU=1.;
  if(!IsDATA){
    w_GenNorm = ev.MCweight()*weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    w_BR      = GetBRWeight();
    w_PU      = GetPileUpWeight(nPileUp, 0);
    weight *= w_GenNorm * w_BR * w_PU;
  }
  FillHist("CutFlow", 0., weight, 20, 0., 20.);

  bool PassTrig=false;
  if(TetraLep or TriLep){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu) or ev.PassTrigger(TrigList_DblEG);
    else{
      if     (MuEG)  PassTrig = ev.PassTrigger(TrigList_MuEG);
      else if(DblMu) PassTrig = (!ev.PassTrigger(TrigList_MuEG)) and ev.PassTrigger(TrigList_DblMu);
      else if(DblEG) PassTrig = (!(ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu))) and ev.PassTrigger(TrigList_DblEG);
    }
  }
  else if(SS2l or TopCR_OS2l){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_DblMu) or ev.PassTrigger(TrigList_DblEG) or ev.PassTrigger(TrigList_MuEG);
    else{
      if     (DblMu) PassTrig = ev.PassTrigger(TrigList_DblMu);
      else if(DblEG) PassTrig = ev.PassTrigger(TrigList_DblEG);
      else if(MuEG)  PassTrig = ev.PassTrigger(TrigList_MuEG);
    }
  }
  if(!PassTrig) return;
  FillHist("CutFlow", 1., weight, 20, 0., 20.);
  if(!PassMETFilter()) return;
  FillHist("CutFlow", 2., weight, 20, 0., 20.);

  bool PreCutPass=false;
  vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);
  if(TriLep and (muonPreColl.size()+electronPreColl.size())>2 ) PreCutPass=true;
  if(TetraLep and (muonPreColl.size()+electronPreColl.size())>3 ) PreCutPass=true;
  if(SS2l){
    int NPreMu=muonPreColl.size(), NPreEl=electronPreColl.size();
    if( (NPreMu+NPreEl)>2 ) PreCutPass=true;
    else if(NPreMu==2 and SumCharge(muonPreColl)!=0) PreCutPass=true;
    else if(NPreEl==2 and SumCharge(electronPreColl)!=0) PreCutPass=true;
  }
  if(((SS2l && IsDATA && FlipRun) or TopCR_OS2l) && (muonPreColl.size()+electronPreColl.size())>1) PreCutPass=true;
  
  if(!PreCutPass) return;


  TString IDSSLabel = "SS", TLabel = FakeRun? "L":"T";
  //TString IDSSLabel = SS2l or TopCR_OS2l? "SS":"";
  vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TopHN17"+TLabel, 10., 2.4);
  vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+TLabel, 15., 2.5);
  vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"L", 15., 2.5);
  vector<Muon>     muonVetoColl      = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronVetoColl  = SelectElectrons(electronPreColl, "TopHN17L", 10., 2.5);


  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  vector<Jet> jetPreColl = GetAllJets();
  sort(jetPreColl.begin(), jetPreColl.end(), PtComparing);
  vector<Jet> jetColl  = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVeto");
  vector<Jet> bjetColl = SelBJets(jetColl, param_jets);


  Particle vMET = ev.GetMETVector();
  Particle vMET_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);
  //Particle vMET_T1xy = GetvMET("T1xyCorr");


  vector<Gen> truthColl;


  bool EventCand = false;
  if(SS2l or TopCR_OS2l){ EventCand = (muonTightColl.size()+electronTightColl.size())==2; }
  if(TriLep){ EventCand = (muonLooseColl.size()+electronLooseColl.size())==3; }
  if(TetraLep){ EventCand = (muonLooseColl.size()+electronLooseColl.size())>3; }

  float w_TopPtRW = 1., w_Prefire = 1., sf_Trig = 1., w_FR=1.;
  float sf_MuTk = 1., sf_MuID = 1., sf_MuIso = 1., sf_ElReco = 1., sf_ElID = 1., sf_BTag = 1.;
//  float w_CF = IsDATA && FlipRun? 0.:1.;
  if((!IsDATA) and EventCand){
    //if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    //w_TopPtRW = mcCorr->GetTopPtReweight(truthColl);
    sf_MuID   = GetMuonSF(muonTightColl, "TopHNID_TkMu", "ID");
    //sf_MuIso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    sf_ElReco = GetElectronSF(electronTightColl, "", "Reco");
    sf_ElID   = GetElectronSF(electronTightColl, "TopHNID"+IDSSLabel, "ID");
    sf_BTag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    TString SFKey_Trig_All = muonTightColl.size()==2? "DiMuIso_HNTopID":electronTightColl.size()==2? "DiElIso_HNTopIDSS":"EMuIso_HNTopIDSS"; 
    sf_Trig   = SS2l or TopCR_OS2l? mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig_All, ""):1.;
    w_Prefire = GetPrefireWeight(0);
    //cout<<" w_PU:"<<w_PU<<" w_Prefire:"<<w_Prefire<<" sf_Trig:"<<sf_Trig<<endl;
    //cout<<"sf_MuTk"<<sf_MuTk<<" sf_MuID:"<<sf_MuID<<" sf_MuIso:"<<sf_MuIso<<" sf_ElReco:"<<sf_ElReco<<" sf_ElID:"<<sf_ElID<<" sf_BTag:"<<sf_BTag<<endl;
  }
  if(IsDATA && FakeRun && EventCand){
    w_FR      = CalcTestFakeWeight(muonLooseColl, electronLooseColl,
                  "TopHN17T", "TopHN17L", "TopHN17"+IDSSLabel+"T", "TopHN17"+IDSSLabel+"L", bjetColl.size(), 0);
  }
//  if(IsDATA && FlipRun && EventCand){
//    for(unsigned int ie=0; ie<electronTightColl.size(); ie++){
//      w_CF += GetCFRSF(electronTightColl.at(ie), "App2Bin2_Fin", "DataEff");
//    }
//  }
  weight *= w_TopPtRW * w_Prefire * sf_Trig * w_FR;
  weight *= sf_MuTk * sf_MuID * sf_MuIso * sf_ElReco * sf_ElID * sf_BTag;
//  weight *= w_CF;

 
  if(SS2l){
    if(SB_SS2L){
      if(!(IsDATA && FlipRun)){
        MakePlotSS2L(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight, "_SB_M"); 
      }
      else{
        for(unsigned ie=0; ie<electronTightColl.size(); ie++){
          double w_CFN = GetCFRSF(electronTightColl.at(ie), "App2Bin2_Fin", "DataEff");
          double PTN   = GetFlipCorrPT(electronTightColl.at(ie), "App2Bin2_Fin", "DataScaleSmear");
          std::vector<Electron> ElTCollN = electronTightColl;
          ElTCollN.at(ie).SetPtEtaPhiM(PTN, electronTightColl.at(ie).Eta(), electronTightColl.at(ie).Phi(), electronTightColl.at(ie).M());
          MakePlotSS2L(muonTightColl, muonLooseColl, muonVetoColl, ElTCollN, electronLooseColl, electronVetoColl,
                          jetColl, bjetColl, vMET_xyCorr, weight*w_CFN, "_SB_M");
        }
      }
    }
    if(TopBSrc){ CheckTopBSrc(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                              jetColl, bjetColl, vMET_xyCorr, weight, ""); }
    if(CFlip){
      if(!(IsDATA && FlipRun)){
        CheckChargeFlip(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight, "");
      }
      else{
        for(unsigned ie=0; ie<electronTightColl.size(); ie++){
          double w_CFN = GetCFRSF(electronTightColl.at(ie), "App2Bin2_Fin", "DataEff");
          double PTN   = GetFlipCorrPT(electronTightColl.at(ie), "App2Bin2_Fin", "DataScaleSmear");
          std::vector<Electron> ElTCollN = electronTightColl;
          ElTCollN.at(ie).SetPtEtaPhiM(PTN, electronTightColl.at(ie).Eta(), electronTightColl.at(ie).Phi(), electronTightColl.at(ie).M());
          CheckChargeFlip(muonTightColl, muonLooseColl, muonVetoColl, ElTCollN, electronLooseColl, electronVetoColl,
                          jetColl, bjetColl, vMET_xyCorr, weight*w_CFN, "");
          if(SystRun){
            CheckChargeFlip(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                            jetColl, bjetColl, vMET_xyCorr, weight*w_CFN, "_SystUp_CFPT");
            CheckChargeFlip(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                            jetColl, bjetColl, vMET_xyCorr, weight*w_CFN, "_SystDown_CFPT");
          }
        }
      }
    }
  }
  if(TopCR_OS2l){
    PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight, "");
  }
  if(TriLep){
    if(ConvCR){
      CheckConvCR(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                  jetColl, bjetColl, vMET_xyCorr, weight, "");
    }
    if(ConvVar){
      vector<Muon>     MuTIDLIsoColl         = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
      vector<Electron> ElTIDTIsoNoConvColl   = SelectElectrons(electronPreColl, "TopHN17SSTIDTIsoNoCvNoMHit", 15., 2.5);
      vector<Electron> ElTIDLIsoNoConvColl   = SelectElectrons(electronPreColl, "TopHN17SSTIDLIsoNoCvNoMHit", 15., 2.5);
      vector<Electron> ElLIDLIsoNoConvColl   = SelectElectrons(electronPreColl, "TopHN17SSLIDLIsoNoCvNoMHit", 15., 2.5);
      vector<Electron> ElLIDLIsoNoConv10Coll = SelectElectrons(electronPreColl, "TopHN17SSLIDLIsoNoCvNoMHit", 10., 2.5);
      vector<Jet> jetNewColl  = SelectJets(jetPreColl, MuTIDLIsoColl, ElLIDLIsoNoConv10Coll, "tight", 25., 2.4, "LVeto");
      vector<Jet> bjetNewColl = SelBJets(jetNewColl, param_jets);

      CheckConvVar(MuTIDLIsoColl, MuTIDLIsoColl, MuTIDLIsoColl, ElTIDTIsoNoConvColl, ElLIDLIsoNoConvColl, ElLIDLIsoNoConv10Coll,
                   jetNewColl, bjetNewColl, vMET_xyCorr, weight, "_LCv");
      CheckConvVar(MuTIDLIsoColl, MuTIDLIsoColl, MuTIDLIsoColl, ElTIDLIsoNoConvColl, ElLIDLIsoNoConvColl, ElLIDLIsoNoConv10Coll,
                   jetNewColl, bjetNewColl, vMET_xyCorr, weight, "_LIsoCv");
      CheckConvVar(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_xyCorr, weight, "_TID");
    }
    if(FkCR3l){
      CheckFkCR3l(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                  jetColl, bjetColl, vMET_xyCorr, weight, "");
    }
  }
  if(SystRun && ((!IsDATA) or FakeRun)){
    vector<Jet> jetJESUpColl    = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJESUp");
    vector<Jet> jetJESDownColl  = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJESDown");
    vector<Jet> jetJERUpColl    = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJERUp");
    vector<Jet> jetJERDownColl  = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJERDown");
    vector<Jet> bjetJESUpColl   = SelBJets(jetJESUpColl, param_jets);
    vector<Jet> bjetJESDownColl = SelBJets(jetJESDownColl, param_jets);
    vector<Jet> bjetJERUpColl   = SelBJets(jetJERUpColl, param_jets);
    vector<Jet> bjetJERDownColl = SelBJets(jetJERDownColl, param_jets);
    Particle vMET_T1xy_JESUp    = GetvMET("T1xyCorr", "SystUpJES");
    Particle vMET_T1xy_JESDown  = GetvMET("T1xyCorr", "SystDownJES");
    Particle vMET_T1xy_JERUp    = GetvMET("T1xyCorr", "SystUpJER");
    Particle vMET_T1xy_JERDown  = GetvMET("T1xyCorr", "SystDownJER");
    Particle vMET_T1xy_UnclUp   = GetvMET("T1xyCorr", "SystUpUncl");
    Particle vMET_T1xy_UnclDown = GetvMET("T1xyCorr", "SystDownUncl");

    float w_PUUp=1., w_PUDown=1., w_PrefireUp=1., w_PrefireDown=1., w_FRUp=1., w_FRDown=1.;
    float sf_ElRecoUp=1., sf_ElRecoDown=1.;
    float sf_BTag_JESUp=1., sf_BTag_JESDown=1., sf_BTag_JERUp=1., sf_BTag_JERDown=1.;
    float sf_BTag_LTagUp=1., sf_BTag_LTagDown=1., sf_BTag_HTagUp=1., sf_BTag_HTagDown=1.;
    if(!IsDATA && EventCand){
      w_PrefireUp   = GetPrefireWeight(1), w_PrefireDown = GetPrefireWeight(-1);
      w_PUUp        = weight!=0.? GetPileUpWeight(nPileUp, 1)/GetPileUpWeight(nPileUp, 0):0.;
      w_PUDown      = weight!=0.? GetPileUpWeight(nPileUp,-1)/GetPileUpWeight(nPileUp, 0):0.;
      sf_ElRecoUp   = GetElectronSF(electronTightColl, "", "RecoSystUp");
      sf_ElRecoDown = GetElectronSF(electronTightColl, "", "RecoSystDown");
      sf_BTag_JESUp    = mcCorr->GetBTaggingReweight_1a(jetJESUpColl, param_jets);
      sf_BTag_JESDown  = mcCorr->GetBTaggingReweight_1a(jetJESDownColl, param_jets);
      sf_BTag_JERUp    = mcCorr->GetBTaggingReweight_1a(jetJERUpColl, param_jets);
      sf_BTag_JERDown  = mcCorr->GetBTaggingReweight_1a(jetJERDownColl, param_jets);
      sf_BTag_LTagUp   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets, "SystUpLTag");
      sf_BTag_LTagDown = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets, "SystDownLTag");
      sf_BTag_HTagUp   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets, "SystUpHTag");
      sf_BTag_HTagDown = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets, "SystDownHTag");
    }
    if(IsDATA && FakeRun && EventCand){
      w_FRUp   = CalcTestFakeWeight(muonLooseColl, electronLooseColl,
                  "TopHN17T", "TopHN17L", "TopHN17"+IDSSLabel+"T", "TopHN17"+IDSSLabel+"L", bjetColl.size(),  1);
      w_FRDown = CalcTestFakeWeight(muonLooseColl, electronLooseColl,
                  "TopHN17T", "TopHN17L", "TopHN17"+IDSSLabel+"T", "TopHN17"+IDSSLabel+"L", bjetColl.size(), -1);
    }

    float w_base = w_GenNorm * w_BR * w_TopPtRW * sf_MuTk * sf_MuID * sf_MuIso * sf_ElID * sf_Trig;
    float weight_PUUp        = w_base * w_PUUp   * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag;
    float weight_PUDown      = w_base * w_PUDown * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag;
    float weight_PrefireUp   = w_base * w_PU     * w_PrefireUp   * w_FR     * sf_ElReco     * sf_BTag;
    float weight_PrefireDown = w_base * w_PU     * w_PrefireDown * w_FR     * sf_ElReco     * sf_BTag;
    float weight_FRUp        = w_base * w_PU     * w_Prefire     * w_FRUp   * sf_ElReco     * sf_BTag;
    float weight_FRDown      = w_base * w_PU     * w_Prefire     * w_FRDown * sf_ElReco     * sf_BTag;
    float weight_ElRecoUp    = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElRecoUp   * sf_BTag;
    float weight_ElRecoDown  = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElRecoDown * sf_BTag;
    float weight_JESUp       = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_JESUp;
    float weight_JESDown     = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_JESDown;
    float weight_JERUp       = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_JERUp;
    float weight_JERDown     = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_JERDown;
    float weight_LTagUp      = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_LTagUp;
    float weight_LTagDown    = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_LTagDown;
    float weight_HTagUp      = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_HTagUp;
    float weight_HTagDown    = w_base * w_PU     * w_Prefire     * w_FR     * sf_ElReco     * sf_BTag_HTagDown;

    if(SB_SS2L){
      if(!IsDATA){
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PUUp, "_SB_M_SystUp_PU");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PUDown, "_SB_M_SystDown_PU");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PrefireUp, "_SB_M_SystUp_Pref");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PrefireDown, "_SB_M_SystDown_Pref");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_ElRecoUp, "_SB_M_SystUp_ElReco");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_ElRecoDown, "_SB_M_SystDown_ElReco");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, weight_JESUp, "_SB_M_SystUp_JES");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, weight_JESDown, "_SB_M_SystDown_JES");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, weight_JERUp, "_SB_M_SystUp_JER");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, weight_JERDown, "_SB_M_SystDown_JER");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_T1xy_UnclUp, weight, "_SB_M_SystUp_Uncl");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_T1xy_UnclDown, weight, "_SB_M_SystDown_Uncl");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_LTagUp, "_SB_M_SystUp_LTag");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_LTagDown, "_SB_M_SystDown_LTag");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_HTagUp, "_SB_M_SystUp_HTag");
        MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_HTagDown, "_SB_M_SystDown_HTag");
      }
      else if(IsDATA && FakeRun){
          MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                       jetColl, bjetColl, vMET_xyCorr, weight_FRUp, "_SB_M_SystUp_FR");
          MakePlotSS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                       jetColl, bjetColl, vMET_xyCorr, weight_FRDown, "_SB_M_SystDown_FR");
      }

    }

    if(TopCR_OS2l){
      if(!IsDATA){
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PUUp, "_SystUp_PU");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PUDown, "_SystDown_PU");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PrefireUp, "_SystUp_Pref");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_PrefireDown, "_SystDown_Pref");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_ElRecoUp, "_SystUp_ElReco");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_ElRecoDown, "_SystDown_ElReco");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, weight_JESUp, "_SystUp_JES");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, weight_JESDown, "_SystDown_JES");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, weight_JERUp, "_SystUp_JER");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, weight_JERDown, "_SystDown_JER");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_T1xy_UnclUp, weight, "_SystUp_Uncl");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_T1xy_UnclDown, weight, "_SystDown_Uncl");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_LTagUp, "_SystUp_LTag");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_LTagDown, "_SystDown_LTag");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_HTagUp, "_SystUp_HTag");
        PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight_HTagDown, "_SystDown_HTag");
      }
      else if(IsDATA && FakeRun){
          PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                       jetColl, bjetColl, vMET_xyCorr, weight_FRUp, "_SystUp_FR");
          PlotTop2LCR_OS2L(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                       jetColl, bjetColl, vMET_xyCorr, weight_FRDown, "_SystDown_FR");
      }
    }

    if(FkCR3l){
      if(!IsDATA){
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PUUp, "_SystUp_PU");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PUDown, "_SystDown_PU");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PrefireUp, "_SystUp_Pref");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PrefireDown, "_SystDown_Pref");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_ElRecoUp, "_SystUp_ElReco");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_ElRecoDown, "_SystDown_ElReco");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, weight_JESUp, "_SystUp_JES");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, weight_JESDown, "_SystDown_JES");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, weight_JERUp, "_SystUp_JER");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, weight_JERDown, "_SystDown_JER");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_T1xy_UnclUp, weight, "_SystUp_Uncl");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_T1xy_UnclDown, weight, "_SystDown_Uncl");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_LTagUp, "_SystUp_LTag");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_LTagDown, "_SystDown_LTag");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_HTagUp, "_SystUp_HTag");
        CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_HTagDown, "_SystDown_HTag");
      }
      else if(IsDATA && FakeRun){
          CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                      jetColl, bjetColl, vMET_xyCorr, weight_FRUp, "_SystUp_FR");
          CheckFkCR3l(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                      jetColl, bjetColl, vMET_xyCorr, weight_FRDown, "_SystDown_FR");
      }
    }
    if(ConvCR){
      if(!IsDATA){
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PUUp, "_SystUp_PU");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PUDown, "_SystDown_PU");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PrefireUp, "_SystUp_Pref");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_PrefireDown, "_SystDown_Pref");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_ElRecoUp, "_SystUp_ElReco");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_ElRecoDown, "_SystDown_ElReco");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, weight_JESUp, "_SystUp_JES");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, weight_JESDown, "_SystDown_JES");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, weight_JERUp, "_SystUp_JER");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, weight_JERDown, "_SystDown_JER");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_T1xy_UnclUp, weight, "_SystUp_Uncl");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_T1xy_UnclDown, weight, "_SystDown_Uncl");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_LTagUp, "_SystUp_LTag");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_LTagDown, "_SystDown_LTag");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_HTagUp, "_SystUp_HTag");
        CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                    jetColl, bjetColl, vMET_xyCorr, weight_HTagDown, "_SystDown_HTag");
      }
      else if(IsDATA && FakeRun){
          CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                      jetColl, bjetColl, vMET_xyCorr, weight_FRUp, "_SystUp_FR");
          CheckConvCR(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                      jetColl, bjetColl, vMET_xyCorr, weight_FRDown, "_SystDown_FR");
      }
    }
    if(CFlip){
      if(!IsDATA){
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_PUUp, "_SystUp_PU");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_PUDown, "_SystDown_PU");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_PrefireUp, "_SystUp_Pref");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_PrefireDown, "_SystDown_Pref");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_ElRecoUp, "_SystUp_ElReco");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_ElRecoDown, "_SystDown_ElReco");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, weight_JESUp, "_SystUp_JES");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, weight_JESDown, "_SystDown_JES");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, weight_JERUp, "_SystUp_JER");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, weight_JERDown, "_SystDown_JER");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_T1xy_UnclUp, weight, "_SystUp_Uncl");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_T1xy_UnclDown, weight, "_SystDown_Uncl");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_LTagUp, "_SystUp_LTag");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_LTagDown, "_SystDown_LTag");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_HTagUp, "_SystUp_HTag");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_HTagDown, "_SystDown_HTag");
      }
      else if(IsDATA && FakeRun){
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_FRUp, "_SystUp_FR");
        CheckChargeFlip(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                        jetColl, bjetColl, vMET_xyCorr, weight_FRDown, "_SystDown_FR");
      }
    }
  }
}


void ControlPlots::CheckChargeFlip(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                   vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                   vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( FakeRun      and weight==0.  ) return; 
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NMuT==2){
    Label="_2M"+Label;
    if(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    float Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(DataYear>2016 && Mll<4) return; 
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      if(ConvRun && GenLepInfo>=100) return;
      if(FlipRun && (GenLepInfo>999 or GenLepInfo<100)) return;
    }
  }
  else if(NElT==2){
    Label="_2E"+Label;
    int aSumQ = abs(SumCharge(ElLColl));
    if(IsDATA && FlipRun){ if(aSumQ!=0) return; }
    else                 { if(aSumQ==0) return; }
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      //int IdxFlipped = GenLepInfo % 1000;
      if(ConvRun && GenLepInfo>=100) return;
      if(FlipRun && (GenLepInfo>999 or GenLepInfo<100)) return;
    }
  }

  int Nj=JetColl.size(), Nb=BJetColl.size();
  FillHist("h2D_Nj_Nb"+Label, Nj, Nb, weight, 10, 0., 10., 5, 0., 5.);
   
  float Mll=-1., PTl1=-1, Etal1=999., PTl2=-1., Etal2=-1., MTW=-1., MTllv=-1., MET=-1., dEtall=-1;
  bool IsOnZ=false, IsBJOrtho=false;
  if(NMuT==2){
    Mll   = (MuTColl.at(0)+MuTColl.at(1)).M();
    PTl1  = MuTColl.at(0).Pt();  Etal1 = MuTColl.at(0).Eta();
    PTl2  = MuTColl.at(1).Pt();  Etal2 = MuTColl.at(1).Eta();
    MTW   = MT(MuTColl.at(0), vMET);
    MTllv = MT((MuTColl.at(0)+MuTColl.at(1)), vMET);
  }
  else if(NElT==2){
    Mll   = (ElTColl.at(0)+ElTColl.at(1)).M();
    PTl1  = ElTColl.at(0).Pt();  Etal1 = ElTColl.at(0).Eta();
    PTl2  = ElTColl.at(1).Pt();  Etal2 = ElTColl.at(1).Eta();
    MTW   = MT(ElTColl.at(0), vMET);
    MTllv = MT((ElTColl.at(0)+ElTColl.at(1)), vMET);
  }
  MET = vMET.Pt(); dEtall = fabs(Etal1-Etal2);
  IsOnZ = fabs(Mll-91.2)<15.; IsBJOrtho = Nb==0 && Nj<3;

  
  if(IsBJOrtho){
    TString SelTag="_0Blt3J";
    FillHist("Mll"+SelTag+Label, Mll, weight, 40, 0., 200.);
    FillHist("PTl1"+SelTag+Label, PTl1, weight, 20, 0., 200.);
    FillHist("PTl2"+SelTag+Label, PTl2, weight, 20, 0., 200.);
    FillHist("Etal1"+SelTag+Label, Etal1, weight, 20, -5., 5.);
    FillHist("Etal2"+SelTag+Label, Etal2, weight, 20, -5., 5.);
    FillHist("dEtall"+SelTag+Label, dEtall, weight, 25, 0., 5.);
    FillHist("MET"+SelTag+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+SelTag+Label, MTW, weight, 20, 0., 200.);
    FillHist("MTllv"+SelTag+Label, MTllv, weight, 20, 0., 200.);
  }
  if(IsOnZ){
    TString SelTag="_OnZ";
    FillHist("NCount"+SelTag+Label, 0., weight, 1, 0., 1.);
    FillHist("Mll"+SelTag+Label, Mll, weight, 30, 60., 120.);
    FillHist("PTl1"+SelTag+Label, PTl1, weight, 20, 0., 200.);
    FillHist("PTl2"+SelTag+Label, PTl2, weight, 20, 0., 200.);
    FillHist("Etal1"+SelTag+Label, Etal1, weight, 20, -5., 5.);
    FillHist("Etal2"+SelTag+Label, Etal2, weight, 20, -5., 5.);
    FillHist("dEtall"+SelTag+Label, dEtall, weight, 25, 0., 5.);
    FillHist("MET"+SelTag+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+SelTag+Label, MTW, weight, 20, 0., 200.);
    FillHist("MTllv"+SelTag+Label, MTllv, weight, 20, 0., 200.);
  }

}


void ControlPlots::CheckFkCR3l(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                               vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                               vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  int NLepT=NMuT+NElT;
  if( FakeRun      and weight==0.  ) return; 
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if( NLepT!=3 ) return;
  bool PassTrigAcc=false;
  if( NMuT>1 && MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 ) PassTrigAcc=true;
  if( NElT>1 && ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15 ) PassTrigAcc=true;
  if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 ) PassTrigAcc=true;
  if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>25 ) PassTrigAcc=true;
  if(!PassTrigAcc) return;
  int Qtot_Mu = SumCharge(MuTColl), Qtot_El = SumCharge(ElTColl), Qtot_Lep = Qtot_Mu+Qtot_El;
  if( abs(Qtot_Lep)!=1 ) return;
  if( NMuT==2 && Qtot_Mu!=0 ) return;
  if( NElT==2 && Qtot_El!=0 ) return;
  if( BJetColl.size()>0 ) return;

  int NJet = JetColl.size();
  float PTl1_Z=0., PTl2_Z=0., PTl_Fk=0., Etal1_Z=0., Etal2_Z=0., Etal_Fk=0., dRlWj1=-1.;
  float MOSSF_Z=0., MTW=0., MOSSF1=0., MOSSF2=0., M3l=0.;
  if(NMuT==2 && NElT==1){
    Label   = "_2M1E"+Label;
    PTl1_Z  = MuTColl.at(0).Pt() , PTl2_Z  = MuTColl.at(1).Pt();
    Etal1_Z = MuTColl.at(0).Eta(), Etal2_Z = MuTColl.at(1).Eta();
    PTl_Fk  = ElTColl.at(0).Pt() , Etal_Fk = ElTColl.at(0).Eta();
    MOSSF_Z = (MuTColl.at(0)+MuTColl.at(1)).M(), MTW = MT(ElTColl.at(0), vMET), MOSSF1=MOSSF_Z; 
    M3l     = (MuTColl.at(0)+MuTColl.at(1)+ElTColl.at(0)).M();
    dRlWj1  = NJet>0? ElTColl.at(0).DeltaR(JetColl.at(0)):-1.;
  }
  else if(NMuT==1 && NElT==2){
    Label   = "_1M2E"+Label;
    PTl1_Z  = ElTColl.at(0).Pt() , PTl2_Z  = ElTColl.at(1).Pt();
    Etal1_Z = ElTColl.at(0).Eta(), Etal2_Z = ElTColl.at(1).Eta();
    PTl_Fk  = MuTColl.at(0).Pt() , Etal_Fk = MuTColl.at(0).Eta();
    MOSSF_Z = (ElTColl.at(0)+ElTColl.at(1)).M(), MTW = MT(MuTColl.at(0), vMET), MOSSF1=MOSSF_Z;
    M3l     = (ElTColl.at(0)+ElTColl.at(1)+MuTColl.at(0)).M();
    dRlWj1  = NJet>0? MuTColl.at(0).DeltaR(JetColl.at(0)):-1.;
  }
  else if(NMuT==3){
    Label = "_3M0E"+Label;
    int IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    MOSSF1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    MOSSF2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    int IdxL1_Z=-1, IdxL2_Z=-1, IdxL_Fk=-1;
    if(fabs(MOSSF1-91.2)<fabs(MOSSF2-91.2)){
      IdxL1_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS1).Pt()? IdxOS:IdxSS1;
      IdxL2_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS1).Pt()? IdxSS1:IdxOS;
      IdxL_Fk = IdxSS2;
    }
    else{
      IdxL1_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS2).Pt()? IdxOS:IdxSS2;
      IdxL2_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS2).Pt()? IdxSS2:IdxOS;
      IdxL_Fk = IdxSS1;
    }
    PTl1_Z  = MuTColl.at(IdxL1_Z).Pt() , PTl2_Z  = MuTColl.at(IdxL2_Z).Pt();
    Etal1_Z = MuTColl.at(IdxL1_Z).Eta(), Etal2_Z = MuTColl.at(IdxL2_Z).Eta();
    PTl_Fk  = MuTColl.at(IdxL_Fk).Pt() , Etal_Fk = MuTColl.at(IdxL_Fk).Eta();
    MOSSF_Z = (MuTColl.at(IdxL1_Z)+MuTColl.at(IdxL2_Z)).M(), MTW = MT(MuTColl.at(IdxL_Fk), vMET);
    M3l     = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M();
    dRlWj1  = NJet>0? MuTColl.at(IdxL_Fk).DeltaR(JetColl.at(0)):-1.;
  }
  else if(NElT==3){
    Label = "_0M3E"+Label;
    int IdxOS  = TriElChargeIndex(ElTColl, "OS");
    int IdxSS1 = TriElChargeIndex(ElTColl, "SS1");
    int IdxSS2 = TriElChargeIndex(ElTColl, "SS2");
    MOSSF1 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS1)).M();
    MOSSF2 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS2)).M();
    int IdxL1_Z=-1, IdxL2_Z=-1, IdxL_Fk=-1;
    if(fabs(MOSSF1-91.2)<fabs(MOSSF2-91.2)){
      IdxL1_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS1).Pt()? IdxOS:IdxSS1;
      IdxL2_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS1).Pt()? IdxSS1:IdxOS;
      IdxL_Fk = IdxSS2;
    }
    else{
      IdxL1_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS2).Pt()? IdxOS:IdxSS2;
      IdxL2_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS2).Pt()? IdxSS2:IdxOS;
      IdxL_Fk = IdxSS1;
    }
    PTl1_Z  = ElTColl.at(IdxL1_Z).Pt() , PTl2_Z  = ElTColl.at(IdxL2_Z).Pt();
    Etal1_Z = ElTColl.at(IdxL1_Z).Eta(), Etal2_Z = ElTColl.at(IdxL2_Z).Eta();
    PTl_Fk  = ElTColl.at(IdxL_Fk).Pt() , Etal_Fk = ElTColl.at(IdxL_Fk).Eta();
    MOSSF_Z = (ElTColl.at(IdxL1_Z)+ElTColl.at(IdxL2_Z)).M(), MTW = MT(ElTColl.at(IdxL_Fk), vMET);
    M3l     = (ElTColl.at(0)+ElTColl.at(1)+ElTColl.at(2)).M();
    dRlWj1  = NJet>0? ElTColl.at(IdxL_Fk).DeltaR(JetColl.at(0)):-1.;
  }
  else return;

  if(MOSSF1<12 or ((NMuT>2 or NElT>2) and MOSSF2<12)) return;
  bool OnZ=fabs(MOSSF_Z-91.2)<10.;
  //if(M3l<91.2+15.) return;

  if( ConvRun && (!IsDATA) ){
    vector<Gen> TruthColl = GetGens();
    int NFk=0, NPt20=0;
    for(unsigned int im=0; im<MuTColl.size(); im++){
      int LepType = GetLeptonType_JH(MuTColl.at(im), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      if(MuTColl.at(im).Pt()>20) NPt20++;
    }
    for(unsigned int ie=0; ie<ElTColl.size(); ie++){
      int LepType = GetLeptonType_JH(ElTColl.at(ie), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      if(ElTColl.at(ie).Pt()>20) NPt20++;
    }
    if(NFk>0) return;
    if(MCSample.Contains("DY") && NPt20>2) return;
    //else if(MCSample.Contains("ZG") && NPt20!=3) return;
  }

  Label = "_3l"+Label;
  FillHist("M3l"+Label, M3l, weight, 50, 0., 500.);
  FillHist("MOSSF1"+Label, MOSSF1, weight, 40, 0., 200.);
  FillHist("MOSSF2"+Label, MOSSF2, weight, 40, 0., 200.);
  FillHist("PTl1_Z"+Label, PTl1_Z, weight, 40, 0., 200.);
  FillHist("PTl2_Z"+Label, PTl2_Z, weight, 40, 0., 200.);
  FillHist("PTl_Fk"+Label, PTl_Fk, weight, 40, 0., 200.);
  FillHist("Etal1_Z"+Label, Etal1_Z, weight, 20, -5., 5.);
  FillHist("Etal2_Z"+Label, Etal2_Z, weight, 20, -5., 5.);
  FillHist("Etal_Fk"+Label, Etal_Fk, weight, 20, -5., 5.);
  FillHist("MOSSF_Z"+Label, MOSSF_Z, weight, 30, 60., 120.);
  FillHist("Nj"+Label, JetColl.size(), weight, 10, 0., 10.);
  FillHist("Nb"+Label, BJetColl.size(), weight, 5, 0., 5.);
  FillHist("MET"+Label, vMET.Pt(), weight, 20, 0., 200.);
  FillHist("MTW"+Label, MTW, weight, 20, 0., 200.); 
  if(NJet>0){
    FillHist("PTj1"+Label, JetColl.at(0).Pt(), weight, 20, 0., 200.);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("dRlj1"+Label, dRlWj1, weight, 25, 0., 5.);
  }
  if(OnZ){
    Label.ReplaceAll("_3l","_3lOnZ");
    FillHist("M3l"+Label, M3l, weight, 50, 0., 500.);
    FillHist("MOSSF1"+Label, MOSSF1, weight, 40, 0., 200.);
    FillHist("MOSSF2"+Label, MOSSF2, weight, 40, 0., 200.);
    FillHist("PTl1_Z"+Label, PTl1_Z, weight, 40, 0., 200.);
    FillHist("PTl2_Z"+Label, PTl2_Z, weight, 40, 0., 200.);
    FillHist("PTl_Fk"+Label, PTl_Fk, weight, 40, 0., 200.);
    FillHist("Etal1_Z"+Label, Etal1_Z, weight, 20, -5., 5.);
    FillHist("Etal2_Z"+Label, Etal2_Z, weight, 20, -5., 5.);
    FillHist("Etal_Fk"+Label, Etal_Fk, weight, 20, -5., 5.);
    FillHist("MOSSF_Z"+Label, MOSSF_Z, weight, 30, 60., 120.);
    FillHist("Nj"+Label, JetColl.size(), weight, 10, 0., 10.);
    FillHist("Nb"+Label, BJetColl.size(), weight, 5, 0., 5.);
    FillHist("MET"+Label, vMET.Pt(), weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 20, 0., 200.); 
    if(NJet>0){
      FillHist("PTj1"+Label, JetColl.at(0).Pt(), weight, 20, 0., 200.);
      FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
      FillHist("dRlj1"+Label, dRlWj1, weight, 25, 0., 5.);
    }
  }

}


void ControlPlots::CheckConvVar(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  int NLepT=NMuT+NElT;
  if( FakeRun      and weight==0.  ) return; 
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if( NLepT!=3 ) return;
  bool PassTrigAcc=false;
  if( NMuT>1 && MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 ) PassTrigAcc=true;
  if( NElT>1 && ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15 ) PassTrigAcc=true;
  if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 ) PassTrigAcc=true;
  if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>25 ) PassTrigAcc=true;
  if(!PassTrigAcc) return;
  int Qtot_Mu = SumCharge(MuTColl), Qtot_El = SumCharge(ElTColl), Qtot_Lep = Qtot_Mu+Qtot_El;
  if( abs(Qtot_Lep)!=1 ) return;
  if( NMuT==2 && Qtot_Mu!=0 ) return;
  if( NElT==2 && Qtot_El!=0 ) return;


  float MOSSF1=0., MOSSF2=0., M3l=0.;
  float NMissHit=-1., PassConvCut=0., MiniRelIso=-1.; 
  int IdxL1_Z=-1, IdxL2_Z=-1, IdxL_Cv=-1, LepType_Cv=0;
  if(NMuT==2 && NElT==1){
    Label   = "_2M1E"+Label;
    MOSSF1  = (MuTColl.at(0)+MuTColl.at(1)).M();
    M3l     = (MuTColl.at(0)+MuTColl.at(1)+ElTColl.at(0)).M();
    NMissHit    = ElTColl.at(0).NMissingHits();
    PassConvCut = ElTColl.at(0).PassConversionVeto()? 1.:0.;
    MiniRelIso  = ElTColl.at(0).MiniRelIso();
    IdxL1_Z=0, IdxL2_Z=1, IdxL_Cv=0;
  }
  else if(NMuT==1 && NElT==2){
    Label   = "_1M2E"+Label;
    MOSSF1     = (ElTColl.at(0)+ElTColl.at(1)).M();
    M3l        = (ElTColl.at(0)+ElTColl.at(1)+MuTColl.at(0)).M();
    MiniRelIso = MuTColl.at(0).MiniRelIso();
    IdxL1_Z=0, IdxL2_Z=1, IdxL_Cv=0;
  }
  else if(NMuT==3 && NElT==0){
    Label   = "_3M0E"+Label;
    int IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    MOSSF1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    MOSSF2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    //int IdxL1_Z=-1, IdxL2_Z=-1, IdxL_Cv=-1;
    if(fabs(MOSSF1-91.2)<fabs(MOSSF2-91.2)){
      IdxL1_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS1).Pt()? IdxOS:IdxSS1;
      IdxL2_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS1).Pt()? IdxSS1:IdxOS;
      IdxL_Cv = IdxSS2;
    }
    else{
      IdxL1_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS2).Pt()? IdxOS:IdxSS2;
      IdxL2_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS2).Pt()? IdxSS2:IdxOS;
      IdxL_Cv = IdxSS1;
    }
    M3l        = (MuTColl.at(IdxL1_Z)+MuTColl.at(IdxL2_Z)+MuTColl.at(IdxL_Cv)).M();
    MiniRelIso = MuTColl.at(IdxL_Cv).MiniRelIso();
  }
  else if(NMuT==0 && NElT==3){
    Label   = "_0M3E"+Label;
    int IdxOS  = TriElChargeIndex(ElTColl, "OS");
    int IdxSS1 = TriElChargeIndex(ElTColl, "SS1");
    int IdxSS2 = TriElChargeIndex(ElTColl, "SS2");
    MOSSF1 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS1)).M();
    MOSSF2 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS2)).M();
    //int IdxL1_Z=-1, IdxL2_Z=-1, IdxL_Cv=-1;
    if(fabs(MOSSF1-91.2)<fabs(MOSSF2-91.2)){
      IdxL1_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS1).Pt()? IdxOS:IdxSS1;
      IdxL2_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS1).Pt()? IdxSS1:IdxOS;
      IdxL_Cv = IdxSS2;
    }
    else{
      IdxL1_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS2).Pt()? IdxOS:IdxSS2;
      IdxL2_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS2).Pt()? IdxSS2:IdxOS;
      IdxL_Cv = IdxSS1;
    }
    M3l         = (ElTColl.at(IdxL1_Z)+ElTColl.at(IdxL2_Z)+ElTColl.at(IdxL_Cv)).M();
    NMissHit    = ElTColl.at(IdxL_Cv).NMissingHits();
    PassConvCut = ElTColl.at(IdxL_Cv).PassConversionVeto()? 1.:0.;
    MiniRelIso  = ElTColl.at(IdxL_Cv).MiniRelIso();
  }

  if( MOSSF1<12 or ((NMuT>2 or NElT>2) and MOSSF2<12) ) return;

  bool PassZGSel=true;
  if( BJetColl.size()>0 ) PassZGSel=false;
  if( fabs(MOSSF1-91.2)<10 or fabs(MOSSF2-91.2)<10 ) PassZGSel=false;
  if( fabs(M3l-91.2)>10. )  PassZGSel=false;
  if(PassZGSel) Label = "_ZGSel"+Label;

  int NFk=0, NCv=0, NPr=0, NElse=0, NPt20=0;
  if( ConvRun && (!IsDATA) ){
    vector<Gen> TruthColl = GetGens();
    for(unsigned int im=0; im<MuTColl.size(); im++){
      int LepType = GetLeptonType_JH(MuTColl.at(im), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0) NPr++;
      else NElse++;
      if(MuTColl.at(im).Pt()>20) NPt20++;
      if((NMuT==1 or NMuT==3) && ((int) im)==IdxL_Cv) LepType_Cv=LepType;
    }
    for(unsigned int ie=0; ie<ElTColl.size(); ie++){
      int LepType = GetLeptonType_JH(ElTColl.at(ie), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0) NPr++;
      else NElse++;
      if(ElTColl.at(ie).Pt()>20) NPt20++;
      if((NElT==1 or NElT==3) && ((int) ie)==IdxL_Cv) LepType_Cv=LepType;
    }
    //if(NFk>0) return;
    //if(MCSample.Contains("DY") && NPt20>2) return;
    //else if(MCSample.Contains("ZG") && NPt20!=3) return;
  }

  if( (!ConvRun) or NFk==0){
    FillHist("NPr"+Label, NPr, weight, 4, 0., 4.);
    FillHist("NFk"+Label, NFk, weight, 4, 0., 4.);
    FillHist("NCv"+Label, NCv, weight, 4, 0., 4.);
  }
  if(NElT!=1) return;
  if(LepType_Cv>0) Label="_Pr"+Label;
  else if(LepType_Cv==-1) Label="_FkMisID"+Label;
  else if(LepType_Cv<-1 && LepType_Cv>-5) Label="_FkHad"+Label;
  else if(LepType_Cv<-4) Label="_Cv"+Label;
  FillHist("MiniRelIso" +Label, MiniRelIso, weight, 20, 0., 0.4);
  FillHist("NMissHit"   +Label, NMissHit, weight, 10, 0., 10.);
  FillHist("PassConvCut"+Label, PassConvCut, weight, 2, 0., 2.);
  FillHist("LepType_Cv"+Label, LepType_Cv, weight, 13, -6., 7.);
  FillHist("h2D_NMissHit_IsCv"+Label, NMissHit, PassConvCut, weight, 5, 0., 5., 2, 0., 2.);

}


void ControlPlots::CheckConvCR(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                               vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                               vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  int NLepT=NMuT+NElT;
  if( FakeRun      and weight==0.  ) return; 
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if( NLepT!=3 ) return;
  bool PassTrigAcc=false;
  if( NMuT>1 && MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 ) PassTrigAcc=true;
  if( NElT>1 && ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15 ) PassTrigAcc=true;
  if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 ) PassTrigAcc=true;
  if( NElT>0 && NMuT>0 && ElTColl.at(0).Pt()>15 && MuTColl.at(0).Pt()>25 ) PassTrigAcc=true;
  if(!PassTrigAcc) return;
  int Qtot_Mu = SumCharge(MuTColl), Qtot_El = SumCharge(ElTColl), Qtot_Lep = Qtot_Mu+Qtot_El;
  if( abs(Qtot_Lep)!=1 ) return;
  if( NMuT==2 && Qtot_Mu!=0 ) return;
  if( NElT==2 && Qtot_El!=0 ) return;
  if( BJetColl.size()>0 ) return;

  //int NJet = JetColl.size();
  float PTl1_Z=0., PTl2_Z=0., PTl_Cv=0., Etal1_Z=0., Etal2_Z=0., Etal_Cv=0.;
  float MTW=0., MOSSF1=0., MOSSF2=0., M3l=0.;
  float dRll_Z=0., dRll_Z1Cv=0., dRll_Z2Cv=0.;
  if(NMuT==2 && NElT==1){
    Label   = "_2M1E"+Label;
    PTl1_Z  = MuTColl.at(0).Pt() , PTl2_Z  = MuTColl.at(1).Pt();
    Etal1_Z = MuTColl.at(0).Eta(), Etal2_Z = MuTColl.at(1).Eta();
    PTl_Cv  = ElTColl.at(0).Pt() , Etal_Cv = ElTColl.at(0).Eta();
    MOSSF1  = (MuTColl.at(0)+MuTColl.at(1)).M(), MTW = MT(ElTColl.at(0), vMET);
    M3l     = (MuTColl.at(0)+MuTColl.at(1)+ElTColl.at(0)).M();
    dRll_Z    = MuTColl.at(0).DeltaR(MuTColl.at(1));
    dRll_Z1Cv = MuTColl.at(0).DeltaR(ElTColl.at(0));
    dRll_Z2Cv = MuTColl.at(1).DeltaR(ElTColl.at(0));
  }
  else if(NMuT==1 && NElT==2){
    Label   = "_1M2E"+Label;
    PTl1_Z  = ElTColl.at(0).Pt() , PTl2_Z  = ElTColl.at(1).Pt();
    Etal1_Z = ElTColl.at(0).Eta(), Etal2_Z = ElTColl.at(1).Eta();
    PTl_Cv  = MuTColl.at(0).Pt() , Etal_Cv = MuTColl.at(0).Eta();
    MOSSF1  = (ElTColl.at(0)+ElTColl.at(1)).M(), MTW = MT(MuTColl.at(0), vMET);
    M3l     = (ElTColl.at(0)+ElTColl.at(1)+MuTColl.at(0)).M();
    dRll_Z    = ElTColl.at(0).DeltaR(ElTColl.at(1));
    dRll_Z1Cv = ElTColl.at(0).DeltaR(MuTColl.at(0));
    dRll_Z2Cv = ElTColl.at(1).DeltaR(MuTColl.at(0));
  }
  else if(NMuT==3 && NElT==0){
    Label   = "_3M0E"+Label;
    int IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    MOSSF1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    MOSSF2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    int IdxL1_Z=-1, IdxL2_Z=-1, IdxL_Cv=-1;
    if(fabs(MOSSF1-91.2)<fabs(MOSSF2-91.2)){
      IdxL1_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS1).Pt()? IdxOS:IdxSS1;
      IdxL2_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS1).Pt()? IdxSS1:IdxOS;
      IdxL_Cv = IdxSS2;
    }
    else{
      IdxL1_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS2).Pt()? IdxOS:IdxSS2;
      IdxL2_Z = MuTColl.at(IdxOS).Pt()>MuTColl.at(IdxSS2).Pt()? IdxSS2:IdxOS;
      IdxL_Cv = IdxSS1;
    }
    PTl1_Z  = MuTColl.at(IdxL1_Z).Pt() , PTl2_Z = MuTColl.at(IdxL2_Z).Pt();
    Etal1_Z = MuTColl.at(IdxL1_Z).Eta(), Etal2_Z = MuTColl.at(IdxL2_Z).Eta();
    PTl_Cv  = MuTColl.at(IdxL_Cv).Pt() , Etal_Cv = MuTColl.at(IdxL_Cv).Eta();
    MTW     = MT(MuTColl.at(IdxL_Cv), vMET);
    M3l     = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M();
    dRll_Z    = MuTColl.at(IdxL1_Z).DeltaR(MuTColl.at(IdxL2_Z));
    dRll_Z1Cv = MuTColl.at(IdxL1_Z).DeltaR(MuTColl.at(IdxL_Cv));
    dRll_Z2Cv = MuTColl.at(IdxL2_Z).DeltaR(MuTColl.at(IdxL_Cv));
  }
  else if(NMuT==0 && NElT==3){
    Label   = "_0M3E"+Label;
    int IdxOS  = TriElChargeIndex(ElTColl, "OS");
    int IdxSS1 = TriElChargeIndex(ElTColl, "SS1");
    int IdxSS2 = TriElChargeIndex(ElTColl, "SS2");
    MOSSF1 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS1)).M();
    MOSSF2 = (ElTColl.at(IdxOS)+ElTColl.at(IdxSS2)).M();
    int IdxL1_Z=-1, IdxL2_Z=-1, IdxL_Cv=-1;
    if(fabs(MOSSF1-91.2)<fabs(MOSSF2-91.2)){
      IdxL1_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS1).Pt()? IdxOS:IdxSS1;
      IdxL2_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS1).Pt()? IdxSS1:IdxOS;
      IdxL_Cv = IdxSS2;
    }
    else{
      IdxL1_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS2).Pt()? IdxOS:IdxSS2;
      IdxL2_Z = ElTColl.at(IdxOS).Pt()>ElTColl.at(IdxSS2).Pt()? IdxSS2:IdxOS;
      IdxL_Cv = IdxSS1;
    }
    PTl1_Z  = ElTColl.at(IdxL1_Z).Pt() , PTl2_Z = ElTColl.at(IdxL2_Z).Pt();
    Etal1_Z = ElTColl.at(IdxL1_Z).Eta(), Etal2_Z = ElTColl.at(IdxL2_Z).Eta();
    PTl_Cv  = ElTColl.at(IdxL_Cv).Pt() , Etal_Cv = ElTColl.at(IdxL_Cv).Eta();
    MTW     = MT(ElTColl.at(IdxL_Cv), vMET);
    M3l     = (ElTColl.at(0)+ElTColl.at(1)+ElTColl.at(2)).M();
    dRll_Z    = ElTColl.at(IdxL1_Z).DeltaR(ElTColl.at(IdxL2_Z));
    dRll_Z1Cv = ElTColl.at(IdxL1_Z).DeltaR(ElTColl.at(IdxL_Cv));
    dRll_Z2Cv = ElTColl.at(IdxL2_Z).DeltaR(ElTColl.at(IdxL_Cv));
  }

  if( MOSSF1<12 or ((NMuT>2 or NElT>2) and MOSSF2<12) ) return;
  if( fabs(MOSSF1-91.2)<10 or fabs(MOSSF2-91.2)<10    ) return;
  if( fabs(M3l-91.2)>10. ) return;

  Label = "_ZGSel"+Label;
  if( ConvRun && (!IsDATA) ){
    vector<Gen> TruthColl = GetGens();
    int NFk=0, NCv=0, NPr=0, NElse=0, NPt20=0;
    for(unsigned int im=0; im<MuTColl.size(); im++){
      int LepType = GetLeptonType_JH(MuTColl.at(im), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0) NPr++;
      else NElse++;
      if(MuTColl.at(im).Pt()>20) NPt20++;
    }
    for(unsigned int ie=0; ie<ElTColl.size(); ie++){
      int LepType = GetLeptonType_JH(ElTColl.at(ie), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0) NPr++;
      else NElse++;
      if(ElTColl.at(ie).Pt()>20) NPt20++;
    }
    FillHist("NPr"+Label, NPr, weight, 4, 0., 4.);
    FillHist("NFk"+Label, NFk, weight, 4, 0., 4.);
    FillHist("NCv"+Label, NCv, weight, 4, 0., 4.);
    FillHist("NElse"+Label, NElse, weight, 4, 0., 4.);
    if(NFk>0) return;
    if(MCSample.Contains("DY") && NPt20>2) return;
    else if(MCSample.Contains("ZG") && NPt20!=3) return;
  }
  int NPt20=0;
  for(unsigned int im=0; im<MuTColl.size(); im++){
    if(MuTColl.at(im).Pt()>20) NPt20++;
  }
  for(unsigned int ie=0; ie<ElTColl.size(); ie++){
    if(ElTColl.at(ie).Pt()>20) NPt20++;
  }

  FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  FillHist("M3l"+Label, M3l, weight, 30, 60., 120.);
  FillHist("MOSSF1"+Label, MOSSF1, weight, 15, 0., 150.);
  FillHist("MOSSF2"+Label, MOSSF2, weight, 15, 0., 150.);
  FillHist("PTl1_Z"+Label, PTl1_Z, weight, 24, 0., 120.);
  FillHist("PTl2_Z"+Label, PTl2_Z, weight, 16, 0., 80.);
  FillHist("PTl_Cv"+Label, PTl_Cv, weight, 16, 0., 80.);
  FillHist("Etal1_Z"+Label, Etal1_Z, weight, 20, -5., 5.);
  FillHist("Etal2_Z"+Label, Etal2_Z, weight, 20, -5., 5.);
  FillHist("Etal_Cv"+Label, Etal_Cv, weight, 20, -5., 5.);
  FillHist("Nj"+Label, JetColl.size(), weight, 10, 0., 10.);
  FillHist("Nb"+Label, BJetColl.size(), weight, 5, 0., 5.);
  FillHist("MET"+Label, vMET.Pt(), weight, 20, 0., 200.);
  FillHist("MTW"+Label, MTW, weight, 20, 0., 200.); 
  FillHist("dRll_Z"+Label, dRll_Z, weight, 25, 0., 5.);
  FillHist("dRll_Z1Cv"+Label, dRll_Z1Cv, weight, 25, 0., 5.);
  FillHist("dRll_Z2Cv"+Label, dRll_Z2Cv, weight, 25, 0., 5.);

  if(NPt20==3){
    Label.ReplaceAll("ZGSel", "ZGSel20");
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
    FillHist("M3l"+Label, M3l, weight, 30, 60., 120.);
    FillHist("MOSSF1"+Label, MOSSF1, weight, 15, 0., 150.);
    FillHist("MOSSF2"+Label, MOSSF2, weight, 15, 0., 150.);
    FillHist("PTl1_Z"+Label, PTl1_Z, weight, 24, 0., 120.);
    FillHist("PTl2_Z"+Label, PTl2_Z, weight, 16, 0., 80.);
    FillHist("PTl_Cv"+Label, PTl_Cv, weight, 16, 0., 80.);
    FillHist("Etal1_Z"+Label, Etal1_Z, weight, 20, -5., 5.);
    FillHist("Etal2_Z"+Label, Etal2_Z, weight, 20, -5., 5.);
    FillHist("Etal_Cv"+Label, Etal_Cv, weight, 20, -5., 5.);
  }
  else{
    Label.ReplaceAll("ZGSel", "ZGSelN20");
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
    FillHist("M3l"+Label, M3l, weight, 30, 60., 120.);
    FillHist("MOSSF1"+Label, MOSSF1, weight, 15, 0., 150.);
    FillHist("MOSSF2"+Label, MOSSF2, weight, 15, 0., 150.);
    FillHist("PTl1_Z"+Label, PTl1_Z, weight, 24, 0., 120.);
    FillHist("PTl2_Z"+Label, PTl2_Z, weight, 16, 0., 80.);
    FillHist("PTl_Cv"+Label, PTl_Cv, weight, 16, 0., 80.);
    FillHist("Etal1_Z"+Label, Etal1_Z, weight, 20, -5., 5.);
    FillHist("Etal2_Z"+Label, Etal2_Z, weight, 20, -5., 5.);
    FillHist("Etal_Cv"+Label, Etal_Cv, weight, 20, -5., 5.);
  }

}


void ControlPlots::PlotTop2LCR_OS2L(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                    vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                    vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) or (NElT==1 and NMuT==1) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NMuT==2){
    if(MuTColl.at(0).Charge()==MuTColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    float Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mll<12 or fabs(Mll-91.2)<15) return; 
    if(BJetColl.size()<1) return;
    if(JetColl.size() <3) return;
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      if(GenLepInfo>999) return;
    }
    InitializeTreeVars();
    SetVarSS2L(MuTColl, MuLColl, MuVColl, ElTColl, ElLColl, ElVColl, JetColl, BJetColl, vMET, weight, "");
    PlotParameters(Label+"_2M");
  }
  if(NElT==2){
    if(ElTColl.at(0).Charge()==ElTColl.at(1).Charge()) return;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    float Mll = (ElTColl.at(0)+ElTColl.at(1)).M();
    if(Mll<12 or fabs(Mll-91.2)<15) return; 
    if(BJetColl.size()<1) return;
    if(JetColl.size() <3) return;
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      if(GenLepInfo>999) return;
    }
    InitializeTreeVars();
    SetVarSS2L(MuTColl, MuLColl, MuVColl, ElTColl, ElLColl, ElVColl, JetColl, BJetColl, vMET, weight, "");
    PlotParameters(Label+"_2E");
  }
  if(NElT==1 && NMuT==1){
    if(ElTColl.at(0).Charge()==MuTColl.at(0).Charge()) return;
    if(!(MuTColl.at(0).Pt()>10 && ElTColl.at(0).Pt()>15)) return; 
    if(!(MuTColl.at(0).Pt()>25 || ElTColl.at(0).Pt()>25)) return; 
    if(BJetColl.size()<1) return;
    if(JetColl.size() <3) return;
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      if(GenLepInfo>999) return;
    }
    Muon TmpMu1 = MuTColl.at(0); Muon TmpMu2 = MuTColl.at(0);
    TmpMu2.SetPtEtaPhiE(ElTColl.at(0).Pt(), ElTColl.at(0).Eta(), ElTColl.at(0).Phi(), ElTColl.at(0).E());
    vector<Muon> TmpMuColl = {TmpMu1, TmpMu2}; vector<Electron> NullColl;
    std::sort(TmpMuColl.begin(), TmpMuColl.end(), PtComparing);
    InitializeTreeVars();
    SetVarSS2L(TmpMuColl, TmpMuColl, TmpMuColl, NullColl, NullColl, NullColl, JetColl, BJetColl, vMET, weight, "");
    PlotParameters(Label+"_EM");

    disc_BDTG = MVAreader_El->EvaluateMVA("BDTG method");
    FillHist("disc_BDTG2"+Label+"_EM", disc_BDTG, w_tot, 40, -1., 1.);
  }

}


void ControlPlots::CheckTopBSrc(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NMuT==2){
    Label+="_2M";
    if(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    float Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(DataYear>2016 && Mll<4) return; 
    if(BJetColl.size()<1) return;
    if(JetColl.size() <3) return;
    vector<Gen> TruthColl = GetGens();
    
    int NFkCv=0, TypePtFake=-1, NPr=0, NInConv=0, NExConv=0, NFk=0, NElse=0;
    for(unsigned int it_m=0; it_m<MuTColl.size(); it_m++){
      int LepType = GetLeptonType(MuTColl.at(it_m), TruthColl);
      if(LepType>0 && LepType<4) NPr++;
      else if(LepType==4 or LepType==5)  NInConv++;
      else if(LepType==-5 or LepType==-6) NExConv++;
      else if(LepType<0 && LepType>-5) NFk++;
      else NElse++;
      if(LepType>0) continue;
      FillHist("LepType"+Label, LepType, weight, 6, -6., 0.);
      NFkCv++;
      if     (it_m==0     ) TypePtFake=1;
      else if(TypePtFake<0) TypePtFake=2;
      else                  TypePtFake=3;

      if(LepType==-2){
        int MatchedIdx  = GenMatchedIdx(MuTColl.at(it_m),TruthColl);
        int MotherIdx   = FirstNonSelfMotherIdx(MatchedIdx,TruthColl);
        int GrMotherIdx = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
        int fMPID       = MotherIdx<0?   -1:abs(TruthColl.at(MotherIdx).PID()); 
        int fGrMPID     = GrMotherIdx<0? -1:abs(TruthColl.at(GrMotherIdx).PID());
        int MPIDType=0, GrMPIDType=0;
        if     (fMPID%1000>550 && fMPID%1000<560 && (fMPID/1000)%10!=5) MPIDType=1;
        else if(fMPID%1000>500 && fMPID%1000<550 && (fMPID/1000)%10!=5) MPIDType=2;
        else if(fMPID%1000>440 && fMPID%1000<450 && (fMPID/1000)%10!=5 && (fMPID/1000)%10!=4) MPIDType=3;
        else if(fMPID%1000>400 && fMPID%1000<440 && (fMPID/1000)%10!=5 && (fMPID/1000)%10!=4) MPIDType=4;
        else if((fMPID/1000)%10==5) MPIDType=5;
        else if((fMPID/1000)%10==4) MPIDType=6;
        else MPIDType=7;

        if     (fGrMPID%1000>550 && fGrMPID%1000<560 && (fGrMPID/1000)%10!=5) GrMPIDType=1;
        else if(fGrMPID%1000>500 && fGrMPID%1000<550 && (fGrMPID/1000)%10!=5) GrMPIDType=2;
        else if(fGrMPID%1000>440 && fGrMPID%1000<450 && (fGrMPID/1000)%10!=5 && (fGrMPID/1000)%10!=4) GrMPIDType=3;
        else if(fGrMPID%1000>400 && fGrMPID%1000<440 && (fGrMPID/1000)%10!=5 && (fGrMPID/1000)%10!=4) GrMPIDType=4;
        else if((fGrMPID/1000)%10==5) GrMPIDType=5;
        else if((fGrMPID/1000)%10==4) GrMPIDType=6;
        else GrMPIDType=7;

        FillHist("MPIDType"+Label, MPIDType, weight, 10, 0., 10.);
        FillHist("GrMPIDType"+Label, GrMPIDType, weight, 10, 0., 10.);
      }
    }
    FillHist("NFkCv"+Label, NFkCv, weight, 3, 0., 3.);
    FillHist("TypePtFake"+Label, TypePtFake, weight, 5, -1., 4.);
    if     (NPr==2                        ) FillHist("LepTypePair"+Label, 0., weight, 10, 0., 10.);
    else if(NPr==1 && NFk==1              ) FillHist("LepTypePair"+Label, 1., weight, 10, 0., 10.);
    else if(NPr==1 && NInConv==1          ) FillHist("LepTypePair"+Label, 2., weight, 10, 0., 10.);
    else if(NPr==1 && NExConv==1          ) FillHist("LepTypePair"+Label, 3., weight, 10, 0., 10.);
    else if(NPr==0 && NFk==2              ) FillHist("LepTypePair"+Label, 4., weight, 10, 0., 10.);
    else if(NPr==0 && NFk==1 && NInConv==1) FillHist("LepTypePair"+Label, 5., weight, 10, 0., 10.);
    else if(NPr==0 && NFk==1 && NExConv==1) FillHist("LepTypePair"+Label, 6., weight, 10, 0., 10.);
    else if(NPr==0 && NInConv==2          ) FillHist("LepTypePair"+Label, 7., weight, 10, 0., 10.);
    else if(NPr==0 && NExConv==2          ) FillHist("LepTypePair"+Label, 8., weight, 10, 0., 10.);

    int NB_b=0, NB_c=0, NB_l=0;
    for(unsigned int it_b=0; it_b<BJetColl.size(); it_b++){
      bool NearB=false, NearCNotB=false;
      for(unsigned int it_t=2; it_t<TruthColl.size(); it_t++){
        int apid  = abs(TruthColl.at(it_t).PID());
        //int ampid = abs(TruthColl.at(TruthColl.at(it_t).MotherIndex()).PID());
        if(!(apid==5 or apid==4)) continue;
        if(BJetColl.at(it_b).DeltaR(TruthColl.at(it_t))>0.5) continue;
        if(TruthColl.at(it_t).Pt()/BJetColl.at(it_b).Pt()<0.5) continue;
        if     (apid==5){ NearB=true; break; }
        else if(apid==4){ NearCNotB=true; break; }
      }
      if     (NearB)     NB_b++; 
      else if(NearCNotB) NB_c++;
      else               NB_l++;
    }
    FillHist("NB_b"+Label, NB_b, weight, 5, 0., 5.);
    FillHist("NB_c"+Label, NB_c, weight, 5, 0., 5.);
    FillHist("NB_l"+Label, NB_l, weight, 5, 0., 5.);
    if(BJetColl.size()==1 && NB_b==0 && NB_c==0) FillHist("NB_Comp"+Label, 0., weight, 19, 0., 19.);
    if(BJetColl.size()==1 && NB_b==0 && NB_c==1) FillHist("NB_Comp"+Label, 1., weight, 19, 0., 19.);
    if(BJetColl.size()==1 && NB_b==1 && NB_c==0) FillHist("NB_Comp"+Label, 2., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==0 && NB_c==0) FillHist("NB_Comp"+Label, 3., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==0 && NB_c==1) FillHist("NB_Comp"+Label, 4., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==0 && NB_c==2) FillHist("NB_Comp"+Label, 5., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==1 && NB_c==0) FillHist("NB_Comp"+Label, 6., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==1 && NB_c==1) FillHist("NB_Comp"+Label, 7., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==2 && NB_c==0) FillHist("NB_Comp"+Label, 8., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==0) FillHist("NB_Comp"+Label, 9., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==1) FillHist("NB_Comp"+Label, 10., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==2) FillHist("NB_Comp"+Label, 11., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==3) FillHist("NB_Comp"+Label, 12., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==1 && NB_c==0) FillHist("NB_Comp"+Label, 13., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==1 && NB_c==1) FillHist("NB_Comp"+Label, 14., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==1 && NB_c==2) FillHist("NB_Comp"+Label, 15., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==2 && NB_c==0) FillHist("NB_Comp"+Label, 16., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==2 && NB_c==1) FillHist("NB_Comp"+Label, 17., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==3 && NB_c==0) FillHist("NB_Comp"+Label, 18., weight, 19, 0., 19.);
  }
  if(NElT==2){
    Label+="_2E";
    if(ElTColl.at(0).Charge()!=ElTColl.at(1).Charge()) return;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    if(BJetColl.size()<1) return;
    if(JetColl.size() <3) return;
    vector<Gen> TruthColl = GetGens();
    
    int NFkCv=0, TypePtFake=-1, NPr=0, NInConv=0, NExConv=0, NFk=0, NElse=0;
    for(unsigned int it_e=0; it_e<ElTColl.size(); it_e++){
      int LepType = GetLeptonType(ElTColl.at(it_e), TruthColl);
      if(LepType>0 && LepType<4) NPr++;
      else if(LepType==4 or LepType==5)  NInConv++;
      else if(LepType==-5 or LepType==-6) NExConv++;
      else if(LepType<0 && LepType>-5) NFk++;
      else NElse++;
      if(LepType>0) continue;
      FillHist("LepType"+Label, LepType, weight, 6, -6., 0.);
      NFkCv++;
      if     (it_e==0     ) TypePtFake=1;
      else if(TypePtFake<0) TypePtFake=2;
      else                  TypePtFake=3;

      if(LepType==-2){
        int MatchedIdx  = GenMatchedIdx(ElTColl.at(it_e),TruthColl);
        int MotherIdx   = FirstNonSelfMotherIdx(MatchedIdx,TruthColl);
        int GrMotherIdx = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
        int fMPID       = MotherIdx<0?   -1:abs(TruthColl.at(MotherIdx).PID()); 
        int fGrMPID     = GrMotherIdx<0? -1:abs(TruthColl.at(GrMotherIdx).PID());
        int MPIDType=0, GrMPIDType=0;
        if     (fMPID%1000>550 && fMPID%1000<560 && (fMPID/1000)%10!=5) MPIDType=1;
        else if(fMPID%1000>500 && fMPID%1000<550 && (fMPID/1000)%10!=5) MPIDType=2;
        else if(fMPID%1000>440 && fMPID%1000<450 && (fMPID/1000)%10!=5 && (fMPID/1000)%10!=4) MPIDType=3;
        else if(fMPID%1000>400 && fMPID%1000<440 && (fMPID/1000)%10!=5 && (fMPID/1000)%10!=4) MPIDType=4;
        else if((fMPID/1000)%10==5) MPIDType=5;
        else if((fMPID/1000)%10==4) MPIDType=6;
        else MPIDType=7;

        if     (fGrMPID%1000>550 && fGrMPID%1000<560 && (fGrMPID/1000)%10!=5) GrMPIDType=1;
        else if(fGrMPID%1000>500 && fGrMPID%1000<550 && (fGrMPID/1000)%10!=5) GrMPIDType=2;
        else if(fGrMPID%1000>440 && fGrMPID%1000<450 && (fGrMPID/1000)%10!=5 && (fGrMPID/1000)%10!=4) GrMPIDType=3;
        else if(fGrMPID%1000>400 && fGrMPID%1000<440 && (fGrMPID/1000)%10!=5 && (fGrMPID/1000)%10!=4) GrMPIDType=4;
        else if((fGrMPID/1000)%10==5) GrMPIDType=5;
        else if((fGrMPID/1000)%10==4) GrMPIDType=6;
        else GrMPIDType=7;

        FillHist("MPIDType"+Label, MPIDType, weight, 10, 0., 10.);
        FillHist("GrMPIDType"+Label, GrMPIDType, weight, 10, 0., 10.);
      }
    }
    FillHist("NFkCv"+Label, NFkCv, weight, 3, 0., 3.);
    FillHist("TypePtFake"+Label, TypePtFake, weight, 5, -1., 4.);
    if     (NPr==2                        ) FillHist("LepTypePair"+Label, 0., weight, 10, 0., 10.);
    else if(NPr==1 && NFk==1              ) FillHist("LepTypePair"+Label, 1., weight, 10, 0., 10.);
    else if(NPr==1 && NInConv==1          ) FillHist("LepTypePair"+Label, 2., weight, 10, 0., 10.);
    else if(NPr==1 && NExConv==1          ) FillHist("LepTypePair"+Label, 3., weight, 10, 0., 10.);
    else if(NPr==0 && NFk==2              ) FillHist("LepTypePair"+Label, 4., weight, 10, 0., 10.);
    else if(NPr==0 && NFk==1 && NInConv==1) FillHist("LepTypePair"+Label, 5., weight, 10, 0., 10.);
    else if(NPr==0 && NFk==1 && NExConv==1) FillHist("LepTypePair"+Label, 6., weight, 10, 0., 10.);
    else if(NPr==0 && NInConv==2          ) FillHist("LepTypePair"+Label, 7., weight, 10, 0., 10.);
    else if(NPr==0 && NExConv==2          ) FillHist("LepTypePair"+Label, 8., weight, 10, 0., 10.);


    int NB_b=0, NB_c=0, NB_l=0;
    for(unsigned int it_b=0; it_b<BJetColl.size(); it_b++){
      bool NearB=false, NearCNotB=false;
      for(unsigned int it_t=2; it_t<TruthColl.size(); it_t++){
        int apid  = abs(TruthColl.at(it_t).PID());
        //int ampid = abs(TruthColl.at(TruthColl.at(it_t).MotherIndex()).PID());
        if(!(apid==5 or apid==4)) continue;
        if(BJetColl.at(it_b).DeltaR(TruthColl.at(it_t))>0.5) continue;
        if(TruthColl.at(it_t).Pt()/BJetColl.at(it_b).Pt()<0.5) continue;
        if     (apid==5){ NearB=true; break; }
        else if(apid==4){ NearCNotB=true; break; }
      }
      if     (NearB)     NB_b++; 
      else if(NearCNotB) NB_c++;
      else               NB_l++;
    }
    FillHist("NB_b"+Label, NB_b, weight, 5, 0., 5.);
    FillHist("NB_c"+Label, NB_c, weight, 5, 0., 5.);
    FillHist("NB_l"+Label, NB_l, weight, 5, 0., 5.);
    if(BJetColl.size()==1 && NB_b==0 && NB_c==0) FillHist("NB_Comp"+Label, 0., weight, 19, 0., 19.);
    if(BJetColl.size()==1 && NB_b==0 && NB_c==1) FillHist("NB_Comp"+Label, 1., weight, 19, 0., 19.);
    if(BJetColl.size()==1 && NB_b==1 && NB_c==0) FillHist("NB_Comp"+Label, 2., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==0 && NB_c==0) FillHist("NB_Comp"+Label, 3., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==0 && NB_c==1) FillHist("NB_Comp"+Label, 4., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==0 && NB_c==2) FillHist("NB_Comp"+Label, 5., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==1 && NB_c==0) FillHist("NB_Comp"+Label, 6., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==1 && NB_c==1) FillHist("NB_Comp"+Label, 7., weight, 19, 0., 19.);
    if(BJetColl.size()==2 && NB_b==2 && NB_c==0) FillHist("NB_Comp"+Label, 8., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==0) FillHist("NB_Comp"+Label, 9., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==1) FillHist("NB_Comp"+Label, 10., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==2) FillHist("NB_Comp"+Label, 11., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==0 && NB_c==3) FillHist("NB_Comp"+Label, 12., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==1 && NB_c==0) FillHist("NB_Comp"+Label, 13., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==1 && NB_c==1) FillHist("NB_Comp"+Label, 14., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==1 && NB_c==2) FillHist("NB_Comp"+Label, 15., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==2 && NB_c==0) FillHist("NB_Comp"+Label, 16., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==2 && NB_c==1) FillHist("NB_Comp"+Label, 17., weight, 19, 0., 19.);
    if(BJetColl.size()==3 && NB_b==3 && NB_c==0) FillHist("NB_Comp"+Label, 18., weight, 19, 0., 19.);
  }

}


void ControlPlots::SetVarSS2L(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                              vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                              vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(BJetColl.size()<1) return;
  if(JetColl.size() <3) return;
  if(NMuT==2){
    InitializeTreeVars();
    Nj      = JetColl.size();
    Nb      = BJetColl.size();
    Ptl1    = MuTColl.at(0).Pt();
    Ptl2    = MuTColl.at(1).Pt();
    Ptj1    = JetColl.at(0).Pt();
    Ptj2    = JetColl.at(1).Pt();
    Ptj3    = JetColl.at(2).Pt();
    Ptb1    = BJetColl.at(0).Pt();
    Ptb2    = BJetColl.size()<2? -1.:BJetColl.at(1).Pt();
    MET     = vMET.Pt();
    dEtall  = abs(MuTColl.at(0).Eta()-MuTColl.at(1).Eta());
    dRll    = MuTColl.at(0).DeltaR(MuTColl.at(1));
    dRjj12  = JetColl.at(0).DeltaR(JetColl.at(1));
    dRjj23  = JetColl.at(1).DeltaR(JetColl.at(2));
    dRjj13  = JetColl.at(0).DeltaR(JetColl.at(2));
    dRlj11  = MuTColl.at(0).DeltaR(JetColl.at(0));
    dRlj12  = MuTColl.at(0).DeltaR(JetColl.at(1));
    dRlj13  = MuTColl.at(0).DeltaR(JetColl.at(2));
    dRlj21  = MuTColl.at(1).DeltaR(JetColl.at(0));
    dRlj22  = MuTColl.at(1).DeltaR(JetColl.at(1));
    dRlj23  = MuTColl.at(1).DeltaR(JetColl.at(2));
    dRlb11  = MuTColl.at(0).DeltaR(BJetColl.at(0));
    dRlb21  = MuTColl.at(1).DeltaR(BJetColl.at(0));
    MSSSF   = (MuTColl.at(0)+MuTColl.at(1)).M();
    Mbl11   = (MuTColl.at(0)+BJetColl.at(0)).M();
    Mbl12   = (MuTColl.at(1)+BJetColl.at(0)).M();
    Mbl21   = BJetColl.size()<2? -1.:(MuTColl.at(0)+BJetColl.at(1)).M();
    Mbl22   = BJetColl.size()<2? -1.:(MuTColl.at(1)+BJetColl.at(1)).M();
    Mlj11   = (MuTColl.at(0)+JetColl.at(0)).M();
    Mlj12   = (MuTColl.at(0)+JetColl.at(1)).M();
    Mlj13   = (MuTColl.at(0)+JetColl.at(2)).M();
    Mlj21   = (MuTColl.at(1)+JetColl.at(0)).M();
    Mlj22   = (MuTColl.at(1)+JetColl.at(1)).M();
    Mlj23   = (MuTColl.at(1)+JetColl.at(2)).M();
    MTvl1   = MT(MuTColl.at(0),vMET);
    MTvl2   = MT(MuTColl.at(1),vMET);
    Mllj1   = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(0)).M();
    Mllj2   = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(1)).M();
    Mllj3   = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(2)).M();
    Mllj4   = JetColl.size()<4? -1.:(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(3)).M();
    Mllb1   = (MuTColl.at(0)+MuTColl.at(1)+BJetColl.at(0)).M();
    Mllb2   = BJetColl.size()<2? -1.:(MuTColl.at(0)+MuTColl.at(1)+BJetColl.at(1)).M();
    Mlljj12 = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(0)+JetColl.at(1)).M();
    Mlljj13 = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(0)+JetColl.at(2)).M();
    Mlljj14 = JetColl.size()<4? -1.:(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(0)+JetColl.at(3)).M();
    Mlljj23 = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(1)+JetColl.at(2)).M();
    Mlljj24 = JetColl.size()<4? -1.:(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(1)+JetColl.at(3)).M();
    Mlljj34 = JetColl.size()<4? -1.:(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(2)+JetColl.at(3)).M();
    Mljj112 = (MuTColl.at(0)+JetColl.at(0)+JetColl.at(1)).M();
    Mljj113 = (MuTColl.at(0)+JetColl.at(0)+JetColl.at(2)).M();
    Mljj114 = JetColl.size()<4? -1.:(MuTColl.at(0)+JetColl.at(0)+JetColl.at(3)).M();
    Mljj123 = (MuTColl.at(0)+JetColl.at(1)+JetColl.at(2)).M();
    Mljj124 = JetColl.size()<4? -1.:(MuTColl.at(0)+JetColl.at(1)+JetColl.at(3)).M();
    Mljj134 = JetColl.size()<4? -1.:(MuTColl.at(0)+JetColl.at(2)+JetColl.at(3)).M();
    Mljj212 = (MuTColl.at(1)+JetColl.at(0)+JetColl.at(1)).M();
    Mljj213 = (MuTColl.at(1)+JetColl.at(0)+JetColl.at(2)).M();
    Mljj214 = JetColl.size()<4? -1.:(MuTColl.at(1)+JetColl.at(0)+JetColl.at(3)).M();
    Mljj223 = (MuTColl.at(1)+JetColl.at(1)+JetColl.at(2)).M();
    Mljj224 = JetColl.size()<4? -1.:(MuTColl.at(1)+JetColl.at(1)+JetColl.at(3)).M();
    Mljj234 = JetColl.size()<4? -1.:(MuTColl.at(1)+JetColl.at(2)+JetColl.at(3)).M();
    Mjj12   = (JetColl.at(0)+JetColl.at(1)).M();
    Mjj13   = (JetColl.at(0)+JetColl.at(2)).M();
    Mjj14   = JetColl.size()<4? -1.:(JetColl.at(0)+JetColl.at(3)).M();
    Mjj23   = (JetColl.at(1)+JetColl.at(2)).M();
    Mjj24   = JetColl.size()<4? -1.:(JetColl.at(1)+JetColl.at(3)).M();
    Mjj34   = JetColl.size()<4? -1.:(JetColl.at(2)+JetColl.at(3)).M();

    //Vars requiring complex algo.
    HT      = 0;
      for(unsigned int itj=0; itj<JetColl.size(); itj++){ HT+=JetColl.at(itj).Pt(); }
    MET2HT  = pow(MET,2.)/HT;
    int Idxj1W_2jL=-1, Idxj2W_2jL=-1; float bestmlljj=-1;
    int Idxj1W_1jL=-1; float bestmllj=-1;
    int Idxj1W1_H=-1, Idxj2W1_H=-1, Idxj1W2_H=-1, Idxj2W2_H=-1; float bestmjj1=-1, bestmjj2=-1;
    for(unsigned int itj1=0; itj1<JetColl.size(); itj1++){
      if(bestmllj<0){ Idxj1W_1jL=itj1; bestmllj=(MuTColl.at(0)+MuTColl.at(1)).M(); }
      else{
        float tmpmljj = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(itj1)).M();
        if(fabs(tmpmljj-80.4)<fabs(bestmllj-80.4)){ bestmllj=tmpmljj; Idxj1W_1jL=itj1; }
      }
      for(unsigned int itj2=itj1+1; itj2<JetColl.size(); itj2++){
        if(bestmlljj<0){ Idxj1W_2jL=itj1; Idxj2W_2jL=itj2; bestmlljj=(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(itj1)+JetColl.at(itj2)).M(); }
        else{
          float tmpmlljj = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(itj1)+JetColl.at(itj2)).M();
          if(fabs(tmpmlljj-80.4)<fabs(bestmlljj-80.4)){ bestmlljj=tmpmlljj; Idxj1W_2jL=itj1, Idxj2W_2jL=itj2; }
        }
        if(bestmjj1<0){ Idxj1W1_H=itj1, Idxj2W1_H=itj2; bestmjj1=(JetColl.at(itj1)+JetColl.at(itj2)).M(); }
        else{
          float tmpmjj = (JetColl.at(itj1)+JetColl.at(itj2)).M();
          if(fabs(tmpmjj-80.4)<fabs(bestmjj1-80.4)){ bestmjj1=tmpmjj; Idxj1W1_H=itj1, Idxj2W1_H=itj2; }
        }
      }
    }
    for(unsigned int itj1=0; itj1<JetColl.size(); itj1++){
      for(unsigned int itj2=itj1+1; itj2<JetColl.size(); itj2++){
        if((int) itj1==Idxj1W1_H or (int) itj1==Idxj2W1_H or (int) itj2==Idxj1W1_H or (int) itj2==Idxj2W1_H) continue;
        if(bestmjj2<0){ Idxj1W2_H=itj1, Idxj2W2_H=itj2; bestmjj2=(JetColl.at(itj1)+JetColl.at(itj2)).M(); }
        else{
          float tmpmjj = (JetColl.at(itj1)+JetColl.at(itj2)).M();
          if(fabs(tmpmjj-80.4)<fabs(bestmjj2-80.4)){ bestmjj2=tmpmjj; Idxj1W2_H=itj1, Idxj2W2_H=itj2; }
        }
      }
    }
    MllW_2jL = bestmlljj;
    MllW_1jL = bestmllj;
    MllW1_H  = (MuTColl.at(0)+MuTColl.at(1)+JetColl.at(Idxj1W1_H)+JetColl.at(Idxj2W1_H)).M();
    MllW2_H  = bestmjj2<0? -1.:(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(Idxj1W2_H)+JetColl.at(Idxj2W2_H)).M();
    MjjW1    = bestmjj1;
    MjjW2    = bestmjj2;
    Ml1W_2jL = (MuTColl.at(0)+JetColl.at(Idxj1W_2jL)+JetColl.at(Idxj2W_2jL)).M();
    Ml1W_1jL = (MuTColl.at(0)+JetColl.at(Idxj1W_1jL)).M();
    Ml2W_2jL = (MuTColl.at(1)+JetColl.at(Idxj1W_2jL)+JetColl.at(Idxj2W_2jL)).M();
    Ml2W_1jL = (MuTColl.at(1)+JetColl.at(Idxj1W_1jL)).M();
    Ml1W1_H  = (MuTColl.at(0)+JetColl.at(Idxj1W1_H)+JetColl.at(Idxj2W1_H)).M();
    Ml1W2_H  = bestmjj2<0? -1.:(MuTColl.at(0)+JetColl.at(Idxj1W2_H)+JetColl.at(Idxj2W2_H)).M();
    Ml2W1_H  = (MuTColl.at(1)+JetColl.at(Idxj1W1_H)+JetColl.at(Idxj2W1_H)).M();
    Ml2W2_H  = bestmjj2<0? -1.:(MuTColl.at(1)+JetColl.at(Idxj1W2_H)+JetColl.at(Idxj2W2_H)).M();
    w_tot    = weight;
    disc_BDTG = MVAreader_Mu->EvaluateMVA("BDTG method");
  }
  if(NElT==2){
    InitializeTreeVars();
    Nj      = JetColl.size();
    Nb      = BJetColl.size();
    Ptl1    = ElTColl.at(0).Pt();
    Ptl2    = ElTColl.at(1).Pt();
    Ptj1    = JetColl.at(0).Pt();
    Ptj2    = JetColl.at(1).Pt();
    Ptj3    = JetColl.at(2).Pt();
    Ptb1    = BJetColl.at(0).Pt();
    Ptb2    = BJetColl.size()<2? -1.:BJetColl.at(1).Pt();
    MET     = vMET.Pt();
    dEtall  = abs(ElTColl.at(0).Eta()-ElTColl.at(1).Eta());
    dRll    = ElTColl.at(0).DeltaR(ElTColl.at(1));
    dRjj12  = JetColl.at(0).DeltaR(JetColl.at(1));
    dRjj23  = JetColl.at(1).DeltaR(JetColl.at(2));
    dRjj13  = JetColl.at(0).DeltaR(JetColl.at(2));
    dRlj11  = ElTColl.at(0).DeltaR(JetColl.at(0));
    dRlj12  = ElTColl.at(0).DeltaR(JetColl.at(1));
    dRlj13  = ElTColl.at(0).DeltaR(JetColl.at(2));
    dRlj21  = ElTColl.at(1).DeltaR(JetColl.at(0));
    dRlj22  = ElTColl.at(1).DeltaR(JetColl.at(1));
    dRlj23  = ElTColl.at(1).DeltaR(JetColl.at(2));
    dRlb11  = ElTColl.at(0).DeltaR(BJetColl.at(0));
    dRlb21  = ElTColl.at(1).DeltaR(BJetColl.at(0));
    MSSSF   = (ElTColl.at(0)+ElTColl.at(1)).M();
    Mbl11   = (ElTColl.at(0)+BJetColl.at(0)).M();
    Mbl12   = (ElTColl.at(1)+BJetColl.at(0)).M();
    Mbl21   = BJetColl.size()<2? -1.:(ElTColl.at(0)+BJetColl.at(1)).M();
    Mbl22   = BJetColl.size()<2? -1.:(ElTColl.at(1)+BJetColl.at(1)).M();
    Mlj11   = (ElTColl.at(0)+JetColl.at(0)).M();
    Mlj12   = (ElTColl.at(0)+JetColl.at(1)).M();
    Mlj13   = (ElTColl.at(0)+JetColl.at(2)).M();
    Mlj21   = (ElTColl.at(1)+JetColl.at(0)).M();
    Mlj22   = (ElTColl.at(1)+JetColl.at(1)).M();
    Mlj23   = (ElTColl.at(1)+JetColl.at(2)).M();
    MTvl1   = MT(ElTColl.at(0),vMET);
    MTvl2   = MT(ElTColl.at(1),vMET);
    Mllj1   = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(0)).M();
    Mllj2   = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(1)).M();
    Mllj3   = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(2)).M();
    Mllj4   = JetColl.size()<4? -1.:(ElTColl.at(0)+ElTColl.at(1)+JetColl.at(3)).M();
    Mllb1   = (ElTColl.at(0)+ElTColl.at(1)+BJetColl.at(0)).M();
    Mllb2   = BJetColl.size()<2? -1.:(ElTColl.at(0)+ElTColl.at(1)+BJetColl.at(1)).M();
    Mlljj12 = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(0)+JetColl.at(1)).M();
    Mlljj13 = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(0)+JetColl.at(2)).M();
    Mlljj14 = JetColl.size()<4? -1.:(ElTColl.at(0)+ElTColl.at(1)+JetColl.at(0)+JetColl.at(3)).M();
    Mlljj23 = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(1)+JetColl.at(2)).M();
    Mlljj24 = JetColl.size()<4? -1.:(ElTColl.at(0)+ElTColl.at(1)+JetColl.at(1)+JetColl.at(3)).M();
    Mlljj34 = JetColl.size()<4? -1.:(ElTColl.at(0)+ElTColl.at(1)+JetColl.at(2)+JetColl.at(3)).M();
    Mljj112 = (ElTColl.at(0)+JetColl.at(0)+JetColl.at(1)).M();
    Mljj113 = (ElTColl.at(0)+JetColl.at(0)+JetColl.at(2)).M();
    Mljj114 = JetColl.size()<4? -1.:(ElTColl.at(0)+JetColl.at(0)+JetColl.at(3)).M();
    Mljj123 = (ElTColl.at(0)+JetColl.at(1)+JetColl.at(2)).M();
    Mljj124 = JetColl.size()<4? -1.:(ElTColl.at(0)+JetColl.at(1)+JetColl.at(3)).M();
    Mljj134 = JetColl.size()<4? -1.:(ElTColl.at(0)+JetColl.at(2)+JetColl.at(3)).M();
    Mljj212 = (ElTColl.at(1)+JetColl.at(0)+JetColl.at(1)).M();
    Mljj213 = (ElTColl.at(1)+JetColl.at(0)+JetColl.at(2)).M();
    Mljj214 = JetColl.size()<4? -1.:(ElTColl.at(1)+JetColl.at(0)+JetColl.at(3)).M();
    Mljj223 = (ElTColl.at(1)+JetColl.at(1)+JetColl.at(2)).M();
    Mljj224 = JetColl.size()<4? -1.:(ElTColl.at(1)+JetColl.at(1)+JetColl.at(3)).M();
    Mljj234 = JetColl.size()<4? -1.:(ElTColl.at(1)+JetColl.at(2)+JetColl.at(3)).M();
    Mjj12   = (JetColl.at(0)+JetColl.at(1)).M();
    Mjj13   = (JetColl.at(0)+JetColl.at(2)).M();
    Mjj14   = JetColl.size()<4? -1.:(JetColl.at(0)+JetColl.at(3)).M();
    Mjj23   = (JetColl.at(1)+JetColl.at(2)).M();
    Mjj24   = JetColl.size()<4? -1.:(JetColl.at(1)+JetColl.at(3)).M();
    Mjj34   = JetColl.size()<4? -1.:(JetColl.at(2)+JetColl.at(3)).M();

    //Vars requiring complex algo.
    HT      = 0;
      for(unsigned int itj=0; itj<JetColl.size(); itj++){ HT+=JetColl.at(itj).Pt(); }
    MET2HT  = pow(MET,2.)/HT;
    int Idxj1W_2jL=-1, Idxj2W_2jL=-1; float bestmlljj=-1;
    int Idxj1W_1jL=-1; float bestmllj=-1;
    int Idxj1W1_H=-1, Idxj2W1_H=-1, Idxj1W2_H=-1, Idxj2W2_H=-1; float bestmjj1=-1, bestmjj2=-1;
    for(unsigned int itj1=0; itj1<JetColl.size(); itj1++){
      if(bestmllj<0){ Idxj1W_1jL=itj1; bestmllj=(ElTColl.at(0)+ElTColl.at(1)).M(); }
      else{
        float tmpmljj = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(itj1)).M();
        if(fabs(tmpmljj-80.4)<fabs(bestmllj-80.4)){ bestmllj=tmpmljj; Idxj1W_1jL=itj1; }
      }
      for(unsigned int itj2=itj1+1; itj2<JetColl.size(); itj2++){
        if(bestmlljj<0){ Idxj1W_2jL=itj1; Idxj2W_2jL=itj2; bestmlljj=(ElTColl.at(0)+ElTColl.at(1)+JetColl.at(itj1)+JetColl.at(itj2)).M(); }
        else{
          float tmpmlljj = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(itj1)+JetColl.at(itj2)).M();
          if(fabs(tmpmlljj-80.4)<fabs(bestmlljj-80.4)){ bestmlljj=tmpmlljj; Idxj1W_2jL=itj1, Idxj2W_2jL=itj2; }
        }
        if(bestmjj1<0){ Idxj1W1_H=itj1, Idxj2W1_H=itj2; bestmjj1=(JetColl.at(itj1)+JetColl.at(itj2)).M(); }
        else{
          float tmpmjj = (JetColl.at(itj1)+JetColl.at(itj2)).M();
          if(fabs(tmpmjj-80.4)<fabs(bestmjj1-80.4)){ bestmjj1=tmpmjj; Idxj1W1_H=itj1, Idxj2W1_H=itj2; }
        }
      }
    }
    for(unsigned int itj1=0; itj1<JetColl.size(); itj1++){
      for(unsigned int itj2=itj1+1; itj2<JetColl.size(); itj2++){
        if((int) itj1==Idxj1W1_H or (int) itj1==Idxj2W1_H or (int) itj2==Idxj1W1_H or (int) itj2==Idxj2W1_H) continue;
        if(bestmjj2<0){ Idxj1W2_H=itj1, Idxj2W2_H=itj2; bestmjj2=(JetColl.at(itj1)+JetColl.at(itj2)).M(); }
        else{
          float tmpmjj = (JetColl.at(itj1)+JetColl.at(itj2)).M();
          if(fabs(tmpmjj-80.4)<fabs(bestmjj2-80.4)){ bestmjj2=tmpmjj; Idxj1W2_H=itj1, Idxj2W2_H=itj2; }
        }
      }
    }
    MllW_2jL = bestmlljj;
    MllW_1jL = bestmllj;
    MllW1_H  = (ElTColl.at(0)+ElTColl.at(1)+JetColl.at(Idxj1W1_H)+JetColl.at(Idxj2W1_H)).M();
    MllW2_H  = bestmjj2<0? -1.:(ElTColl.at(0)+ElTColl.at(1)+JetColl.at(Idxj1W2_H)+JetColl.at(Idxj2W2_H)).M();
    MjjW1    = bestmjj1;
    MjjW2    = bestmjj2;
    Ml1W_2jL = (ElTColl.at(0)+JetColl.at(Idxj1W_2jL)+JetColl.at(Idxj2W_2jL)).M();
    Ml1W_1jL = (ElTColl.at(0)+JetColl.at(Idxj1W_1jL)).M();
    Ml2W_2jL = (ElTColl.at(1)+JetColl.at(Idxj1W_2jL)+JetColl.at(Idxj2W_2jL)).M();
    Ml2W_1jL = (ElTColl.at(1)+JetColl.at(Idxj1W_1jL)).M();
    Ml1W1_H  = (ElTColl.at(0)+JetColl.at(Idxj1W1_H)+JetColl.at(Idxj2W1_H)).M();
    Ml1W2_H  = bestmjj2<0? -1.:(ElTColl.at(0)+JetColl.at(Idxj1W2_H)+JetColl.at(Idxj2W2_H)).M();
    Ml2W1_H  = (ElTColl.at(1)+JetColl.at(Idxj1W1_H)+JetColl.at(Idxj2W1_H)).M();
    Ml2W2_H  = bestmjj2<0? -1.:(ElTColl.at(1)+JetColl.at(Idxj1W2_H)+JetColl.at(Idxj2W2_H)).M();
    w_tot    = weight;
    disc_BDTG = MVAreader_Mu->EvaluateMVA("BDTG method");
  }

}



void ControlPlots::MakePlotSS2L(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  bool CheckSB_M = Label.Contains("SB_M"), CheckSB_Iso = Label.Contains("SB_Iso");
  if(!CheckSB_Iso){
    if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
    if( !(NMuT==NMuL and NElT==NElL) ) return;
  }
  else{
    if( !( (NMuL==2 and NElL==0) or (NElL==2 and NMuL==0) ) ) return;
    if( (NMuT==NMuL and NElT==NElL) ) return;
  }
  if( !(NMuL==NMuV and NElL==NElV) ) return;
  if( FakeRun      and weight==0.  ) return; 
  if(NMuL==2){
    if(MuLColl.at(0).Charge()!=MuLColl.at(1).Charge()) return;
    if(!(MuLColl.at(0).Pt()>20 && MuLColl.at(1).Pt()>10)) return;
    float Mll = (MuLColl.at(0)+MuLColl.at(1)).M();
    if(DataYear>2016 && Mll<4) return; 
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      if(ConvRun && GenLepInfo>=100) return;
      if(FlipRun && (GenLepInfo>999 or GenLepInfo<100)) return;
    }

    if(BJetColl.size()<1) return;
      FillHist("NJPreCut"+Label+"_2M", JetColl.size(), weight, 10, 0., 10.);
    if(JetColl.size() <3) return;
    if(CheckSB_M && Mll<80) return;
    InitializeTreeVars();
    if(!CheckSB_Iso) SetVarSS2L(MuTColl, MuLColl, MuVColl, ElTColl, ElLColl, ElVColl, JetColl, BJetColl, vMET, weight, "");
    else             SetVarSS2L(MuLColl, MuLColl, MuVColl, ElLColl, ElLColl, ElVColl, JetColl, BJetColl, vMET, weight, "");
    PlotParametersSS(Label+"_2M");
  }
  if(NElL==2){
    int aSumQ = abs(SumCharge(ElLColl));
    if(IsDATA && FlipRun){ if(aSumQ!=0) return; }
    else                 { if(aSumQ==0) return; }
    //if(ElLColl.at(0).Charge()!=ElLColl.at(1).Charge()) return;
    if(!(ElLColl.at(0).Pt()>25 && ElLColl.at(1).Pt()>15)) return;
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      int IdxFlipped = GenLepInfo % 1000;
      if(ConvRun && GenLepInfo>=100) return;
      if(FlipRun && (GenLepInfo>999 or GenLepInfo<100)) return;
      if(FlipRun && IdxFlipped<2){ weight *= GetCFRSF(ElTColl.at(IdxFlipped), "App2Bin1_Fin"); }
    }
    float Mll = (ElLColl.at(0)+ElLColl.at(1)).M();

    if(BJetColl.size()<1) return;
      FillHist("NJPreCut"+Label+"_2E", JetColl.size(), weight, 10, 0., 10.);
    if(JetColl.size() <3) return;
    if(CheckSB_M && (Mll<101.2)) return;
    InitializeTreeVars();
    if(!CheckSB_Iso) SetVarSS2L(MuTColl, MuLColl, MuVColl, ElTColl, ElLColl, ElVColl, JetColl, BJetColl, vMET, weight, "");
    else             SetVarSS2L(MuLColl, MuLColl, MuVColl, ElLColl, ElLColl, ElVColl, JetColl, BJetColl, vMET, weight, "");
    PlotParametersSS(Label+"_2E");
  }

}


void ControlPlots::PlotParameters(TString Label){

  FillHist("Nj"+Label, Nj, w_tot, 10, 0., 10.);
  FillHist("Nb"+Label, Nb, w_tot, 5, 0., 5.);
  FillHist("Ptl1"+Label, Ptl1, w_tot, 25, 0., 250.);
  FillHist("Ptl2"+Label, Ptl2, w_tot, 30, 0., 150.);
  FillHist("Ptj1"+Label, Ptj1, w_tot, 25, 0., 500.);
  FillHist("Ptj2"+Label, Ptj2, w_tot, 30, 0., 300.);
  FillHist("Ptj3"+Label, Ptj3, w_tot, 20, 0., 200.);
  FillHist("Ptb1"+Label, Ptb1, w_tot, 35, 0., 350.);
  FillHist("Ptb2"+Label, max((Float_t)0.,Ptb2), w_tot, 15, 0., 150.);
  FillHist("MET"+Label, MET, w_tot, 30, 0., 300.);
  FillHist("HT"+Label, HT, w_tot, 25, 0., 1000.);
  FillHist("MET2HT"+Label, MET2HT, w_tot, 30, 0., 150.);
  FillHist("dEtall"+Label, dEtall, w_tot, 25, 0., 5.);
  FillHist("dRll"+Label, dRll, w_tot, 25, 0., 5.);
  FillHist("dRjj12"+Label, dRjj12, w_tot, 25, 0., 5.);
  FillHist("dRjj23"+Label, dRjj23, w_tot, 25, 0., 5.);
  FillHist("dRjj13"+Label, dRjj13, w_tot, 25, 0., 5.);
  FillHist("dRlj11"+Label, dRlj11, w_tot, 25, 0., 5.);
  FillHist("dRlj12"+Label, dRlj12, w_tot, 25, 0., 5.);
  FillHist("dRlj13"+Label, dRlj13, w_tot, 25, 0., 5.);
  FillHist("dRlj21"+Label, dRlj21, w_tot, 25, 0., 5.);
  FillHist("dRlj22"+Label, dRlj22, w_tot, 25, 0., 5.);
  FillHist("dRlj23"+Label, dRlj23, w_tot, 25, 0., 5.);
  FillHist("dRlb11"+Label, dRlb11, w_tot, 25, 0., 5.);
  FillHist("dRlb21"+Label, dRlb21, w_tot, 25, 0., 5.);
  FillHist("MSSSF"+Label, MSSSF, w_tot, 40, 0., 400.);
  FillHist("Mbl11"+Label, Mbl11, w_tot, 25, 0., 500.);
  FillHist("Mbl12"+Label, Mbl12, w_tot, 35, 0., 350.);
  FillHist("Mbl21"+Label, max((Float_t)0.,Mbl21), w_tot, 20, 0., 200.);
  FillHist("Mbl22"+Label, max((Float_t)0.,Mbl22), w_tot, 20, 0., 200.);
  FillHist("Mlj11"+Label, Mlj11, w_tot, 20, 0., 800.);
  FillHist("Mlj12"+Label, Mlj12, w_tot, 25, 0., 500.);
  FillHist("Mlj13"+Label, Mlj13, w_tot, 20, 0., 400.);
  FillHist("Mlj21"+Label, Mlj21, w_tot, 20, 0., 400.);
  FillHist("Mlj22"+Label, Mlj22, w_tot, 30, 0., 300.);
  FillHist("Mlj23"+Label, Mlj23, w_tot, 30, 0., 300.);
  FillHist("MTvl1"+Label, MTvl1, w_tot, 30, 0., 300.);
  FillHist("MTvl2"+Label, MTvl2, w_tot, 20, 0., 200.);
  FillHist("Mllj1"+Label, Mllj1, w_tot, 20, 0., 800.);
  FillHist("Mllj2"+Label, Mllj2, w_tot, 24, 0., 600.);
  FillHist("Mllj3"+Label, Mllj3, w_tot, 24, 0., 600.);
  FillHist("Mllj4"+Label, std::max((Float_t)0.,Mllj4), w_tot, 20, 0., 400.);
  FillHist("Mllb1"+Label, Mllb1, w_tot, 25, 0., 500.);
  FillHist("Mllb2"+Label, std::max((Float_t)0.,Mllb2), w_tot, 30, 0., 300.);
  FillHist("Mlljj12"+Label, Mlljj12, w_tot, 24, 0., 1200.);
  FillHist("Mlljj13"+Label, Mlljj13, w_tot, 25, 0., 1000.);
  FillHist("Mlljj14"+Label, std::max((Float_t)0.,Mlljj14), w_tot, 16, 0., 800.);
  FillHist("Mlljj23"+Label, Mlljj23, w_tot, 25, 0., 1000.);
  FillHist("Mlljj24"+Label, std::max((Float_t)0.,Mlljj24), w_tot, 15, 0., 600.);
  FillHist("Mlljj34"+Label, std::max((Float_t)0.,Mlljj34), w_tot, 15, 0., 600.);
  FillHist("Mljj112"+Label, Mljj112, w_tot, 20, 0., 1000.);
  FillHist("Mljj113"+Label, Mljj113, w_tot, 25, 0., 1000.);
  FillHist("Mljj114"+Label, std::max((Float_t)0.,Mljj114), w_tot, 20, 0., 800.);
  FillHist("Mljj123"+Label, Mljj123, w_tot, 20, 0., 800.);
  FillHist("Mljj124"+Label, std::max((Float_t)0.,Mljj124), w_tot, 24, 0., 600.);
  FillHist("Mljj134"+Label, std::max((Float_t)0.,Mljj134), w_tot, 24, 0., 600.);
  FillHist("Mljj212"+Label, Mljj212, w_tot, 25, 0., 1000.);
  FillHist("Mljj213"+Label, Mljj213, w_tot, 20, 0., 800.);
  FillHist("Mljj214"+Label, std::max((Float_t)0.,Mljj214), w_tot, 15, 0., 600.);
  FillHist("Mljj223"+Label, Mljj223, w_tot, 24, 0., 600.);
  FillHist("Mljj224"+Label, std::max((Float_t)0.,Mljj224), w_tot, 20, 0., 500.);
  FillHist("Mljj234"+Label, std::max((Float_t)0.,Mljj234), w_tot, 16, 0., 400.);
  FillHist("Mjj12"+Label, Mjj12, w_tot, 20, 0., 800.);
  FillHist("Mjj13"+Label, Mjj13, w_tot, 24, 0., 600.);
  FillHist("Mjj14"+Label, std::max((Float_t)0.,Mjj14), w_tot, 20, 0., 500.);
  FillHist("Mjj23"+Label, Mjj23, w_tot, 24, 0., 600.);
  FillHist("Mjj24"+Label, std::max((Float_t)0.,Mjj24), w_tot, 20, 0., 400.);
  FillHist("Mjj34"+Label, std::max((Float_t)0.,Mjj34), w_tot, 15, 0., 300.);
  FillHist("MllW_2jL"+Label, MllW_2jL, w_tot, 25, 0., 500.);
  FillHist("MllW_1jL"+Label, MllW_1jL, w_tot, 25, 0., 250.);
  FillHist("MllW1_H"+Label, MllW1_H, w_tot, 24, 0., 600.);
  FillHist("MllW2_H"+Label, max((Float_t)0.,MllW2_H), w_tot, 28, 0., 700.);
  FillHist("MjjW1"+Label, MjjW1, w_tot, 40, 0., 200.);
  FillHist("MjjW2"+Label, std::max((Float_t)0.,MjjW2), w_tot, 20, 0., 400.);
  FillHist("Ml1W_2jL"+Label, Ml1W_2jL, w_tot, 25, 0., 500.);
  FillHist("Ml1W_1jL"+Label, Ml1W_1jL, w_tot, 25, 0., 500.);
  FillHist("Ml2W_2jL"+Label, Ml2W_2jL, w_tot, 40, 0., 400.);
  FillHist("Ml2W_1jL"+Label, Ml2W_1jL, w_tot, 30, 0., 300.);
  FillHist("Ml1W1_H"+Label, Ml1W1_H, w_tot, 30, 0., 600.);
  FillHist("Ml1W2_H"+Label, std::max((Float_t)0.,Ml1W2_H), w_tot, 15, 0., 600.);
  FillHist("Ml2W1_H"+Label, Ml2W1_H, w_tot, 40, 0., 400.);
  FillHist("Ml2W2_H"+Label, std::max((Float_t)0.,Ml2W2_H), w_tot, 20, 0., 500.);
  FillHist("disc_BDTG"+Label, disc_BDTG, w_tot, 40, -1., 1.);

}


void ControlPlots::PlotParametersSS(TString Label){

  FillHist("Nj"+Label, Nj, w_tot, 10, 0., 10.);
  FillHist("Nb"+Label, Nb, w_tot, 5, 0., 5.);
  FillHist("Ptl1"+Label, Ptl1, w_tot, 15, 0., 300.);
  FillHist("Ptl2"+Label, Ptl2, w_tot, 15, 0., 150.);
  FillHist("Ptj1"+Label, Ptj1, w_tot, 10, 0., 500.);
  FillHist("Ptj2"+Label, Ptj2, w_tot, 12, 0., 300.);
  FillHist("Ptj3"+Label, Ptj3, w_tot,  8, 0., 200.);
  FillHist("Ptb1"+Label, Ptb1, w_tot, 12, 0., 300.);
  FillHist("Ptb2"+Label, max((Float_t)0.,Ptb2), w_tot, 6, 0., 150.);
  FillHist("MET"+Label, MET, w_tot, 12, 0., 300.);
  FillHist("HT"+Label, HT, w_tot, 12, 0., 900.);
  FillHist("MET2HT"+Label, MET2HT, w_tot, 15, 0., 150.);
  FillHist("dEtall"+Label, dEtall, w_tot, 12, 0., 4.8);
  FillHist("dRll"+Label, dRll, w_tot, 15, 0., 6.);
  FillHist("dRjj12"+Label, dRjj12, w_tot, 15, 0., 6.);
  FillHist("dRjj23"+Label, dRjj23, w_tot, 15, 0., 6.);
  FillHist("dRjj13"+Label, dRjj13, w_tot, 15, 0., 6.);
  FillHist("dRlj11"+Label, dRlj11, w_tot, 15, 0., 6.);
  FillHist("dRlj12"+Label, dRlj12, w_tot, 15, 0., 6.);
  FillHist("dRlj13"+Label, dRlj13, w_tot, 15, 0., 6.);
  FillHist("dRlj21"+Label, dRlj21, w_tot, 15, 0., 6.);
  FillHist("dRlj22"+Label, dRlj22, w_tot, 15, 0., 6.);
  FillHist("dRlj23"+Label, dRlj23, w_tot, 15, 0., 6.);
  FillHist("dRlb11"+Label, dRlb11, w_tot, 15, 0., 6.);
  FillHist("dRlb21"+Label, dRlb21, w_tot, 15, 0., 6.);
  FillHist("MSSSF"+Label, MSSSF, w_tot, 40, 0., 400.);
  FillHist("Mbl11"+Label, Mbl11, w_tot, 10, 0., 500.);
  FillHist("Mbl12"+Label, Mbl12, w_tot, 16, 0., 400.);
  FillHist("Mbl21"+Label, max((Float_t)0.,Mbl21), w_tot, 12, 0., 300.);
  FillHist("Mbl22"+Label, max((Float_t)0.,Mbl22), w_tot, 12, 0., 300.);
  FillHist("Mlj11"+Label, Mlj11, w_tot, 10, 0., 800.);
  FillHist("Mlj12"+Label, Mlj12, w_tot, 10, 0., 500.);
  FillHist("Mlj13"+Label, Mlj13, w_tot, 10, 0., 500.);
  FillHist("Mlj21"+Label, Mlj21, w_tot, 10, 0., 500.);
  FillHist("Mlj22"+Label, Mlj22, w_tot, 15, 0., 450.);
  FillHist("Mlj23"+Label, Mlj23, w_tot, 15, 0., 450.);
  FillHist("MTvl1"+Label, MTvl1, w_tot, 16, 0., 400.);
  FillHist("MTvl2"+Label, MTvl2, w_tot, 12, 0., 300.);
  FillHist("Mllj1"+Label, Mllj1, w_tot, 20, 0., 1000.);
  FillHist("Mllj2"+Label, Mllj2, w_tot, 16, 0., 800.);
  FillHist("Mllj3"+Label, Mllj3, w_tot, 16, 0., 800.);	
  FillHist("Mllj4"+Label, std::max((Float_t)0.,Mllj4), w_tot, 12, 0., 600.);
  FillHist("Mllb1"+Label, Mllb1, w_tot, 12, 0., 600.);
  FillHist("Mllb2"+Label, std::max((Float_t)0.,Mllb2), w_tot, 8, 0., 400.);
  FillHist("Mlljj12"+Label, Mlljj12, w_tot, 15, 0., 1500.);
  FillHist("Mlljj13"+Label, Mlljj13, w_tot, 24, 0., 1200.);
  FillHist("Mlljj14"+Label, std::max((Float_t)0.,Mlljj14), w_tot, 10, 0., 1000.);
  FillHist("Mlljj23"+Label, Mlljj23, w_tot, 20, 0., 1000.);
  FillHist("Mlljj24"+Label, std::max((Float_t)0.,Mlljj24), w_tot, 16, 0., 800.);
  FillHist("Mlljj34"+Label, std::max((Float_t)0.,Mlljj34), w_tot, 16, 0., 800.);
  FillHist("Mljj112"+Label, Mljj112, w_tot, 12, 0., 1200.);
  FillHist("Mljj113"+Label, Mljj113, w_tot, 12, 0., 1200.);
  FillHist("Mljj114"+Label, std::max((Float_t)0.,Mljj114), w_tot, 10, 0., 1000.);
  FillHist("Mljj123"+Label, Mljj123, w_tot, 10, 0., 1000.);
  FillHist("Mljj124"+Label, std::max((Float_t)0.,Mljj124), w_tot, 16, 0., 800.);
  FillHist("Mljj134"+Label, std::max((Float_t)0.,Mljj134), w_tot, 16, 0., 800.);
  FillHist("Mljj212"+Label, Mljj212, w_tot, 20, 0., 1000.);
  FillHist("Mljj213"+Label, Mljj213, w_tot, 20, 0., 1000.);
  FillHist("Mljj214"+Label, std::max((Float_t)0.,Mljj214), w_tot, 16, 0., 800.);
  FillHist("Mljj223"+Label, Mljj223, w_tot, 16, 0., 800.);
  FillHist("Mljj224"+Label, std::max((Float_t)0.,Mljj224), w_tot, 12, 0., 600.);
  FillHist("Mljj234"+Label, std::max((Float_t)0.,Mljj234), w_tot, 12, 0., 600.);
  FillHist("Mjj12"+Label, Mjj12, w_tot, 20, 0., 1000.);
  FillHist("Mjj13"+Label, Mjj13, w_tot, 16, 0., 800.);
  FillHist("Mjj14"+Label, std::max((Float_t)0.,Mjj14), w_tot, 10, 0., 500.);
  FillHist("Mjj23"+Label, Mjj23, w_tot, 12, 0., 600.);
  FillHist("Mjj24"+Label, std::max((Float_t)0.,Mjj24), w_tot, 8, 0., 400.);
  FillHist("Mjj34"+Label, std::max((Float_t)0.,Mjj34), w_tot, 6, 0., 300.);
  FillHist("MllW_2jL"+Label, MllW_2jL, w_tot, 12, 0., 600.);
  FillHist("MllW_1jL"+Label, MllW_1jL, w_tot, 12, 0., 300.);
  FillHist("MllW1_H"+Label, MllW1_H, w_tot, 16, 0., 800.);
  FillHist("MllW2_H"+Label, max((Float_t)0.,MllW2_H), w_tot, 16, 0., 800.);
  FillHist("MjjW1"+Label, MjjW1, w_tot, 20, 0., 200.);
  FillHist("MjjW2"+Label, std::max((Float_t)0.,MjjW2), w_tot, 8, 0., 400.);
  FillHist("Ml1W_2jL"+Label, Ml1W_2jL, w_tot, 12, 0., 600.);
  FillHist("Ml1W_1jL"+Label, Ml1W_1jL, w_tot, 12, 0., 600.);
  FillHist("Ml2W_2jL"+Label, Ml2W_2jL, w_tot, 16, 0., 400.);
  FillHist("Ml2W_1jL"+Label, Ml2W_1jL, w_tot, 10, 0., 400.);
  FillHist("Ml1W1_H"+Label, Ml1W1_H, w_tot, 16, 0., 800.);
  FillHist("Ml1W2_H"+Label, std::max((Float_t)0.,Ml1W2_H), w_tot, 16, 0., 800.);
  FillHist("Ml2W1_H"+Label, Ml2W1_H, w_tot, 12, 0., 600.);
  FillHist("Ml2W2_H"+Label, std::max((Float_t)0.,Ml2W2_H), w_tot, 12, 0., 600.);
  FillHist("disc_BDTG"+Label, disc_BDTG, w_tot, 40, -1., 1.);

}



void ControlPlots::executeEventFromParameter(AnalyzerParameter param){
}


void ControlPlots::InitializeTreeVars(){

  Nj=-1, Nb=-1;
  Ptl1=-1, Ptl2=-1, Ptj1=-1, Ptj2=-1, Ptj3=-1, Ptb1=-1, Ptb2=-1, MET=-1, HT=-1, MET2HT=-1;
  dEtall=-1, dRll=-1, dRjj12=-1, dRjj23=-1, dRjj13=-1;
  dRlj11=-1, dRlj12=-1, dRlj13=-1, dRlj21=-1, dRlj22=-1, dRlj23=-1;
  dRlb11=-1, dRlb21=-1;
  MSSSF=-1, Mbl11=-1, Mbl12=-1, Mbl21=-1, Mbl22=-1, Mlj11=-1, Mlj12=-1, Mlj13=-1, Mlj21=-1, Mlj22=-1, Mlj23=-1;
  MTvl1=-1, MTvl2=-1, Mllj1=-1, Mllj2=-1, Mllj3=-1, Mllj4=-1, Mllb1=-1, Mllb2=-1;
  Mlljj12=-1, Mlljj13=-1, Mlljj14=-1, Mlljj23=-1, Mlljj24=-1, Mlljj34=-1;
  Mljj112=-1, Mljj113=-1, Mljj114=-1, Mljj123=-1, Mljj124=-1, Mljj134=-1;
  Mljj212=-1, Mljj213=-1, Mljj214=-1, Mljj223=-1, Mljj224=-1, Mljj234=-1;
  Mjj12=-1, Mjj13=-1, Mjj14=-1, Mjj23=-1, Mjj24=-1, Mjj34=-1;
  MllW_2jL=-1, MllW_1jL=-1, MllW1_H=-1, MllW2_H=-1, MjjW1=-1, MjjW2=-1;
  Ml1W_2jL=-1, Ml1W_1jL=-1, Ml2W_2jL=-1, Ml2W_1jL=-1, Ml1W1_H=-1, Ml1W2_H=-1, Ml2W1_H=-1, Ml2W2_H=-1;
  w_tot=-1, disc_BDTG=-999.;

}



float ControlPlots::CalcTestFakeWeight(vector<Muon>& MuColl, vector<Electron>& ElColl, TString MuIDT, TString MuIDL, TString ElIDT, TString ElIDL, int NBJet, int SystDir){

  float w_FR=-1.; int NLepLNotT=0;
  for(unsigned int im=0; im<MuColl.size(); im++){
    if(MuColl.at(im).PassID(MuIDT)) continue;
    float FR = GetTestMuFR(MuColl.at(im), MuIDT+"_"+MuIDL+"_Incl", SystDir);
    w_FR*=-FR/(1.-FR);
    NLepLNotT++;
  }

  for(unsigned int ie=0; ie<ElColl.size(); ie++){
    if(ElColl.at(ie).PassID(ElIDT)) continue;
    TString SelCond = NBJet==0? "Incl":"HasB";
    float FR = GetTestElFR(ElColl.at(ie), ElIDT+"_"+ElIDL+"_"+SelCond, SystDir);
    w_FR*=-FR/(1.-FR);
    NLepLNotT++;
  }
  if(NLepLNotT==0) w_FR=0.;
  
  return w_FR;  
}


float ControlPlots::GetTestElFR(Electron& El, TString Key, int SystDir){

  float FR=0., PTCorr=0., fEta=fabs(El.Eta());
  if(Key=="TopHN17SST_TopHN17SSL_Incl"){
    float TightIso=0.1;
    PTCorr = El.CalcPtCone(El.MiniRelIso(), TightIso);
    if(fEta<0.8){
      if     (PTCorr<15) FR=0.;
      else if(PTCorr<20) FR=0.379217*(1.+(SystDir==0? 0.:SystDir>0? 0.354774 :-0.0106956));
      else if(PTCorr<25) FR=0.224299*(1.+(SystDir==0? 0.:SystDir>0? 0.281698 :-0.0359124));
      else if(PTCorr<35) FR=0.154058*(1.+(SystDir==0? 0.:SystDir>0? 0.22729  :-0.0616697));
      else               FR=0.159764*(1.+(SystDir==0? 0.:SystDir>0? 0.0699418:-0.209737 ));
    }
    else if(fEta<1.479){
      if     (PTCorr<15) FR=0.;
      else if(PTCorr<20) FR=0.3647  *(1.+(SystDir==0? 0.:SystDir>0? 0.351692:-0.027291 ));
      else if(PTCorr<25) FR=0.236944*(1.+(SystDir==0? 0.:SystDir>0? 0.466218:-0.0097618));
      else if(PTCorr<35) FR=0.183473*(1.+(SystDir==0? 0.:SystDir>0? 0.647966:-0.0151806));
      else               FR=0.167406*(1.+(SystDir==0? 0.:SystDir>0? 0.192329:-0.0774894));
    }
    else{
      if     (PTCorr<15) FR=0.;
      else if(PTCorr<20) FR=0.43818*(1.+(SystDir==0? 0.:SystDir>0? 0.267022:-0.0116451));
      else if(PTCorr<25) FR=0.30755*(1.+(SystDir==0? 0.:SystDir>0? 0.463056:-0.0087176));
      else if(PTCorr<35) FR=0.24232*(1.+(SystDir==0? 0.:SystDir>0? 0.264124:-0.0071143));
      else               FR=0.18749*(1.+(SystDir==0? 0.:SystDir>0? 0.559945:-0.107224 ));
    }
  }
  else if(Key=="TopHN17SST_TopHN17SSL_HasB"){
    float TightIso=0.1;
    PTCorr = El.CalcPtCone(El.MiniRelIso(), TightIso);
    if(fEta<0.8){
      if     (PTCorr<15) FR=0.;
      else if(PTCorr<20) FR=0.379217*(1.353263+(SystDir==0? 0.:SystDir>0? 0.0106956:-0.354774 ));
      else if(PTCorr<25) FR=0.224299*(1.249031+(SystDir==0? 0.:SystDir>0? 0.0359124:-0.281698 ));
      else if(PTCorr<35) FR=0.154058*(1.222836+(SystDir==0? 0.:SystDir>0? 0.0616697:-0.22729  ));
      else               FR=0.159764*(1.      +(SystDir==0? 0.:SystDir>0? 0.209737 :-0.0699418));
    }
    else if(fEta<1.479){
      if     (PTCorr<15) FR=0.;
      else if(PTCorr<20) FR=0.3647  *(1.351683+(SystDir==0? 0.:SystDir>0? 0.027291 :-0.351692));
      else if(PTCorr<25) FR=0.236944*(1.465871+(SystDir==0? 0.:SystDir>0? 0.0097618:-0.466218));
      else if(PTCorr<35) FR=0.183473*(1.64359 +(SystDir==0? 0.:SystDir>0? 0.0151806:-0.647966));
      else               FR=0.167406*(1.184503+(SystDir==0? 0.:SystDir>0? 0.0774894:-0.192329));
    }
    else{
      if     (PTCorr<15) FR=0.;
      else if(PTCorr<20) FR=0.43818 *(1.266234+(SystDir==0? 0.:SystDir>0? 0.0116451:-0.267022));
      else if(PTCorr<25) FR=0.307551*(1.460599+(SystDir==0? 0.:SystDir>0? 0.0087176:-0.463056));
      else if(PTCorr<35) FR=0.24232 *(1.263883+(SystDir==0? 0.:SystDir>0? 0.0071143:-0.264124));
      else               FR=0.18749 *(1.559198+(SystDir==0? 0.:SystDir>0? 0.107224 :-0.559945));
    }
  }

  return FR;
}


float ControlPlots::GetTestMuFR(Muon& Mu, TString Key, int SystDir){

  float FR=0., PTCorr=0., fEta=fabs(Mu.Eta());
  if(Key=="TopHN17T_TopHN17L_Incl"){
    float TightIso=0.1;
    PTCorr = Mu.CalcPtCone(Mu.MiniRelIso(), TightIso);
    if(fEta<0.9){
      if     (PTCorr<10) FR=0.;
      else if(PTCorr<15) FR=0.468723*(1.+(SystDir==0? 0.:SystDir>0? 0.0996232:-0.00039495));
      else if(PTCorr<20) FR=0.258502*(1.+(SystDir==0? 0.:SystDir>0? 0.100746 :-0.128182  ));
      else if(PTCorr<30) FR=0.170125*(1.+(SystDir==0? 0.:SystDir>0? 0.247411 :-0.0809608 ));
      else               FR=0.10884 *(1.+(SystDir==0? 0.:SystDir>0? 0.128926 :-0.166817  ));
    }
    else if(fEta<1.6){
      if     (PTCorr<10) FR=0.;
      else if(PTCorr<15) FR=0.507634*(1.+(SystDir==0? 0.:SystDir>0? 0.0767083:-0.0461033));
      else if(PTCorr<20) FR=0.302245*(1.+(SystDir==0? 0.:SystDir>0? 0.0849491:-0.0551178));
      else if(PTCorr<30) FR=0.223215*(1.+(SystDir==0? 0.:SystDir>0? 0.121284 :-0.144019 ));
      else               FR=0.150397*(1.+(SystDir==0? 0.:SystDir>0? 0.0901934:-0.212049 ));
    }
    else{
      if     (PTCorr<10) FR=0.;
      else if(PTCorr<15) FR=0.583916*(1.+(SystDir==0? 0.:SystDir>0? 0.126021 :-0.00025308));
      else if(PTCorr<20) FR=0.330751*(1.+(SystDir==0? 0.:SystDir>0? 0.153935 :-0.0504944 ));
      else if(PTCorr<30) FR=0.276111*(1.+(SystDir==0? 0.:SystDir>0? 0.17514  :-0.00782237));
      else               FR=0.211517*(1.+(SystDir==0? 0.:SystDir>0? 0.0826472:-0.0732807 ));
    }
  }

  return FR;
}




float ControlPlots::GetCFRSF(Electron& El, TString Tag, TString Option){

  if     (IsDATA  && Option=="") return 1.;
  else if(IsDATA  && Option!="DataEff") return 0.;
  else if(!IsDATA && Option=="DataEff") return 0.;

  bool ReturnMCEff=Option.Contains("MCEff"), ReturnDataEff=Option.Contains("DataEff");
  float PT=El.Pt(), fEta=fabs(El.Eta());
  int EtaRegIndex1 = fEta<1.47? 0:1;
  int EtaRegIndex2 = fEta<0.8? 0:fEta<1.47? 1: fEta<2.? 2:3;
  int PtRegIndex1  = PT<35? 0:PT<50? 1: 2;
  int PtRegIndex2  = PT<35? 0:PT<50? 1:PT<100? 2:3;
  int BinIndex1 = EtaRegIndex1*3+PtRegIndex1+1; 
  int BinIndex2 = EtaRegIndex2+1; 
  int BinIndex3 = PtRegIndex2 +1; 

  float SF=1., MCEff=1., DataEff=1., ReturnVal=1.;
  if(Tag.Contains("App2Bin1_")){
    int Idx = BinIndex1;
    if(Tag.EndsWith("It1")) SF= Idx==1? 1.49:Idx==2? 1.50:Idx==3? 1.82:Idx==4? 1.51:Idx==5? 1.41:Idx==6? 1.28:1.;
    if(Tag.EndsWith("It2")) SF= Idx==1? 1.03:Idx==2? 1.02:Idx==3? 1.15:Idx==4? 1.02:Idx==5? 1.01:Idx==6? 0.97:1.;
    if(Tag.EndsWith("It3")) SF= Idx==1? 1.01:Idx==2? 1.01:Idx==3? 1.08:Idx==4? 1.00:Idx==5? 1.00:Idx==6? 0.98:1.;
    if(Tag.EndsWith("It4")) SF= Idx==1? 1.01:Idx==2? 1.00:Idx==3? 1.05:Idx==4? 1.00:Idx==5? 1.00:Idx==6? 0.99:1.;
    if(Tag.EndsWith("Fin")){
       SF      = Idx==1? 1.58:Idx==2? 1.54:Idx==3? 2.45:Idx==4? 1.55:Idx==5? 1.41:Idx==6? 1.21:1.;
       DataEff = Idx==1? 7.98E-5:Idx==2? 8.11E-5:Idx==3? 1.50E-4:Idx==4? 1.48E-3:Idx==5? 1.48E-3:Idx==6? 1.67E-3:1.;
       MCEff   = Idx==1? 5.06E-5:Idx==2? 5.25E-5:Idx==3? 6.15E-5:Idx==4? 9.58E-4:Idx==5? 1.04E-3:Idx==6? 1.39E-3:1.;
    }
  }
  else if(Tag.Contains("App2Bin2_")){
    int Idx = BinIndex2;
    if(Tag.EndsWith("It1")) SF= Idx==1? 1.48:Idx==2? 1.56:Idx==3? 1.36:Idx==4? 1.39:1.;
    if(Tag.EndsWith("It2")) SF= Idx==1? 1.05:Idx==2? 1.08:Idx==3? 0.99:Idx==4? 1.00:1.;
    if(Tag.EndsWith("It3")) SF= Idx==1? 1.04:Idx==2? 1.05:Idx==3? 0.99:Idx==4? 1.00:1.;
    if(Tag.EndsWith("It4")) SF= Idx==1? 1.03:Idx==2? 1.03:Idx==3? 1.00:Idx==4? 1.00:1.;
    if(Tag.EndsWith("Fin")){
      SF      = Idx==1? 1.70:Idx==2? 1.87:Idx==3? 1.34:Idx==4? 1.36:1.;
      DataEff = Idx==1? 3.07E-5:Idx==2? 2.03E-4:Idx==3? 1.11E-3:Idx==4? 2.00E-3:1.;
      MCEff   = Idx==1? 1.81E-5:Idx==2? 1.09E-4:Idx==3? 8.31E-4:Idx==4? 1.45E-3:1.;
    }
  }
  else if(Tag.Contains("App2Bin3_")){
    int Idx = BinIndex3;
    if(Tag.EndsWith("It1")) SF= Idx==1? 1.40:Idx==2? 1.44:Idx==3? 1.42:Idx==4? 1.03:1.;
    if(Tag.EndsWith("It2")) SF= Idx==1? 1.00:Idx==2? 1.00:Idx==3? 1.00:Idx==4? 1.13:1.;
    if(Tag.EndsWith("It3")) SF= Idx==1? 1.00:Idx==2? 1.00:Idx==3? 1.00:Idx==4? 1.05:1.;
    if(Tag.EndsWith("It4")) SF= Idx==1? 1.00:Idx==2? 1.00:Idx==3? 1.00:Idx==4? 1.01:1.;
    if(Tag.EndsWith("Fin")){
      SF      = Idx==1? 1.39:Idx==2? 1.46:Idx==3? 1.43:Idx==4? 1.24:1.;
      DataEff = Idx==1? 4.22E-4:Idx==2? 4.11E-4:Idx==3? 5.15E-4:Idx==4? 6.60E-4:1.;
      MCEff   = Idx==1? 3.02E-4:Idx==2? 2.84E-4:Idx==3? 3.62E-4:Idx==4? 5.33E-4:1.;
    }
  }

  if     (ReturnMCEff  ) ReturnVal = MCEff;
  else if(ReturnDataEff) ReturnVal = DataEff;
  else                   ReturnVal = SF; 

  return ReturnVal;
}


float ControlPlots::GetFlipCorrPT(Electron& El, TString Tag, TString Option){

  if(!(IsDATA  && Option.Contains("Data"))) return 0.;

  bool DoScale = Option.Contains("Scale"), DoSmear = Option.Contains("Smear");
  float PT=El.Pt(), fEta=fabs(El.Eta());
  int EtaRegIndex1 = fEta<1.47? 0:1;
  int EtaRegIndex2 = fEta<0.8? 0:fEta<1.47? 1: fEta<2.? 2:3;
  int PtRegIndex1  = PT<35? 0:PT<50? 1: 2;
  int PtRegIndex2  = PT<35? 0:PT<50? 1:PT<100? 2:3;
  int BinIndex1 = EtaRegIndex1*3+PtRegIndex1+1; 
  int BinIndex2 = EtaRegIndex2+1; 
  int BinIndex3 = PtRegIndex2 +1; 

  float PTScaleCorr=0., PTResCorr=0., ReturnPT=PT, RelRes=0.;
  if(Tag.Contains("App2Bin1_")){
    int Idx = BinIndex1;
    PTScaleCorr = Idx==1? 0.978:Idx==2? 0.981:Idx==3? 0.988:Idx==4? 0.972:Idx==5? 0.985:Idx==6? 0.990:0.;
    PTResCorr   = Idx==1?  1.38:Idx==2?  1.47:Idx==3?  1.31:Idx==4?  1.15:Idx==5?  1.08:Idx==6? 0.981:0.;
    RelRes      = Idx==1? 0.0322:Idx==2? 0.0276:Idx==3? 0.0264:Idx==4? 0.0388:Idx==5? 0.0375:Idx==6? 0.0382:0.;
  }
  else if(Tag.Contains("App2Bin2_")){//POG res.+const. scale res. matching mee dist.
    int Idx = BinIndex2;
    PTScaleCorr = 0.985;
    PTResCorr   = 1.75;
    RelRes      = Idx==1? 0.01:Idx==2? 0.0153:Idx==3? 0.0217:Idx==4? 0.0208:0.;
  }
  else if(Tag.Contains("App2Bin3_")){
    int Idx = BinIndex3;
    PTScaleCorr = Idx==1? 0.976:Idx==2? 0.985:Idx==3? 0.992:Idx==4? 0.999:0.;
    PTResCorr   = Idx==1?  1.29:Idx==2?  1.41:Idx==3?  1.32:Idx==4?  1.18:0.;
    RelRes      = Idx==1? 0.0343:Idx==2? 0.0287:Idx==3? 0.0285:Idx==4? 0.0297:0.;
  }
  PTScaleCorr = 1.-2.*(1.-PTScaleCorr);
  PTResCorr   = 1.-2.*(1.-PTResCorr);

  if(DoScale) ReturnPT = PT*PTScaleCorr;
  if(DoSmear) ReturnPT = ReturnPT*( 1.+gRandom->Gaus(0,RelRes)*sqrt(max(pow(PTResCorr,2.)-1,0.)) );

  return ReturnPT;
}




int ControlPlots::GetGenLepInfo(vector<Electron>& ElColl, vector<Muon>& MuColl, TString Option){

  vector<Gen> TruthColl = GetGens();
  int NFk=0, NFlip=0, NCv=0, IdxFlipped=9;
  for(unsigned int im=0; im<MuColl.size(); im++){
    int LepType = GetLeptonType_JH(MuColl.at(im), TruthColl);
    if(LepType<0 && LepType>-5) NFk++;
  }
  for(unsigned int ie=0; ie<ElColl.size(); ie++){
    int LepType = GetLeptonType_JH(ElColl.at(ie), TruthColl);
    if(LepType<0 && LepType>-5) NFk++;
    else if(LepType<-4) NCv++;
    else if(LepType>0){
      int Idx_Closest    = GenMatchedIdx(ElColl.at(ie),TruthColl);
      int IdxType_NearEl = LepType>3? GetPrElType_InSameSCRange(Idx_Closest, TruthColl, "IdxType"):Idx_Closest;
      int Idx_NearEl     = LepType>3? IdxType_NearEl/10:Idx_Closest;
      if(ElColl.at(ie).Charge()*TruthColl.at(Idx_NearEl).PID()>0){ NFlip++; IdxFlipped=ie; }
    }
  }

  int ReturnVal = NFk*1000+NFlip*100+NCv*10+IdxFlipped;

  return ReturnVal;

}



ControlPlots::ControlPlots(){

  TMVA::Tools::Instance();
  MVAreader_Mu = new TMVA::Reader();
  MVAreader_El = new TMVA::Reader();

}


ControlPlots::~ControlPlots(){

  delete MVAreader_Mu;
  delete MVAreader_El;

}


void ControlPlots::InitializeReader(TMVA::Reader *MVAreader, TString FileName){

  MVAreader->AddVariable("Nj"      , &Nj      );
  MVAreader->AddVariable("Nb"      , &Nb      );
  MVAreader->AddVariable("Ptl1"    , &Ptl1    );
  MVAreader->AddVariable("Ptl2"    , &Ptl2    );
  MVAreader->AddVariable("Ptj1"    , &Ptj1    );
  MVAreader->AddVariable("Ptj2"    , &Ptj2    );
  MVAreader->AddVariable("Ptj3"    , &Ptj3    );
  MVAreader->AddVariable("Ptb1"    , &Ptb1    );
  MVAreader->AddVariable("dEtall"  , &dEtall  );
  MVAreader->AddVariable("dRll"    , &dRll    );
  MVAreader->AddVariable("dRjj12"  , &dRjj12  );
  MVAreader->AddVariable("dRjj23"  , &dRjj23  );
  MVAreader->AddVariable("dRjj13"  , &dRjj13  );
  MVAreader->AddVariable("dRlj11"  , &dRlj11  );
  MVAreader->AddVariable("dRlj12"  , &dRlj12  );
  MVAreader->AddVariable("dRlj13"  , &dRlj13  );
  MVAreader->AddVariable("dRlj21"  , &dRlj21  );
  MVAreader->AddVariable("dRlj22"  , &dRlj22  );
  MVAreader->AddVariable("dRlj23"  , &dRlj23  );
  MVAreader->AddVariable("dRlb11"  , &dRlb11  );
  MVAreader->AddVariable("dRlb21"  , &dRlb21  );
  MVAreader->AddVariable("MSSSF"   , &MSSSF   );
  MVAreader->AddVariable("Mbl11"   , &Mbl11   );
  MVAreader->AddVariable("Mbl12"   , &Mbl12   );
  MVAreader->AddVariable("Mllb1"   , &Mllb1   );
  MVAreader->AddVariable("MTvl1"   , &MTvl1   );
  MVAreader->AddVariable("MTvl2"   , &MTvl2   );
  MVAreader->AddVariable("Mlj11"   , &Mlj11   );
  MVAreader->AddVariable("Mlj12"   , &Mlj12   );
  MVAreader->AddVariable("Mlj13"   , &Mlj13   );
  MVAreader->AddVariable("Mlj21"   , &Mlj21   );
  MVAreader->AddVariable("Mlj22"   , &Mlj22   );
  MVAreader->AddVariable("Mlj23"   , &Mlj23   );
  MVAreader->AddVariable("Mllj1"   , &Mllj1   );
  MVAreader->AddVariable("Mllj2"   , &Mllj2   );
  MVAreader->AddVariable("Mllj3"   , &Mllj3   );
  MVAreader->AddVariable("Mlljj12" , &Mlljj12 );
  MVAreader->AddVariable("Mlljj13" , &Mlljj13 );
  MVAreader->AddVariable("Mlljj23" , &Mlljj23 );
  MVAreader->AddVariable("Mljj112" , &Mljj112 );
  MVAreader->AddVariable("Mljj113" , &Mljj113 );
  MVAreader->AddVariable("Mljj123" , &Mljj123 );
  MVAreader->AddVariable("Mljj212" , &Mljj212 );
  MVAreader->AddVariable("Mljj213" , &Mljj213 );
  MVAreader->AddVariable("Mljj223" , &Mljj223 );
  MVAreader->AddVariable("Mjj12"   , &Mjj12   );
  MVAreader->AddVariable("Mjj13"   , &Mjj13   );
  MVAreader->AddVariable("Mjj23"   , &Mjj23   );
  MVAreader->AddVariable("HT"      , &HT      );
  MVAreader->AddVariable("MET2HT"  , &MET2HT  );
  MVAreader->AddVariable("MllW_2jL", &MllW_2jL);
  MVAreader->AddVariable("MllW_1jL", &MllW_1jL);
  MVAreader->AddVariable("MllW1_H" , &MllW1_H );
  MVAreader->AddVariable("Ml1W_2jL", &Ml1W_2jL);
  MVAreader->AddVariable("Ml1W_1jL", &Ml1W_1jL);
  MVAreader->AddVariable("Ml2W_2jL", &Ml2W_2jL);
  MVAreader->AddVariable("Ml2W_1jL", &Ml2W_1jL);
  MVAreader->AddVariable("Ml1W1_H" , &Ml1W1_H );
  MVAreader->AddVariable("Ml1W2_H" , &Ml1W2_H );
  MVAreader->AddVariable("Ml2W1_H" , &Ml2W1_H );
  MVAreader->AddVariable("Ml2W2_H" , &Ml2W2_H );
  MVAreader->AddVariable("MjjW1"   , &MjjW1   );
  MVAreader->AddVariable("MjjW2"   , &MjjW2   );
  MVAreader->AddVariable("Ptb2"    , &Ptb2    );
  MVAreader->AddVariable("Mbl21"   , &Mbl21   );
  MVAreader->AddVariable("Mbl22"   , &Mbl22   );
  MVAreader->AddVariable("Mllb2"   , &Mllb2   );
  MVAreader->AddVariable("Mllj4"   , &Mllj4   );
  MVAreader->AddVariable("Mlljj14" , &Mlljj14 );
  MVAreader->AddVariable("Mlljj24" , &Mlljj24 );
  MVAreader->AddVariable("Mlljj34" , &Mlljj34 );
  MVAreader->AddVariable("Mljj114" , &Mljj114 );
  MVAreader->AddVariable("Mljj124" , &Mljj124 );
  MVAreader->AddVariable("Mljj134" , &Mljj134 );
  MVAreader->AddVariable("Mljj214" , &Mljj214 );
  MVAreader->AddVariable("Mljj224" , &Mljj224 );
  MVAreader->AddVariable("Mljj234" , &Mljj234 );
  MVAreader->AddVariable("Mjj14"   , &Mjj14   );
  MVAreader->AddVariable("Mjj24"   , &Mjj24   );
  MVAreader->AddVariable("Mjj34"   , &Mjj34   );
  MVAreader->AddVariable("MllW2_H" , &MllW2_H );

  TString AnalyzerPath=std::getenv("SKFlat_WD"), MVAPath="/data/Run2Legacy_v4/2017/MVAClassifier/";
  MVAreader->BookMVA("BDTG method", AnalyzerPath+MVAPath+FileName);

}

