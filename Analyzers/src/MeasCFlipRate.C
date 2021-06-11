#include "MeasCFlipRate.h"

void MeasCFlipRate::initializeAnalyzer(){

  TriLep=false, TetraLep=false, SS2l=false, OS2l=false;
  CFlip=false, MCCFRate=false, CFMCClos=false, MDists=false;
  FakeRun=false, ConvRun=false, FlipRun=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("OS2l"))       OS2l       = true;
    if(Userflags.at(i).Contains("SS2l"))       SS2l       = true;
    if(Userflags.at(i).Contains("TriLep"))     TriLep     = true;
    if(Userflags.at(i).Contains("TetraLep"))   TetraLep   = true; 
    if(Userflags.at(i).Contains("CFlip"))      CFlip      = true;
    if(Userflags.at(i).Contains("CFMCClos"))     CFMCClos     = true;
    if(Userflags.at(i).Contains("MDists"))     MDists     = true;
    if(Userflags.at(i).Contains("MCCFRate"))   MCCFRate   = true;
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

  outfile->cd();

}


void MeasCFlipRate::executeEvent(){


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
  else if(SS2l or OS2l){
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
  int NPreMu=muonPreColl.size(), NPreEl=electronPreColl.size();
  if(TriLep and (NPreMu+NPreEl)>2 ) PreCutPass=true;
  if(TetraLep and (NPreMu+NPreEl)>3 ) PreCutPass=true;
  if(SS2l){
    if( (NPreMu+NPreEl)>2 ) PreCutPass=true;
    else if(NPreMu==2 and SumCharge(muonPreColl)!=0) PreCutPass=true;
    else if(NPreEl==2 and SumCharge(electronPreColl)!=0) PreCutPass=true;
  }
  if(OS2l && (NPreMu+NPreEl)>1) PreCutPass=true;
  if(OS2l && (MCCFRate or CFMCClos or MDists) && (NPreMu+NPreEl)>1 && NPreEl>0) PreCutPass=true;
  
  if(!PreCutPass) return;


  TString IDSSLabel = "SS", TLabel = FakeRun? "L":"T";
  //TString IDSSLabel = SS2l or OS2l? "SS":"";
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
  if(SS2l or OS2l){
    EventCand = (muonTightColl.size()+electronTightColl.size())==2;
    if(MCCFRate or CFMCClos or MDists) EventCand = EventCand && electronTightColl.size()>0;
  }
  if(TriLep){ EventCand = (muonLooseColl.size()+electronLooseColl.size())==3; }
  if(TetraLep){ EventCand = (muonLooseColl.size()+electronLooseColl.size())>3; }

  float w_TopPtRW = 1., w_Prefire = 1., sf_Trig = 1., w_FR=1.;
  float sf_MuTk = 1., sf_MuID = 1., sf_MuIso = 1., sf_ElReco = 1., sf_ElID = 1., sf_BTag = 1.;
  if((!IsDATA) and EventCand){
    //if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    //w_TopPtRW = mcCorr->GetTopPtReweight(truthColl);
    sf_MuID   = GetMuonSF(muonTightColl, "TopHNID_TkMu", "ID");
    //sf_MuIso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    sf_ElReco = GetElectronSF(electronTightColl, "", "Reco");
    sf_ElID   = GetElectronSF(electronTightColl, "TopHNID"+IDSSLabel, "ID");
    sf_BTag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    TString SFKey_Trig_All = muonTightColl.size()==2? "DiMuIso_HNTopID":electronTightColl.size()==2? "DiElIso_HNTopIDSS":"EMuIso_HNTopIDSS"; 
    sf_Trig   = SS2l or OS2l? mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig_All, ""):1.;
    w_Prefire = GetPrefireWeight(0);
    //cout<<" w_PU:"<<w_PU<<" w_Prefire:"<<w_Prefire<<" sf_Trig:"<<sf_Trig<<endl;
    //cout<<"sf_MuTk"<<sf_MuTk<<" sf_MuID:"<<sf_MuID<<" sf_MuIso:"<<sf_MuIso<<" sf_ElReco:"<<sf_ElReco<<" sf_ElID:"<<sf_ElID<<" sf_BTag:"<<sf_BTag<<endl;
  }
  if(IsDATA && FakeRun && EventCand){
    w_FR      = CalcTestFakeWeight(muonLooseColl, electronLooseColl,
                  "TopHN17T", "TopHN17L", "TopHN17"+IDSSLabel+"T", "TopHN17"+IDSSLabel+"L", bjetColl.size(), 0);
  }
  weight *= w_TopPtRW * w_Prefire * sf_Trig * w_FR;
  weight *= sf_MuTk * sf_MuID * sf_MuIso * sf_ElReco * sf_ElID * sf_BTag;

 
  if(OS2l){
    if(MCCFRate){ MeasMCCFRate(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                               jetColl, bjetColl, vMET_xyCorr, weight, ""); }
    if(CFMCClos){ CheckMCCFClosure(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                                 jetColl, bjetColl, vMET_xyCorr, weight, ""); }
    if(MDists){ GetMassDists(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                             jetColl, bjetColl, vMET_xyCorr, weight, ""); }
  }
  if(SS2l){
    if(CFlip){ CheckChargeFlip(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                               jetColl, bjetColl, vMET_xyCorr, weight, ""); }
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


void MeasCFlipRate::GetMassDists(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                 vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                 vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  bool IsSS2l=false;
  double Mll=-1., PT1=-1., Eta1=999., PT2=-1., Eta2=999., fEta1=999., fEta2=999.;
  if( FakeRun      and weight==0.  ) return; 
  if( !( NElT==2 and NMuT==0 ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NElT==2){
    if( ElTColl.at(0).Charge()==ElTColl.at(1).Charge()  ) IsSS2l=true;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    Mll = (ElTColl.at(0)+ElTColl.at(1)).M();
    PT1 = ElTColl.at(0).Pt(), Eta1 = ElTColl.at(0).Eta(), fEta1=fabs(Eta1);
    PT2 = ElTColl.at(1).Pt(), Eta2 = ElTColl.at(1).Eta(), fEta2=fabs(Eta2);
    if( !(Mll>60. && Mll<120.) ) return;
  }
  else return;

  if(IsSS2l) Label="_SS"+Label;
  else       Label="_OS"+Label;

  int IdxFlipped=-1;
  float CFSF1=1., CFSF2=1., CFSF3=1.;
  if( IsSS2l && (!IsDATA) ){
    vector<Gen> TruthColl = GetGens();
    int NFk=0, NFlip=0, NCv=0;
    for(unsigned int im=0; im<MuTColl.size(); im++){
      int LepType = GetLeptonType_JH(MuTColl.at(im), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
    }
    for(unsigned int ie=0; ie<ElTColl.size(); ie++){
      int LepType = GetLeptonType_JH(ElTColl.at(ie), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0){
        int Idx_Closest    = GenMatchedIdx(ElTColl.at(ie),TruthColl); 
        int IdxType_NearEl = LepType>3? GetPrElType_InSameSCRange(Idx_Closest, TruthColl, "IdxType"):Idx_Closest;
        int Idx_NearEl     = LepType>3? IdxType_NearEl/10:Idx_Closest;
        if(ElTColl.at(ie).Charge()*TruthColl.at(Idx_NearEl).PID()>0){
          NFlip++; IdxFlipped=ie;
//          CFSF1 = GetCFRSF(ElTColl.at(ie), "App2Bin1_It1")*GetCFRSF(ElTColl.at(ie), "App2Bin1_It2")*GetCFRSF(ElTColl.at(ie), "App2Bin1_It3")*GetCFRSF(ElTColl.at(ie), "App2Bin1_It4");
//          CFSF2 = GetCFRSF(ElTColl.at(ie), "App2Bin2_It1")*GetCFRSF(ElTColl.at(ie), "App2Bin2_It2")*GetCFRSF(ElTColl.at(ie), "App2Bin2_It3")*GetCFRSF(ElTColl.at(ie), "App2Bin2_It4");
//          CFSF3 = GetCFRSF(ElTColl.at(ie), "App2Bin3_It1")*GetCFRSF(ElTColl.at(ie), "App2Bin3_It2")*GetCFRSF(ElTColl.at(ie), "App2Bin3_It3")*GetCFRSF(ElTColl.at(ie), "App2Bin3_It4");
        }
      } 
    }
    if( !(NFk==0 && NFlip>0) ) return;
  }


  //Approach0
  int EtaRegIndex0=0, BinIndex0=0;
  if(fEta1<1.47){ EtaRegIndex0 = fEta2<1.47? 1:2;}
  else{ EtaRegIndex0 = fEta2<1.47? -2:3; }
  BinIndex0 = abs(EtaRegIndex0);
  TString BinIdxStr0 = TString::Itoa(BinIndex0,10);
  FillHist("Mll_Bin0_"+BinIdxStr0+Label, Mll, weight, 60, 60., 120.);
  if(IsSS2l && !IsDATA){
    FillHist("NCnt_Bin0"+Label, BinIndex0, weight, 4, 0., 4.);
    if((EtaRegIndex0>0 && IdxFlipped==1) or (EtaRegIndex0<0 && IdxFlipped==0) ) FillHist("NCntCF_Bin0"+Label, BinIndex0, weight, 4, 0., 4.);
  }

  //Approach1
  int BinIndex1=0, BinIndex2=0;
  int EtaRegIndex1=0;
  if     (fEta1<0.8 ){ EtaRegIndex1=fEta2<0.8?  1:fEta2<1.47?  2:fEta2<2.?  3: 4; }
  else if(fEta1<1.47){ EtaRegIndex1=fEta2<0.8? -2:fEta2<1.47?  5:fEta2<2.?  6: 7; }
  else if(fEta1<2.  ){ EtaRegIndex1=fEta2<0.8? -3:fEta2<1.47? -6:fEta2<2.?  8: 9; }
  else               { EtaRegIndex1=fEta2<0.8? -4:fEta2<1.47? -7:fEta2<2.? -9:10; }

  if     (EtaRegIndex1==1){ BinIndex1=1; BinIndex2=1; }
  else if(EtaRegIndex1==5){ BinIndex1=2; BinIndex2=2; }
  else if(abs(EtaRegIndex1)==3){
    if     (EtaRegIndex1>0) BinIndex1= PT2<35.? 3:PT2<50.? 4: 5;
    else if(EtaRegIndex1<0) BinIndex1= PT1<35.? 3:PT1<50.? 4: 5;
    BinIndex2=3;
  }
  else if(abs(EtaRegIndex1)==4){
    if     (EtaRegIndex1>0) BinIndex1= PT2<35.? 6:PT2<50.? 7: 8;
    else if(EtaRegIndex1<0) BinIndex1= PT1<35.? 6:PT1<50.? 7: 8;
    BinIndex2=4;
  }
  TString BinIdxStr1 = TString::Itoa(BinIndex1,10);
  TString BinIdxStr2 = TString::Itoa(BinIndex2,10);
  FillHist("Mll_Bin1_"+BinIdxStr1+Label, Mll, weight, 60, 60., 120.);
  FillHist("Mll_Bin2_"+BinIdxStr2+Label, Mll, weight, 60, 60., 120.);

  //Approach2
  for(unsigned int ie=0; ie<ElTColl.size(); ie++){
    float PT=ElTColl.at(ie).Pt(), fEta=fabs(ElTColl.at(ie).Eta());
    int TmpEtaRegIndex1 = fEta<1.47? 0:1;
    int TmpEtaRegIndex2 = fEta<0.8? 0:fEta<1.47? 1: fEta<2.? 2:3;
    int TmpPtRegIndex1  = PT<35? 0:PT<50? 1: 2;
    int TmpPtRegIndex2  = PT<35? 0:PT<50? 1:PT<100? 2:3;
    int TmpBinIndex1 = TmpEtaRegIndex1*3+TmpPtRegIndex1+1; 
    int TmpBinIndex2 = TmpEtaRegIndex2+1; 
    int TmpBinIndex3 = TmpPtRegIndex2 +1; 
    TString TmpBinIdxStr1 = TString::Itoa(TmpBinIndex1,10);
    TString TmpBinIdxStr2 = TString::Itoa(TmpBinIndex2,10);
    TString TmpBinIdxStr3 = TString::Itoa(TmpBinIndex3,10);
    FillHist("Mll_App2Bin1_"+TmpBinIdxStr1+Label, Mll, weight*CFSF1, 60, 60., 120.);
    FillHist("Mll_App2Bin2_"+TmpBinIdxStr2+Label, Mll, weight*CFSF2, 60, 60., 120.);
    FillHist("Mll_App2Bin3_"+TmpBinIdxStr3+Label, Mll, weight*CFSF3, 60, 60., 120.);
    if(!IsDATA){
    //if(IsSS2l && !IsDATA){
      FillHist("NCnt_App2Bin1"+Label, TmpBinIndex1, weight*CFSF1, 7, 0., 7.);
      if(IdxFlipped==(int) ie) FillHist("NCntCF_App2Bin1"+Label, TmpBinIndex1, weight*CFSF1, 7, 0., 7.);
      FillHist("NCnt_App2Bin2"+Label, TmpBinIndex2, weight*CFSF2, 5, 0., 5.);
      if(IdxFlipped==(int) ie) FillHist("NCntCF_App2Bin2"+Label, TmpBinIndex2, weight*CFSF2, 5, 0., 5.);
      FillHist("NCnt_App2Bin3"+Label, TmpBinIndex3, weight*CFSF3, 5, 0., 5.);
      if(IdxFlipped==(int) ie) FillHist("NCntCF_App2Bin3"+Label, TmpBinIndex3, weight*CFSF3, 5, 0., 5.);
    }
  }

}


void MeasCFlipRate::CheckMCCFClosure(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                     vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                     vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)

{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  bool IsSS2l=false, IsOnZ=false;
  double Mll=-1., PT1=-1., Eta1=999., PT2=-1., Eta2=999., dRll=-1.;
  if( FakeRun      and weight==0.  ) return; 
  if( !( (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NElT==2){
    if( ElTColl.at(0).Charge()==ElTColl.at(1).Charge()  ) IsSS2l=true;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    Mll  = (ElTColl.at(0)+ElTColl.at(1)).M();
    dRll = ElTColl.at(0).DeltaR(ElTColl.at(1));
    PT1 = ElTColl.at(0).Pt(), Eta1 = ElTColl.at(0).Eta();
    PT2 = ElTColl.at(1).Pt(), Eta2 = ElTColl.at(1).Eta();
  }
  else return;
//  else if(NMuT==1 && NElT==1){
//    if( ElTColl.at(0).Charge()==MuTColl.at(0).Charge()  ) IsSS2l=true;
//    if(!(MuTColl.at(0).Pt()>10 && ElTColl.at(0).Pt()>15)) return; 
//    if(!(MuTColl.at(0).Pt()>25 || ElTColl.at(0).Pt()>25)) return; 
//    PT1 = ElTColl.at(0).Pt(), Eta1 = ElTColl.at(0).Eta();
//  }

  if(IsSS2l) Label="_SS"+Label;
  else       Label="_OS"+Label;
  if(fabs(Mll-91.2)<15.) IsOnZ=true;

  int IdxFlipped=-1;
  float CFSF1=1., CFSF2=1., CFSF3=1.;
  if( IsSS2l && (!IsDATA) ){
    vector<Gen> TruthColl = GetGens();
    int NFk=0, NFlip=0, NCv=0;
    for(unsigned int im=0; im<MuTColl.size(); im++){
      int LepType = GetLeptonType_JH(MuTColl.at(im), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
    }
    for(unsigned int ie=0; ie<ElTColl.size(); ie++){
      int LepType = GetLeptonType_JH(ElTColl.at(ie), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0){
        int Idx_Closest    = GenMatchedIdx(ElTColl.at(ie),TruthColl); 
        int IdxType_NearEl = LepType>3? GetPrElType_InSameSCRange(Idx_Closest, TruthColl, "IdxType"):Idx_Closest;
        int Idx_NearEl     = LepType>3? IdxType_NearEl/10:Idx_Closest;
        if(ElTColl.at(ie).Charge()*TruthColl.at(Idx_NearEl).PID()>0){
          NFlip++; IdxFlipped=ie; 
          CFSF1 = GetCFRSF(ElTColl.at(ie), "App2Bin1_Fin");
          CFSF2 = GetCFRSF(ElTColl.at(ie), "App2Bin2_Fin");
          CFSF3 = GetCFRSF(ElTColl.at(ie), "App2Bin3_Fin");
        }
      } 
    }
    if( !(NFk==0 && NFlip>0) ) return;
  }


//  float CFR1=0., CFR2=0., CFR3=0.;
  float CFR1_1=0., CFR2_1=0., CFR1_2=0., CFR2_2=0., CFR1_3=0., CFR2_3=0.;
  float PT1_1=0., PT2_1=0., PT1_2=0., PT2_2=0., PT1_3=0., PT2_3=0.;
//  for(unsigned int ie=0; ie<ElTColl.size(); ie++){
//    if(IsSS2l) break;
//    TString TypeStr=IsDATA? "Data":"MC";
//    CFR1 += GetCFRSF(ElTColl.at(ie), "App2Bin1_Fin", TypeStr+"Eff");
//    CFR2 += GetCFRSF(ElTColl.at(ie), "App2Bin2_Fin", TypeStr+"Eff");
//    CFR3 += GetCFRSF(ElTColl.at(ie), "App2Bin3_Fin", TypeStr+"Eff");
//  }
  if(!IsSS2l){
    TString TypeStr=IsDATA? "Data":"MC";
    CFR1_1=GetCFRSF(ElTColl.at(0), "App2Bin1_Fin", TypeStr+"Eff");
    CFR1_2=GetCFRSF(ElTColl.at(0), "App2Bin2_Fin", TypeStr+"Eff");
    CFR1_3=GetCFRSF(ElTColl.at(0), "App2Bin3_Fin", TypeStr+"Eff");
    CFR2_1=GetCFRSF(ElTColl.at(1), "App2Bin1_Fin", TypeStr+"Eff");
    CFR2_2=GetCFRSF(ElTColl.at(1), "App2Bin2_Fin", TypeStr+"Eff");
    CFR2_3=GetCFRSF(ElTColl.at(1), "App2Bin3_Fin", TypeStr+"Eff");
    PT1_1=GetFlipCorrPT(ElTColl.at(0), "App2Bin1_Fin", TypeStr+"ScaleSmear");
    PT1_2=GetFlipCorrPT(ElTColl.at(0), "App2Bin2_Fin", TypeStr+"ScaleSmear");
    PT1_3=GetFlipCorrPT(ElTColl.at(0), "App2Bin3_Fin", TypeStr+"ScaleSmear");
    PT2_1=GetFlipCorrPT(ElTColl.at(1), "App2Bin1_Fin", TypeStr+"ScaleSmear");
    PT2_2=GetFlipCorrPT(ElTColl.at(1), "App2Bin2_Fin", TypeStr+"ScaleSmear");
    PT2_3=GetFlipCorrPT(ElTColl.at(1), "App2Bin3_Fin", TypeStr+"ScaleSmear");
  }


  //Exp from CFR
  if((!IsSS2l) && (!FakeRun)){
    for(int it=1; it<=3; it++){
      for(unsigned int ie=0; ie<ElTColl.size(); ie++){
        TString TagStr="", TypeStr=IsDATA? "_ExpRD":"_ExpMC"; float Prob=1., CorrPT1=0., CorrPT2=0.;
        if     (it==1){ TagStr="_App2Bin1"; Prob = ie==0? CFR1_1:CFR2_1; CorrPT1=PT1_1; CorrPT2=PT2_1;}
        else if(it==2){ TagStr="_App2Bin2"; Prob = ie==0? CFR1_2:CFR2_2; CorrPT1=PT1_2; CorrPT2=PT2_2;}
        else if(it==3){ TagStr="_App2Bin3"; Prob = ie==0? CFR1_3:CFR2_3; CorrPT1=PT1_3; CorrPT2=PT2_3;}
        if(IsDATA){
          TLorentzVector e1Corr, e2Corr;
          if(ie==0){
            e1Corr.SetPtEtaPhiM(CorrPT1, ElTColl.at(0).Eta(), ElTColl.at(0).Phi(), ElTColl.at(0).M());
            e2Corr.SetPtEtaPhiM(ElTColl.at(1).Pt(), ElTColl.at(1).Eta(), ElTColl.at(1).Phi(), ElTColl.at(1).M());
          }
          else{
            e1Corr.SetPtEtaPhiM(ElTColl.at(0).Pt(), ElTColl.at(0).Eta(), ElTColl.at(0).Phi(), ElTColl.at(0).M());
            e2Corr.SetPtEtaPhiM(CorrPT2, ElTColl.at(1).Eta(), ElTColl.at(1).Phi(), ElTColl.at(1).M());
          }
          vector<TLorentzVector> NewElTColl = {e1Corr, e2Corr};
          sort(NewElTColl.begin(), NewElTColl.end(), PtComparing);
          PT1=e1Corr.Pt(), PT2=e2Corr.Pt(), Mll=(e1Corr+e2Corr).M(); IsOnZ=fabs(Mll-91.2)<15.;
//          if(it==2){
//            printf("PTe1:%.3e PTe2:%.3e\n", ElTColl.at(0).Pt(), ElTColl.at(1).Pt());
//            printf("PTe1N:%.3e, PTe1:%.3e\n", NewElTColl.at(0).Pt(), ElTColl.at(0).Pt());
//            printf("PTe2N:%.3e, PTe2:%.3e\n", NewElTColl.at(1).Pt(), ElTColl.at(1).Pt());
//          }
        }

        //printf("%d %d %.2e %.2e\n", it, ie, PT1, Mll);

        FillHist("PTl1"+Label+TagStr+TypeStr, PT1, weight*Prob, 30, 0., 300.); 
        FillHist("PTl2"+Label+TagStr+TypeStr, PT2, weight*Prob, 20, 0., 200.); 
        FillHist("Etal1"+Label+TagStr+TypeStr, Eta1, weight*Prob, 20, -5., 5.); 
        FillHist("Etal2"+Label+TagStr+TypeStr, Eta2, weight*Prob, 20, -5., 5.); 
        FillHist("Mll"+Label+TagStr+TypeStr, Mll, weight*Prob, 30, 0., 300.);
        FillHist("NCnt"+Label+TagStr+TypeStr, 0., weight*Prob, 1, 0., 1.);
        if(IsOnZ){
          FillHist("PTl1_OnZ"+Label+TagStr+TypeStr, PT1, weight*Prob, 30, 0., 300.); 
          FillHist("PTl2_OnZ"+Label+TagStr+TypeStr, PT2, weight*Prob, 20, 0., 200.); 
          FillHist("Etal1_OnZ"+Label+TagStr+TypeStr, Eta1, weight*Prob, 20, -5., 5.); 
          FillHist("Etal2_OnZ"+Label+TagStr+TypeStr, Eta2, weight*Prob, 20, -5., 5.); 
          FillHist("Mll_OnZ"+Label+TagStr+TypeStr, Mll, weight*Prob, 30, 60., 120.);
          FillHist("dRll_OnZ"+Label+TagStr+TypeStr, dRll, weight*Prob, 25, 0., 5.);
          FillHist("NCnt_OnZ"+Label+TagStr+TypeStr, 0., weight*Prob, 1, 0., 1.);
        }
//      }
        if(!IsDATA){
      //for(unsigned int ie=0; ie<ElTColl.size() && (!IsDATA); ie++){
        if(!(Mll>60 && Mll<120)) continue;
        float PT=ElTColl.at(ie).Pt(), fEta = fabs(ElTColl.at(ie).Eta());
        int EtaRegIndex1 = fEta<1.47? 0:1;
        int EtaRegIndex2 = fEta<0.8? 0:fEta<1.47? 1: fEta<2.? 2:3;
        int PtRegIndex1  = PT<35? 0:PT<50? 1: 2;
        int PtRegIndex2  = PT<35? 0:PT<50? 1:PT<100? 2:3;
        int BinIndex1 = EtaRegIndex1*3+PtRegIndex1+1; 
        int BinIndex2 = EtaRegIndex2+1; 
        int BinIndex3 = PtRegIndex2 +1; 
        int BinIndex = it==1? BinIndex1:it==2? BinIndex2:it==3? BinIndex3:0;

        float TmpCFR1=0., TmpCFR2=0., TmpCFR3=0.;
        TmpCFR1 = GetCFRSF(ElTColl.at(ie), "App2Bin1_Fin", "MCEff");
        TmpCFR2 = GetCFRSF(ElTColl.at(ie), "App2Bin2_Fin", "MCEff");
        TmpCFR3 = GetCFRSF(ElTColl.at(ie), "App2Bin3_Fin", "MCEff");
        float TmpProb = it==1? TmpCFR1:it==2? TmpCFR2:it==3? TmpCFR3:0.;
        FillHist("BinCnt"+Label+TagStr+"_MC", BinIndex, weight, 7, 0., 7.);
        FillHist("BinCnt"+Label+TagStr+TypeStr, BinIndex, weight*Prob, 7, 0., 7.);
        FillHist("BinCntCF"+Label+TagStr+TypeStr, BinIndex, weight*TmpProb, 7, 0., 7.);
        }
      }
    }
  }
  else{
    //Obs & Exp from SF*MC
    for(int im=1; im<=2; im++){
      for(int it=1; it<=3; it++){
        if(IsDATA && im!=1) break;
        if(im==1 && it>1) continue;
        TString TypeStr=IsDATA? "_ObsRD":im==1? "_ObsMC":"_ExpSFMC";
        TString TagStr=""; float SF=1.;
        if(im==2 && (!IsDATA)){
          if     (it==1){ TagStr="_App2Bin1"; SF = CFSF1; }
          else if(it==2){ TagStr="_App2Bin2"; SF = CFSF2; }
          else if(it==3){ TagStr="_App2Bin3"; SF = CFSF3; }
        }
        FillHist("PTl1"+Label+TagStr+TypeStr, PT1, weight*SF, 30, 0., 300.); 
        FillHist("PTl2"+Label+TagStr+TypeStr, PT2, weight*SF, 20, 0., 200.); 
        FillHist("Etal1"+Label+TagStr+TypeStr, Eta1, weight*SF, 20, -5., 5.); 
        FillHist("Etal2"+Label+TagStr+TypeStr, Eta2, weight*SF, 20, -5., 5.); 
        FillHist("Mll"+Label+TagStr+TypeStr, Mll, weight*SF, 30, 0., 300.);
        FillHist("NCnt"+Label+TagStr+TypeStr, 0., weight*SF, 1, 0., 1.);
        if(IsOnZ){
          FillHist("PTl1_OnZ"+Label+TagStr+TypeStr, PT1, weight*SF, 30, 0., 300.); 
          FillHist("PTl2_OnZ"+Label+TagStr+TypeStr, PT2, weight*SF, 20, 0., 200.); 
          FillHist("Etal1_OnZ"+Label+TagStr+TypeStr, Eta1, weight*SF, 20, -5., 5.); 
          FillHist("Etal2_OnZ"+Label+TagStr+TypeStr, Eta2, weight*SF, 20, -5., 5.); 
          FillHist("dRll_OnZ"+Label+TagStr+TypeStr, dRll, weight*SF, 25, 0., 5.);
          FillHist("Mll_OnZ"+Label+TagStr+TypeStr, Mll, weight*SF, 30, 60., 120.);
          FillHist("NCnt_OnZ"+Label+TagStr+TypeStr, 0., weight*SF, 1, 0., 1.);
        }
      }
    }

    for(int it=1; it<=3 && (!IsDATA); it++){
      for(unsigned int ie=0; ie<ElTColl.size(); ie++){
        if(!(Mll>60 && Mll<120)) continue;
        TString TagStr=it==1? "_App2Bin1":it==2? "_App2Bin2":it==3? "_App2Bin3":"";
        TString TypeStr="_ObsMC";
        float PT=ElTColl.at(ie).Pt(), fEta = fabs(ElTColl.at(ie).Eta());
        int EtaRegIndex1 = fEta<1.47? 0:1;
        int EtaRegIndex2 = fEta<0.8? 0:fEta<1.47? 1: fEta<2.? 2:3;
        int PtRegIndex1  = PT<35? 0:PT<50? 1: 2;
        int PtRegIndex2  = PT<35? 0:PT<50? 1:PT<100? 2:3;
        int BinIndex1 = EtaRegIndex1*3+PtRegIndex1+1; 
        int BinIndex2 = EtaRegIndex2+1; 
        int BinIndex3 = PtRegIndex2 +1; 
        int BinIndex = it==1? BinIndex1:it==2? BinIndex2:it==3? BinIndex3:0;
        FillHist("BinCnt"+Label+TagStr+TypeStr, BinIndex, weight, 7, 0., 7.);
        if(IdxFlipped==(int) ie) FillHist("BinCntCF"+Label+TagStr+TypeStr, BinIndex, weight, 7, 0., 7.);
      }
    }
  }//SS2l ends

}


float MeasCFlipRate::GetTmpCFRate(float PT, float fEta){

  float CFR=0.;
  if(fEta<0) fEta=fabs(fEta);
  if     (fEta<0.8 ) CFR=1.74e-05;
  else if(fEta<1.47) CFR=9.38e-05;
  else if(fEta<2.  ) CFR=PT<35.? 8.13e-04:PT<50.? 6.61e-04:1.14e-03;
  else               CFR=PT<35.? 1.01e-03:PT<50.? 1.06e-03:1.75e-03;

  return CFR;
}


void MeasCFlipRate::MeasMCCFRate(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                 vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                 vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  //int Nj=JetColl.size(), Nb=BJetColl.size();
  bool IsSS2l=false, IsOnZ=false;
  double Mll=-1., PT1=-1., Eta1=999., PT2=-1., Eta2=999., fEta1=999., fEta2=999.;
  if( FakeRun      and weight==0.  ) return; 
  if( !( (NMuT==1 and NElT==1) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NMuT==2){
    Label="_2M0E"+Label;
    if( MuTColl.at(0).Charge()==MuTColl.at(1).Charge()  ) IsSS2l=true;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(DataYear>2016 && Mll<4) return; 
  }
  else if(NElT==2){
    //Label="_0M2E"+Label;
    if( ElTColl.at(0).Charge()==ElTColl.at(1).Charge()  ) IsSS2l=true;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    Mll = (ElTColl.at(0)+ElTColl.at(1)).M();
    PT1 = ElTColl.at(0).Pt(), Eta1 = ElTColl.at(0).Eta(), fEta1=fabs(Eta1);
    PT2 = ElTColl.at(1).Pt(), Eta2 = ElTColl.at(1).Eta(), fEta2=fabs(Eta2);
  }
  else if(NMuT==1 && NElT==1){
    //Label="_1M1E"+Label;
    if( ElTColl.at(0).Charge()==MuTColl.at(0).Charge()  ) IsSS2l=true;
    if(!(MuTColl.at(0).Pt()>10 && ElTColl.at(0).Pt()>15)) return; 
    if(!(MuTColl.at(0).Pt()>25 || ElTColl.at(0).Pt()>25)) return; 
    PT1 = ElTColl.at(0).Pt(), Eta1 = ElTColl.at(0).Eta(), fEta1=fabs(Eta1);
  }

  if(IsSS2l) Label="_SS"+Label;
  else       Label="_OS"+Label;
  if(fabs(Mll-91.2)<15.) IsOnZ=true;

  int IdxFlipped=-1;
  if( IsSS2l && (!IsDATA) ){
    vector<Gen> TruthColl = GetGens();
    int NFk=0, NFlip=0, NCv=0;
    for(unsigned int im=0; im<MuTColl.size(); im++){
      int LepType = GetLeptonType_JH(MuTColl.at(im), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
    }
    for(unsigned int ie=0; ie<ElTColl.size(); ie++){
      int LepType = GetLeptonType_JH(ElTColl.at(ie), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0){
        int Idx_Closest    = GenMatchedIdx(ElTColl.at(ie),TruthColl); 
        int IdxType_NearEl = LepType>3? GetPrElType_InSameSCRange(Idx_Closest, TruthColl, "IdxType"):Idx_Closest;
        int Idx_NearEl     = LepType>3? IdxType_NearEl/10:Idx_Closest;
        if(ElTColl.at(ie).Charge()*TruthColl.at(Idx_NearEl).PID()>0){ NFlip++; IdxFlipped=ie; }
      } 
    }
    if( !(NFk==0 && NFlip>0) ) return;
  }

  const int NPtBinEdges=7, NEtaBinEdges=9, NfEtaBinEdges=5, NEtaBinEdges2=11, NfEtaBinEdges2=6;
  double PtBinEdges[NPtBinEdges] = {0., 15., 25., 35., 50., 100., 200.};
  double EtaBinEdges[NEtaBinEdges] = {-2.5, -2., -1.47, -0.8, 0., 0.8, 1.47, 2., 2.5};
  double fEtaBinEdges[NfEtaBinEdges] = {0., 0.8, 1.47, 2., 2.5};
  double EtaBinEdges2[NEtaBinEdges2] = {-2.5, -2., -1.47, -1., -0.5, 0., 0.5, 1., 1.47, 2., 2.5};
  double fEtaBinEdges2[NfEtaBinEdges2] = {0., 0.5, 1., 1.47, 2., 2.5};

  double PTcf=IsSS2l? ElTColl.at(IdxFlipped).Pt():-1., Etacf=IsSS2l? ElTColl.at(IdxFlipped).Eta():999., fEtacf=fabs(Etacf);
  PT1=min(max(PT1,PtBinEdges[1]),PtBinEdges[NPtBinEdges-1]); 
  PT2=min(max(PT2,PtBinEdges[1]),PtBinEdges[NPtBinEdges-1]); 
  PTcf=min(max(PTcf,PtBinEdges[1]),PtBinEdges[NPtBinEdges-1]); 

  //1:3-eta, 2:4-eta, 3: 4-eta+2-pt 4: 3-eta2
  int BinIndex1=0, BinIndex2=0, BinIndex3=0, BinIndex4=0, BinIndex5=0, BinIndex6=0, BinIndex7=0;
  int EtaRegIndex1=0;
  if(NElT==2){
    if     (fEta1<0.8 ){ EtaRegIndex1=fEta2<0.8?  1:fEta2<1.47?  2:fEta2<2.?  3: 4; }
    else if(fEta1<1.47){ EtaRegIndex1=fEta2<0.8? -2:fEta2<1.47?  5:fEta2<2.?  6: 7; }
    else if(fEta1<2.  ){ EtaRegIndex1=fEta2<0.8? -3:fEta2<1.47? -6:fEta2<2.?  8: 9; }
    else               { EtaRegIndex1=fEta2<0.8? -4:fEta2<1.47? -7:fEta2<2.? -9:10; }

    if     (fEta1<0.8 ){ BinIndex1=fEta2<0.8? 1:fEta2<1.47? 2:fEta2<2.? 3:4; }
    else if(fEta1<1.47){ BinIndex1=fEta2<0.8? 2:fEta2<1.47? 5:fEta2<2.? 6:7; }
    else if(fEta1<2.  ){ BinIndex1=fEta2<0.8? 3:fEta2<1.47? 6:fEta2<2.? 8:9; }
    else               { BinIndex1=fEta2<0.8? 4:fEta2<1.47? 7:fEta2<2.? 9:10; }
    if     (fEta1<0.8 ){ BinIndex2=fEta2<0.8? 1:fEta2<1.47? 2:3; }
    else if(fEta1<1.47){ BinIndex2=fEta2<0.8? 2:fEta2<1.47? 4:5; }
    else               { BinIndex2=fEta2<0.8? 3:fEta2<1.47? 5:6; }
    if     (fEta1<1.47){ BinIndex3=fEta2<1.47? 1:fEta2<2.? 2:3; }
    else if(fEta1<2.  ){ BinIndex3=fEta2<1.47? 2:fEta2<2.? 4:5; }
    else               { BinIndex3=fEta2<1.47? 3:fEta2<2.? 5:6; }

    if     (EtaRegIndex1==1){ BinIndex7=1; }
    else if(EtaRegIndex1==5){ BinIndex7=2; }
    else if(abs(EtaRegIndex1)==3){
      if     (EtaRegIndex1>0) BinIndex7= PT2<35.? 3:PT2<50.? 4: 5;
      else if(EtaRegIndex1<0) BinIndex7= PT1<35.? 3:PT1<50.? 4: 5;
    }
    else if(abs(EtaRegIndex1)==4){
      if     (EtaRegIndex1>0) BinIndex7= PT2<35.? 6:PT2<50.? 7: 8;
      else if(EtaRegIndex1<0) BinIndex7= PT1<35.? 6:PT1<50.? 7: 8;
    }


    if(fEta1<1.47 && PT1<50){
      if(fEta2<1.47) BinIndex4=PT2<50? 1:2;
      else           BinIndex4=PT2<50? 3:4;
    }
    else if(fEta1<1.47 && PT1>=50){
      if(fEta2<1.47) BinIndex4=PT2<50? 2:5;
      else           BinIndex4=PT2<50? 6:7;
    }
    else if(PT1<50.){
      if(fEta2<1.47) BinIndex4=PT2<50? 3:6;
      else           BinIndex4=PT2<50? 8:9;
    }
    else{
      if(fEta2<1.47) BinIndex4=PT2<50? 4:7;
      else           BinIndex4=PT2<50? 9:10;
    }


    if(fEta1<1.47){
      if(fEta2<1.47) BinIndex5=1;
      else           BinIndex5=PT2<50.? 2:3;
    }
    else if(PT1<50.){
      if(fEta2<1.47) BinIndex5=2;
      else           BinIndex5=PT2<50.? 4:5;
    }
    else{
      if(fEta2<1.47) BinIndex5=3;
      else           BinIndex5=PT2<50.? 5:6;
    }


    if(fEta1<1.47){
      if(fEta2<1.47) BinIndex6=1;
      else           BinIndex6=PT2<35? 2:PT2<50? 3:4;
    }
    else if(PT1<35.){
      if(fEta2<1.47) BinIndex6=2;
      else           BinIndex6=PT2<35? 5:PT2<50? 6:7;
    }
    else if(PT1<50.){
      if(fEta2<1.47) BinIndex6=3;
      else           BinIndex6=PT2<35? 6:PT2<50? 8:9;
    }
    else{
      if(fEta2<1.47) BinIndex6=4;
      else           BinIndex6=PT2<35? 7:PT2<50? 9:10;
    }
  }
  TString EtaRegStr = TString::Itoa(abs(EtaRegIndex1),10);

  if(IsSS2l){
    FillHist("IdxFlipped"+Label, IdxFlipped, weight, 2, 0., 2.);
    FillHist("h1D_PT"+Label, PTcf, weight, NPtBinEdges-1, PtBinEdges);
    FillHist("h1D_Eta"+Label, Etacf, weight, NEtaBinEdges-1, EtaBinEdges);
    FillHist("h2D_PTEta"+Label, PTcf, Etacf, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges-1, EtaBinEdges);
    FillHist("h2D_PTfEta"+Label, PTcf, fEtacf, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
    FillHist("h1D_Eta2"+Label, Etacf, weight, NEtaBinEdges2-1, EtaBinEdges2);
    FillHist("h2D_PTEta2"+Label, PTcf, Etacf, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges2-1, EtaBinEdges2);
    FillHist("h2D_PTfEta2"+Label, PTcf, fEtacf, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges2-1, fEtaBinEdges2);
  }
  else{
    FillHist("h1D_PT"+Label, PT1, weight, NPtBinEdges-1, PtBinEdges);
    FillHist("h1D_Eta"+Label, Eta1, weight, NEtaBinEdges-1, EtaBinEdges);
    FillHist("h2D_PTEta"+Label, PT1, Eta1, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges-1, EtaBinEdges);
    FillHist("h2D_PTfEta"+Label, PT1, fEta1, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
    FillHist("h1D_Eta2"+Label, Eta1, weight, NEtaBinEdges2-1, EtaBinEdges2);
    FillHist("h2D_PTEta2"+Label, PT1, Eta1, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges2-1, EtaBinEdges2);
    FillHist("h2D_PTfEta2"+Label, PT1, fEta1, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges2-1, fEtaBinEdges2);
    if(NElT==2){
      FillHist("h1D_PT"+Label, PT2, weight, NPtBinEdges-1, PtBinEdges);
      FillHist("h1D_Eta"+Label, Eta2, weight, NEtaBinEdges-1, EtaBinEdges);
      FillHist("h2D_PTEta"+Label, PT2, Eta2, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges-1, EtaBinEdges);
      FillHist("h2D_PTfEta"+Label, PT2, fEta2, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
      FillHist("h1D_Eta2"+Label, Eta2, weight, NEtaBinEdges2-1, EtaBinEdges2);
      FillHist("h2D_PTEta2"+Label, PT2, Eta2, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges2-1, EtaBinEdges2);
      FillHist("h2D_PTfEta2"+Label, PT2, fEta2, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges2-1, fEtaBinEdges2);
    }
  }
  if(IsOnZ && NElT==2){
    Label="_OnZ"+Label;
    if(IsSS2l){
      FillHist("IdxFlipped"+Label, IdxFlipped, weight, 2, 0., 2.);
      FillHist("h1D_PT"+Label, PTcf, weight, NPtBinEdges-1, PtBinEdges);
      FillHist("h1D_Eta"+Label, Etacf, weight, NEtaBinEdges-1, EtaBinEdges);
      FillHist("h2D_PTEta"+Label, PTcf, Etacf, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges-1, EtaBinEdges);
      FillHist("h2D_PTfEta"+Label, PTcf, fEtacf, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
      FillHist("h1D_Eta2"+Label, Etacf, weight, NEtaBinEdges2-1, EtaBinEdges2);
      FillHist("h2D_PTEta2"+Label, PTcf, Etacf, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges2-1, EtaBinEdges2);
      FillHist("h2D_PTfEta2"+Label, PTcf, fEtacf, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges2-1, fEtaBinEdges2);
      FillHist("NCount_Div1"+Label, BinIndex1, weight, 11, 0., 11.);
      FillHist("NCount_Div2"+Label, BinIndex2, weight, 7, 0., 7.);
      FillHist("NCount_Div3"+Label, BinIndex3, weight, 7, 0., 7.);
      FillHist("NCount_Div4"+Label, BinIndex4, weight, 11, 0., 11.);
      FillHist("NCount_Div5"+Label, BinIndex5, weight, 7, 0., 7.);
      FillHist("NCount_Div6"+Label, BinIndex6, weight, 11, 0., 11.);
      FillHist("NCount_Div7"+Label, BinIndex7, weight, 9, 0., 9.);
      FillHist("h2D_PT1PT2_BinIdx"+EtaRegStr+Label, PT1, PT2, weight, NPtBinEdges-1, PtBinEdges, NPtBinEdges-1, PtBinEdges);
      if(EtaRegIndex1>0) FillHist("h2D_PT1PT2_EtaReg"+EtaRegStr+Label, PT1, PT2, weight, NPtBinEdges-1, PtBinEdges, NPtBinEdges-1, PtBinEdges);
      if(EtaRegIndex1<0) FillHist("h2D_PT1PT2_EtaReg"+EtaRegStr+Label, PT2, PT1, weight, NPtBinEdges-1, PtBinEdges, NPtBinEdges-1, PtBinEdges);
    }
    else{
      FillHist("h1D_PT"+Label, PT1, weight, NPtBinEdges-1, PtBinEdges);
      FillHist("h1D_Eta"+Label, Eta1, weight, NEtaBinEdges-1, EtaBinEdges);
      FillHist("h2D_PTEta"+Label, PT1, Eta1, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges-1, EtaBinEdges);
      FillHist("h2D_PTfEta"+Label, PT1, fEta1, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
      FillHist("h1D_Eta2"+Label, Eta1, weight, NEtaBinEdges2-1, EtaBinEdges2);
      FillHist("h2D_PTEta2"+Label, PT1, Eta1, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges2-1, EtaBinEdges2);
      FillHist("h2D_PTfEta2"+Label, PT1, fEta1, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges2-1, fEtaBinEdges2);
      if(NElT==2){
        FillHist("h1D_PT"+Label, PT2, weight, NPtBinEdges-1, PtBinEdges);
        FillHist("h1D_Eta"+Label, Eta2, weight, NEtaBinEdges-1, EtaBinEdges);
        FillHist("h2D_PTEta"+Label, PT2, Eta2, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges-1, EtaBinEdges);
        FillHist("h2D_PTfEta"+Label, PT2, fEta2, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges-1, fEtaBinEdges);
        FillHist("h1D_Eta2"+Label, Eta2, weight, NEtaBinEdges2-1, EtaBinEdges2);
        FillHist("h2D_PTEta2"+Label, PT2, Eta2, weight, NPtBinEdges-1, PtBinEdges, NEtaBinEdges2-1, EtaBinEdges2);
        FillHist("h2D_PTfEta2"+Label, PT2, fEta2, weight, NPtBinEdges-1, PtBinEdges, NfEtaBinEdges2-1, fEtaBinEdges2);
        FillHist("NCount_Div1"+Label, BinIndex1, weight, 11, 0., 11.);
        FillHist("NCount_Div2"+Label, BinIndex2, weight, 7, 0., 7.);
        FillHist("NCount_Div3"+Label, BinIndex3, weight, 7, 0., 7.);
        FillHist("NCount_Div4"+Label, BinIndex4, weight, 11, 0., 11.);
        FillHist("NCount_Div5"+Label, BinIndex5, weight, 7, 0., 7.);
        FillHist("NCount_Div6"+Label, BinIndex6, weight, 11, 0., 11.);
        FillHist("NCount_Div7"+Label, BinIndex7, weight, 9, 0., 9.);
        FillHist("h2D_PT1PT2_BinIdx"+EtaRegStr+Label, PT1, PT2, weight, NPtBinEdges-1, PtBinEdges, NPtBinEdges-1, PtBinEdges);
        if(EtaRegIndex1>0) FillHist("h2D_PT1PT2_EtaReg"+EtaRegStr+Label, PT1, PT2, weight, NPtBinEdges-1, PtBinEdges, NPtBinEdges-1, PtBinEdges);
        if(EtaRegIndex1<0) FillHist("h2D_PT1PT2_EtaReg"+EtaRegStr+Label, PT2, PT1, weight, NPtBinEdges-1, PtBinEdges, NPtBinEdges-1, PtBinEdges);
      }
    }
  }

}


void MeasCFlipRate::CheckChargeFlip(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
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
  }
  else if(NElT==2){
    Label="_2E"+Label;
    if(ElTColl.at(0).Charge()!=ElTColl.at(1).Charge()) return;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
  }

  int Nj=JetColl.size(), Nb=BJetColl.size();
  FillHist("h2D_Nj_Nb"+Label, Nj, Nb, weight, 10, 0., 10., 5, 0., 5.);
  if( !(Nb==0 && Nj<3) ) return;
   
  float Mll=-1., PTl1=-1, Etal1=999., PTl2=-1., Etal2=-1., MTW=-1., MTllv=-1.;
  float MET=-1.;
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
  MET = vMET.Pt();
  
  if( (ConvRun or FlipRun) && (!IsDATA) ){
    vector<Gen> TruthColl = GetGens();
    int NFk=0, NFlip=0, NCv=0;
    for(unsigned int im=0; im<MuTColl.size(); im++){
      int LepType = GetLeptonType_JH(MuTColl.at(im), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
    }
    for(unsigned int ie=0; ie<ElTColl.size(); ie++){
      int LepType = GetLeptonType_JH(ElTColl.at(ie), TruthColl);
      if(LepType<0 && LepType>-5) NFk++;
      else if(LepType<-4) NCv++;
      else if(LepType>0){
        int Idx_Closest    = GenMatchedIdx(ElTColl.at(ie),TruthColl); 
        int IdxType_NearEl = LepType>3? GetPrElType_InSameSCRange(Idx_Closest, TruthColl, "IdxType"):Idx_Closest;
        int Idx_NearEl     = LepType>3? IdxType_NearEl/10:Idx_Closest;
        if(ElTColl.at(ie).Charge()*TruthColl.at(Idx_NearEl).PID()>0) NFlip++;
      } 
    }
    if( ConvRun && (!(NFk==0 && NFlip==0)) ) return;
    if( FlipRun && (!(NFk==0 && NFlip>0)) ) return;
  }

  FillHist("Mll"+Label, Mll, weight, 40, 0., 200.);
  FillHist("PTl1"+Label, PTl1, weight, 20, 0., 200.);
  FillHist("PTl2"+Label, PTl2, weight, 20, 0., 200.);
  FillHist("Etal1"+Label, Etal1, weight, 20, -5., 5.);
  FillHist("Etal2"+Label, Etal2, weight, 20, -5., 5.);
  FillHist("MET"+Label, MET, weight, 20, 0., 200.);
  FillHist("MTW"+Label, MTW, weight, 20, 0., 200.);
  FillHist("MTllv"+Label, MTllv, weight, 20, 0., 200.);
  if(fabs(Mll-91.2)<15){
    Label="_OnZ"+Label;
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
    FillHist("Mll"+Label, Mll, weight, 30, 60., 120.);
    FillHist("PTl1"+Label, PTl1, weight, 20, 0., 200.);
    FillHist("PTl2"+Label, PTl2, weight, 20, 0., 200.);
    FillHist("Etal1"+Label, Etal1, weight, 20, -5., 5.);
    FillHist("Etal2"+Label, Etal2, weight, 20, -5., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
  }

}


float MeasCFlipRate::CalcTestFakeWeight(vector<Muon>& MuColl, vector<Electron>& ElColl, TString MuIDT, TString MuIDL, TString ElIDT, TString ElIDL, int NBJet, int SystDir){

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


float MeasCFlipRate::GetTestElFR(Electron& El, TString Key, int SystDir){

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


float MeasCFlipRate::GetTestMuFR(Muon& Mu, TString Key, int SystDir){

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



float MeasCFlipRate::GetFlipCorrPT(Electron& El, TString Tag, TString Option){

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
  else if(Tag.Contains("App2Bin2_")){
    int Idx = BinIndex2;
//    PTScaleCorr = Idx==1? 0.978:Idx==2? 0.982:Idx==3? 0.981:Idx==4? 0.988:0.;
//    PTResCorr   = Idx==1?  1.62:Idx==2?  1.30:Idx==3?  1.14:Idx==4?  1.04:0.;
//    RelRes      = Idx==1? 0.0265:Idx==2? 0.032:Idx==3? 0.0375:Idx==4? 0.0397:0.;
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


float MeasCFlipRate::GetCFRSF(Electron& El, TString Tag, TString Option){

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



MeasCFlipRate::MeasCFlipRate(){

}


MeasCFlipRate::~MeasCFlipRate(){

}

