#include "FakeRateMeas.h"

void FakeRateMeas::initializeAnalyzer(){

  ElFR=false, MuFR=false, MeasFR=false, MeasPU=false, PrVal=false, METMTWCut=false;
  SS2l=false, FkVal=false,
  SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("ElFR"))      ElFR      = true; 
    if(Userflags.at(i).Contains("MuFR"))      MuFR      = true; 
    if(Userflags.at(i).Contains("MeasFR"))    MeasFR    = true; 
    if(Userflags.at(i).Contains("MeasPU"))    MeasPU    = true; 
    if(Userflags.at(i).Contains("METMTWCut")) METMTWCut = true; 
    if(Userflags.at(i).Contains("PrVal"))     PrVal     = true; 
    if(Userflags.at(i).Contains("FkVal"))     FkVal     = true; 
    if(Userflags.at(i).Contains("SystRun"))   SystRun   = true; 
  }

  if(DataYear==2016){
  }
  if(DataYear==2017){
    TrigList_MuFR.push_back("HLT_Mu17_TrkIsoVVL_v");
    TrigList_MuFR.push_back("HLT_Mu8_TrkIsoVVL_v");

    TrigList_ElFR.push_back("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
    TrigList_ElFR.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
    TrigList_ElFR.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  }
  if(DataYear==2018){
  }

  //Set up the tagger map only for the param settings to be used.
  std::vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void FakeRateMeas::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;
  if(!IsDATA){
    weight*=ev.MCweight()*weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    weight*=GetBRWeight();
    weight*=GetPileUpWeight(nPileUp, 0);
  }
  FillHist("CutFlow", 0., weight, 20, 0., 20.);

  bool PassTrig=false;
  if     (MuFR){ PassTrig = ev.PassTrigger(TrigList_MuFR); }
  else if(ElFR){ PassTrig = ev.PassTrigger(TrigList_ElFR); }
  if(!PassTrig) return;
  FillHist("CutFlow", 1., weight, 20, 0., 20.);
  if(!PassMETFilter()) return;
  FillHist("CutFlow", 2., weight, 20, 0., 20.);

  bool PreCutPass=false;
  vector<Muon>     muonPreColl     = GetMuons("NOCUT", 5., 2.4);
  vector<Electron> electronPreColl = GetElectrons("NOCUT", 5., 2.5);
  sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);
  if     (MuFR and muonPreColl.size()>0    ) PreCutPass=true;
  else if(ElFR and electronPreColl.size()>0) PreCutPass=true;
  if(!PreCutPass) return;


  TString IDSSLabel = "SS";
  vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TopHN17T", 10., 2.4);
  vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"T", 10., 2.5);
  vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"L", 10., 2.5);
  vector<Electron> electronVetoColl  = SelectElectrons(electronPreColl, "TopHN17L", 10., 2.5);


  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  vector<Jet> jetPreColl = GetAllJets();
  sort(jetPreColl.begin(), jetPreColl.end(), PtComparing);
  vector<Jet> jetColl  = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVeto");
  vector<Jet> bjetColl = SelBJets(jetColl, param_jets);


  Particle vMET = ev.GetMETVector();
  Particle vMET_T1xy(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);


  std::vector<Gen> truthColl;


  bool EventCand = false;
  if     (MuFR){ EventCand = muonLooseColl.size()==1; }
  else if(ElFR){ EventCand = electronLooseColl.size()==1; }

  float w_topptrw = 1., w_prefire = 1., sf_trig = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
  float w_Norm_Mu17=1., w_Norm_Mu8=1., w_Norm_El23=1., w_Norm_El12=1.;
  if((!IsDATA) and EventCand){
    //if(MCSample.Contains("TT") and MCSample.Contains("powheg")) truthColl = GetGens();
    //w_topptrw = mcCorr->GetTopPtReweight(truthColl);
    //sf_muid   = GetMuonSF(muonTightColl, "TopHNID_TkMu", "ID");
    //sf_muiso  = GetMuonSF(muonTightColl, "POGTIso_POGTID", "Iso");
    sf_elreco = GetElectronSF(electronLooseColl, "", "Reco");
    //sf_elid   = GetElectronSF(electronTightColl, "TopHNID"+IDSSLabel, "ID");
    sf_btag   = mcCorr->GetBTaggingReweight_1a(jetColl, param_jets);
    TString SFKey_Trig_All = muonTightColl.size()==2? "DiMuIso_HNTopID":electronTightColl.size()==2? "DiElIso_HNTopIDSS":"EMuIso_HNTopIDSS"; 
    //sf_trig   = SS2l? mcCorr->GetTriggerSF(electronTightColl, muonTightColl, SFKey_Trig_All, ""):1.;
    w_prefire = GetPrefireWeight(0);
    w_Norm_Mu17 = GetResidualNormSF("TrigMu17"), w_Norm_Mu8 = GetResidualNormSF("TrigMu8");
    w_Norm_El23 = GetResidualNormSF("TrigEl23"), w_Norm_El12 = GetResidualNormSF("TrigEl12");
    //cout<<"w_gen:"<<w_gen<<" w_lumi:"<<w_lumi<<" w_PU:"<<w_PU<<" w_prefire:"<<w_prefire<<" sf_trig:"<<sf_trig<<endl;
    //cout<<"sf_mutk"<<sf_mutk<<" sf_muid:"<<sf_muid<<" sf_muiso:"<<sf_muiso<<" sf_elreco:"<<sf_elreco<<" sf_elid:"<<sf_elid<<" sf_btag:"<<sf_btag<<endl;
  }
  weight *= w_topptrw * w_prefire * sf_trig;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

 
  if(MuFR){
    if(MeasPU){
      MeasPileUp(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigMu17");
      MeasPileUp(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigMu8");
    }
    if(PrVal){
      CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigMu17");
      CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigMu8");
    }
    if(MeasFR){
      MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu17, "TrigMu17_Cent");
      MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu8, "TrigMu8_Cent");
    }
    if(METMTWCut){
      OptimizeMETMTWCuts(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                         jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigMu8");
    }
  }
  if(ElFR){
    if(MeasPU){
      MeasPileUp(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigEl23_Pt15");
      MeasPileUp(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigEl12_Pt15");
    }
    if(PrVal){
      CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigEl23_Pt15");
      CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                 jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigEl12_Pt15");
    }
    if(MeasFR){
      MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El23, "TrigEl23_Pt15_Cent");
      MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El12, "TrigEl12_Pt15_Cent");
    }
    if(METMTWCut){
      OptimizeMETMTWCuts(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                         jetColl, bjetColl, vMET_T1xy, ev, weight, "TrigEl12_Pt15");
    }
  }
  if(SystRun && ((!IsDATA) or MeasFR)){
    float nvtx_reweight_Mu17 = NvtxReweight("TrigMu17");
    float nvtx_reweight_Mu8  = NvtxReweight("TrigMu8");
    float nvtx_reweight_El23 = NvtxReweight("TrigEl23_Pt15");
    float nvtx_reweight_El12 = NvtxReweight("TrigEl12_Pt15");
    float w_NormUp_Mu17=1.1, w_NormDown_Mu17=0.9, w_NormUp_Mu8=1.1, w_NormDown_Mu8=0.9;
    float w_NormUp_El23=1.1, w_NormDown_El23=0.9, w_NormUp_El12=1.1, w_NormDown_El12=0.9;
    if(IsDATA){ w_NormUp_Mu17=1., w_NormDown_Mu17=1., w_NormUp_Mu8=1., w_NormDown_Mu8=1.;
                w_NormUp_El23=1., w_NormDown_El23=1., w_NormUp_El12=1., w_NormDown_El12=1.; }
    float w_PUUp   = weight!=0.? GetPileUpWeight(nPileUp,  1)/GetPileUpWeight(nPileUp, 0):0.;
    float w_PUDown = weight!=0.? GetPileUpWeight(nPileUp, -1)/GetPileUpWeight(nPileUp, 0):0.;
    vector<Jet> jet20Coll       = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 20., 2.4, "LVeto");
    vector<Jet> jetJESUpColl    = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJESUp");
    vector<Jet> jetJESDownColl  = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJESDown");
    vector<Jet> jetJERUpColl    = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJERUp");
    vector<Jet> jetJERDownColl  = SelectJets(jetPreColl, muonLooseColl, electronVetoColl, "tight", 25., 2.4, "LVetoSystJERDown");
    vector<Jet> bjetJESUpColl   = SelBJets(jetJESUpColl, param_jets);
    vector<Jet> bjetJESDownColl = SelBJets(jetJESDownColl, param_jets);
    vector<Jet> bjetJERUpColl   = SelBJets(jetJERUpColl, param_jets);
    vector<Jet> bjetJERDownColl = SelBJets(jetJERDownColl, param_jets);
    Particle vMET_T1xy_JESUp    (pfMET_Type1_PhiCor_pt_shifts->at(2)*TMath::Cos(pfMET_Type1_PhiCor_phi_shifts->at(2)), pfMET_Type1_PhiCor_pt_shifts->at(2)*TMath::Sin(pfMET_Type1_PhiCor_phi_shifts->at(2)), 0., pfMET_Type1_PhiCor_pt_shifts->at(2));
    Particle vMET_T1xy_JESDown  (pfMET_Type1_PhiCor_pt_shifts->at(3)*TMath::Cos(pfMET_Type1_PhiCor_phi_shifts->at(3)), pfMET_Type1_PhiCor_pt_shifts->at(3)*TMath::Sin(pfMET_Type1_PhiCor_phi_shifts->at(3)), 0., pfMET_Type1_PhiCor_pt_shifts->at(3));
    Particle vMET_T1xy_JERUp    (pfMET_Type1_PhiCor_pt_shifts->at(0)*TMath::Cos(pfMET_Type1_PhiCor_phi_shifts->at(0)), pfMET_Type1_PhiCor_pt_shifts->at(0)*TMath::Sin(pfMET_Type1_PhiCor_phi_shifts->at(0)), 0., pfMET_Type1_PhiCor_pt_shifts->at(0));
    Particle vMET_T1xy_JERDown  (pfMET_Type1_PhiCor_pt_shifts->at(1)*TMath::Cos(pfMET_Type1_PhiCor_phi_shifts->at(1)), pfMET_Type1_PhiCor_pt_shifts->at(1)*TMath::Sin(pfMET_Type1_PhiCor_phi_shifts->at(1)), 0., pfMET_Type1_PhiCor_pt_shifts->at(1));
    Particle vMET_T1xy_UnclUp   (pfMET_Type1_PhiCor_pt_shifts->at(10)*TMath::Cos(pfMET_Type1_PhiCor_phi_shifts->at(10)), pfMET_Type1_PhiCor_pt_shifts->at(10)*TMath::Sin(pfMET_Type1_PhiCor_phi_shifts->at(10)), 0., pfMET_Type1_PhiCor_pt_shifts->at(10));
    Particle vMET_T1xy_UnclDown (pfMET_Type1_PhiCor_pt_shifts->at(11)*TMath::Cos(pfMET_Type1_PhiCor_phi_shifts->at(11)), pfMET_Type1_PhiCor_pt_shifts->at(11)*TMath::Sin(pfMET_Type1_PhiCor_phi_shifts->at(11)), 0., pfMET_Type1_PhiCor_pt_shifts->at(11));

    if(MuFR){
      if(MeasFR){
        //TrigMu17
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu17*w_NormUp_Mu17, "TrigMu17_PrUp");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu17*w_NormDown_Mu17, "TrigMu17_PrDown");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu17, "TrigMu17_HasB");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jet20Coll, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu17, "TrigMu17_JetPt20");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu17, "TrigMu17_JetPt60");


        //TrigMu8
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu8*w_NormUp_Mu8, "TrigMu8_PrUp");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu8*w_NormDown_Mu8, "TrigMu8_PrDown");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu8, "TrigMu8_HasB");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jet20Coll, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu8, "TrigMu8_JetPt20");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_Mu8, "TrigMu8_JetPt60");
      }
      if(PrVal){
        //TrigMu17
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_Mu17, "TrigMu17_SystUp_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_Mu17, "TrigMu17_SystDown_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUUp, "TrigMu17_SystUp_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUDown, "TrigMu17_SystDown_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, ev, weight, "TrigMu17_SystUp_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, ev, weight, "TrigMu17_SystDown_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, ev, weight, "TrigMu17_SystUp_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, ev, weight, "TrigMu17_SystDown_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclUp, ev, weight, "TrigMu17_SystUp_Uncl");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclDown, ev, weight, "TrigMu17_SystDown_Uncl");

        //TrigMu8
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_Mu8, "TrigMu8_SystUp_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_Mu8, "TrigMu8_SystDown_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUUp, "TrigMu8_SystUp_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUDown, "TrigMu8_SystDown_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, ev, weight, "TrigMu8_SystUp_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, ev, weight, "TrigMu8_SystDown_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, ev, weight, "TrigMu8_SystUp_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, ev, weight, "TrigMu8_SystDown_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclUp, ev, weight, "TrigMu8_SystUp_Uncl");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclDown, ev, weight, "TrigMu8_SystDown_Uncl");
      }//End of PrVal
    }//End of MuFR
    if(ElFR){
      if(MeasFR){
        //TrigEl23_Pt15
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El23*w_NormUp_El23, "TrigEl23_Pt15_PrUp");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El23*w_NormDown_El23, "TrigEl23_Pt15_PrDown");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El23, "TrigEl23_Pt15_HasB");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El23, "TrigEl23_Pt15_JetPt30");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El23, "TrigEl23_Pt15_JetPt60");


        //TrigEl12_Pt15
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El12*w_NormUp_El12, "TrigEl12_Pt15_PrUp");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El12*w_NormDown_El12, "TrigEl12_Pt15_PrDown");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El12, "TrigEl12_Pt15_HasB");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El12, "TrigEl12_Pt15_JetPt30");
        MeasFakeRate(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_Norm_El12, "TrigEl12_Pt15_JetPt60");
      }
      if(PrVal){
        //TrigEl23_Pt15
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_El23, "TrigEl23_Pt15_SystUp_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_El23, "TrigEl23_Pt15_SystDown_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUUp, "TrigEl23_Pt15_SystUp_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUDown, "TrigEl23_Pt15_SystDown_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, ev, weight, "TrigEl23_Pt15_SystUp_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, ev, weight, "TrigEl23_Pt15_SystDown_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, ev, weight, "TrigEl23_Pt15_SystUp_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, ev, weight, "TrigEl23_Pt15_SystDown_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclUp, ev, weight, "TrigEl23_Pt15_SystUp_Uncl");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclDown, ev, weight, "TrigEl23_Pt15_SystDown_Uncl");

        //TrigEl12_Pt15
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_El12, "TrigEl12_Pt15_SystUp_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*nvtx_reweight_El12, "TrigEl12_Pt15_SystDown_Nvtx");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUUp, "TrigEl12_Pt15_SystUp_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy, ev, weight*w_PUDown, "TrigEl12_Pt15_SystDown_PU");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESUpColl, bjetJESUpColl, vMET_T1xy_JESUp, ev, weight, "TrigEl12_Pt15_SystUp_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJESDownColl, bjetJESDownColl, vMET_T1xy_JESDown, ev, weight, "TrigEl12_Pt15_SystDown_JES");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERUpColl, bjetJERUpColl, vMET_T1xy_JERUp, ev, weight, "TrigEl12_Pt15_SystUp_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetJERDownColl, bjetJERDownColl, vMET_T1xy_JERDown, ev, weight, "TrigEl12_Pt15_SystDown_JER");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclUp, ev, weight, "TrigEl12_Pt15_SystUp_Uncl");
        CheckPromptValid(muonTightColl, muonLooseColl, muonLooseColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_T1xy_UnclDown, ev, weight, "TrigEl12_Pt15_SystDown_Uncl");
      }//End of PrVal
    }//End of ElFR
  }//End of SystRun

}


void FakeRateMeas::CheckPromptValid(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                    vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                    vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuL==1 && NElL==0) or (NElL==1 && NMuL==0) or (NMuL==2 && NElL==0) or (NMuL==0 && NElL==2) ) ) return;
  if( !(NMuL==NMuV && NElL==NElV) ) return;
  if( !(NMuL==NMuT && NElL==NElT) ) return;
  if(NMuT==1){
    float PT = MuTColl.at(0).Pt(), MTW = MT(MuTColl.at(0),vMET), MET=vMET.Pt();
    bool TrigMu17 = Label.Contains("TrigMu17"), TrigMu8 = Label.Contains("TrigMu8"), TrigSel=false, PassJetReq=false;
    if     (TrigMu17 && Ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v") && PT>20) TrigSel=true;
    else if(TrigMu8  && Ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v")  && PT>10) TrigSel=true;
    else return;
    if(!TrigSel) return;
    Label="_"+Label;
    FillHist("Rho"+Label, Rho, weight, 35, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl"+Label, PT, weight, 40, 0., 200.);
    FillHist("Etal"+Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);

    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(MuTColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    Label="_1j"+Label;
    FillHist("Rho"+Label, Rho, weight, 35, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl"+Label, PT, weight, 40, 0., 200.);
    FillHist("Etal" +Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("Ptj1" +Label, JetColl.at(0).Pt() , weight, 20, 0., 200.);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("dRlj1"+Label, MuTColl.at(0).DeltaR(JetColl.at(0)), weight, 25, 0., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);

    if( !(MET>50) ) return;
    Label.ReplaceAll("_1j","_1jMET");
    FillHist("Rho"+Label, Rho, weight, 35, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl"+Label, PT, weight, 40, 0., 200.);
    FillHist("Etal" +Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("Ptj1" +Label, JetColl.at(0).Pt() , weight, 20, 0., 200.);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("dRlj1"+Label, MuTColl.at(0).DeltaR(JetColl.at(0)), weight, 25, 0., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);

    if( !(MTW>90) ) return;
    Label.ReplaceAll("_1jMET","_1jMETMTW");
    FillHist("Rho"+Label, Rho, weight, 35, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl"+Label, PT, weight, 40, 0., 200.);
    FillHist("Etal" +Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("Ptj1" +Label, JetColl.at(0).Pt() , weight, 20, 0., 200.);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("dRlj1"+Label, MuTColl.at(0).DeltaR(JetColl.at(0)), weight, 25, 0., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  }
  else if(NElT==1){
    float PT = ElTColl.at(0).Pt(), MTW = MT(ElTColl.at(0),vMET), MET=vMET.Pt();
    bool TrigEl23 = Label.Contains("TrigEl23"), TrigEl12 = Label.Contains("TrigEl12"), TrigEl8 = Label.Contains("TrigEl8");
    bool TrigSel=false, PassJetReq=false;
    if     (TrigEl23 && Ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>25) TrigSel=true;
    else if(TrigEl12 && Ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>15) TrigSel=true;
    else if(TrigEl8  && Ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")  && PT>10) TrigSel=true;
    if(!TrigSel) return;
    Label="_"+Label;
    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(ElTColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    Label="_1j"+Label;
    FillHist("Rho"+Label, Rho, weight, 35, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl"+Label, PT, weight, 40, 0., 200.);
    FillHist("Etal" +Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("Ptj1" +Label, JetColl.at(0).Pt() , weight, 20, 0., 200.);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("dRlj1"+Label, ElTColl.at(0).DeltaR(JetColl.at(0)), weight, 25, 0., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);

    if( !(MET>50) ) return;
    Label.ReplaceAll("_1j","_1jMET");
    FillHist("Rho"+Label, Rho, weight, 35, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl"+Label, PT, weight, 40, 0., 200.);
    FillHist("Etal" +Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("Ptj1" +Label, JetColl.at(0).Pt() , weight, 20, 0., 200.);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("dRlj1"+Label, ElTColl.at(0).DeltaR(JetColl.at(0)), weight, 25, 0., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);

    if( !(MTW>90) ) return;
    Label.ReplaceAll("_1jMET","_1jMETMTW");
    FillHist("Rho"+Label, Rho, weight, 35, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl"+Label, PT, weight, 40, 0., 200.);
    FillHist("Etal" +Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("Ptj1" +Label, JetColl.at(0).Pt() , weight, 20, 0., 200.);
    FillHist("Etaj1"+Label, JetColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("dRlj1"+Label, ElTColl.at(0).DeltaR(JetColl.at(0)), weight, 25, 0., 5.);
    FillHist("MET"+Label, MET, weight, 20, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  }
  else if(NMuT==2){
    float PT1 = MuTColl.at(0).Pt(), PT2 = MuTColl.at(1).Pt(); 
    bool TrigMu17 = Label.Contains("TrigMu17"), TrigMu8 = Label.Contains("TrigMu8"), TrigSel=false, PassJetReq=false;
    if     (TrigMu17 && Ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v") && PT1>20 && PT2>10) TrigSel=true;
    else if(TrigMu8  && Ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v")  && PT1>10 && PT2>10) TrigSel=true;
    else return;

    if(!TrigSel) return;
    if(MuTColl.at(0).Charge()==MuTColl.at(1).Charge()) return;

    float Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(fabs(Mll-91.2)>15.) return;
    Label="_2lMQ_"+Label;
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl1"+Label, PT1, weight, 24, 0., 120.);
    FillHist("Etal1"+Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("PTl2"+Label, PT2, weight, 24, 0., 120.);
    FillHist("Etal2"+Label, MuTColl.at(1).Eta(), weight, 20, -5., 5.);
    FillHist("Mll"+Label, Mll, weight, 30, 60., 120.); 
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);

    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(MuTColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    Label.ReplaceAll("_2lMQ","_1j2lMQ");
    FillHist("PTl1"+Label, PT1, weight, 24, 0., 120.);
    FillHist("Etal1"+Label, MuTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("PTl2"+Label, PT2, weight, 24, 0., 120.);
    FillHist("Etal2"+Label, MuTColl.at(1).Eta(), weight, 20, -5., 5.);
    FillHist("Mll"+Label, Mll, weight, 30, 60., 120.); 
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  }
  else if(NElT==2){
    float PT1 = ElTColl.at(0).Pt(), PT2 = ElTColl.at(1).Pt(); 
    bool TrigEl23 = Label.Contains("TrigEl23"), TrigEl12 = Label.Contains("TrigEl12"), TrigEl8 = Label.Contains("TrigEl8");
    bool TrigSel=false, PassJetReq=false;
    if     (TrigEl23 && Ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT1>25) TrigSel=true;
    else if(TrigEl12 && Ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT1>15) TrigSel=true;
    else if(TrigEl8  && Ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")  && PT1>10) TrigSel=true;
    if(!TrigSel) return;
    if(Label.Contains("Pt15") && PT2<15) return;
    else if(Label.Contains("PT10") && PT2<10) return;
    if(ElTColl.at(0).Charge()==ElTColl.at(1).Charge()) return;

    float Mll = (ElTColl.at(0)+ElTColl.at(1)).M();
    if(fabs(Mll-91.2)>15.) return;
    Label="_2lMQ_"+Label;
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);
    FillHist("PTl1"+Label, PT1, weight, 24, 0., 120.);
    FillHist("Etal1"+Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("PTl2"+Label, PT2, weight, 24, 0., 120.);
    FillHist("Etal2"+Label, ElTColl.at(1).Eta(), weight, 20, -5., 5.);
    FillHist("Mll"+Label, Mll, weight, 30, 60., 120.); 
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);

    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(ElTColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    Label.ReplaceAll("_2lMQ","_1j2lMQ");
    FillHist("PTl1"+Label, PT1, weight, 24, 0., 120.);
    FillHist("Etal1"+Label, ElTColl.at(0).Eta(), weight, 20, -5., 5.);
    FillHist("PTl2"+Label, PT2, weight, 24, 0., 120.);
    FillHist("Etal2"+Label, ElTColl.at(1).Eta(), weight, 20, -5., 5.);
    FillHist("Mll"+Label, Mll, weight, 30, 60., 120.); 
    FillHist("NCount"+Label, 0., weight, 1, 0., 1.);
  }

}


void FakeRateMeas::MeasPileUp(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                              vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                              vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuL==1 and NElL==0) or (NElL==1 and NMuL==0) ) ) return;
  if( !(NMuL==NMuV and NElL==NElV) ) return;
  if( !(NMuL==NMuT and NElL==NElT) ) return;
  if(NMuT==1){
    float PT = MuTColl.at(0).Pt(), MTW = MT(MuTColl.at(0),vMET), MET=vMET.Pt();
    bool TrigMu17 = Label.Contains("TrigMu17"), TrigMu8 = Label.Contains("TrigMu8"), TrigSel=false, PassJetReq=false;
    if     (TrigMu17 && Ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v") && PT>20) TrigSel=true;
    else if(TrigMu8  && Ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v")  && PT>10) TrigSel=true;
    else return;
    if(!TrigSel) return;
    FillHist("Rho"+Label, Rho, weight, 70, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);

    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(MuTColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    FillHist("Rho_1j"+Label, Rho, weight, 70, 0., 70.);
    FillHist("NPV_1j"+Label, nPV, weight, 80, 0., 80.);

    if( !(MET>50 && MTW>80) ) return;
    FillHist("Rho_1jMETMTW"+Label, Rho, weight, 70, 0., 70.);
    FillHist("NPV_1jMETMTW"+Label, nPV, weight, 80, 0., 80.);
  }
  else if(NElT==1){
    float PT = ElTColl.at(0).Pt(), MTW = MT(ElTColl.at(0),vMET), MET=vMET.Pt();
    bool TrigEl23 = Label.Contains("TrigEl23"), TrigEl12 = Label.Contains("TrigEl12"), TrigEl8 = Label.Contains("TrigEl8");
    bool TrigSel=false, PassJetReq=false;
    if     (TrigEl23 && Ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>25) TrigSel=true;
    else if(TrigEl12 && Ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>15) TrigSel=true;
    else if(TrigEl8  && Ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")  && PT>10) TrigSel=true;
    if(!TrigSel) return;
    FillHist("Rho"+Label, Rho, weight, 70, 0., 70.);
    FillHist("NPV"+Label, nPV, weight, 80, 0., 80.);

    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(ElTColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    FillHist("Rho_1j"+Label, Rho, weight, 70, 0., 70.);
    FillHist("NPV_1j"+Label, nPV, weight, 80, 0., 80.);

    if( !(MET>50 && MTW>80) ) return;
    FillHist("Rho_1jMETMTW"+Label, Rho, weight, 70, 0., 70.);
    FillHist("NPV_1jMETMTW"+Label, nPV, weight, 80, 0., 80.);
  }
}



void FakeRateMeas::OptimizeMETMTWCuts(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                      vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                      vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label)
{
  int NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  //int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuL==1 and NElL==0) or (NElL==1 and NMuL==0) ) ) return;
  if( !(NMuL==NMuV and NElL==NElV) ) return;
  if(NMuL==1){
    float PT = MuLColl.at(0).Pt(), MTW = MT(MuLColl.at(0),vMET), MET=vMET.Pt();
    bool TrigMu8 = Label.Contains("TrigMu8"), TrigSel=false, PassJetReq=false;
    if(TrigMu8 && Ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v") && PT>10) TrigSel=true;
    else return;
    if(!TrigSel) return;

    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(MuLColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    Label="_"+Label;

    FillHist("SumW"+Label, 0., weight, 1, 0., 1.);
    FillHist("MET"+Label, MET, weight, 40, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);
    for(float METCut=0.; METCut<100.; METCut+=5.){
      for(float MTWCut=0.; MTWCut<100.; MTWCut+=5.){
        if(MET<METCut && MTW<MTWCut){ FillHist("SumWInMETMTW"+Label, METCut-1E-4, MTWCut-1E-4, weight, 20, 0., 100., 20, 0., 100.); }
      }
    }
  }
  else if(NElL==1){
    float PT = ElLColl.at(0).Pt(), MTW = MT(ElLColl.at(0),vMET), MET=vMET.Pt();
    bool TrigEl12 = Label.Contains("TrigEl12"), TrigEl8 = Label.Contains("TrigEl8");
    bool PTCut15 = Label.Contains("Pt15"), PTCut10 = Label.Contains("Pt10");
    bool TrigSel=false, PassJetReq=false;
    if     (TrigEl12 && Ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) TrigSel=true;
    else if(TrigEl8  && Ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) TrigSel=true;
    if(!TrigSel) return;
    if(PTCut15 && PT<15) return;
    else if(PTCut10 && PT<10) return; 

    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<40) continue;
      if(ElLColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    Label="_"+Label;

    FillHist("SumW"+Label, 0., weight, 1, 0., 1.);
    FillHist("MET"+Label, MET, weight, 40, 0., 200.);
    FillHist("MTW"+Label, MTW, weight, 40, 0., 200.);
    for(float METCut=0.; METCut<100.; METCut+=5.){
      for(float MTWCut=0.; MTWCut<100.; MTWCut+=5.){
        if(MET<METCut && MTW<MTWCut){ FillHist("SumWInMETMTW"+Label, METCut-1E-4, MTWCut-1E-4, weight, 20, 0., 100., 20, 0., 100.); }
      }
    }
  }
}


void FakeRateMeas::MeasFakeRate(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, Event& Ev, float weight, TString Label)
{
  //vector<Gen> TruthColl; if(CheckComposition){ TruthColl = GetGens(); }
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuL==1 and NElL==0) or (NElL==1 and NMuL==0) ) ) return;
  if( !(NMuL==NMuV and NElL==NElV) ) return;
  if(NMuL==1){
    float TightIso=0.1, RelIso=MuLColl.at(0).MiniRelIso(), MET=vMET.Pt(), MTW=MT(MuLColl.at(0),vMET);
    float PT=MuLColl.at(0).Pt(), PTCorr=MuLColl.at(0).CalcPtCone(RelIso, TightIso), Eta=MuLColl.at(0).Eta(), fEta=fabs(Eta);
    const int NPTEdges1D   = 10; double PTEdges1D[NPTEdges1D]     = {0., 10., 15., 20., 25., 30., 40., 50., 70., 100.};
    const int NPTEdges2D   = 7 ; double PTEdges2D[NPTEdges2D]     = {0., 10., 15., 20., 30., 50., 100.};
    const int NEtaEdges1D  = 7 ; double EtaEdges1D[NEtaEdges1D]   = {-2.4, -1.6, -0.9, 0., 0.9, 1.6, 2.4};
    const int NfEtaEdges2D = 4 ; double fEtaEdges2D[NfEtaEdges2D] = {0., 0.9, 1.6, 2.4};
    bool TrigMu17 = Label.Contains("TrigMu17"), TrigMu8 = Label.Contains("TrigMu8"), TrigSel=false, PassJetReq=false;
    if     (TrigMu17 && Ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v") && PT>20. && PTCorr>=30.){ TrigSel=true; Label.ReplaceAll("TrigMu17",""); }
    else if(TrigMu8  && Ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v")  && PT>10. && PTCorr <30.){ TrigSel=true; Label.ReplaceAll("TrigMu8", ""); }
    if(!TrigSel) return;

    float JetPtCut = Label.Contains("JetPt20")? 20.: Label.Contains("JetPt60")? 60.: 40.;
    bool PassNBCut = Label.Contains("HasB")? BJetColl.size()>0:true;
    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<JetPtCut) continue;
      if(MuLColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    if(!PassNBCut ) return;
    if( !(MET<25 && MTW<25) ) return;

    FillHist("SumWL_PT1D"+Label, PTCorr, weight, NPTEdges1D-1, PTEdges1D);
    FillHist("SumWL_Eta1D"+Label, Eta, weight, NEtaEdges1D-1, EtaEdges1D);
    FillHist("SumWL_PTfEta2D"+Label, PTCorr, fEta, weight, NPTEdges2D-1, PTEdges2D, NfEtaEdges2D-1, fEtaEdges2D);
    if(NMuT>0){
      FillHist("SumWT_PT1D"+Label, PTCorr, weight, NPTEdges1D-1, PTEdges1D);
      FillHist("SumWT_Eta1D"+Label, Eta, weight, NEtaEdges1D-1, EtaEdges1D);
      FillHist("SumWT_PTfEta2D"+Label, PTCorr, fEta, weight, NPTEdges2D-1, PTEdges2D, NfEtaEdges2D-1, fEtaEdges2D);
    }

    //For old codes
    FillHist("MuAllSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    if(NMuT>0) FillHist("MuAllIDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    if(fEta<0.9){
      FillHist("MuBSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
      if(NMuT>0) FillHist("MuBIDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    }
    else if(fEta<1.6){
      FillHist("MuBESumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
      if(NMuT>0) FillHist("MuBEIDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    }
    else{
      FillHist("MuESumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
      if(NMuT>0) FillHist("MuEIDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    }
  }
  else if(NElL==1){
    float TightIso=0.1, RelIso=ElLColl.at(0).MiniRelIso(), MET=vMET.Pt(), MTW=MT(ElLColl.at(0),vMET);
    float PT=ElLColl.at(0).Pt(), PTCorr=ElLColl.at(0).CalcPtCone(RelIso, TightIso), Eta=ElLColl.at(0).Eta(), fEta=fabs(Eta);
    const int NPTEdges1D   = 11; double PTEdges1D[NPTEdges1D]     = {0., 10., 15., 20., 25., 30., 35., 40., 50., 70., 100.};
    const int NPTEdges2D   = 8 ; double PTEdges2D[NPTEdges2D]     = {0., 10., 15., 20., 25., 35., 50., 100.};
    const int NEtaEdges1D  = 7 ; double EtaEdges1D[NEtaEdges1D]   = {-2.5, -1.5, -0.8, 0., 0., 1.5, 2.5};
    const int NfEtaEdges2D = 4 ; double fEtaEdges2D[NfEtaEdges2D] = {0., 0.8, 1.5, 2.5};
    bool TrigEl23 = Label.Contains("TrigEl23"), TrigEl12 = Label.Contains("TrigEl12"), TrigEl8 = Label.Contains("TrigEl8");
    bool PTCut10 = Label.Contains("_Pt10"), PTCut15 = Label.Contains("_Pt15"), TrigSel=false, PassJetReq=false;
    if(PTCut15){
      if     (TrigEl23 && Ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>25. && PTCorr>=35.){
        TrigSel=true; Label.ReplaceAll("TrigEl23_Pt15",""); }
      else if(TrigEl12 && Ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>15. && PTCorr <35.){
        TrigSel=true; Label.ReplaceAll("TrigEl12_Pt15",""); }
    }
    else if(PTCut10){
      if     (TrigEl23 && Ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>25. && PTCorr>=35.){
        TrigSel=true; Label.ReplaceAll("TrigEl23_Pt10",""); }
      else if(TrigEl12 && Ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && PT>15. && PTCorr<35. && PTCorr>=20.){
        TrigSel=true; Label.ReplaceAll("TrigEl12_Pt10",""); }
      else if(TrigEl8  && Ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")  && PT>10. && PTCorr<20.){
        TrigSel=true; Label.ReplaceAll("TrigEl8_Pt10",""); }
    }
    if(!TrigSel) return;

    float JetPtCut = Label.Contains("JetPt30")? 30.: Label.Contains("JetPt60")? 60.: 40.;
    bool PassNBCut = Label.Contains("HasB")? BJetColl.size()>0:true;
    for(unsigned int i=0; i<JetColl.size(); i++){
      if(JetColl.at(i).Pt()<JetPtCut) continue;
      if(ElLColl.at(0).DeltaR(JetColl.at(i))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    if(!PassNBCut ) return;
    if( !(MET<25 && MTW<25) ) return;

    FillHist("SumWL_PT1D"+Label, PTCorr, weight, NPTEdges1D-1, PTEdges1D);
    FillHist("SumWL_Eta1D"+Label, Eta, weight, NEtaEdges1D-1, EtaEdges1D);
    FillHist("SumWL_PTfEta2D"+Label, PTCorr, fEta, weight, NPTEdges2D-1, PTEdges2D, NfEtaEdges2D-1, fEtaEdges2D);
    if(NElT>0){
      FillHist("SumWT_PT1D"+Label, PTCorr, weight, NPTEdges1D-1, PTEdges1D);
      FillHist("SumWT_Eta1D"+Label, Eta, weight, NEtaEdges1D-1, EtaEdges1D);
      FillHist("SumWT_PTfEta2D"+Label, PTCorr, fEta, weight, NPTEdges2D-1, PTEdges2D, NfEtaEdges2D-1, fEtaEdges2D);
    }

    //Old Codes
    FillHist("EleAllSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    if(NElT>0) FillHist("EleAllIDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    if(fEta<0.8){
      FillHist("EleB1SumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
      if(NElT>0) FillHist("EleB1IDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    }
    else if(fEta<1.479){
      FillHist("EleB2SumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
      if(NElT>0) FillHist("EleB2IDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    }
    else{
      FillHist("EleESumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
      if(NElT>0) FillHist("EleEIDSumW_PT_FR1D"+Label, PTCorr, weight, NPTEdges2D-1, PTEdges2D);
    }
  }
}



void FakeRateMeas::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

FakeRateMeas::FakeRateMeas(){

}

FakeRateMeas::~FakeRateMeas(){

}


float FakeRateMeas::GetResidualNormSF(TString Key){

  float SF=1.;
  if     (Key=="TrigMu17") SF=1.586e-03;
  else if(Key=="TrigMu8" ) SF=9.514e-05;
  else if(Key=="TrigEl23") SF=8.397e-04;
  else if(Key=="TrigEl12") SF=6.043e-04;

  return SF;
}


float FakeRateMeas::NvtxReweight(TString Key){
  
  if(IsDATA) return 1.;

  float SF=1.;
  if(Key=="TrigMu17"){
    if     (nPV==1) SF=3.29;
    else if(nPV==2) SF=2.96;
    else if(nPV==3) SF=1.95;
    else if(nPV==4) SF=1.80;
    else if(nPV==5) SF=1.40;
    else if(nPV==6) SF=1.45;
    else if(nPV==7) SF=1.43;
    else if(nPV==8) SF=1.42;
    else if(nPV==9) SF=1.29;
    else if(nPV==10) SF=1.26;
    else if(nPV==11) SF=1.22;
    else if(nPV==12) SF=1.19;
    else if(nPV==13) SF=1.16;
    else if(nPV==14) SF=1.14;
    else if(nPV==15) SF=1.12;
    else if(nPV==16) SF=1.12;
    else if(nPV==17) SF=1.10;
    else if(nPV==18) SF=1.09;
    else if(nPV==19) SF=1.05;
    else if(nPV==20) SF=1.05;
    else if(nPV==21) SF=1.05;
    else if(nPV==22) SF=1.04;
    else if(nPV==23) SF=1.02;
    else if(nPV==24) SF=1.03;
    else if(nPV==25) SF=1.02;
    else if(nPV==26) SF=1.00;
    else if(nPV==27) SF=1.00;
    else if(nPV==28) SF=0.98;
    else if(nPV==29) SF=0.99;
    else if(nPV==30) SF=0.98;
    else if(nPV==31) SF=0.95;
    else if(nPV==32) SF=0.93;
    else if(nPV==33) SF=0.92;
    else if(nPV==34) SF=0.90;
    else if(nPV==35) SF=0.88;
    else if(nPV==36) SF=0.88;
    else if(nPV==37) SF=0.85;
    else if(nPV==38) SF=0.83;
    else if(nPV==39) SF=0.82;
    else if(nPV==40) SF=0.83;
    else if(nPV==41) SF=0.80;
    else if(nPV==42) SF=0.81;
    else if(nPV==43) SF=0.81;
    else if(nPV==44) SF=0.80;
    else if(nPV==45) SF=0.81;
    else if(nPV==46) SF=0.85;
    else if(nPV==47) SF=0.82;
    else if(nPV==48) SF=0.88;
    else if(nPV==49) SF=0.90;
    else if(nPV==50) SF=0.90;
    else if(nPV==51) SF=0.98;
    else if(nPV==52) SF=1.00;
    else if(nPV==53) SF=1.07;
    else if(nPV==54) SF=1.07;
    else if(nPV==55) SF=1.11;
    else if(nPV==56) SF=1.19;
    else if(nPV==57) SF=1.33;
    else if(nPV==58) SF=1.25;
    else if(nPV==59) SF=1.37;
    else if(nPV==60) SF=1.49;
    else if(nPV==61) SF=1.65;
    else if(nPV==62) SF=1.71;
    else if(nPV==63) SF=1.64;
    else if(nPV==64) SF=1.76;
    else if(nPV==65) SF=1.83;
    else if(nPV==66) SF=2.09;
    else if(nPV==67) SF=2.32;
    else if(nPV==68) SF=1.94;
    else if(nPV==69) SF=2.05;
    else if(nPV==70) SF=2.53;
    else if(nPV==71) SF=2.24;
    else if(nPV==72) SF=2.71;
    else if(nPV==73) SF=3.33;
    else if(nPV==74) SF=2.31;
    else if(nPV==75) SF=3.02;
    else if(nPV==76) SF=2.45;
    else if(nPV==77) SF=2.85;
    else if(nPV==78) SF=3.32;
    else if(nPV>=79) SF=3.69;
  }
  else if(Key=="TrigMu8"){
    if     (nPV==1) SF=2.97;
    else if(nPV==2) SF=3.21;
    else if(nPV==3) SF=2.48;
    else if(nPV==4) SF=2.08;
    else if(nPV==5) SF=1.63;
    else if(nPV==6) SF=1.72;
    else if(nPV==7) SF=1.46;
    else if(nPV==8) SF=1.58;
    else if(nPV==9) SF=1.46;
    else if(nPV==10) SF=1.47;
    else if(nPV==11) SF=1.37;
    else if(nPV==12) SF=1.34;
    else if(nPV==13) SF=1.30;
    else if(nPV==14) SF=1.27;
    else if(nPV==15) SF=1.21;
    else if(nPV==16) SF=1.19;
    else if(nPV==17) SF=1.16;
    else if(nPV==18) SF=1.14;
    else if(nPV==19) SF=1.11;
    else if(nPV==20) SF=1.08;
    else if(nPV==21) SF=1.07;
    else if(nPV==22) SF=1.05;
    else if(nPV==23) SF=1.03;
    else if(nPV==24) SF=1.01;
    else if(nPV==25) SF=1.00;
    else if(nPV==26) SF=0.99;
    else if(nPV==27) SF=0.97;
    else if(nPV==28) SF=0.96;
    else if(nPV==29) SF=0.95;
    else if(nPV==30) SF=0.92;
    else if(nPV==31) SF=0.91;
    else if(nPV==32) SF=0.88;
    else if(nPV==33) SF=0.88;
    else if(nPV==34) SF=0.85;
    else if(nPV==35) SF=0.83;
    else if(nPV==36) SF=0.83;
    else if(nPV==37) SF=0.80;
    else if(nPV==38) SF=0.77;
    else if(nPV==39) SF=0.79;
    else if(nPV==40) SF=0.77;
    else if(nPV==41) SF=0.76;
    else if(nPV==42) SF=0.78;
    else if(nPV==43) SF=0.76;
    else if(nPV==44) SF=0.75;
    else if(nPV==45) SF=0.77;
    else if(nPV==46) SF=0.76;
    else if(nPV==47) SF=0.79;
    else if(nPV==48) SF=0.82;
    else if(nPV==49) SF=0.81;
    else if(nPV==50) SF=0.84;
    else if(nPV==51) SF=0.85;
    else if(nPV==52) SF=0.90;
    else if(nPV==53) SF=0.96;
    else if(nPV==54) SF=0.91;
    else if(nPV==55) SF=0.99;
    else if(nPV==56) SF=1.02;
    else if(nPV==57) SF=1.18;
    else if(nPV==58) SF=1.24;
    else if(nPV==59) SF=1.22;
    else if(nPV==60) SF=1.38;
    else if(nPV==61) SF=1.42;
    else if(nPV==62) SF=1.56;
    else if(nPV==63) SF=1.60;
    else if(nPV==64) SF=1.86;
    else if(nPV==65) SF=1.83;
    else if(nPV==66) SF=1.69;
    else if(nPV==67) SF=1.81;
    else if(nPV==68) SF=2.14;
    else if(nPV==69) SF=1.93;
    else if(nPV==70) SF=1.97;
    else if(nPV==71) SF=2.22;
    else if(nPV==72) SF=2.31;
    else if(nPV==73) SF=1.93;
    else if(nPV==74) SF=2.21;
    else if(nPV==75) SF=2.25;
    else if(nPV==76) SF=2.41;
    else if(nPV==77) SF=2.92;
    else if(nPV==78) SF=2.12;
    else if(nPV>=79) SF=2.92;
  }
  else if(Key=="TrigEl23_Pt15"){
    if     (nPV== 1) SF=9.33;
    else if(nPV== 2) SF=6.10;
    else if(nPV== 3) SF=3.09;
    else if(nPV== 4) SF=3.92;
    else if(nPV== 5) SF=3.58;
    else if(nPV== 6) SF=2.92;
    else if(nPV== 7) SF=2.53;
    else if(nPV== 8) SF=2.56;
    else if(nPV== 9) SF=2.37;
    else if(nPV==10) SF=2.02;
    else if(nPV==11) SF=2.01;
    else if(nPV==12) SF=1.73;
    else if(nPV==13) SF=1.64;
    else if(nPV==14) SF=1.48;
    else if(nPV==15) SF=1.36;
    else if(nPV==16) SF=1.29;
    else if(nPV==17) SF=1.26;
    else if(nPV==18) SF=1.21;
    else if(nPV==19) SF=1.14;
    else if(nPV==20) SF=1.14;
    else if(nPV==21) SF=1.13;
    else if(nPV==22) SF=1.09;
    else if(nPV==23) SF=1.03;
    else if(nPV==24) SF=1.04;
    else if(nPV==25) SF=1.01;
    else if(nPV==26) SF=0.98;
    else if(nPV==27) SF=0.98;
    else if(nPV==28) SF=0.91;
    else if(nPV==29) SF=0.93;
    else if(nPV==30) SF=0.89;
    else if(nPV==31) SF=0.91;
    else if(nPV==32) SF=0.84;
    else if(nPV==33) SF=0.81;
    else if(nPV==34) SF=0.82;
    else if(nPV==35) SF=0.77;
    else if(nPV==36) SF=0.73;
    else if(nPV==37) SF=0.70;
    else if(nPV==38) SF=0.66;
    else if(nPV==39) SF=0.64;
    else if(nPV==40) SF=0.65;
    else if(nPV==41) SF=0.57;
    else if(nPV==42) SF=0.61;
    else if(nPV==43) SF=0.54;
    else if(nPV==44) SF=0.55;
    else if(nPV==45) SF=0.60;
    else if(nPV==46) SF=0.55;
    else if(nPV==47) SF=0.59;
    else if(nPV==48) SF=0.57;
    else if(nPV==49) SF=0.61;
    else if(nPV==50) SF=0.52;
    else if(nPV==51) SF=0.58;
    else if(nPV==52) SF=0.59;
    else if(nPV==53) SF=0.66;
    else if(nPV==54) SF=0.72;
    else if(nPV==55) SF=0.69;
    else if(nPV==56) SF=0.73;
    else if(nPV==57) SF=0.81;
    else if(nPV==58) SF=0.73;
    else if(nPV==59) SF=0.80;
    else if(nPV==60) SF=0.75;
    else if(nPV==61) SF=1.01;
    else if(nPV==62) SF=0.99;
    else if(nPV==63) SF=0.79;
    else if(nPV==64) SF=0.82;
    else if(nPV==65) SF=1.56;
    else if(nPV==66) SF=1.20;
    else if(nPV==67) SF=0.98;
    else if(nPV==68) SF=1.60;
    else if(nPV==69) SF=1.87;
    else if(nPV==70) SF=1.44;
    else if(nPV==71) SF=1.13;
    else if(nPV==72) SF=1.12;
    else if(nPV==73) SF=1.82;
    else if(nPV==74) SF=0.58;
    else if(nPV==75) SF=1.88;
    else if(nPV==76) SF=1.80;
    else if(nPV==77) SF=1.50;
    else if(nPV==78) SF=2.02;
    else if(nPV>=79) SF=1.39;
  }
  else if(Key=="TrigEl12_Pt15"){
    if     (nPV== 1) SF=4.87;
    else if(nPV== 2) SF=4.41;
    else if(nPV== 3) SF=2.88;
    else if(nPV== 4) SF=2.31;
    else if(nPV== 5) SF=1.78;
    else if(nPV== 6) SF=1.81;
    else if(nPV== 7) SF=1.79;
    else if(nPV== 8) SF=1.58;
    else if(nPV== 9) SF=1.56;
    else if(nPV==10) SF=1.45;
    else if(nPV==11) SF=1.47;
    else if(nPV==12) SF=1.35;
    else if(nPV==13) SF=1.27;
    else if(nPV==14) SF=1.25;
    else if(nPV==15) SF=1.18;
    else if(nPV==16) SF=1.13;
    else if(nPV==17) SF=1.14;
    else if(nPV==18) SF=1.10;
    else if(nPV==19) SF=1.10;
    else if(nPV==20) SF=1.12;
    else if(nPV==21) SF=1.10;
    else if(nPV==22) SF=1.05;
    else if(nPV==23) SF=1.03;
    else if(nPV==24) SF=1.02;
    else if(nPV==25) SF=1.02;
    else if(nPV==26) SF=0.98;
    else if(nPV==27) SF=0.99;
    else if(nPV==28) SF=0.95;
    else if(nPV==29) SF=0.95;
    else if(nPV==30) SF=0.94;
    else if(nPV==31) SF=0.94;
    else if(nPV==32) SF=0.90;
    else if(nPV==33) SF=0.90;
    else if(nPV==34) SF=0.86;
    else if(nPV==35) SF=0.88;
    else if(nPV==36) SF=0.84;
    else if(nPV==37) SF=0.82;
    else if(nPV==38) SF=0.81;
    else if(nPV==39) SF=0.77;
    else if(nPV==40) SF=0.81;
    else if(nPV==41) SF=0.72;
    else if(nPV==42) SF=0.76;
    else if(nPV==43) SF=0.75;
    else if(nPV==44) SF=0.73;
    else if(nPV==45) SF=0.79;
    else if(nPV==46) SF=0.75;
    else if(nPV==47) SF=0.78;
    else if(nPV==48) SF=0.78;
    else if(nPV==49) SF=0.81;
    else if(nPV==50) SF=0.77;
    else if(nPV==51) SF=0.83;
    else if(nPV==52) SF=0.87;
    else if(nPV==53) SF=0.90;
    else if(nPV==54) SF=0.94;
    else if(nPV==55) SF=1.08;
    else if(nPV==56) SF=1.09;
    else if(nPV==57) SF=1.18;
    else if(nPV==58) SF=1.16;
    else if(nPV==59) SF=1.33;
    else if(nPV==60) SF=1.10;
    else if(nPV==61) SF=1.52;
    else if(nPV==62) SF=1.43;
    else if(nPV==63) SF=1.13;
    else if(nPV==64) SF=1.10;
    else if(nPV==65) SF=2.36;
    else if(nPV==66) SF=1.54;
    else if(nPV==67) SF=1.38;
    else if(nPV==68) SF=1.83;
    else if(nPV==69) SF=2.13;
    else if(nPV==70) SF=2.11;
    else if(nPV==71) SF=1.39;
    else if(nPV==72) SF=1.85;
    else if(nPV==73) SF=2.61;
    else if(nPV==74) SF=1.01;
    else if(nPV==75) SF=2.18;
    else if(nPV==76) SF=3.55;
    else if(nPV==77) SF=2.86;
    else if(nPV==78) SF=2.94;
    else if(nPV>=79) SF=2.11;
  }

  return SF;
}



bool FakeRateMeas::IsNearBJet(Lepton& Lepton, vector<Jet>& BJetColl){

  bool FoundNearB=false;
  for(unsigned int it_b=0; it_b<BJetColl.size(); it_b++){
    if(Lepton.DeltaR(BJetColl.at(it_b))<0.4){ FoundNearB=true; break; }
  }

  return FoundNearB;
}
