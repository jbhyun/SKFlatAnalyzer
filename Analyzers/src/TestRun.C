#include "TestRun.h"

void TestRun::initializeAnalyzer(){

  TriLep=false, TetraLep=false, SS2l=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("SS2l"))   SS2l=true;
    if(Userflags.at(i).Contains("TriLep")) TriLep=true;
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
  vector<JetTagging::Parameters> jtps;
  jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
  mcCorr->SetJetTaggingParameters(jtps);

}


void TestRun::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;
  if(!IsDATA){
    weight*= ev.MCweight()*weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    weight*= GetBRWeight();
    //weight*= GetPileUpWeight(nPileUp, 0);
  }
  FillHist("CutFlow", 0., weight, 20, 0., 20.);

  bool PassTrig=false;
  if(TetraLep or TriLep){
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
  if(!PreCutPass) return;


  TString IDSSLabel = SS2l? "SS":"";
  vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TopHN17T", 10., 2.4);
  vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"T", 10., 2.5);
  vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"L", 10., 2.5);
  vector<Muon>     muonVetoColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronVetoColl = SelectElectrons(electronPreColl, "TopHN17L", 10., 2.5);



  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  vector<Jet> jetNoVetoColl  = GetJets("tight", 25., 2.4);
  sort(jetNoVetoColl.begin(), jetNoVetoColl.end(), PtComparing);
  vector<Jet> bjetNoVetoColl = SelBJets(jetNoVetoColl, param_jets);
  vector<Jet> jetColl  = JetsVetoLeptonInside(jetNoVetoColl, electronVetoColl, muonVetoColl, 0.4);
  vector<Jet> bjetColl = SelBJets(jetColl, param_jets);


  Particle vMET = ev.GetMETVector();
  Particle vMET_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);


  vector<Gen> truthColl;


  bool EventCand = false;
  if(SS2l){ EventCand = muonTightColl.size()==2 or electronTightColl.size()==2 or (muonTightColl.size()==1 && electronTightColl.size()==1); }
  if(TriLep){ EventCand = (muonTightColl.size()+electronTightColl.size())==3; }
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

 
  if(TriLep){
      TestThis(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
               jetColl, bjetColl, vMET_xyCorr, weight, "");
  }

}

void TestRun::TestThis(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                       vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                       vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuT==2 and NElT==1) or (NElT==0 and NMuT==3) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NMuT==2){
    if(!(MuTColl.at(1).Pt()>10 && ElTColl.at(0).Pt()>10)) return;
    vector<Gen> TruthColl=GetGens(); 
    int LepType_Mu1 = GetLeptonType_JH(MuTColl.at(0), TruthColl);
    int LepType_Mu2 = GetLeptonType_JH(MuTColl.at(1), TruthColl);
    int LepType_El1 = GetLeptonType_JH(ElTColl.at(0), TruthColl);

    FillHist("Count_1E2M", 0., weight, 10, 0., 10.);
    if(LepType_Mu1==2 && LepType_Mu2==2) FillHist("Count_1E2M", 1., weight, 10, 0., 10.);
    if(LepType_Mu1==2 && LepType_Mu2==2 && LepType_El1==1) FillHist("Count_1E2M", 2., weight, 10, 0., 10.);
    if(LepType_Mu1==2 && LepType_Mu2==2 && LepType_El1==2) FillHist("Count_1E2M", 3., weight, 10, 0., 10.);
  }
  if(NMuT==3){
    if(!(MuTColl.at(2).Pt()>10)) return;
    vector<Gen> TruthColl=GetGens(); 
    int LepType_Mu1 = GetLeptonType_JH(MuTColl.at(0), TruthColl);
    int LepType_Mu2 = GetLeptonType_JH(MuTColl.at(1), TruthColl);
    int LepType_Mu3 = GetLeptonType_JH(MuTColl.at(2), TruthColl);

    FillHist("Count_3M", 0., weight, 10, 0., 10.);
    if(LepType_Mu1==2 && LepType_Mu2==2) FillHist("Count_3M", 1., weight, 10, 0., 10.);
    if(LepType_Mu1==2 && LepType_Mu2==2 && LepType_Mu3==1) FillHist("Count_3M", 2., weight, 10, 0., 10.);
    if(LepType_Mu1==2 && LepType_Mu2==2 && LepType_Mu3==2) FillHist("Count_3M", 3., weight, 10, 0., 10.);
  }


}


void TestRun::executeEventFromParameter(AnalyzerParameter param){
}



TestRun::TestRun(){

}


TestRun::~TestRun(){

}

