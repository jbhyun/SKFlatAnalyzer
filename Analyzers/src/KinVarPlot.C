#include "KinVarPlot.h"

void KinVarPlot::initializeAnalyzer(){

  TriLep=false, TetraLep=false, SS2l=false, SystRun=false; 
  FakeRun=false, ConvRun=false, FlipRun=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("SS2l"))   SS2l=true;
    if(Userflags.at(i).Contains("TriLep")) TriLep=true;
    if(Userflags.at(i).Contains("TetraLep")) TetraLep=true; 
    if(Userflags.at(i).Contains("SystRun")) SystRun=true; 
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

  MVAreader_Mu->AddVariable("Nj"      , &Nj      );  MVAreader_El->AddVariable("Nj"      , &Nj      );
  MVAreader_Mu->AddVariable("Nb"      , &Nb      );  MVAreader_El->AddVariable("Nb"      , &Nb      );
  MVAreader_Mu->AddVariable("Ptl1"    , &Ptl1    );  MVAreader_El->AddVariable("Ptl1"    , &Ptl1    );
  MVAreader_Mu->AddVariable("Ptl2"    , &Ptl2    );  MVAreader_El->AddVariable("Ptl2"    , &Ptl2    );
  MVAreader_Mu->AddVariable("Ptj1"    , &Ptj1    );  MVAreader_El->AddVariable("Ptj1"    , &Ptj1    );
  MVAreader_Mu->AddVariable("Ptj2"    , &Ptj2    );  MVAreader_El->AddVariable("Ptj2"    , &Ptj2    );
  MVAreader_Mu->AddVariable("Ptj3"    , &Ptj3    );  MVAreader_El->AddVariable("Ptj3"    , &Ptj3    );
  MVAreader_Mu->AddVariable("Ptb1"    , &Ptb1    );  MVAreader_El->AddVariable("Ptb1"    , &Ptb1    );
  MVAreader_Mu->AddVariable("dEtall"  , &dEtall  );  MVAreader_El->AddVariable("dEtall"  , &dEtall  );
  MVAreader_Mu->AddVariable("dRll"    , &dRll    );  MVAreader_El->AddVariable("dRll"    , &dRll    );
  MVAreader_Mu->AddVariable("dRjj12"  , &dRjj12  );  MVAreader_El->AddVariable("dRjj12"  , &dRjj12  );
  MVAreader_Mu->AddVariable("dRjj23"  , &dRjj23  );  MVAreader_El->AddVariable("dRjj23"  , &dRjj23  );
  MVAreader_Mu->AddVariable("dRjj13"  , &dRjj13  );  MVAreader_El->AddVariable("dRjj13"  , &dRjj13  );
  MVAreader_Mu->AddVariable("dRlj11"  , &dRlj11  );  MVAreader_El->AddVariable("dRlj11"  , &dRlj11  );
  MVAreader_Mu->AddVariable("dRlj12"  , &dRlj12  );  MVAreader_El->AddVariable("dRlj12"  , &dRlj12  );
  MVAreader_Mu->AddVariable("dRlj13"  , &dRlj13  );  MVAreader_El->AddVariable("dRlj13"  , &dRlj13  );
  MVAreader_Mu->AddVariable("dRlj21"  , &dRlj21  );  MVAreader_El->AddVariable("dRlj21"  , &dRlj21  );
  MVAreader_Mu->AddVariable("dRlj22"  , &dRlj22  );  MVAreader_El->AddVariable("dRlj22"  , &dRlj22  );
  MVAreader_Mu->AddVariable("dRlj23"  , &dRlj23  );  MVAreader_El->AddVariable("dRlj23"  , &dRlj23  );
  MVAreader_Mu->AddVariable("dRlb11"  , &dRlb11  );  MVAreader_El->AddVariable("dRlb11"  , &dRlb11  );
  MVAreader_Mu->AddVariable("dRlb21"  , &dRlb21  );  MVAreader_El->AddVariable("dRlb21"  , &dRlb21  );
  MVAreader_Mu->AddVariable("MSSSF"   , &MSSSF   );  MVAreader_El->AddVariable("MSSSF"   , &MSSSF   );
  MVAreader_Mu->AddVariable("Mbl11"   , &Mbl11   );  MVAreader_El->AddVariable("Mbl11"   , &Mbl11   );
  MVAreader_Mu->AddVariable("Mbl12"   , &Mbl12   );  MVAreader_El->AddVariable("Mbl12"   , &Mbl12   );
  MVAreader_Mu->AddVariable("Mllb1"   , &Mllb1   );  MVAreader_El->AddVariable("Mllb1"   , &Mllb1   );
  MVAreader_Mu->AddVariable("MTvl1"   , &MTvl1   );  MVAreader_El->AddVariable("MTvl1"   , &MTvl1   );
  MVAreader_Mu->AddVariable("MTvl2"   , &MTvl2   );  MVAreader_El->AddVariable("MTvl2"   , &MTvl2   );
  MVAreader_Mu->AddVariable("Mlj11"   , &Mlj11   );  MVAreader_El->AddVariable("Mlj11"   , &Mlj11   );
  MVAreader_Mu->AddVariable("Mlj12"   , &Mlj12   );  MVAreader_El->AddVariable("Mlj12"   , &Mlj12   );
  MVAreader_Mu->AddVariable("Mlj13"   , &Mlj13   );  MVAreader_El->AddVariable("Mlj13"   , &Mlj13   );
  MVAreader_Mu->AddVariable("Mlj21"   , &Mlj21   );  MVAreader_El->AddVariable("Mlj21"   , &Mlj21   );
  MVAreader_Mu->AddVariable("Mlj22"   , &Mlj22   );  MVAreader_El->AddVariable("Mlj22"   , &Mlj22   );
  MVAreader_Mu->AddVariable("Mlj23"   , &Mlj23   );  MVAreader_El->AddVariable("Mlj23"   , &Mlj23   );
  MVAreader_Mu->AddVariable("Mllj1"   , &Mllj1   );  MVAreader_El->AddVariable("Mllj1"   , &Mllj1   );
  MVAreader_Mu->AddVariable("Mllj2"   , &Mllj2   );  MVAreader_El->AddVariable("Mllj2"   , &Mllj2   );
  MVAreader_Mu->AddVariable("Mllj3"   , &Mllj3   );  MVAreader_El->AddVariable("Mllj3"   , &Mllj3   );
  MVAreader_Mu->AddVariable("Mlljj12" , &Mlljj12 );  MVAreader_El->AddVariable("Mlljj12" , &Mlljj12 );
  MVAreader_Mu->AddVariable("Mlljj13" , &Mlljj13 );  MVAreader_El->AddVariable("Mlljj13" , &Mlljj13 );
  MVAreader_Mu->AddVariable("Mlljj23" , &Mlljj23 );  MVAreader_El->AddVariable("Mlljj23" , &Mlljj23 );
  MVAreader_Mu->AddVariable("Mljj112" , &Mljj112 );  MVAreader_El->AddVariable("Mljj112" , &Mljj112 );
  MVAreader_Mu->AddVariable("Mljj113" , &Mljj113 );  MVAreader_El->AddVariable("Mljj113" , &Mljj113 );
  MVAreader_Mu->AddVariable("Mljj123" , &Mljj123 );  MVAreader_El->AddVariable("Mljj123" , &Mljj123 );
  MVAreader_Mu->AddVariable("Mljj212" , &Mljj212 );  MVAreader_El->AddVariable("Mljj212" , &Mljj212 );
  MVAreader_Mu->AddVariable("Mljj213" , &Mljj213 );  MVAreader_El->AddVariable("Mljj213" , &Mljj213 );
  MVAreader_Mu->AddVariable("Mljj223" , &Mljj223 );  MVAreader_El->AddVariable("Mljj223" , &Mljj223 );
  MVAreader_Mu->AddVariable("Mjj12"   , &Mjj12   );  MVAreader_El->AddVariable("Mjj12"   , &Mjj12   );
  MVAreader_Mu->AddVariable("Mjj13"   , &Mjj13   );  MVAreader_El->AddVariable("Mjj13"   , &Mjj13   );
  MVAreader_Mu->AddVariable("Mjj23"   , &Mjj23   );  MVAreader_El->AddVariable("Mjj23"   , &Mjj23   );
  MVAreader_Mu->AddVariable("HT"      , &HT      );  MVAreader_El->AddVariable("HT"      , &HT      );
  MVAreader_Mu->AddVariable("MET2HT"  , &MET2HT  );  MVAreader_El->AddVariable("MET2HT"  , &MET2HT  );
  MVAreader_Mu->AddVariable("MllW_2jL", &MllW_2jL);  MVAreader_El->AddVariable("MllW_2jL", &MllW_2jL);
  MVAreader_Mu->AddVariable("MllW_1jL", &MllW_1jL);  MVAreader_El->AddVariable("MllW_1jL", &MllW_1jL);
  MVAreader_Mu->AddVariable("MllW1_H" , &MllW1_H );  MVAreader_El->AddVariable("MllW1_H" , &MllW1_H );
  MVAreader_Mu->AddVariable("Ml1W_2jL", &Ml1W_2jL);  MVAreader_El->AddVariable("Ml1W_2jL", &Ml1W_2jL);
  MVAreader_Mu->AddVariable("Ml1W_1jL", &Ml1W_1jL);  MVAreader_El->AddVariable("Ml1W_1jL", &Ml1W_1jL);
  MVAreader_Mu->AddVariable("Ml2W_2jL", &Ml2W_2jL);  MVAreader_El->AddVariable("Ml2W_2jL", &Ml2W_2jL);
  MVAreader_Mu->AddVariable("Ml2W_1jL", &Ml2W_1jL);  MVAreader_El->AddVariable("Ml2W_1jL", &Ml2W_1jL);
  MVAreader_Mu->AddVariable("Ml1W1_H" , &Ml1W1_H );  MVAreader_El->AddVariable("Ml1W1_H" , &Ml1W1_H );
  MVAreader_Mu->AddVariable("Ml1W2_H" , &Ml1W2_H );  MVAreader_El->AddVariable("Ml1W2_H" , &Ml1W2_H );
  MVAreader_Mu->AddVariable("Ml2W1_H" , &Ml2W1_H );  MVAreader_El->AddVariable("Ml2W1_H" , &Ml2W1_H );
  MVAreader_Mu->AddVariable("Ml2W2_H" , &Ml2W2_H );  MVAreader_El->AddVariable("Ml2W2_H" , &Ml2W2_H );
  MVAreader_Mu->AddVariable("MjjW1"   , &MjjW1   );  MVAreader_El->AddVariable("MjjW1"   , &MjjW1   );
  MVAreader_Mu->AddVariable("MjjW2"   , &MjjW2   );  MVAreader_El->AddVariable("MjjW2"   , &MjjW2   );
  MVAreader_Mu->AddVariable("Ptb2"    , &Ptb2    );  MVAreader_El->AddVariable("Ptb2"    , &Ptb2    );
  MVAreader_Mu->AddVariable("Mbl21"   , &Mbl21   );  MVAreader_El->AddVariable("Mbl21"   , &Mbl21   );
  MVAreader_Mu->AddVariable("Mbl22"   , &Mbl22   );  MVAreader_El->AddVariable("Mbl22"   , &Mbl22   );
  MVAreader_Mu->AddVariable("Mllb2"   , &Mllb2   );  MVAreader_El->AddVariable("Mllb2"   , &Mllb2   );
  MVAreader_Mu->AddVariable("Mllj4"   , &Mllj4   );  MVAreader_El->AddVariable("Mllj4"   , &Mllj4   );
  MVAreader_Mu->AddVariable("Mlljj14" , &Mlljj14 );  MVAreader_El->AddVariable("Mlljj14" , &Mlljj14 );
  MVAreader_Mu->AddVariable("Mlljj24" , &Mlljj24 );  MVAreader_El->AddVariable("Mlljj24" , &Mlljj24 );
  MVAreader_Mu->AddVariable("Mlljj34" , &Mlljj34 );  MVAreader_El->AddVariable("Mlljj34" , &Mlljj34 );
  MVAreader_Mu->AddVariable("Mljj114" , &Mljj114 );  MVAreader_El->AddVariable("Mljj114" , &Mljj114 );
  MVAreader_Mu->AddVariable("Mljj124" , &Mljj124 );  MVAreader_El->AddVariable("Mljj124" , &Mljj124 );
  MVAreader_Mu->AddVariable("Mljj134" , &Mljj134 );  MVAreader_El->AddVariable("Mljj134" , &Mljj134 );
  MVAreader_Mu->AddVariable("Mljj214" , &Mljj214 );  MVAreader_El->AddVariable("Mljj214" , &Mljj214 );
  MVAreader_Mu->AddVariable("Mljj224" , &Mljj224 );  MVAreader_El->AddVariable("Mljj224" , &Mljj224 );
  MVAreader_Mu->AddVariable("Mljj234" , &Mljj234 );  MVAreader_El->AddVariable("Mljj234" , &Mljj234 );
  MVAreader_Mu->AddVariable("Mjj14"   , &Mjj14   );  MVAreader_El->AddVariable("Mjj14"   , &Mjj14   );
  MVAreader_Mu->AddVariable("Mjj24"   , &Mjj24   );  MVAreader_El->AddVariable("Mjj24"   , &Mjj24   );
  MVAreader_Mu->AddVariable("Mjj34"   , &Mjj34   );  MVAreader_El->AddVariable("Mjj34"   , &Mjj34   );
  MVAreader_Mu->AddVariable("MllW2_H" , &MllW2_H );  MVAreader_El->AddVariable("MllW2_H" , &MllW2_H );


  TString AnalyzerPath=std::getenv("SKFlat_WD"), MVAPath="/data/Run2Legacy_v4/2017/MVAClassifier/";
  TString FileName_Mu="TMVAClassification_MN20_Mu_BDTG.weights.xml"; 
  TString FileName_El="TMVAClassification_MN20_El_BDTG.weights.xml"; 
  MVAreader_Mu->BookMVA("BDTG method", AnalyzerPath+MVAPath+FileName_Mu);
  MVAreader_El->BookMVA("BDTG method", AnalyzerPath+MVAPath+FileName_El);

  outfile->cd();

}


void KinVarPlot::executeEvent(){


  Event ev = GetEvent();
  float weight = 1.;
  if(!IsDATA){
    weight*= ev.MCweight()*weight_norm_1invpb*GetKFactor()*ev.GetTriggerLumi("Full");
    weight*= GetBRWeight();
    weight*= GetPileUpWeight(nPileUp, 0);
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
    if(FlipRun && IsDATA && (NPreMu+NPreEl)>1) PreCutPass=true;
  }
  if(!PreCutPass) return;


  TString IDSSLabel = SS2l? "SS":""; TString TLabel = FakeRun? "L":"T";
  vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TopHN17"+TLabel, 10., 2.4);
  vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+TLabel, 15., 2.5);
  vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TopHN17"+IDSSLabel+"L", 15., 2.5);
  vector<Muon>     muonVetoColl      = SelectMuons(muonPreColl, "TopHN17L", 10., 2.4);
  vector<Electron> electronVetoColl  = SelectElectrons(electronPreColl, "TopHN17L", 10., 2.5);


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

  float w_topptrw = 1., w_prefire = 1., sf_trig = 1., w_FR = 1.;
  float sf_mutk = 1., sf_muid = 1., sf_muiso = 1., sf_elreco = 1., sf_elid = 1., sf_btag = 1.;
//  float w_CF = IsDATA && FlipRun? 0.:1.;
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
  if(IsDATA && FakeRun && EventCand){
    w_FR      = CalcTestFakeWeight(muonLooseColl, electronLooseColl,
                  "TopHN17T", "TopHN17L", "TopHN17"+IDSSLabel+"T", "TopHN17"+IDSSLabel+"L", bjetColl.size(), 0);
  }
//  if(IsDATA && FlipRun && EventCand){
//    for(unsigned int ie=0; ie<electronTightColl.size(); ie++){
//      w_CF += GetCFRSF(electronTightColl.at(ie), "App2Bin2_Fin", "DataEff");
//    }
//  }
  weight *= w_topptrw * w_prefire * sf_trig * w_FR;
  weight *= sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;
//  weight *= w_CF;

 
  if(SS2l){
    if(!(IsDATA && FlipRun)){
      MakePlotSS2L(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_xyCorr, weight, "");
    }
    else{
      for(unsigned ie=0; ie<electronTightColl.size(); ie++){
        double w_CFN = GetCFRSF(electronTightColl.at(ie), "App2Bin2_Fin", "DataEff");
        double PTN   = GetFlipCorrPT(electronTightColl.at(ie), "App2Bin2_Fin", "DataScaleSmear");
        std::vector<Electron> ElTCollN = electronTightColl;
        ElTCollN.at(ie).SetPtEtaPhiM(PTN, electronTightColl.at(ie).Eta(), electronTightColl.at(ie).Phi(), electronTightColl.at(ie).M());
        MakePlotSS2L(muonTightColl, muonLooseColl, muonVetoColl, ElTCollN, electronLooseColl, electronVetoColl,
                     jetColl, bjetColl, vMET_xyCorr, weight*w_CFN, "");
      }
    }
  }

}

void KinVarPlot::MakePlotSS2L(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                              vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                              vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if( FakeRun      and weight==0.  ) return; 
  if(NMuT==2){
    if(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    double Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(DataYear>2016 && Mll<4) return; 
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      if(ConvRun && GenLepInfo>=100) return;
      if(FlipRun && (GenLepInfo>999 or GenLepInfo<100)) return;
    }
    InitializeTreeVars();
    FillHist("NbPre_2M", BJetColl.size(), weight, 5, 0., 5.);
    if(BJetColl.size()<1) return;
    FillHist("NjPre_2M", JetColl.size(), weight, 10, 0., 10.);
    if(JetColl.size() <3) return;
    FillHist("MSSSFPre_2M", Mll, weight, 25, 0., 250.);

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
    PlotParameters("_2M");
  }
  if(NElT==2){
    int aSumQ = abs(SumCharge(ElLColl));
    if(IsDATA && FlipRun){ if(aSumQ!=0) return; }
    else                 { if(aSumQ==0) return; }
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    if(!IsDATA){
      int GenLepInfo = GetGenLepInfo(ElTColl, MuTColl);
      int IdxFlipped = GenLepInfo % 1000;
      if(ConvRun && GenLepInfo>=100) return;
      if(FlipRun && (GenLepInfo>999 or GenLepInfo<100)) return;
      if(FlipRun && IdxFlipped<2){ weight *= GetCFRSF(ElTColl.at(IdxFlipped), "App2Bin1_Fin"); }
    }
    InitializeTreeVars();
    double Mll = (ElTColl.at(0)+ElTColl.at(1)).M();
    FillHist("NbPre_2E", BJetColl.size(), weight, 5, 0., 5.);
    if(BJetColl.size()<1) return;
    FillHist("NjPre_2E", JetColl.size(), weight, 10, 0., 10.);
    if(JetColl.size() <3) return;
    FillHist("MSSSFPre_2E", Mll, weight, 25, 0., 250.);
    if(fabs(Mll-91.2)<10.) return;

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
    disc_BDTG = MVAreader_El->EvaluateMVA("BDTG method");
    PlotParameters("_2E");
  }

}


void KinVarPlot::PlotParameters(TString Label){

  FillHist("NEvt"+Label, 0., w_tot, 1, 0., 1.);
  FillHist("disc_BDTG"+Label, disc_BDTG, w_tot, 40, -1., 1.);
  FillHist("Nj"+Label, Nj, w_tot, 10, 0., 10.);
  FillHist("Nb"+Label, Nb, w_tot, 5, 0., 5.);
  FillHist("Ptl1"+Label, Ptl1, w_tot, 25, 0., 250.);
  FillHist("Ptl2"+Label, Ptl2, w_tot, 20, 0., 100.);
  FillHist("Ptj1"+Label, Ptj1, w_tot, 25, 0., 500.);
  FillHist("Ptj2"+Label, Ptj2, w_tot, 30, 0., 300.);
  FillHist("Ptj3"+Label, Ptj3, w_tot, 20, 0., 200.);
  FillHist("Ptb1"+Label, Ptb1, w_tot, 35, 0., 350.);
  FillHist("Ptb2"+Label, max((Float_t)0.,Ptb2), w_tot, 15, 0., 150.);
  FillHist("MET"+Label, MET, w_tot, 25, 0., 250.);
  FillHist("HT"+Label, HT, w_tot, 25, 0., 1000.);
  FillHist("MET2HT"+Label, MET2HT, w_tot, 20, 0., 100.);
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
  FillHist("MSSSF"+Label, MSSSF, w_tot, 25, 0., 250.);
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
  FillHist("MTvll"+Label, MTvll, w_tot, 35, 0., 350.);
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

}


void KinVarPlot::executeEventFromParameter(AnalyzerParameter param){
}


void KinVarPlot::InitializeTreeVars(){

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

float KinVarPlot::CalcTestFakeWeight(vector<Muon>& MuColl, vector<Electron>& ElColl, TString MuIDT, TString MuIDL, TString ElIDT, TString ElIDL, int NBJet, int SystDir){

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


float KinVarPlot::GetTestElFR(Electron& El, TString Key, int SystDir){

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


float KinVarPlot::GetTestMuFR(Muon& Mu, TString Key, int SystDir){

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




float KinVarPlot::GetCFRSF(Electron& El, TString Tag, TString Option){

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



float KinVarPlot::GetFlipCorrPT(Electron& El, TString Tag, TString Option){

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



int KinVarPlot::GetGenLepInfo(vector<Electron>& ElColl, vector<Muon>& MuColl, TString Option){

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



KinVarPlot::KinVarPlot(){

  TMVA::Tools::Instance();
  MVAreader_Mu = new TMVA::Reader();
  MVAreader_El = new TMVA::Reader();

}


KinVarPlot::~KinVarPlot(){

  delete MVAreader_Mu;
  delete MVAreader_El;

}
