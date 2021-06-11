#include "KinVarSearch.h"

void KinVarSearch::initializeAnalyzer(){

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

  InitializeTreeVars();
  tree_mm->Branch("Nj", &Nj, "Nj/I");                 tree_ee->Branch("Nj", &Nj, "Nj/I");
  tree_mm->Branch("Nb", &Nb, "Nb/I");                 tree_ee->Branch("Nb", &Nb, "Nb/I");
  tree_mm->Branch("Ptl1", &Ptl1, "Ptl1/F");           tree_ee->Branch("Ptl1", &Ptl1, "Ptl1/F");
  tree_mm->Branch("Ptl2", &Ptl2, "Ptl2/F");           tree_ee->Branch("Ptl2", &Ptl2, "Ptl2/F");
  tree_mm->Branch("Ptj1", &Ptj1, "Ptj1/F");           tree_ee->Branch("Ptj1", &Ptj1, "Ptj1/F");
  tree_mm->Branch("Ptj2", &Ptj2, "Ptj2/F");           tree_ee->Branch("Ptj2", &Ptj2, "Ptj2/F");
  tree_mm->Branch("Ptj3", &Ptj3, "Ptj3/F");           tree_ee->Branch("Ptj3", &Ptj3, "Ptj3/F");
  tree_mm->Branch("Ptb1", &Ptb1, "Ptb1/F");           tree_ee->Branch("Ptb1", &Ptb1, "Ptb1/F");
  tree_mm->Branch("Ptb2", &Ptb2, "Ptb2/F");           tree_ee->Branch("Ptb2", &Ptb2, "Ptb2/F");
  tree_mm->Branch("MET", &MET, "MET/F");              tree_ee->Branch("MET", &MET, "MET/F");
  tree_mm->Branch("HT", &HT, "HT/F");                 tree_ee->Branch("HT", &HT, "HT/F");
  tree_mm->Branch("MET2HT", &MET2HT, "MET2HT/F");     tree_ee->Branch("MET2HT", &MET2HT, "MET2HT/F");
  tree_mm->Branch("dEtall", &dEtall, "dEtall/F");     tree_ee->Branch("dEtall", &dEtall, "dEtall/F");
  tree_mm->Branch("dRll", &dRll, "dRll/F");           tree_ee->Branch("dRll", &dRll, "dRll/F");
  tree_mm->Branch("dRjj12", &dRjj12, "dRjj12/F");     tree_ee->Branch("dRjj12", &dRjj12, "dRjj12/F");
  tree_mm->Branch("dRjj23", &dRjj23, "dRjj23/F");     tree_ee->Branch("dRjj23", &dRjj23, "dRjj23/F");
  tree_mm->Branch("dRjj13", &dRjj13, "dRjj13/F");     tree_ee->Branch("dRjj13", &dRjj13, "dRjj13/F");
  tree_mm->Branch("dRlj11", &dRlj11, "dRlj11/F");     tree_ee->Branch("dRlj11", &dRlj11, "dRlj11/F");
  tree_mm->Branch("dRlj12", &dRlj12, "dRlj12/F");     tree_ee->Branch("dRlj12", &dRlj12, "dRlj12/F");
  tree_mm->Branch("dRlj13", &dRlj13, "dRlj13/F");     tree_ee->Branch("dRlj13", &dRlj13, "dRlj13/F");
  tree_mm->Branch("dRlj21", &dRlj21, "dRlj21/F");     tree_ee->Branch("dRlj21", &dRlj21, "dRlj21/F");
  tree_mm->Branch("dRlj22", &dRlj22, "dRlj22/F");     tree_ee->Branch("dRlj22", &dRlj22, "dRlj22/F");
  tree_mm->Branch("dRlj23", &dRlj23, "dRlj23/F");     tree_ee->Branch("dRlj23", &dRlj23, "dRlj23/F");
  tree_mm->Branch("dRlb11", &dRlb11, "dRlb11/F");     tree_ee->Branch("dRlb11", &dRlb11, "dRlb11/F");
  tree_mm->Branch("dRlb21", &dRlb21, "dRlb21/F");     tree_ee->Branch("dRlb21", &dRlb21, "dRlb21/F");
  tree_mm->Branch("MSSSF", &MSSSF, "MSSSF/F");        tree_ee->Branch("MSSSF", &MSSSF, "MSSSF/F");
  tree_mm->Branch("Mbl11", &Mbl11, "Mbl11/F");        tree_ee->Branch("Mbl11", &Mbl11, "Mbl11/F");
  tree_mm->Branch("Mbl12", &Mbl12, "Mbl12/F");        tree_ee->Branch("Mbl12", &Mbl12, "Mbl12/F");
  tree_mm->Branch("Mbl21", &Mbl21, "Mbl21/F");        tree_ee->Branch("Mbl21", &Mbl21, "Mbl21/F");
  tree_mm->Branch("Mbl22", &Mbl22, "Mbl22/F");        tree_ee->Branch("Mbl22", &Mbl22, "Mbl22/F");
  tree_mm->Branch("Mlj11", &Mlj11, "Mlj11/F");        tree_ee->Branch("Mlj11", &Mlj11, "Mlj11/F");
  tree_mm->Branch("Mlj12", &Mlj12, "Mlj12/F");        tree_ee->Branch("Mlj12", &Mlj12, "Mlj12/F");
  tree_mm->Branch("Mlj13", &Mlj13, "Mlj13/F");        tree_ee->Branch("Mlj13", &Mlj13, "Mlj13/F");
  tree_mm->Branch("Mlj21", &Mlj21, "Mlj21/F");        tree_ee->Branch("Mlj21", &Mlj21, "Mlj21/F");
  tree_mm->Branch("Mlj22", &Mlj22, "Mlj22/F");        tree_ee->Branch("Mlj22", &Mlj22, "Mlj22/F");
  tree_mm->Branch("Mlj23", &Mlj23, "Mlj23/F");        tree_ee->Branch("Mlj23", &Mlj23, "Mlj23/F");
  tree_mm->Branch("MTvl1", &MTvl1, "MTvl1/F");        tree_ee->Branch("MTvl1", &MTvl1, "MTvl1/F");
  tree_mm->Branch("MTvl2", &MTvl2, "MTvl2/F");        tree_ee->Branch("MTvl2", &MTvl2, "MTvl2/F");
  tree_mm->Branch("Mllj1", &Mllj1, "Mllj1/F");        tree_ee->Branch("Mllj1", &Mllj1, "Mllj1/F");
  tree_mm->Branch("Mllj2", &Mllj2, "Mllj2/F");        tree_ee->Branch("Mllj2", &Mllj2, "Mllj2/F");
  tree_mm->Branch("Mllj3", &Mllj3, "Mllj3/F");        tree_ee->Branch("Mllj3", &Mllj3, "Mllj3/F");
  tree_mm->Branch("Mllj4", &Mllj4, "Mllj4/F");        tree_ee->Branch("Mllj4", &Mllj4, "Mllj4/F");
  tree_mm->Branch("Mllb1", &Mllb1, "Mllb1/F");        tree_ee->Branch("Mllb1", &Mllb1, "Mllb1/F");
  tree_mm->Branch("Mllb2", &Mllb2, "Mllb2/F");        tree_ee->Branch("Mllb2", &Mllb2, "Mllb2/F");
  tree_mm->Branch("Mlljj12", &Mlljj12, "Mlljj12/F");  tree_ee->Branch("Mlljj12", &Mlljj12, "Mlljj12/F");
  tree_mm->Branch("Mlljj13", &Mlljj13, "Mlljj13/F");  tree_ee->Branch("Mlljj13", &Mlljj13, "Mlljj13/F");
  tree_mm->Branch("Mlljj14", &Mlljj14, "Mlljj14/F");  tree_ee->Branch("Mlljj14", &Mlljj14, "Mlljj14/F");
  tree_mm->Branch("Mlljj23", &Mlljj23, "Mlljj23/F");  tree_ee->Branch("Mlljj23", &Mlljj23, "Mlljj23/F");
  tree_mm->Branch("Mlljj24", &Mlljj24, "Mlljj24/F");  tree_ee->Branch("Mlljj24", &Mlljj24, "Mlljj24/F");
  tree_mm->Branch("Mlljj34", &Mlljj34, "Mlljj34/F");  tree_ee->Branch("Mlljj34", &Mlljj34, "Mlljj34/F");
  tree_mm->Branch("Mljj112", &Mljj112, "Mljj112/F");  tree_ee->Branch("Mljj112", &Mljj112, "Mljj112/F");
  tree_mm->Branch("Mljj113", &Mljj113, "Mljj113/F");  tree_ee->Branch("Mljj113", &Mljj113, "Mljj113/F");
  tree_mm->Branch("Mljj114", &Mljj114, "Mljj114/F");  tree_ee->Branch("Mljj114", &Mljj114, "Mljj114/F");
  tree_mm->Branch("Mljj123", &Mljj123, "Mljj123/F");  tree_ee->Branch("Mljj123", &Mljj123, "Mljj123/F");
  tree_mm->Branch("Mljj124", &Mljj124, "Mljj124/F");  tree_ee->Branch("Mljj124", &Mljj124, "Mljj124/F");
  tree_mm->Branch("Mljj134", &Mljj134, "Mljj134/F");  tree_ee->Branch("Mljj134", &Mljj134, "Mljj134/F");
  tree_mm->Branch("Mljj212", &Mljj212, "Mljj212/F");  tree_ee->Branch("Mljj212", &Mljj212, "Mljj212/F");
  tree_mm->Branch("Mljj213", &Mljj213, "Mljj213/F");  tree_ee->Branch("Mljj213", &Mljj213, "Mljj213/F");
  tree_mm->Branch("Mljj214", &Mljj214, "Mljj214/F");  tree_ee->Branch("Mljj214", &Mljj214, "Mljj214/F");
  tree_mm->Branch("Mljj223", &Mljj223, "Mljj223/F");  tree_ee->Branch("Mljj223", &Mljj223, "Mljj223/F");
  tree_mm->Branch("Mljj224", &Mljj224, "Mljj224/F");  tree_ee->Branch("Mljj224", &Mljj224, "Mljj224/F");
  tree_mm->Branch("Mljj234", &Mljj234, "Mljj234/F");  tree_ee->Branch("Mljj234", &Mljj234, "Mljj234/F");
  tree_mm->Branch("Mjj12", &Mjj12, "Mjj12/F");        tree_ee->Branch("Mjj12", &Mjj12, "Mjj12/F");
  tree_mm->Branch("Mjj13", &Mjj13, "Mjj13/F");        tree_ee->Branch("Mjj13", &Mjj13, "Mjj13/F");
  tree_mm->Branch("Mjj14", &Mjj14, "Mjj14/F");        tree_ee->Branch("Mjj14", &Mjj14, "Mjj14/F");
  tree_mm->Branch("Mjj23", &Mjj23, "Mjj23/F");        tree_ee->Branch("Mjj23", &Mjj23, "Mjj23/F");
  tree_mm->Branch("Mjj24", &Mjj24, "Mjj24/F");        tree_ee->Branch("Mjj24", &Mjj24, "Mjj24/F");
  tree_mm->Branch("Mjj34", &Mjj34, "Mjj34/F");        tree_ee->Branch("Mjj34", &Mjj34, "Mjj34/F");
  tree_mm->Branch("MjjW1", &MjjW1, "MjjW1/F");          tree_ee->Branch("MjjW1", &MjjW1, "MjjW1/F");
  tree_mm->Branch("MjjW2", &MjjW2, "MjjW2/F");          tree_ee->Branch("MjjW2", &MjjW2, "MjjW2/F");
  tree_mm->Branch("MllW_2jL", &MllW_2jL, "MllW_2jL/F"); tree_ee->Branch("MllW_2jL", &MllW_2jL, "MllW_2jL/F");
  tree_mm->Branch("MllW_1jL", &MllW_1jL, "MllW_1jL/F"); tree_ee->Branch("MllW_1jL", &MllW_1jL, "MllW_1jL/F");
  tree_mm->Branch("MllW1_H", &MllW1_H, "MllW1_H/F");    tree_ee->Branch("MllW1_H", &MllW1_H, "MllW1_H/F");
  tree_mm->Branch("MllW2_H", &MllW2_H, "MllW2_H/F");    tree_ee->Branch("MllW2_H", &MllW2_H, "MllW2_H/F");
  tree_mm->Branch("Ml1W_2jL", &Ml1W_2jL, "Ml1W_2jL/F"); tree_ee->Branch("Ml1W_2jL", &Ml1W_2jL, "Ml1W_2jL/F");
  tree_mm->Branch("Ml1W_1jL", &Ml1W_1jL, "Ml1W_1jL/F"); tree_ee->Branch("Ml1W_1jL", &Ml1W_1jL, "Ml1W_1jL/F");
  tree_mm->Branch("Ml2W_2jL", &Ml2W_2jL, "Ml2W_2jL/F"); tree_ee->Branch("Ml2W_2jL", &Ml2W_2jL, "Ml2W_2jL/F");
  tree_mm->Branch("Ml2W_1jL", &Ml2W_1jL, "Ml2W_1jL/F"); tree_ee->Branch("Ml2W_1jL", &Ml2W_1jL, "Ml2W_1jL/F");
  tree_mm->Branch("Ml1W1_H", &Ml1W1_H, "Ml1W1_H/F");    tree_ee->Branch("Ml1W1_H", &Ml1W1_H, "Ml1W1_H/F");
  tree_mm->Branch("Ml1W2_H", &Ml1W2_H, "Ml1W2_H/F");    tree_ee->Branch("Ml1W2_H", &Ml1W2_H, "Ml1W2_H/F");
  tree_mm->Branch("Ml2W1_H", &Ml2W1_H, "Ml2W1_H/F");    tree_ee->Branch("Ml2W1_H", &Ml2W1_H, "Ml2W1_H/F");
  tree_mm->Branch("Ml2W2_H", &Ml2W2_H, "Ml2W2_H/F");    tree_ee->Branch("Ml2W2_H", &Ml2W2_H, "Ml2W2_H/F");
  tree_mm->Branch("w_tot", &w_tot, "w_tot/F");        tree_ee->Branch("w_tot", &w_tot, "w_tot/F");

  outfile->cd();

}


void KinVarSearch::executeEvent(){


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

 
  if(SS2l){
//    AnalyzeSSDiLepton(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
//                      jetColl, bjetColl, vMET_xyCorr, weight, "");
      MakeTreeSS2L(muonTightColl, muonLooseColl, muonVetoColl, electronTightColl, electronLooseColl, electronVetoColl,
                   jetColl, bjetColl, vMET_xyCorr, weight, "");
  }

}

void KinVarSearch::MakeTreeSS2L(vector<Muon>& MuTColl, vector<Muon>& MuLColl, vector<Muon>& MuVColl,
                                vector<Electron>& ElTColl, vector<Electron>& ElLColl, vector<Electron>& ElVColl,
                                vector<Jet>& JetColl, vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{

  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size(), NMuV=MuVColl.size(), NElV=ElVColl.size();
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if( !(NMuT==NMuV and NElT==NElV) ) return;
  if(NMuT==2){
    if(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) return;
    if(!(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10)) return;
    float Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(DataYear>2016 && Mll<4) return; 
    if(BJetColl.size()<1) return;
    if(JetColl.size() <3) return;
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
    w_tot    = !IsDATA? weight:-1.;
    tree_mm->Fill();
  }
  if(NElT==2){
    if(ElTColl.at(0).Charge()!=ElTColl.at(1).Charge()) return;
    if(!(ElTColl.at(0).Pt()>25 && ElTColl.at(1).Pt()>15)) return;
    float Mll = (MuTColl.at(0)+MuTColl.at(1)).M();
    if(BJetColl.size()<1) return;
    if(JetColl.size() <3) return;
    if(fabs(Mll-91.2)<10.) return;
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
    w_tot    = !IsDATA? weight:-1.;
    tree_ee->Fill();
  }

}



void KinVarSearch::executeEventFromParameter(AnalyzerParameter param){
}


void KinVarSearch::InitializeTreeVars(){

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
  w_tot=-1;

}



KinVarSearch::KinVarSearch(){

  tree_mm = new TTree("Tree_mm", "Tree_mm");
  tree_ee = new TTree("Tree_ee", "Tree_ee");

}


KinVarSearch::~KinVarSearch(){

  delete tree_mm;
  delete tree_ee;

}


void KinVarSearch::WriteHist(){

  outfile->cd();
  tree_mm->Write();
  tree_ee->Write();
  outfile->cd();

}
