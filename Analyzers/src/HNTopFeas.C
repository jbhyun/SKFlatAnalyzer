#include "HNTopFeas.h"

void HNTopFeas::initializeAnalyzer(){

  ElMuMu=false, MuMuMu=false, TetraLep=false, SS2l=false, SystRun=false; 
  for(unsigned int i=0; i<Userflags.size(); i++){
    if(Userflags.at(i).Contains("ElMuMu")) ElMuMu=true; 
    if(Userflags.at(i).Contains("MuMuMu")) MuMuMu=true; 
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
    TrigList_DblMu.push_back("HLT_TripleMu_12_10_5_v");

    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  if(DataYear==2017){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    TrigList_DblMu.push_back("HLT_TripleMu_10_5_5_DZ_v");
    TrigList_DblMu.push_back("HLT_TripleMu_12_10_5_v");

    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
    TrigList_MuEG.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v");

    TrigList_DblEG.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  }
  if(DataYear==2018){
    TrigList_DblMu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    TrigList_DblMu.push_back("HLT_TripleMu_10_5_5_DZ_v");
    TrigList_DblMu.push_back("HLT_TripleMu_12_10_5_v");

    TrigList_MuEG.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    TrigList_MuEG.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v");

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

  bool PassTrig=false;
  if     (MuMuMu){ PassTrig = ev.PassTrigger(TrigList_DblMu); }
  else if(ElMuMu){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu);
    else{
      if     (MuEG)  PassTrig = ev.PassTrigger(TrigList_MuEG);
      else if(DblMu) PassTrig = (!ev.PassTrigger(TrigList_MuEG)) and ev.PassTrigger(TrigList_DblMu);
    }
  }
  else if(TetraLep){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu) or ev.PassTrigger(TrigList_DblEG);
    else{
      if     (MuEG)  PassTrig = ev.PassTrigger(TrigList_MuEG);
      else if(DblMu) PassTrig = (!ev.PassTrigger(TrigList_MuEG)) and ev.PassTrigger(TrigList_DblMu);
      else if(DblEG) PassTrig = (!(ev.PassTrigger(TrigList_MuEG) or ev.PassTrigger(TrigList_DblMu))) and ev.PassTrigger(DblEG);
    }
  }
  else if(SS2l){
    if(!IsDATA) PassTrig = ev.PassTrigger(TrigList_DblMu) or ev.PassTrigger(DblEG);
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
  if(TetraLep and (muonPreColl.size()+electronPreColl.size())>3 ) PreCutPass=true;
  if(SS2l){
    int NPreMu=muonPreColl.size(), NPreEl=electronPreColl.size();
    if( (NPreMu+NPreEl)>2 ) PreCutPass=true;
    else if(NPreMu==2 and SumCharge(muonPreColl)!=0) PreCutPass=true;
    else if(NPreEl==2 and SumCharge(electronPreColl)!=0) PreCutPass=true;
  }
  if(!PreCutPass) return;
  FillHist("CutFlow", 2., weight*TmpW, 10, 0., 10.);


  std::vector<Muon>     muonTightColl     = SelectMuons(muonPreColl, "TESTT", 10., 2.4);
  std::vector<Electron> electronTightColl = SelectElectrons(electronPreColl, "TESTT", 10., 2.5);
  std::vector<Muon>     muonLooseColl     = SelectMuons(muonPreColl, "TESTL", 10., 2.4);
  std::vector<Electron> electronLooseColl = SelectElectrons(electronPreColl, "TESTL", 10., 2.5);

//  std::vector<Muon>     muonTightColl2     = SelectMuons(muonPreColl, "TEST2T", 10., 2.4);
//  std::vector<Electron> electronTightColl2 = SelectElectrons(electronPreColl, "TEST2T", 10., 2.5);
//  std::vector<Muon>     muonLooseColl2     = SelectMuons(muonPreColl, "TEST2L", 10., 2.4);
//  std::vector<Electron> electronLooseColl2 = SelectElectrons(electronPreColl, "TEST2L", 10., 2.5);


  JetTagging::Parameters param_jets = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
  std::vector<Jet> jetNoVetoColl  = GetJets("tight", 25., 2.4);
  std::sort(jetNoVetoColl.begin(), jetNoVetoColl.end(), PtComparing);
  std::vector<Jet> bjetNoVetoColl = SelBJets(jetNoVetoColl, param_jets);
  std::vector<Jet> jetColl  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl, muonLooseColl, 0.4);
  std::vector<Jet> bjetColl = SelBJets(jetColl, param_jets);
//  std::vector<Jet> jetColl2  = JetsVetoLeptonInside(jetNoVetoColl, electronLooseColl2, muonLooseColl2, 0.4);
//  std::vector<Jet> bjetColl2 = SelBJets(jetColl2, param_jets);


  Particle vMET = ev.GetMETVector();
  Particle vMET_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);


  std::vector<Gen> truthColl;


  bool EventCand = false;
  if(SS2l){ EventCand = muonLooseColl.size()>1 or electronLooseColl.size()>1; }
  if(MuMuMu){ EventCand = muonLooseColl.size()>2; }
  if(ElMuMu){ EventCand = electronLooseColl.size()>0 && muonLooseColl.size()>1; }
  if(TetraLep){ EventCand = (muonLooseColl.size()+electronLooseColl.size())>3; }
//  if(TetraLep){ EventCand = (muonLooseColl.size()+electronLooseColl.size())>3
//                           or (muonLooseColl2.size()+electronLooseColl2.size())>3; }

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
  if(TetraLep){
    AnalyzeTetraLepton(muonTightColl, muonLooseColl, electronTightColl, electronLooseColl,
                       jetColl, bjetColl, vMET_xyCorr, weight, "_Iso04");
//    AnalyzeTetraLepton(muonTightColl2, muonLooseColl2, electronTightColl2, electronLooseColl2,
//                       jetColl2, bjetColl2, vMET_xyCorr, weight, "_MiniIso");
  }


}


void HNTopFeas::AnalyzeSSDiLepton(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                  std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NMuT=MuTColl.size(), NElT=ElTColl.size(), NMuL=MuLColl.size(), NElL=ElLColl.size();
  if( !( (NMuT==2 and NElT==0) or (NElT==2 and NMuT==0) ) ) return;
  if( !(NMuT==NMuL and NElT==NElL) ) return;
  if(NMuT==2){
    if(NMuT==2 and MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) return;
    if(NMuT==2 and MuTColl.at(0).Pt()<20) return;
    FillHist("CutFlow_2Mu"+Label, 3., weight, 10, 0., 10.);
    if(BJetColl.size()<2) return;
    FillHist("CutFlow_2Mu"+Label, 4., weight, 10, 0., 10.);
  
    if(JetColl.size()<3) return;
    FillHist("CutFlow_2Mu"+Label, 5., weight, 10, 0., 10.);
  
    float MSSSF = NMuT==2? (MuTColl.at(0)+MuTColl.at(1)).M():NElT==2? (ElTColl.at(0)+ElTColl.at(1)).M():-1;
    FillHist("MSSSF_2Mu"+Label, MSSSF, weight, 30, 0., 300.);
  
    if(MSSSF>80) return;
    FillHist("CutFlow_2Mu"+Label, 6., weight, 10, 0., 10.);
  }
  else if(NElT==2){
    if(NElT==2 and ElTColl.at(0).Charge()!=ElTColl.at(1).Charge()) return;
    if(NElT==2 and !(ElTColl.at(0).Pt()>25 and ElTColl.at(1).Pt()>15)) return;
    FillHist("CutFlow_2El"+Label, 3., weight, 10, 0., 10.);
    if(BJetColl.size()<2) return;
    FillHist("CutFlow_2El"+Label, 4., weight, 10, 0., 10.);
  
    if(JetColl.size()<3) return;
    FillHist("CutFlow_2El"+Label, 5., weight, 10, 0., 10.);
  
    float MSSSF = NMuT==2? (MuTColl.at(0)+MuTColl.at(1)).M():NElT==2? (ElTColl.at(0)+ElTColl.at(1)).M():-1;
    FillHist("MSSSF_2El"+Label, MSSSF, weight, 30, 0., 300.);
  
    if(MSSSF>80) return;
    FillHist("CutFlow_2El"+Label, 6., weight, 10, 0., 10.);
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


void HNTopFeas::AnalyzeTetraLepton(std::vector<Muon>& MuTColl, std::vector<Muon>& MuLColl, std::vector<Electron>& ElTColl, std::vector<Electron>& ElLColl,
                                   std::vector<Jet>& JetColl, std::vector<Jet>& BJetColl, Particle& vMET, float weight, TString Label)
{
  int NelT = ElTColl.size(), NmuT=MuTColl.size(), NelL = ElLColl.size(), NmuL=MuLColl.size();
  if( (NelT+NmuT)!=4 ) return;
  if( (NelL+NmuL)!=4 ) return;
  if( NmuT<2 ) return;//VmN test

  bool PassTrigPt=false;
  if     (NelT==0){ PassTrigPt = MuTColl.at(0).Pt()>20; }
  else if(NelT==1){ PassTrigPt = MuTColl.at(0).Pt()>20 or ElTColl.at(0).Pt()>25; }
  else if(NelT==2){ PassTrigPt = MuTColl.at(0).Pt()>20 or ElTColl.at(0).Pt()>25; }
  if(!PassTrigPt) return;
  FillHist("CutFlow"+Label, 3., weight, 10, 0., 10.);
  
  vector<int> IdxMupColl, IdxMumColl, IdxElpColl, IdxElmColl;
  for(unsigned int i=0; i<MuTColl.size(); i++){
    if(MuTColl.at(i).Charge()>0){IdxMupColl.push_back(i);} 
    else {IdxMumColl.push_back(i);}
  }
  for(unsigned int i=0; i<ElTColl.size(); i++){
    if(ElTColl.at(i).Charge()>0){IdxElpColl.push_back(i);} 
    else {IdxElmColl.push_back(i);}
  }

  bool IsZlike=false, IsQCDlike=false, IsZTo4l=false, OSSF=false, SSSF=false, OddF=false;
  float M4l=0;
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
  if(IsQCDlike) FillHist("MOSSFType"+Label, 0., weight, 10, 0., 10.);
  if(IsZlike)   FillHist("MOSSFType"+Label, 1., weight, 10, 0., 10.);
  if(IsQCDlike) return;
  FillHist("CutFlow"+Label, 4., weight, 10, 0., 10.);

  if     (NelT==0){ M4l = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+MuTColl.at(3)).M(); }
  else if(NelT==1){ M4l = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+ElTColl.at(0)).M(); }
  else if(NelT==2){ M4l = (MuTColl.at(0)+MuTColl.at(1)+ElTColl.at(0)+ElTColl.at(1)).M(); }

  SSSF = (IdxMupColl.size()==2 and IdxElmColl.size()==2) or (IdxMumColl.size()==2 and IdxElpColl.size()==2);
  OSSF = (IdxMupColl.size()==IdxMumColl.size()) and (IdxElpColl.size()==IdxElmColl.size());
  OddF = NelT==1 and NmuT==3;
  IsZTo4l = OSSF and fabs(M4l-91.2)<10;
  if(IsZTo4l)   FillHist("MOSSFType"+Label, 2., weight, 10, 0., 10.);
  
  if     (OSSF and NelT==2) FillHist("Composition_4l"+Label, 0., weight, 10, 0., 10.);
  else if(OSSF and NelT==0) FillHist("Composition_4l"+Label, 1., weight, 10, 0., 10.);
  else if(SSSF and NelT==2) FillHist("Composition_4l"+Label, 2., weight, 10, 0., 10.);
  else if(OddF and NelT==1) FillHist("Composition_4l"+Label, 3., weight, 10, 0., 10.);
  else if(SSSF and NelT==0) FillHist("Composition_4l"+Label, 4., weight, 10, 0., 10.);
  if(!IsZlike){
    if     (OSSF and NelT==2) FillHist("Composition_4lNoZ2l"+Label, 0., weight, 10, 0., 10.);
    else if(OSSF and NelT==0) FillHist("Composition_4lNoZ2l"+Label, 1., weight, 10, 0., 10.);
    else if(SSSF and NelT==2) FillHist("Composition_4lNoZ2l"+Label, 2., weight, 10, 0., 10.);
    else if(OddF and NelT==1) FillHist("Composition_4lNoZ2l"+Label, 3., weight, 10, 0., 10.);
    else if(SSSF and NelT==0) FillHist("Composition_4lNoZ2l"+Label, 4., weight, 10, 0., 10.);
    if(!IsZTo4l){
      if     (OSSF and NelT==2) FillHist("Composition_4lNoZ2l4l"+Label, 0., weight, 10, 0., 10.);
      else if(OSSF and NelT==0) FillHist("Composition_4lNoZ2l4l"+Label, 1., weight, 10, 0., 10.);
      else if(SSSF and NelT==2) FillHist("Composition_4lNoZ2l4l"+Label, 2., weight, 10, 0., 10.);
      else if(OddF and NelT==1) FillHist("Composition_4lNoZ2l4l"+Label, 3., weight, 10, 0., 10.);
      else if(SSSF and NelT==0) FillHist("Composition_4lNoZ2l4l"+Label, 4., weight, 10, 0., 10.);
    }
  }
  
  if(BJetColl.size()==0) return;
  FillHist("CutFlow"+Label, 5., weight, 10, 0., 10.);

  if     (OSSF and NelT==2) FillHist("Composition_4l1b"+Label, 0., weight, 10, 0., 10.);
  else if(OSSF and NelT==0) FillHist("Composition_4l1b"+Label, 1., weight, 10, 0., 10.);
  else if(SSSF and NelT==2) FillHist("Composition_4l1b"+Label, 2., weight, 10, 0., 10.);
  else if(OddF and NelT==1) FillHist("Composition_4l1b"+Label, 3., weight, 10, 0., 10.);
  else if(SSSF and NelT==0) FillHist("Composition_4l1b"+Label, 4., weight, 10, 0., 10.);
  if(!IsZlike){
    if     (OSSF and NelT==2) FillHist("Composition_4l1bNoZ2l"+Label, 0., weight, 10, 0., 10.);
    else if(OSSF and NelT==0) FillHist("Composition_4l1bNoZ2l"+Label, 1., weight, 10, 0., 10.);
    else if(SSSF and NelT==2) FillHist("Composition_4l1bNoZ2l"+Label, 2., weight, 10, 0., 10.);
    else if(OddF and NelT==1) FillHist("Composition_4l1bNoZ2l"+Label, 3., weight, 10, 0., 10.);
    else if(SSSF and NelT==0) FillHist("Composition_4l1bNoZ2l"+Label, 4., weight, 10, 0., 10.);
    if(!IsZTo4l){
      if     (OSSF and NelT==2) FillHist("Composition_4l1bNoZ2l4l"+Label, 0., weight, 10, 0., 10.);
      else if(OSSF and NelT==0) FillHist("Composition_4l1bNoZ2l4l"+Label, 1., weight, 10, 0., 10.);
      else if(SSSF and NelT==2) FillHist("Composition_4l1bNoZ2l4l"+Label, 2., weight, 10, 0., 10.);
      else if(OddF and NelT==1) FillHist("Composition_4l1bNoZ2l4l"+Label, 3., weight, 10, 0., 10.);
      else if(SSSF and NelT==0) FillHist("Composition_4l1bNoZ2l4l"+Label, 4., weight, 10, 0., 10.);
    }
  }

  if(JetColl.size()<2) return;
  FillHist("CutFlow"+Label, 6., weight, 10, 0., 10.);

  if     (OSSF and NelT==2) FillHist("Composition_4l1b2j"+Label, 0., weight, 10, 0., 10.);
  else if(OSSF and NelT==0) FillHist("Composition_4l1b2j"+Label, 1., weight, 10, 0., 10.);
  else if(SSSF and NelT==2) FillHist("Composition_4l1b2j"+Label, 2., weight, 10, 0., 10.);
  else if(OddF and NelT==1) FillHist("Composition_4l1b2j"+Label, 3., weight, 10, 0., 10.);
  else if(SSSF and NelT==0) FillHist("Composition_4l1b2j"+Label, 4., weight, 10, 0., 10.);
  if(!IsZlike){
    if     (OSSF and NelT==2) FillHist("Composition_4l1b2jNoZ2l"+Label, 0., weight, 10, 0., 10.);
    else if(OSSF and NelT==0) FillHist("Composition_4l1b2jNoZ2l"+Label, 1., weight, 10, 0., 10.);
    else if(SSSF and NelT==2) FillHist("Composition_4l1b2jNoZ2l"+Label, 2., weight, 10, 0., 10.);
    else if(OddF and NelT==1) FillHist("Composition_4l1b2jNoZ2l"+Label, 3., weight, 10, 0., 10.);
    else if(SSSF and NelT==0) FillHist("Composition_4l1b2jNoZ2l"+Label, 4., weight, 10, 0., 10.);
    if(!IsZTo4l){
      if     (OSSF and NelT==2) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 0., weight, 10, 0., 10.);
      else if(OSSF and NelT==0) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 1., weight, 10, 0., 10.);
      else if(SSSF and NelT==2) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 2., weight, 10, 0., 10.);
      else if(OddF and NelT==1) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 3., weight, 10, 0., 10.);
      else if(SSSF and NelT==0) FillHist("Composition_4l1b2jNoZ2l4l"+Label, 4., weight, 10, 0., 10.);
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


