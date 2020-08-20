float GetNentries(TString filename){

  TFile* file = new TFile(filename);

  TTree* tree;
  file->GetObject("recoTree/SKFlat", tree);
  Long64_t nentries = tree->GetEntries();

  delete file;

  return nentries;
}


float GetSumW(TString filename){

  TFile* file = new TFile(filename);

  TTree* tree;
  file->GetObject("recoTree/SKFlat", tree);

  double gen_Weight=0;
  TBranch        *b_gen_Weight;

  tree->SetBranchAddress("gen_weight", &gen_Weight, &b_gen_Weight);

  double sum_weight=0.;
  Long64_t nentries = tree->GetEntries();

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    tree->LoadTree(jentry);
    tree->GetEntry(jentry);
    sum_weight += gen_Weight;
  }

  delete file;

  return sum_weight;
}

void LumiWCalc(std::string filepath){

  std::vector<TString> Vec_subfile;
  ifstream filelist(filepath);
  string line;
  while (getline(filelist,line)){
    TString currentfile(line);
    Vec_subfile.push_back(currentfile);
  }

  float Nevent=0; float SumW=0.;

  for(unsigned int it=0; it<Vec_subfile.size(); it++){
    Nevent+=GetNentries(Vec_subfile.at(it));
    SumW+=GetSumW(Vec_subfile.at(it));
  }

  cout<<Nevent<<" "<<SumW<<endl;
   
}
