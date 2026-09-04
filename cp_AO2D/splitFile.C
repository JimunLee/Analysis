#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <map>

using namespace std;

void split(){
  TFile *fin = TFile::Open("/Users/jimun/Analysis/kstar-in-jets/Results_rootfile/al_skimmed.root");
  
  map<int, TString> name = {
    {48527, "JetChHighPt100_al"},
    {45532, "JetChHighPt150_al"}
  };
  
  for (auto& entry:name) {
    int id = entry.first;
    TString outName = entry.second;

    TDirectory *dir = (TDirectory*)fin->Get(Form("kstar-in-o-o_id%d", id));

    if (!dir)
      continue;

    TFile *fout = new TFile(outName + ".root", "RECREATE");
    fout->cd();

    TDirectory *newdir = fout->mkdir("kstar-in-o-o");
    newdir->cd();
    
    TIter nextkey(dir->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();

      newdir->cd();
      obj->Write();
    }
    
    fout->Write();
    fout->Close();
  }

  fin->Close();
};
