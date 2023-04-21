#ifndef SkimmerPisa_h
#define SkimmerPisa_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <string>
#include <TLorentzVector.h>

#include "./BParkBase.h"

using namespace std;

class SkimmerPisa : public BParkBase{
public:

  //! constructor
  SkimmerPisa(TChain *tree=0, int isMC=0 );
  //! destructor
  virtual ~SkimmerPisa();
  //! loop over events
  void Loop();
  void PrepareOutputs(std::string filename);             
  float DeltaR(float x, float y,float z,float w);
  float DeltaPhi(float x, float y);

private:
  
  // Analysis methods
  bool isMcB( int myB );
  int isMcBOtto( int myB );
  bool isMcEleFromJPsi (int myEle );
  void bookOutputTree();
  void bookOutputHistos();

  // settings
  int sampleID;

  // ---- outputs
  TFile* outFile_;
  TTree* outTree_;
  TH1F* h_entries;
  TH1F* h_selection;

  // dataset name
  std::string _datasetName;      

  //---output tree branches variables
  int   theRun;
  int   theEvent;
  int   nvtx;
  int   theSampleID;
  float rho;
  int   hlt7;
  int   hlt9;
  float trg_muon_pt;
  // 
  int selectedBSize;
  //
  vector <float> tag_pt={};   
  vector <float> tag_eta={};
  vector <float> tag_phi={};
  vector <bool> tag_isPF={};
  vector <bool> tag_isPFOverlap={};
  vector <bool> tag_isLowPt={};
  vector <float> tag_id={};
  //
  vector <bool>  tag_matchMcFromJPsi={};
  vector <bool>  tag_matchMc={};
  vector <float> tag_ptMc={};
  //
  vector <float> probe_pt={}; 
  vector <float> probe_eta ={};
  vector <float> probe_phi={};
  vector <bool> probe_isPF ={};
  vector <bool> probe_isPFOverlap ={};
  vector <bool> probe_isLowPt ={};
  vector <float> probe_id={};        
  //
  vector <bool>  probe_matchMcFromJPsi={};
  vector <bool>  probe_matchMc={};
  vector <float> probe_ptMc={};
  //
  vector <float> K_pt={}; 
  vector <float> K_eta={}; 
  vector <float> K_phi={}; 
  //
  vector <float> mll_fullfit={};
  vector <float> mll_raw={};
  //
  vector <float> fit_Bmass={};   
  vector <float> fit_Bpt={};   
  vector <float> fit_Bcos2D={};   
  vector <float> fit_Bcos3D={};   
  vector <float> fit_Bsvprob={};
  vector <float> fit_Bxysig={};
  vector <float> fit_BxysigWrtPV={};
  vector <int> bmatchMC={};   
  //
  vector <float> k_svip2d_vec={};
  vector <float> k_svip3d_vec={};
  vector <float> k_dxy_sig_vec={}; 
  vector <float> tag_dxy_sig_vec={};
  vector <float> probe_dxy_sig_vec={};

  vector <float> tag_iso04_vec={};
  vector <float> tag_iso04_dca_vec={};
  vector <float> tag_iso04_dca_tight_vec={};
  vector <float> tag_3diso04_vec={};
  vector <float> tag_3diso04_dca_vec={};
  vector <float> tag_3diso04_dca_tight_vec={};
  vector <int> tag_ntrk_iso04_vec={};
  vector <int> tag_ntrk_iso04_dca_vec={};
  vector <int> tag_ntrk_iso04_dca_tight_vec={};

  vector <float> probe_iso04_vec={};
  vector <float> probe_iso04_dca_vec={};
  vector <float> probe_iso04_dca_tight_vec={};
  vector <float> probe_3diso04_vec={};
  vector <float> probe_3diso04_dca_vec={};
  vector <float> probe_3diso04_dca_tight_vec={};
  vector <int> probe_ntrk_iso04_vec={};
  vector <int> probe_ntrk_iso04_dca_vec={};
  vector <int> probe_ntrk_iso04_dca_tight_vec={};

  vector <float> k_iso04_vec={};
  vector <float> k_iso04_dca_vec={};
  vector <float> k_iso04_dca_tight_vec={};
  vector <float> k_3diso04_vec={};
  vector <float> k_3diso04_dca_vec={};
  vector <float> k_3diso04_dca_tight_vec={};
  vector <int> k_ntrk_iso04_vec={};
  vector <int> k_ntrk_iso04_dca_vec={};
  vector <int> k_ntrk_iso04_dca_tight_vec={};

  vector <float> b_iso04_vec={};
  vector <float> b_iso04_dca_vec={};
  vector <float> b_iso04_dca_tight_vec={};
  vector <float> b_3diso04_vec={};
  vector <float> b_3diso04_dca_vec={};
  vector <float> b_3diso04_dca_tight_vec={};
  vector <int> b_ntrk_iso04_vec={};
  vector <int> b_ntrk_iso04_dca_vec={};
  vector <int> b_ntrk_iso04_dca_tight_vec={};

  int selectedPairsSize;
};

#endif
