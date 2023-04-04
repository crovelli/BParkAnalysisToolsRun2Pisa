#include <TROOT.h>
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include <TFile.h>
#include <TClass.h>

#include <iostream>
#include <cmath>

#include "TLorentzVector.h"

#include "../include/SkimmerPisa.hh"    

using namespace std;

float SkimmerPisa::DeltaR(float eta1,float  phi1, float eta2, float phi2){
  float PI=3.1415972;
  float deta=eta1-eta2;
  float dphi=phi1-phi2;
  if(dphi>PI){
    dphi-=2.0*PI;
  } else if(dphi<=-PI){
    dphi+=2.0*PI;
  }
  return TMath::Sqrt(deta*deta+dphi*dphi);
}

float SkimmerPisa::DeltaPhi(float  phi1, float phi2){
  float PI=3.1415972;
  float dphi=phi1-phi2;
  if(dphi>PI){
    dphi-=2.0*PI;
  } else if(dphi<=-PI){
    dphi+=2.0*PI;
  }
  return dphi;
}

SkimmerPisa::SkimmerPisa(TChain *tree, int isMC )     
  : BParkBase(tree) {        

  sampleID = isMC;         
}

SkimmerPisa::~SkimmerPisa() {

  // output
  outFile_->cd();
  h_entries -> Write();
  h_selection -> Write();
  outTree_->Write();
  outFile_->Close();
}     

void SkimmerPisa::Loop() {

  if (fChain == 0) return;

  // Loop over events
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  cout << "entries : " <<  nentries << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%500==0) cout << jentry << endl;

    // To keep track of the total number of events 
    h_entries->Fill(5);
    
    // Event info     
    theRun   = run;
    theEvent = event;
    theSampleID = sampleID;

    // # Vertices
    nvtx = PV_npvs;
    
    // Energy density
    rho = fixedGridRhoFastjetAll;
    
    // other weights for the dataset
    float perEveW = 1.;
    if (sampleID>0) perEveW = Generator_weight;

    // Events breakdown  
    h_selection->Fill(0.,perEveW);
    

    // ----------------------------------------------------

    // PV must be there
    bool goodPV = false;
    if (PV_x>-999) goodPV = true;
    if (!goodPV) continue;
    h_selection->Fill(1.,perEveW);

    // Trigger 
    int iHLT_Mu7_IP4 = (int)HLT_Mu7_IP4;
    int iHLT_Mu9_IP6 = (int)HLT_Mu9_IP6;
    hlt7 = iHLT_Mu7_IP4;
    hlt9 = iHLT_Mu9_IP6;
    
    // Triggering muons
    bool nTriggerMuon=0;
    int idTrgMu=-1;
    float ptTrgMu=0.;
    trg_muon_pt=0.;
    if(nMuon>0){
      for (u_int iM=0; iM<nMuon; iM++) {
	if(Muon_isTriggering[iM]) {
	  nTriggerMuon=nTriggerMuon+1;
	  if(Muon_pt[iM]>ptTrgMu) {
	    ptTrgMu=Muon_pt[iM];
	    idTrgMu=iM;
	  }
	}
      }
    }
    trg_muon_pt=ptTrgMu;
    

    // B candidates
    if (nBToKEE<=0) continue;
    h_selection->Fill(2.,perEveW);


    // Minimal Bcandidate requirements
    vector<int> goodBs;
    vector<int> goodOrderBs;
    int goodHighestPt = -1;
    float ptGoodHighestPt = -1.;

    for (u_int iB=0; iB<nBToKEE; iB++) {

      // preparing variables for selection
      int ele1_idx = BToKEE_l1Idx[iB];
      int ele2_idx = BToKEE_l2Idx[iB];
      int k_idx    = BToKEE_kIdx[iB];
      
      float ele1_pt = BToKEE_fit_l1_pt[iB];
      float ele2_pt = BToKEE_fit_l2_pt[iB];
      float k_pt    = BToKEE_fit_k_pt[iB];
      
      float ele1_eta = BToKEE_fit_l1_eta[iB];     
      float ele2_eta = BToKEE_fit_l2_eta[iB];    
      float k_eta    = BToKEE_fit_k_eta[iB];    

      float ele1_phi = BToKEE_fit_l1_phi[iB];     
      float ele2_phi = BToKEE_fit_l2_phi[iB];    
      float k_phi    = BToKEE_fit_k_phi[iB];    

      // Electron ID
      float L1id=-1000.;
      float L2id=-1000.;
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==1 ) { // PFPF
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_pfmvaId[ele2_idx];
      }
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==0 ) { // PFLP
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_mvaId[ele2_idx];
      }
      if( Electron_isPF[ele1_idx]==0 && Electron_isPF[ele2_idx]==1 ) { // PFLP
	L2id=Electron_pfmvaId[ele2_idx];
	L1id=Electron_mvaId[ele1_idx];
      }

      float theXySig = BToKEE_l_xy[iB]/BToKEE_l_xy_unc[iB];        

      // Preselection: for PF-LPT remove B candidates which are also PF-PF
      if (Electron_isLowPt[ele1_idx] && Electron_isPFoverlap[ele1_idx]) continue;
      if (Electron_isLowPt[ele2_idx] && Electron_isPFoverlap[ele2_idx]) continue;

      // Preselection proposed by George & Otto on top of nanoAOD presel
      bool commonPresel = false;
      bool etaPresel = fabs(ele1_eta)<2.4 && fabs(ele2_eta)<2.4 && fabs(k_eta)<2.4;
      if (Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==1) { // PFPF
	commonPresel = (BToKEE_svprob[iB]>0.000001) && (BToKEE_fit_cos2D[iB]>0.85) && (theXySig>3.) && (k_pt>0.5) && (ele1_pt)>2.0 && (ele2_pt)>2.0 && (L1id)>-1.5 && (L2id)>-3;
      } else {
	commonPresel = (BToKEE_svprob[iB]>0.000001) && (BToKEE_fit_cos2D[iB]>0.95) && (fabs(BToKEE_k_svip3d[iB])<0.06) && (k_pt>0.5) && (ele1_pt)>2.0 && (ele2_pt>1.0) && (L1id>-2.0) && (L2id>0);
      }
      if (!etaPresel)    continue;
      if (!commonPresel) continue;

      // Order Bs based on pT
      if (BToKEE_fit_pt[iB]>ptGoodHighestPt) {
	ptGoodHighestPt = BToKEE_fit_pt[iB];
	goodHighestPt = iB;	
      }

      goodBs.push_back(iB);
    }
    
    if (goodBs.size()>0) h_selection->Fill(3.,perEveW);    
    selectedBSize = goodBs.size();
    
    // At least one good B candidate
    if (goodBs.size()<=0) continue;
    h_selection->Fill(4.,perEveW);



    // --------------------------------------------------------------------
    // Pickup only the highest pT B
    goodOrderBs.push_back(goodHighestPt);
    if (goodOrderBs.size()<=0 || goodOrderBs.size()>1) cout << "PROBLEM! " << endl;

    for (u_int iB=0; iB<goodOrderBs.size(); iB++) {

      // B-infos for further cuts offline
      int thisB = goodOrderBs[iB];

      // MC truth
      bool isThisAMcB = -1;
      if (sampleID>0) isThisAMcB = isMcBOtto(thisB);            

      float thisBmass   = BToKEE_fit_mass[thisB];
      float thisBpt     = BToKEE_fit_pt[thisB];
      float thisBcos    = BToKEE_fit_cos2D[thisB];
      float thisBcos3d  = BToKEE_fit_cos3D[thisB]; 
      float thisBsvprob = BToKEE_svprob[thisB];
      float thisBxysig  = BToKEE_l_xy[thisB]/BToKEE_l_xy_unc[thisB];
      float thisBxysigWrtPV = BToKEE_l_xy_pv[thisB]/BToKEE_l_xy_unc_pv[thisB]; 

      // Index
      int k_idx    = BToKEE_kIdx[thisB];  
      int ele1_idx = BToKEE_l1Idx[thisB];
      int ele2_idx = BToKEE_l2Idx[thisB];

      // Kine
      float k_pt     = BToKEE_fit_k_pt[thisB];   
      float ele1_pt  = BToKEE_fit_l1_pt[thisB];
      float ele2_pt  = BToKEE_fit_l2_pt[thisB];
      float k_eta    = BToKEE_fit_k_eta[thisB];   
      float ele1_eta = BToKEE_fit_l1_eta[thisB];
      float ele2_eta = BToKEE_fit_l2_eta[thisB];
      float k_phi    = BToKEE_fit_k_phi[thisB];   
      float ele1_phi = BToKEE_fit_l1_phi[thisB];
      float ele2_phi = BToKEE_fit_l2_phi[thisB];

      // Ele-ID
      float L1id=-1000.;
      float L2id=-1000.;
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==1 ) { 
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_pfmvaId[ele2_idx];
      }
      if( Electron_isPF[ele1_idx]==1 && Electron_isPF[ele2_idx]==0 ) { 
	L1id=Electron_pfmvaId[ele1_idx];
	L2id=Electron_mvaId[ele2_idx];
      }
      if( Electron_isPF[ele1_idx]==0 && Electron_isPF[ele2_idx]==1 ) { 
	L2id=Electron_pfmvaId[ele2_idx];
	L1id=Electron_mvaId[ele1_idx];
      }
      
      // Vectors
      TLorentzVector kTLV(0,0,0,0);
      kTLV.SetPtEtaPhiM(k_pt,k_eta,k_phi,0.493);
      TLorentzVector ele1TLV(0,0,0,0);
      ele1TLV.SetPtEtaPhiM(ele1_pt,ele1_eta,ele1_phi,0.000511);
      TLorentzVector ele2TLV(0,0,0,0);
      ele2TLV.SetPtEtaPhiM(ele2_pt,ele2_eta,ele2_phi,0.000511);

      // Vertices
      float kvx = ProbeTracks_vx[k_idx];
      float kvy = ProbeTracks_vy[k_idx];
      float kvz = ProbeTracks_vz[k_idx];
      
      float l1vx = Electron_vx[ele1_idx];
      float l1vy = Electron_vy[ele1_idx];
      float l1vz = Electron_vz[ele1_idx];
      
      float l2vx = Electron_vx[ele2_idx];
      float l2vy = Electron_vy[ele2_idx];
      float l2vz = Electron_vz[ele2_idx];
      
      fit_Bmass.push_back(thisBmass);
      fit_Bpt.push_back(thisBpt);
      fit_Bcos2D.push_back(thisBcos);
      fit_Bcos3D.push_back(thisBcos3d);
      fit_Bsvprob.push_back(thisBsvprob);
      fit_Bxysig.push_back(thisBxysig);
      fit_BxysigWrtPV.push_back(thisBxysigWrtPV);
      if(sampleID>0) bmatchMC.push_back(isThisAMcB);
      // 
      tag_pt.push_back(ele1_pt);
      tag_eta.push_back(ele1_eta);
      tag_phi.push_back(ele1_phi);
      tag_isPF.push_back(Electron_isPF[ele1_idx]);           
      tag_isPFOverlap.push_back(Electron_isPFoverlap[ele1_idx]); 
      tag_isLowPt.push_back(Electron_isLowPt[ele1_idx]);
      tag_id.push_back(L1id);
      tag_iso04_vec.push_back(BToKEE_l1_iso04[thisB]);
      tag_3diso04_vec.push_back(BToKEE_l1_3diso04[thisB]);
      //
      probe_pt.push_back(ele2_pt);
      probe_eta.push_back(ele2_eta);
      probe_phi.push_back(ele2_phi);
      probe_isPF.push_back(Electron_isPF[ele2_idx]);
      probe_isPFOverlap.push_back(Electron_isPFoverlap[ele2_idx]); 
      probe_isLowPt.push_back(Electron_isLowPt[ele2_idx]); 
      probe_id.push_back(L2id);
      probe_iso04_vec.push_back(BToKEE_l2_iso04[thisB]);
      probe_3diso04_vec.push_back(BToKEE_l2_3diso04[thisB]);
      // 
      K_pt.push_back(k_pt);
      K_eta.push_back(k_eta);
      K_phi.push_back(k_phi);
      k_iso04_vec.push_back(BToKEE_k_iso04[thisB]);
      k_3diso04_vec.push_back(BToKEE_k_3diso04[thisB]);
      // 
      mll_fullfit.push_back(BToKEE_mll_fullfit[thisB]);
      mll_raw.push_back(BToKEE_mll_raw[thisB]);
      //
      // b isolation
      b_iso04_vec.push_back(BToKEE_b_iso04[thisB]);
      b_3diso04_vec.push_back(BToKEE_b_3diso04[thisB]);

      // I.P.
      float BToKEE_l1_dxy_sig=(Electron_dxy[ele1_idx]) /Electron_dxyErr[ele1_idx];
      float BToKEE_l2_dxy_sig=(Electron_dxy[ele2_idx]) /Electron_dxyErr[ele2_idx];
      float BToKEE_k_dxy_sig=ProbeTracks_dxyS[k_idx];
      k_svip2d_vec.push_back(BToKEE_k_svip2d[thisB]);
      k_svip3d_vec.push_back(BToKEE_k_svip3d[thisB]);
      k_dxy_sig_vec.push_back(BToKEE_k_dxy_sig);
      tag_dxy_sig_vec.push_back(BToKEE_l1_dxy_sig);
      probe_dxy_sig_vec.push_back(BToKEE_l2_dxy_sig);
      
      if (sampleID>0) {     // MC      
	bool isTagAMcEleFromJPsi   = isMcEleFromJPsi(ele1_idx);
	bool isProbeAMcEleFromJPsi = isMcEleFromJPsi(ele2_idx);
	bool isTagAMcEle   = (Electron_genPartIdx[ele1_idx]>-0.5);
	bool isProbeAMcEle = (Electron_genPartIdx[ele2_idx]>-0.5);

	if(isTagAMcEle){
	  tag_ptMc.push_back(GenPart_pt[Electron_genPartIdx[ele1_idx]]);
	}else{
	  tag_ptMc.push_back(0.);
	}
	if(isProbeAMcEle){
	  probe_ptMc.push_back(GenPart_pt[Electron_genPartIdx[ele2_idx]]);
	}else{
	  probe_ptMc.push_back(0.);
	}
	tag_matchMc.push_back(isTagAMcEle);
	probe_matchMc.push_back(isProbeAMcEle);
	tag_matchMcFromJPsi.push_back(isTagAMcEleFromJPsi);
	probe_matchMcFromJPsi.push_back(isProbeAMcEleFromJPsi);
	
      } else {
	tag_matchMcFromJPsi.push_back(0);  
	probe_matchMcFromJPsi.push_back(0);  
	tag_matchMc.push_back(0);  
	probe_matchMc.push_back(0);  
	tag_ptMc.push_back(0.);
	probe_ptMc.push_back(0.);
      }
      
    } // Loop over good Bs
      

    // At least one tag and one probe
    selectedPairsSize = tag_pt.size();
    if (selectedPairsSize<=0) continue;
    h_selection->Fill(5.,perEveW);

    
    // Filling the output tree
    outTree_->Fill();
    
    // Cleaning all vectors used for the selection
    goodBs.clear();
    goodOrderBs.clear();
    
    // Cleaning all vectors used for the output tree, ready for a new entry
    tag_pt.clear();  
    tag_eta.clear();  
    tag_phi.clear();  
    tag_isPF.clear();  
    tag_isPFOverlap.clear();
    tag_isLowPt.clear();  
    tag_id.clear();
    //
    probe_pt.clear();  
    probe_eta.clear();  
    probe_phi.clear();  
    probe_isPF.clear();  
    probe_isPFOverlap.clear();
    probe_isLowPt.clear();  
    probe_id.clear();  
    //
    if(sampleID>0){
      tag_ptMc.clear();
      probe_ptMc.clear();
      tag_matchMcFromJPsi.clear();  
      tag_matchMc.clear();
      bmatchMC.clear();
      probe_matchMcFromJPsi.clear();  
      probe_matchMc.clear();  
    }  

    mll_fullfit.clear();
    mll_raw.clear();

    fit_Bmass.clear();
    fit_Bpt.clear();
    fit_Bcos2D.clear();
    fit_Bcos3D.clear();
    fit_Bsvprob.clear();    
    fit_Bxysig.clear();    
    fit_BxysigWrtPV.clear();    

    K_pt.clear();
    K_eta.clear();
    K_phi.clear();

    k_svip2d_vec.clear();
    k_svip3d_vec.clear();
    k_dxy_sig_vec.clear();
    tag_dxy_sig_vec.clear();
    probe_dxy_sig_vec.clear();

    k_iso04_vec.clear();
    k_3diso04_vec.clear();
    b_iso04_vec.clear();
    b_3diso04_vec.clear();
    tag_iso04_vec.clear();
    tag_3diso04_vec.clear();
    probe_iso04_vec.clear();
    probe_3diso04_vec.clear();
  }
}

bool SkimmerPisa::isMcB( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx      = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx    = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genGMotherIdx   = GenPart_genPartIdxMother[k_genMotherIdx];
  int k_genPdgId        = GenPart_pdgId[k_genPartIdx];
  int k_genMotherPdgId  = GenPart_pdgId[k_genMotherIdx];
  int k_genGMotherPdgId = GenPart_pdgId[k_genGMotherIdx];

  int ele1_genPartIdx      = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx    = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genGMotherIdx   = GenPart_genPartIdxMother[ele1_genMotherIdx];
  int ele1_genPdgId        = GenPart_pdgId[ele1_genPartIdx];
  int ele1_genMotherPdgId  = GenPart_pdgId[ele1_genMotherIdx];
  int ele1_genGMotherPdgId = GenPart_pdgId[ele1_genGMotherIdx];

  int ele2_genPartIdx      = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx    = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genGMotherIdx   = GenPart_genPartIdxMother[ele2_genMotherIdx];
  int ele2_genPdgId        = GenPart_pdgId[ele2_genPartIdx];
  int ele2_genMotherPdgId  = GenPart_pdgId[ele2_genMotherIdx];
  int ele2_genGMotherPdgId = GenPart_pdgId[ele2_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  bool RK_res1 = abs(ele1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genGMotherPdgId) && (k_genMotherPdgId==ele2_genGMotherPdgId);

  bool RK_res3 = abs(ele1_genMotherPdgId)==521 && abs(k_genMotherPdgId)==521;
  bool RK_res4 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genMotherPdgId) && (k_genMotherPdgId==ele2_genMotherPdgId);

  bool RK_res = okMatch  && (( RK_res1 && RK_res2) || (RK_res3 &&RK_res4)) ;

  return RK_res;
}

bool SkimmerPisa::isMcBOtto( int theB ) {
  
  // taking index
  int ele1_idx = BToKEE_l1Idx[theB];
  int ele2_idx = BToKEE_l2Idx[theB];
  int k_idx    = BToKEE_kIdx[theB];

  // Gen tree
  int k_genPartIdx      = ProbeTracks_genPartIdx[k_idx];  
  int k_genMotherIdx    = GenPart_genPartIdxMother[k_genPartIdx];
  int k_genGMotherIdx   = GenPart_genPartIdxMother[k_genMotherIdx];
  int k_genPdgId        = GenPart_pdgId[k_genPartIdx];
  int k_genMotherPdgId  = GenPart_pdgId[k_genMotherIdx];
  int k_genGMotherPdgId = GenPart_pdgId[k_genGMotherIdx];

  int ele1_genPartIdx      = Electron_genPartIdx[ele1_idx];  
  int ele1_genMotherIdx    = GenPart_genPartIdxMother[ele1_genPartIdx];
  int ele1_genGMotherIdx   = GenPart_genPartIdxMother[ele1_genMotherIdx];
  int ele1_genPdgId        = GenPart_pdgId[ele1_genPartIdx];
  int ele1_genMotherPdgId  = GenPart_pdgId[ele1_genMotherIdx];
  int ele1_genGMotherPdgId = GenPart_pdgId[ele1_genGMotherIdx];

  int ele2_genPartIdx      = Electron_genPartIdx[ele2_idx];  
  int ele2_genMotherIdx    = GenPart_genPartIdxMother[ele2_genPartIdx];
  int ele2_genGMotherIdx   = GenPart_genPartIdxMother[ele2_genMotherIdx];
  int ele2_genPdgId        = GenPart_pdgId[ele2_genPartIdx];
  int ele2_genMotherPdgId  = GenPart_pdgId[ele2_genMotherIdx];
  int ele2_genGMotherPdgId = GenPart_pdgId[ele2_genGMotherIdx];

  // Reco TVector3
  float reco_ele1_pt  = BToKEE_fit_l1_pt[theB];
  float reco_ele2_pt  = BToKEE_fit_l2_pt[theB];
  float reco_k_pt     = BToKEE_fit_k_pt[theB];
  float reco_ele1_eta = BToKEE_fit_l1_eta[theB];
  float reco_ele2_eta = BToKEE_fit_l2_eta[theB];
  float reco_k_eta    = BToKEE_fit_k_eta[theB];
  float reco_ele1_phi = BToKEE_fit_l1_phi[theB];
  float reco_ele2_phi = BToKEE_fit_l2_phi[theB];
  float reco_k_phi    = BToKEE_fit_k_phi[theB];
  // 
  TVector3 reco_ele1(0,0,0);
  TVector3 reco_ele2(0,0,0);
  TVector3 reco_k(0,0,0);
  reco_ele1.SetPtEtaPhi(reco_ele1_pt, reco_ele1_eta, reco_ele1_phi);
  reco_ele2.SetPtEtaPhi(reco_ele2_pt, reco_ele2_eta, reco_ele2_phi);
  reco_k.SetPtEtaPhi(reco_k_pt,reco_k_eta,reco_k_phi);

  // Gen TVector3  
  float gen_ele1_pt  = GenPart_pt[ele1_genPartIdx];
  float gen_ele2_pt  = GenPart_pt[ele2_genPartIdx];      
  float gen_k_pt     = GenPart_pt[k_genPartIdx];      
  float gen_ele1_eta = GenPart_eta[ele1_genPartIdx];
  float gen_ele2_eta = GenPart_eta[ele2_genPartIdx];      
  float gen_k_eta    = GenPart_eta[k_genPartIdx];      
  float gen_ele1_phi = GenPart_phi[ele1_genPartIdx];
  float gen_ele2_phi = GenPart_phi[ele2_genPartIdx];      
  float gen_k_phi    = GenPart_phi[k_genPartIdx];      
  // 
  TVector3 gen_ele1(0,0,0);
  TVector3 gen_ele2(0,0,0);
  TVector3 gen_k(0,0,0);
  gen_ele1.SetPtEtaPhi(gen_ele1_pt, gen_ele1_eta, gen_ele1_phi);
  gen_ele2.SetPtEtaPhi(gen_ele2_pt, gen_ele2_eta, gen_ele2_phi);
  gen_k.SetPtEtaPhi(gen_k_pt,gen_k_eta,gen_k_phi);

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele1_genPartIdx>-0.5 && ele2_genPartIdx>-0.5 && k_genPartIdx>-0.5);
  if (!okMatch) return okMatch;

  float dR1 = gen_ele1.DeltaR(reco_ele1);
  float dR2 = gen_ele2.DeltaR(reco_ele2);
  float dRK = gen_k.DeltaR(reco_k);
  bool RK_dR = (dR1<0.3) && (dR2<0.3) && (dRK<0.3);

  bool RK_res1 = abs(ele1_genMotherPdgId)==443 && abs(k_genMotherPdgId)==521;
  bool RK_res2 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genGMotherPdgId) && (k_genMotherPdgId==ele2_genGMotherPdgId);

  bool RK_res3 = abs(ele1_genMotherPdgId)==521 && abs(k_genMotherPdgId)==521;
  bool RK_res4 = (ele1_genMotherPdgId==ele2_genMotherPdgId) && (k_genMotherPdgId==ele1_genMotherPdgId) && (k_genMotherPdgId==ele2_genMotherPdgId);

  bool RK_res = okMatch && RK_dR && (( RK_res1 && RK_res2) || (RK_res3 &&RK_res4)) ;

  return RK_res;
}

bool SkimmerPisa::isMcEleFromJPsi( int ele_idx ) {

  // Gen tree
  int ele_genPartIdx      = Electron_genPartIdx[ele_idx];  
  int ele_genMotherIdx    = GenPart_genPartIdxMother[ele_genPartIdx];
  int ele_genGMotherIdx   = GenPart_genPartIdxMother[ele_genMotherIdx];
  int ele_genPdgId        = GenPart_pdgId[ele_genPartIdx];
  int ele_genMotherPdgId  = GenPart_pdgId[ele_genMotherIdx];
  int ele_genGMotherPdgId = GenPart_pdgId[ele_genGMotherIdx];

  // B -> K J/psi(ll) at gen level
  bool okMatch = (ele_genPartIdx>-0.5) && (abs(ele_genMotherPdgId)==443);

  return okMatch;
}

void SkimmerPisa::PrepareOutputs(std::string filename) 
{
  _datasetName=filename;

  std::string outname = _datasetName+".root";    
  cout << "output: " << outname << endl;
  outFile_ = new TFile(outname.c_str(),"RECREATE");

  bookOutputTree();
  bookOutputHistos();
};


void SkimmerPisa::bookOutputTree() 
{
  outTree_ = new TTree("TaPtree", "TaPtree");
  
  cout << "Booking tree" << endl;
  
  outTree_->Branch("theRun", &theRun, "theRun/I");    
  outTree_->Branch("theEvent", &theEvent, "theEvent/I");    
  outTree_->Branch("nvtx", &nvtx, "nvtx/I");    
  outTree_->Branch("sampleID", &sampleID, "sampleID/I");    
  outTree_->Branch("rho", &rho, "rho/F");    

  outTree_->Branch("hlt7", &hlt7, "hlt7/I");
  outTree_->Branch("hlt9", &hlt9, "hlt9/I");
  outTree_->Branch("trg_muon_pt", &trg_muon_pt, "trg_muon_pt/F");

  outTree_->Branch("selectedBSize",  &selectedBSize,  "selectedBSize/I");   
  
  outTree_->Branch("tag_pt", "std::vector<float>", &tag_pt);  
  outTree_->Branch("tag_eta", "std::vector<float>", &tag_eta);  
  outTree_->Branch("tag_phi", "std::vector<float>", &tag_phi);  
  outTree_->Branch("tag_isPF", "std::vector<bool>", &tag_isPF);  
  outTree_->Branch("tag_isPFOverlap", "std::vector<bool>", &tag_isPFOverlap);  
  outTree_->Branch("tag_isLowPt", "std::vector<bool>", &tag_isLowPt);  
  outTree_->Branch("tag_id", "std::vector<float>", &tag_id);  

  outTree_->Branch("probe_pt", "std::vector<float>", &probe_pt);  
  outTree_->Branch("probe_eta", "std::vector<float>", &probe_eta);  
  outTree_->Branch("probe_phi", "std::vector<float>", &probe_phi);  
  outTree_->Branch("probe_isPF", "std::vector<bool>", &probe_isPF);  
  outTree_->Branch("probe_isPFOverlap", "std::vector<bool>", &probe_isPFOverlap);  
  outTree_->Branch("probe_isLowPt", "std::vector<bool>", &probe_isLowPt);  
  outTree_->Branch("probe_id", "std::vector<float>", &probe_id);  

  if(sampleID>0){
    outTree_->Branch("tag_ptMc", "std::vector<float>", &tag_ptMc);
    outTree_->Branch("tag_matchMcFromJPsi", "std::vector<bool>", &tag_matchMcFromJPsi);  
    outTree_->Branch("tag_matchMc", "std::vector<bool>", &tag_matchMc);  
  }

  if(sampleID>0){
    outTree_->Branch("probe_ptMc", "std::vector<float>", &probe_ptMc);
    outTree_->Branch("probe_matchMcFromJPsi", "std::vector<bool>", &probe_matchMcFromJPsi);  
    outTree_->Branch("probe_matchMc", "std::vector<bool>", &probe_matchMc);  
  }

  outTree_->Branch("K_pt", "std::vector<float>", &K_pt);  
  outTree_->Branch("K_eta", "std::vector<float>", &K_eta);  
  outTree_->Branch("K_phi", "std::vector<float>", &K_phi);  

  outTree_->Branch("mll_fullfit", "std::vector<float>", &mll_fullfit);  
  outTree_->Branch("mll_raw", "std::vector<float>", &mll_raw);  

  outTree_->Branch("fit_Bmass",   "std::vector<float>", &fit_Bmass);  
  outTree_->Branch("fit_Bpt",     "std::vector<float>", &fit_Bpt);  
  outTree_->Branch("fit_Bcos2D",  "std::vector<float>", &fit_Bcos2D);  
  outTree_->Branch("fit_Bcos3D",  "std::vector<float>", &fit_Bcos3D);  
  outTree_->Branch("fit_Bsvprob", "std::vector<float>", &fit_Bsvprob);  
  outTree_->Branch("fit_Bxysig",  "std::vector<float>", &fit_Bxysig);  
  outTree_->Branch("fit_BxysigWrtPV",  "std::vector<float>", &fit_BxysigWrtPV);  

  outTree_->Branch("bmatchMC", "std::vector<bool>", &bmatchMC);  

  outTree_->Branch("k_svip2d", "std::vector<float>", &k_svip2d_vec);
  outTree_->Branch("k_svip3d", "std::vector<float>", &k_svip3d_vec);
  outTree_->Branch("k_dxy_sig", "std::vector<float>", &k_dxy_sig_vec);
  outTree_->Branch("tag_dxy_sig", "std::vector<float>", &tag_dxy_sig_vec);
  outTree_->Branch("probe_dxy_sig", "std::vector<float>", &probe_dxy_sig_vec);

  outTree_->Branch("b_iso04",     "std::vector<float>", &b_iso04_vec);
  outTree_->Branch("b_3diso04",   "std::vector<float>", &b_3diso04_vec);
  outTree_->Branch("k_iso04",     "std::vector<float>", &k_iso04_vec);
  outTree_->Branch("k_3diso04",   "std::vector<float>", &k_3diso04_vec);
  outTree_->Branch("tag_iso04",   "std::vector<float>", &tag_iso04_vec);
  outTree_->Branch("tag_3diso04", "std::vector<float>", &tag_3diso04_vec);
  outTree_->Branch("probe_iso04", "std::vector<float>", &probe_iso04_vec);
  outTree_->Branch("probe_3diso04", "std::vector<float>", &probe_3diso04_vec);

  outTree_->Branch("selectedPairsSize",  &selectedPairsSize,  "selectedPairsSize/I");   
}

void SkimmerPisa::bookOutputHistos() 
{
  cout << "Booking histos" << endl;
  //
  h_entries   = new TH1F("h_entries",  "Number of entries",   3,  3.5, 6.5);
  h_selection = new TH1F("h_selection","Selection breakdown", 6, -0.5, 5.5);
}

