//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#ifndef EVENTUTILS_H_
#define EVENTUTILS_H_

#include "TLorentzVector.h"
#include "initOniaTree.C"

namespace HI {
  enum TRIGGERBIT {
    HLT_HIL1DoubleMu0_v1 = 0,
    HLT_HIL1DoubleMu0_2HF_v1 = 1,
    HLT_HIL1DoubleMu0_2HF0_v1 = 2,
    HLT_HIL1DoubleMu10_v1 = 3,
    HLT_HIL2DoubleMu0_NHitQ_v2 = 4,
    HLT_HIL2DoubleMu0_NHitQ_2HF_v1 = 5,
    HLT_HIL2DoubleMu0_NHitQ_2HF0_v1 = 6,
    HLT_HIL1DoubleMu0_2HF_Cent30100_v1 = 7,
    HLT_HIL1DoubleMu0_2HF0_Cent30100_v1 = 8,
    HLT_HIL2DoubleMu0_2HF_Cent30100_NHitQ_v1 = 9,
    HLT_HIL1DoubleMu0_Cent30_v1 = 10,
    HLT_HIL2DoubleMu0_2HF0_Cent30100_NHitQ_v1 = 11,
    HLT_HIL2DoubleMu0_Cent30_NHitQ_v1 = 12,
    HLT_HIL2DoubleMu0_Cent30_OS_NHitQ_v1 = 13,
    HLT_HIL3DoubleMu0_Cent30_v1 = 14,
    HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1 = 15,
    HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v1 = 16,
    HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1 = 17,
    HLT_HIL3DoubleMu0_OS_m7to14_v1 = 18,
    HLT_HIL2Mu3_NHitQ10_2HF_v1 = 19,
    HLT_HIL2Mu3_NHitQ10_2HF0_v1 = 20,
    HLT_HIL3Mu3_NHitQ15_2HF_v1 = 21,
    HLT_HIL3Mu3_NHitQ15_2HF0_v1 = 22,
    HLT_HIL2Mu5_NHitQ10_2HF_v1 = 23,
    HLT_HIL2Mu5_NHitQ10_2HF0_v1 = 24,
    HLT_HIL3Mu5_NHitQ15_2HF_v1 = 25,
    HLT_HIL3Mu5_NHitQ15_2HF0_v1 = 26,
    HLT_HIL2Mu7_NHitQ10_2HF_v1 = 27,
    HLT_HIL2Mu7_NHitQ10_2HF0_v1 = 28,
    HLT_HIL3Mu7_NHitQ15_2HF_v1 = 29,
    HLT_HIL3Mu7_NHitQ15_2HF0_v1 = 30,
    HLT_HIL2Mu15_v2 = 31,
    HLT_HIL2Mu15_2HF_v1 = 32,
    HLT_HIL2Mu15_2HF0_v1 = 33,
    HLT_HIL3Mu15_v1 = 34,
    HLT_HIL3Mu15_2HF_v1 = 35,
    HLT_HIL3Mu15_2HF0_v1 = 36,
    HLT_HIL2Mu20_v1 = 37,
    HLT_HIL2Mu20_2HF_v1 = 38,
    HLT_HIL2Mu20_2HF0_v1 = 39,
    HLT_HIL3Mu20_v1 = 40,
    HLT_HIL3Mu20_2HF_v1 = 41,
    HLT_HIL3Mu20_2HF0_v1 = 42
  };

  Double_t findNcoll(int hiBin) {
    const int nbins = 200;
    const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
    return Ncoll[hiBin];
  };
  
  Double_t findNcollAverage(int hiBinLow, int hiBinHigh) {
    Double_t w=0;
    const int nbins = 200;
    const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
    for(int i=hiBinLow; i<hiBinHigh; i++)  w+=Ncoll[i]/(hiBinHigh-hiBinLow);
    return w;
  };
};

namespace PP {
  enum TRIGGERBIT {
    HLT_HIL1DoubleMu0_v1 = 0,
    HLT_HIL1DoubleMu10_v1 = 1,
    HLT_HIL2DoubleMu0_NHitQ_v1 = 2,
    HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1 = 3,
    HLT_HIL3DoubleMu0_OS_m7to14_v1 = 4,
    HLT_HIL2Mu3_NHitQ10_v1 = 5,
    HLT_HIL3Mu3_NHitQ15_v1 = 6,
    HLT_HIL2Mu5_NHitQ10_v1 = 7,
    HLT_HIL3Mu5_NHitQ15_v1 = 8,
    HLT_HIL2Mu7_NHitQ10_v1 = 9,
    HLT_HIL3Mu7_NHitQ15_v1 = 10,
    HLT_HIL2Mu15_v1 = 11,
    HLT_HIL3Mu15_v1 = 12,
    HLT_HIL2Mu20_v1 = 13,
    HLT_HIL3Mu20_v1 = 14
  };
};

namespace RecoQQ {
  void iniBranches(TChain* fChain)
  {
    fChain->SetBranchStatus("HLTriggers",1); 
    fChain->SetBranchStatus("Reco_QQ_trig",1);
    fChain->SetBranchStatus("Reco_QQ_VtxProb",1); 
    fChain->SetBranchStatus("Reco_QQ_mupl_isGoodMuon",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_isGoodMuon",1);
    fChain->SetBranchStatus("Reco_QQ_mupl_nTrkWMea",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_nTrkWMea",1);
    fChain->SetBranchStatus("Reco_QQ_mupl_nPixWMea",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_nPixWMea",1);
    fChain->SetBranchStatus("Reco_QQ_mupl_dxy",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_dxy",1);
    fChain->SetBranchStatus("Reco_QQ_mupl_dz",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_dz",1);

    fChain->SetBranchStatus("Reco_QQ_mupl_nTrkHits",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_nTrkHits",1);
    fChain->SetBranchStatus("Reco_QQ_mupl_normChi2_global",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_normChi2_global",1);
    fChain->SetBranchStatus("Reco_QQ_mupl_normChi2_inner",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_normChi2_inner",1);
    fChain->SetBranchStatus("Reco_QQ_mupl_TrkMuArb",1); 
    fChain->SetBranchStatus("Reco_QQ_mumi_TrkMuArb",1);

  }
    
  Bool_t isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit) 
  {
    Bool_t cond = true;
    cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) ); 
    cond = cond && ( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
    return cond;
  };
  
  Bool_t isGlobalMuonInAccept2011 (TLorentzVector* Muon) 
  {
    return (fabs(Muon->Eta()) < 2.4 &&
	    ((fabs(Muon->Eta()) < 1.0 && Muon->Pt() >= 3.4) ||
	     (1.0 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 1.5 && Muon->Pt() >= 5.8-2.4*fabs(Muon->Eta())) ||
	     (1.5 <= fabs(Muon->Eta()) && Muon->Pt() >= 3.3667-7.0/9.0*fabs(Muon->Eta()))));
  };

  Bool_t areMuonsInAcceptance2011 (Int_t iRecoQQ)
  {
    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
    return ( isGlobalMuonInAccept2011(RecoQQmupl) && isGlobalMuonInAccept2011(RecoQQmumi) );
  };
  
  Bool_t isGlobalMuonInAccept2015 (TLorentzVector* Muon) 
  {
  return (fabs(Muon->Eta()) < 2.4 &&
          ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
           (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
           (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
  };

  Bool_t areMuonsInAcceptance2015 (Int_t iRecoQQ)
  {
    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
    return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
  };  
  
  Bool_t passQualityCuts2015 (Int_t iRecoQQ) 
  {
    Bool_t cond = true;
    
    // cond = cond && (Reco_QQ_mumi_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mumi_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mumi_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 20.);
    
    // cond = cond && (Reco_QQ_mupl_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mupl_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mupl_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 20.);
    
    cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);
    
    return cond;
  }; 

 Bool_t passQualityCuts2011 (Int_t iRecoQQ) 
  {
    Bool_t cond = true;
    
    cond = cond && (Reco_QQ_mumi_nTrkHits[iRecoQQ] > 10);
    cond = cond && (Reco_QQ_mumi_normChi2_global[iRecoQQ] < 20.0);
    cond = cond && (Reco_QQ_mumi_normChi2_inner[iRecoQQ] < 4.0);
    cond = cond && (Reco_QQ_mumi_TrkMuArb[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ]>0);
    cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 3.0);
    cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 15.0);

    cond = cond && (Reco_QQ_mupl_nTrkHits[iRecoQQ] > 10);
    cond = cond && (Reco_QQ_mupl_normChi2_global[iRecoQQ] < 20.0);
    cond = cond && (Reco_QQ_mupl_normChi2_inner[iRecoQQ] < 4.0);
    cond = cond && (Reco_QQ_mupl_TrkMuArb[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ]>0);
    cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 3.0);
    cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 15.0);
    
    cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);
    
    return cond;
  }; 

}

#endif /* EVENTUTILS_H_ */
