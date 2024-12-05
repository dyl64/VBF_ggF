using namespace std;

bool Nan4MomVec(vector<float> fatjetM){
  bool NanIn4Mom = false;
  for (int i = 0; i < fatjetM.size(); i++){
    if (fatjetM.at(i) != fatjetM.at(i)) {
      NanIn4Mom = true;
      return NanIn4Mom;
    }
  }
  return NanIn4Mom;
}

// Selection which requires events in which:
// - 2 fatjet
// - t1 OR t2
// - pT_lead > pT1 & pT_sublead > pT2
bool FirstEventSelection_oldConf(vector<int> list_reordered, bool t1, bool t2, vector<float> fatjet_pt){
  bool IsGood = false;
  if (t1 || t2){
    if (fatjet_pt.size()>=2) {
      if ((fatjet_pt.at(list_reordered.at(0))>480.) && (fatjet_pt.at(list_reordered.at(1))>250.)){
        IsGood = true;
      }
    }
  }
  return IsGood;
}


int FirstEventSelection(vector<int> list_reordered, bool t1_pt, bool t2_ptm, vector<float> fatjet_pt, vector<float> fatjet_m){
  int IsGood = 0;
  if (fatjet_pt.size()>=2){
    if (t1_pt){
      for (unsigned int idx = 0; idx < fatjet_pt.size(); idx++){
        if (fatjet_pt.at(list_reordered.at(idx))>450.){
          IsGood = 1;
        }
      }
    }
    if (t2_ptm){
      for (unsigned int idx = 0; idx < fatjet_pt.size(); idx++){
        if ((fatjet_pt.at(list_reordered.at(idx))>410) && (fatjet_m.at(list_reordered.at(idx))>40)){
          IsGood = 2;
        }
      }
    }
  } 
  return IsGood;
}

int FirstEventSelection_ST(vector<int> list_reordered, bool trigger, vector<float> fatjet_pt, vector<float> fatjet_m, float p_th, float m_th){
  int IsGood = 0;
  //  cout << list_reordered.at(0) << " " << trigger << " " << fatjet_pt.at(0) << " " << fatjet_m.at(0) << " " << p_th << " " << m_th << endl;
  
  if (fatjet_pt.size()>=2 && list_reordered.size()){
    if (trigger || 1){
      for (unsigned int idx = 0; idx < fatjet_pt.size(); idx++){
        if (fatjet_pt.at(list_reordered.at(idx))>p_th && fatjet_m.at(list_reordered.at(idx))>m_th){
          IsGood = 1;
        }
      }
    }
  }
  return IsGood;
}

// "Boosted" selection.
// Is the fatjet boosted= 2m/pT<1.
bool IsBoosted(float fatjet_pt, float fatjet_m){
  bool IsGood = false;
  if (2*fatjet_m/fatjet_pt<1.) {
    IsGood = true;
  }
  return IsGood;
}

std::vector<int> ListSortedFJ(vector<float> fatjet_pt){
  struct sort_pred { 
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
  };
  std::vector<std::pair<int, double> > vec_idxpt;
  for (unsigned int i = 0; i<fatjet_pt.size(); i++) {
    std::pair <int,double> p;
    p.first = i;
    p.second = fatjet_pt.at(i);
    vec_idxpt.push_back(p);
  }
  std::sort(vec_idxpt.begin(), vec_idxpt.end(), sort_pred());

  std::vector<int> list_ordered_fj;
  for (unsigned int j = 0; j<vec_idxpt.size(); j++) {
    list_ordered_fj.push_back(vec_idxpt.at(j).first);
  } 
  return list_ordered_fj;
}

// Vector of trackjet ordered in pt. It's a vector of pairs [idx_i,pt_i] sorted in pt_i.
std::vector<std::pair<int, double> > GetSortVRjets(int fj_idx, vector<int> indices ,vector<float> trackjet_pt, double thr){
  std::vector<std::pair<int, double> > vec_idxpt;
  struct sort_pred { 
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
  };
  for (unsigned int idx = 0; idx < indices.size(); idx++){
    if (indices.at(idx) == fj_idx){
      if (trackjet_pt.at(idx)<thr) continue;
      std::pair <int,double> p;
      p.first = idx;
      p.second = trackjet_pt.at(idx);
      
      vec_idxpt.push_back(p);
    }
  }
  std::sort(vec_idxpt.begin(), vec_idxpt.end(), sort_pred());
  return vec_idxpt;
}

// Vector of trackjet ordered in FTAG. It's a vector of pairs [idx_i,FTAG] sorted in FTAG.
std::vector<std::pair<int, double> > GetSortVRjets_FTAG(int fj_idx, vector<int> indices ,vector<float> trackjet_pt ,vector<int> trackjet_FTAG, double thr){
  std::vector<std::pair<int, double> > vec_idxFTAG;
  struct sort_pred { 
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
  };
  for (unsigned int idx = 0; idx < indices.size(); idx++){
    if (indices.at(idx) == fj_idx){
      if (trackjet_pt.at(idx)<thr) continue;
      std::pair <int,double> p;
      p.first = idx;
      p.second = trackjet_FTAG.at(idx);
      
      vec_idxFTAG.push_back(p);
    }
  }
  std::sort(vec_idxFTAG.begin(), vec_idxFTAG.end(), sort_pred());
  return vec_idxFTAG;
}

// "Resolved" selection.
// Is the fatjet such that:
// - n_trackjet >=2
// - DR(lead,sublead)/R>1

bool HasTwoTracksandIsResolved_eachcouple(int fj_idx, float fatjet_pt, vector<int> indices ,vector<float> trackjet_px, vector<float> trackjet_py, vector<float> trackjet_pz, vector<float> trackjet_pt){
  bool HasTwoTracksandIsResolved = true;

  std::vector<std::pair<int, double> > vec_idxpt_5gev;
  // list of trackjets in the fatjet with pT > 5 GeV, pT ordered
  vec_idxpt_5gev = GetSortVRjets(fj_idx, indices, trackjet_pt, 5.);

  std::vector<std::pair<int, double> > vec_idxpt_10gev;
  // list of trackjets in the fatjet with pT > 10 GeV, pT ordered
  vec_idxpt_10gev = GetSortVRjets(fj_idx, indices, trackjet_pt, 10.);

  if (vec_idxpt_10gev.size()<2){
    HasTwoTracksandIsResolved = false;
    return HasTwoTracksandIsResolved;
  }
  
  // just for leading and subleading of the list more then 10 GeV (the only ones inspected for b-tagging)

  for (unsigned int i = 0; i < 2; i++){
    // all the list more then 5 GeV
    for (unsigned int j = 0; j < vec_idxpt_5gev.size(); j++){
      // check they are not the same
      if (vec_idxpt_10gev.at(i).first == vec_idxpt_5gev.at(j).first) continue;
      TVector3 trackjet_i;
      TVector3 trackjet_j;
      trackjet_i.SetXYZ(trackjet_px.at(vec_idxpt_10gev.at(i).first), trackjet_py.at(vec_idxpt_10gev.at(i).first), trackjet_pz.at(vec_idxpt_10gev.at(i).first));
      trackjet_j.SetXYZ(trackjet_px.at(vec_idxpt_5gev.at(j).first), trackjet_py.at(vec_idxpt_5gev.at(j).first), trackjet_pz.at(vec_idxpt_5gev.at(j).first));
      float R_L = min(max(0.02,min(0.4,30./vec_idxpt_10gev.at(i).second)),max(0.02,min(0.4,30./vec_idxpt_5gev.at(j).second)));
      float DeltaR = trackjet_i.DeltaR(trackjet_j);
      if (DeltaR/R_L<1.){
        HasTwoTracksandIsResolved = false;
        return HasTwoTracksandIsResolved;
      }
    }
  }
  return HasTwoTracksandIsResolved;
}

// Combines the two previous selection in one. Takes a list of fatjet indices and returns the list of indices which pass the "Boosted" and "Resolved" selection
std::vector<int> fjIdxAftSel(vector<int> list_reordered, vector<float> fatJetPt, vector<float> fatJetM, vector<int> vrJetIdFatJet, vector<float> vrJetPx, vector<float> vrJetPy, vector<float> vrJetPz, vector<float> vrJetPt){

  std::vector<int> fatjet_aft_presel;
  for (unsigned int nfj_i = 0; nfj_i < list_reordered.size(); nfj_i++){

    if (!IsBoosted(fatJetPt.at(list_reordered.at(nfj_i)), fatJetM.at(list_reordered.at(nfj_i)))) continue;
    if (!HasTwoTracksandIsResolved_eachcouple(list_reordered.at(nfj_i),fatJetPt.at(list_reordered.at(nfj_i)),vrJetIdFatJet, vrJetPx, vrJetPy, vrJetPz, vrJetPt)) continue;
    std::vector<std::pair<int, double> > vec_idxpt_10gev;
    vec_idxpt_10gev = GetSortVRjets(list_reordered.at(nfj_i), vrJetIdFatJet, vrJetPt, 10.);
    if (vec_idxpt_10gev.size()<2) continue;
    fatjet_aft_presel.push_back(list_reordered.at(nfj_i));
  }
  return fatjet_aft_presel;
}

// Combines the two previous selection in one. Takes a list of fatjet indices and returns the list of indices which pass the "Boosted" and "Resolved" selection
std::vector<int> fjIdxAftSel2(vector<int> list_reordered, vector<float> fatJetPt, vector<float> fatJetM){

  std::vector<int> fatjet_aft_presel;
  for (unsigned int nfj_i = 0; nfj_i < list_reordered.size(); nfj_i++){

    if (!IsBoosted(fatJetPt.at(list_reordered.at(nfj_i)), fatJetM.at(list_reordered.at(nfj_i)))) continue;
    fatjet_aft_presel.push_back(list_reordered.at(nfj_i));
  }
  return fatjet_aft_presel;
}

vector<vector<int>> idxFJcandVRjets_newSR_ST(bool trigger, vector<float> fatJetPt, vector<float> fatJetM, vector<int> vrJetIdFatJet, vector<float> vrJetPx, vector<float> vrJetPy, vector<float> vrJetPz,vector<float> vrJetPt, float pt, float m){
  
  // first component is the FJ index, second and third component are the VRjet indices which need to be inspected for btagging
  
  vector<vector<int>> idx_final;
  
  // pT reordered largeR indices (possible effects of late calibration cancelled)
  vector<int> list_reorderd_fj;
  list_reorderd_fj = ListSortedFJ(fatJetPt);

  // trigger fired in this event: 0 None, 1 trigger_pt440, 2 trigger_pt390_m30, 2 trigger_pt440 & trigger_pt390_m30
  // >= 2 fatjet required and offline threshold applied

  int trigger_fired = FirstEventSelection_ST(list_reorderd_fj, trigger, fatJetPt, fatJetM, pt, m);
  
  //cout << "trigger_fired:\t" << trigger_fired << endl;
  if (trigger_fired == 0) return idx_final;

  // pT reordered list of largeR jet with have two tracks, 2m/pT<1, DR/R<1
  vector<int> fatjet_aft_presel = fjIdxAftSel(list_reorderd_fj, fatJetPt, fatJetM, vrJetIdFatJet, vrJetPx, vrJetPy, vrJetPz, vrJetPt);

  if (fatjet_aft_presel.size()<1) return idx_final;
 
  for (int idx = 0; idx < fatjet_aft_presel.size(); idx ++){
    //only leading or subleading
    if (idx > 1) continue;
    
    // get sorted VRjet belonging to this largeRjet
    vector<pair<int,double>> VRjet = GetSortVRjets(fatjet_aft_presel.at(idx), vrJetIdFatJet, vrJetPt, 10.);

    // create a vector (idx_fj,idx_VR_lead,idx_VR_sublead)
    vector<int> idx_final_1component;
    idx_final_1component.push_back(fatjet_aft_presel.at(idx));
    idx_final_1component.push_back(VRjet.at(0).first);
    idx_final_1component.push_back(VRjet.at(1).first);
    idx_final.push_back(idx_final_1component);
  }
  return idx_final;
}

vector<vector<int>> idxFJcand_newSR_ST(bool trigger, vector<float> fatJetPt, vector<float> fatJetM, float pt, float m){
  
  // first component is the FJ index, second and third component are the VRjet indices which need to be inspected for btagging
  
  vector<vector<int>> idx_final;
  
  // pT reordered largeR indices (possible effects of late calibration cancelled)
  vector<int> list_reorderd_fj;
  list_reorderd_fj = ListSortedFJ(fatJetPt);

  // trigger fired in this event: 0 None, 1 trigger_pt440, 2 trigger_pt390_m30, 2 trigger_pt440 & trigger_pt390_m30
  // >= 2 fatjet required and offline threshold applied

  int trigger_fired = FirstEventSelection_ST(list_reorderd_fj, trigger, fatJetPt, fatJetM, pt, m);
  
  //cout << "trigger_fired:\t" << trigger_fired << endl;
  if (trigger_fired == 0) return idx_final;

  // pT reordered list of largeR jet with have two tracks, 2m/pT<1, DR/R<1
  vector<int> fatjet_aft_presel = fjIdxAftSel2(list_reorderd_fj, fatJetPt, fatJetM);

  if (fatjet_aft_presel.size()<1) return idx_final;
 
  for (int idx = 0; idx < fatjet_aft_presel.size(); idx ++){
    //only leading or subleading
    if (idx > 1) continue;
    
    // create a vector (idx_fj)
    vector<int> idx_final_1component;
    idx_final_1component.push_back(fatjet_aft_presel.at(idx));
    idx_final.push_back(idx_final_1component);
  }
  return idx_final;
}

vector<vector<int>> idxFJcandVRjets_newSR(bool HLT_j440_a10t_lcw_jes_L1J100MC,bool HLT_j390_a10t_lcw_jes_30smcINF_L1J100MC, vector<float> fatJetPt, vector<float> fatJetM, vector<int> vrJetIdFatJet, vector<float> vrJetPx, vector<float> vrJetPy, vector<float> vrJetPz,vector<float> vrJetPt){
  // first component is the FJ index, second and third component are the VRjet indices which need to be inspected for btagging
  vector<vector<int>> idx_final;
  // pT reordered largeR indices (possible effects of late calibration cancelled)
  vector<int> list_reorderd_fj;
  list_reorderd_fj = ListSortedFJ(fatJetPt);
  
  // trigger fired in this event: 0 None, 1 trigger_pt440, 2 trigger_pt390_m30, 2 trigger_pt440 & trigger_pt390_m30
  // >= 2 fatjet required and offline threshold applied
  int trigger_fired = FirstEventSelection(list_reorderd_fj, HLT_j440_a10t_lcw_jes_L1J100MC, HLT_j390_a10t_lcw_jes_30smcINF_L1J100MC, fatJetPt, fatJetM);
  if (trigger_fired == 0) return idx_final;

  // pT reordered list of largeR jet with have two tracks, 2m/pT<1, DR/R<1
  vector<int> fatjet_aft_presel = fjIdxAftSel(list_reorderd_fj, fatJetPt, fatJetM, vrJetIdFatJet, vrJetPx, vrJetPy, vrJetPz, vrJetPt);
  
  if (fatjet_aft_presel.size()<1) return idx_final;

  for (int idx = 0; idx < fatjet_aft_presel.size(); idx ++){
    if (idx > 1) continue; 
    if (fatJetPt.at(fatjet_aft_presel.at(idx))<250. || fatJetM.at(fatjet_aft_presel.at(idx))<40.) continue;
    
    // get sorted VRjet belonging to this largeRjet
    vector<pair<int,double>> VRjet = GetSortVRjets(fatjet_aft_presel.at(idx), vrJetIdFatJet, vrJetPt, 10.);

    // create a vector (idx_fj,idx_VR_lead,idx_VR_sublead)
    vector<int> idx_final_1component;
    idx_final_1component.push_back(fatjet_aft_presel.at(idx));
    idx_final_1component.push_back(VRjet.at(0).first);
    idx_final_1component.push_back(VRjet.at(1).first);
    idx_final.push_back(idx_final_1component);
  }
  return idx_final;
}

pair<int,double> newSR_fatjet(vector<vector<int>> indices, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  double SF = 1.;
  pair<int,double> p;
  for (int i = 0; i < indices.size(); i++ ){
    if ((btag_IS77.at(indices.at(i).at(1))==1) && (btag_IS77.at(indices.at(i).at(2))==1)){
      SF = SF*btag_SF77.at(indices.at(i).at(1))*btag_SF77.at(indices.at(i).at(2));
      idx_btag = i;
      p.first = indices.at(i).at(0);
      p.second = SF;
      return p;
    }
    else{
      SF = SF*btag_SF77.at(indices.at(i).at(1))*btag_SF77.at(indices.at(i).at(2));
    }
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> newCR_fatjet(vector<vector<int>> indices, double SF, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  SF = 1.;
  pair<int,double> p;
  for (int i = 0; i < indices.size(); i++ ){
    if ((btag_IS77.at(indices.at(i).at(1))==0) && (btag_IS77.at(indices.at(i).at(2))==0)){
      SF = SF*btag_SF77.at(indices.at(i).at(1))*btag_SF77.at(indices.at(i).at(2));
      p.first = indices.at(i).at(0);
      p.second = SF;
      return p;
    }
    else{
      SF = SF*btag_SF77.at(indices.at(i).at(1))*btag_SF77.at(indices.at(i).at(2));
    }
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> newCR_fatjet_mixWP(vector<vector<int>> indices, float SF, vector<int> btag_IS77, vector<float> btag_SF77, vector<int> btag_IS85, vector<float> btag_SF85){
  int idx_btag = -1;
  //SF = 1.;
  pair<int,double> p;
  for (int i = 0; i < indices.size(); i++ ){
    if (((btag_IS77.at(indices.at(i).at(1))==0) && (btag_IS77.at(indices.at(i).at(2))==0)) && 
      ((btag_IS85.at(indices.at(i).at(1))==1) || (btag_IS85.at(indices.at(i).at(2))==1))){
      SF = SF*btag_SF77.at(indices.at(i).at(1))*btag_SF77.at(indices.at(i).at(2))*btag_SF85.at(indices.at(i).at(1))*btag_SF85.at(indices.at(i).at(2));
      p.first = indices.at(i).at(0);
      p.second = SF;
      return p;
    }
    else{
      SF = SF*btag_SF77.at(indices.at(i).at(1))*btag_SF77.at(indices.at(i).at(2))*btag_SF85.at(indices.at(i).at(1))*btag_SF85.at(indices.at(i).at(2));
    }
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> newSR_fatjet_DATA(vector<vector<int>> indices, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  float SF = 1.;
  pair<int,double> p;
  for (int i = 0; i < indices.size(); i++ ){
    if ((btag_IS77.at(indices.at(i).at(1))==1) && (btag_IS77.at(indices.at(i).at(2))==1)){
      p.first = indices.at(i).at(0);
      p.second = SF;
      return p;
    }
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> newCR_fatjet_DATA(vector<vector<int>> indices, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  float SF = 1.;
  pair<int,double> p;
  for (int i = 0; i < indices.size(); i++ ){
    if ((btag_IS77.at(indices.at(i).at(1))==0) && (btag_IS77.at(indices.at(i).at(2))==0)){
      p.first = indices.at(i).at(0);
      p.second = SF;
      return p;
    }
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> newCR_fatjet_mixWP_DATA(vector<vector<int>> indices, float SF, vector<int> btag_IS77, vector<float> btag_SF77, vector<int> btag_IS85, vector<float> btag_SF85){
  int idx_btag = -1;
  pair<int,double> p;
  for (int i = 0; i < indices.size(); i++ ){
    if (((btag_IS77.at(indices.at(i).at(1))==0) && (btag_IS77.at(indices.at(i).at(2))==0)) && 
      ((btag_IS85.at(indices.at(i).at(1))==1) || (btag_IS85.at(indices.at(i).at(2))==1))){
      p.first = indices.at(i).at(0);
      p.second = SF;
      return p;
    }
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

vector<vector<int>> idxFJcandVRjets_newtrig_oldSR(bool HLT_j440_a10t_lcw_jes_L1J100MC,bool HLT_j390_a10t_lcw_jes_30smcINF_L1J100MC, vector<float> fatJetPt, vector<float> fatJetM, vector<int> vrJetIdFatJet, vector<float> vrJetPx, vector<float> vrJetPy, vector<float> vrJetPz,vector<float> vrJetPt){
  // first component is the FJ index
  // second and third component are the VRjet indices which need to be inspected for btagging
  vector<vector<int>> idx_final;
  // pT reordered largeR indices (possible effects of late calibration cancelled)
  vector<int> list_reorderd_fj;
  list_reorderd_fj = ListSortedFJ(fatJetPt);
  // trigger fired in this event: 0 None, 1 trigger_pt440, 2 trigger_pt390_m30, 2 trigger_pt440 & trigger_pt390_m30
  // >= 2 fatjet required and offline threshold applied
  int trigger_fired = FirstEventSelection(list_reorderd_fj, HLT_j440_a10t_lcw_jes_L1J100MC, HLT_j390_a10t_lcw_jes_30smcINF_L1J100MC, fatJetPt, fatJetM);
  if (trigger_fired == 0) return idx_final;

  float mass_cut_cand = 0.;
  float pt_cut_cand = 0.;
  if(trigger_fired==1){
    pt_cut_cand = 460.;
    mass_cut_cand = 0.;
  }
  else if(trigger_fired==2){
    pt_cut_cand = 410.;
    mass_cut_cand = 40.;
  }
  // pT reordered list of largeR jet with have two tracks, 2m/pT<1, DR/R<1
  vector<int> fatjet_aft_presel = fjIdxAftSel(list_reorderd_fj, fatJetPt, fatJetM, vrJetIdFatJet, vrJetPx, vrJetPy, vrJetPz, vrJetPt);
  if (fatjet_aft_presel.size()<1) return idx_final;

  for (int idx = 0; idx < fatjet_aft_presel.size(); idx ++){
    if (idx > 1) continue; 
    if ((fatJetPt.at(fatjet_aft_presel.at(idx))<pt_cut_cand) || (fatJetM.at(fatjet_aft_presel.at(idx))<mass_cut_cand)) continue;
    // get sorted VRjet belonging to this largeRjet
    //vector<pair<int,double>> VRjet = GetSortVRjets(fatjet_aft_presel.at(idx), vrJetIdFatJet, vrJetPt, 0.1*fatJetPt.at(fatjet_aft_presel.at(idx)));
    vector<pair<int,double>> VRjet = GetSortVRjets(fatjet_aft_presel.at(idx), vrJetIdFatJet, vrJetPt, 10.);
    // create a vector (idx_fj,idx_VR_lead,idx_VR_sublead)
    vector<int> idx_final_1component;
    idx_final_1component.push_back(fatjet_aft_presel.at(idx));
    idx_final_1component.push_back(VRjet.at(0).first);
    idx_final_1component.push_back(VRjet.at(1).first);
    idx_final.push_back(idx_final_1component);
  }
  return idx_final;
}

pair<int,double> oldSR_fatjet(vector<vector<int>> indices, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  float SF = 1.;
  pair<int,double> p;
  if ((btag_IS77.at(indices.at(0).at(1))==1) && (btag_IS77.at(indices.at(0).at(2))==1)){
    SF = SF*btag_SF77.at(indices.at(0).at(1))*btag_SF77.at(indices.at(0).at(2));
    p.first = indices.at(0).at(0);
    p.second = SF;
    return p;
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> oldCR_fatjet(vector<vector<int>> indices, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  float SF = 1.;
  pair<int,double> p;
  if ((btag_IS77.at(indices.at(0).at(1))==0) && (btag_IS77.at(indices.at(0).at(2))==0)){
    SF = SF*btag_SF77.at(indices.at(0).at(1))*btag_SF77.at(indices.at(0).at(2));
    p.first = indices.at(0).at(0);
    p.second = SF;
    return p;
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> oldSR_fatjet_DATA(vector<vector<int>> indices, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  float SF = 1.;
  pair<int,double> p;
  if ((btag_IS77.at(indices.at(0).at(1))==1) && (btag_IS77.at(indices.at(0).at(2))==1)){
    p.first = indices.at(0).at(0);
    p.second = SF;
    return p;
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

pair<int,double> oldCR_fatjet_DATA(vector<vector<int>> indices, vector<int> btag_IS77, vector<float> btag_SF77){
  int idx_btag = -1;
  float SF = 1.;
  pair<int,double> p;
  if ((btag_IS77.at(indices.at(0).at(1))==0) && (btag_IS77.at(indices.at(0).at(2))==0)){
    p.first = indices.at(0).at(0);
    p.second = SF;
    return p;
  }
  p.first = idx_btag;
  p.second = SF;
  return p;
}

vector<vector<int>> idxFJcandVRjets_oldSR(bool HLT_j460_a10t_lcw_jes_L1J100MC,bool HLT_j420_a10_lcw_L1J100MC, vector<float> fatJetPt, vector<float> fatJetM, vector<int> vrJetIdFatJet, vector<float> vrJetPx, vector<float> vrJetPy, vector<float> vrJetPz,vector<float> vrJetPt){
  // first component is the FJ index
  // second and third component are the VRjet indices which need to be inspected for btagging
  vector<vector<int>> idx_final;
  // pT reordered largeR indices (possible effects of late calibration cancelled)
  vector<int> list_reorderd_fj;
  list_reorderd_fj = ListSortedFJ(fatJetPt);
  // trigger fired
  // >= 2 fatjet required and offline threshold applied
  if (!FirstEventSelection_oldConf(list_reorderd_fj, HLT_j460_a10t_lcw_jes_L1J100MC, HLT_j420_a10_lcw_L1J100MC, fatJetPt)) return idx_final;
  float mass_cut_cand = 40.;
  float pt_cut_cand = 480.;

  // pT reordered list of largeR jet with have two tracks, 2m/pT<1, DR/R<1
  vector<int> fatjet_aft_presel = fjIdxAftSel(list_reorderd_fj, fatJetPt, fatJetM, vrJetIdFatJet, vrJetPx, vrJetPy, vrJetPz, vrJetPt);
  if (fatjet_aft_presel.size()<1) return idx_final;

  for (int idx = 0; idx < fatjet_aft_presel.size(); idx ++){
    if (idx > 1) continue; 
    if ((fatJetPt.at(fatjet_aft_presel.at(idx))<pt_cut_cand) || (fatJetM.at(fatjet_aft_presel.at(idx))<mass_cut_cand)) continue;
    // get sorted VRjet belonging to this largeRjet
    //vector<pair<int,double>> VRjet = GetSortVRjets(fatjet_aft_presel.at(idx), vrJetIdFatJet, vrJetPt, 0.1*fatJetPt.at(fatjet_aft_presel.at(idx)));
    vector<pair<int,double>> VRjet = GetSortVRjets(fatjet_aft_presel.at(idx), vrJetIdFatJet, vrJetPt, 10.);
    // create a vector (idx_fj,idx_VR_lead,idx_VR_sublead)
    vector<int> idx_final_1component;
    idx_final_1component.push_back(fatjet_aft_presel.at(idx));
    idx_final_1component.push_back(VRjet.at(0).first);
    idx_final_1component.push_back(VRjet.at(1).first);
    idx_final.push_back(idx_final_1component);
  }
  return idx_final;
}

// Takes a fatjet vector and chack if it's matched with the Higgs vector
bool IsParticle(vector<TVector3> particle,TVector3 jet){
  unsigned int IsH = false;
  TVector3 v_def;
  if (particle.size()!=0){
    for(int particle_idx = 0; particle_idx<particle.size(); particle_idx++){
      if ((particle.at(particle_idx)).DeltaR(jet) < 0.4) {
        IsH = true;
      }
    }
  }
  return IsH;
}

vector<TVector3> VecHiggs(vector<float> bpx, vector<float> bpy, vector<float> bpz, vector<int> bpdgId, vector<int> status, double pdg_match){
  vector<TVector3> list_boson;
  TVector3 boson;
  for (unsigned int bos = 0; bos < bpx.size(); bos++){
    if (status.at(bos) == 62) {
      if (bpdgId.at(bos) == pdg_match){
        boson.SetXYZ(bpx.at(bos),bpy.at(bos),bpz.at(bos));
        list_boson.push_back(boson);
      }
    }
  }
  return list_boson;
}

