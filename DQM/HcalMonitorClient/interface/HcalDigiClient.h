#ifndef HcalDigiClient_H
#define HcalDigiClient_H

#include "DQM/HcalMonitorClient/interface/HcalBaseClient.h"
#include "DQMServices/Core/interface/DQMStore.h"


class HcalDigiClient : public HcalBaseClient {

 public:
  
  /// Constructor
  HcalDigiClient();
  
  /// Destructor
  ~HcalDigiClient();

  void init(const edm::ParameterSet& ps, DQMStore* dbe, string clientName);    

  /// Analyze
  void analyze(void);
  
  /// BeginJob
  void beginJob(void);
  
  /// EndJob
  void endJob(void);
  
  /// BeginRun
  void beginRun(void);
  
  /// EndRun
  void endRun(void);
  
  /// Setup
  void setup(void);
  
  /// Cleanup
  void cleanup(void);
  
  
  ///process report
  void report();
  
  /// WriteDB
  void htmlOutput(int run, string htmlDir, string htmlName);
  void htmlExpertOutput(int run, string htmlDir, string htmlName);
  void getHistograms();
  void loadHistograms(TFile* f);
  
  void resetAllME();
  void createTests();

private:

  TH2F* gl_occ_geo_[4];
  TH2F* gl_occ_elec_[3];
  TH1F* gl_occ_eta_;
  TH1F* gl_occ_phi_;
  TH2F* gl_err_geo_;
  TH2F* gl_err_elec_[3];

  TH1F* gl_num_digi_;
  TH1F* gl_num_bqdigi_;
  TH1F* gl_bqdigi_frac_;
  TH1F* gl_capid_t0_;

  TH2F* sub_occ_geo_[4][4];
  TH2F* sub_occ_elec_[4][3];
  TH1F* sub_occ_eta_[4];
  TH1F* sub_occ_phi_[4];

  TH2F* sub_err_geo_[4];
  TH2F* sub_err_elec_[4][3];

  TH1F* sub_num_bqdigi_[4];
  TH1F* sub_bqdigi_frac_[4];
  TH1F* sub_capid_t0_[4];
  TH1F* sub_digi_shape_[4];
  TH1F* sub_digi_size_[4];

  TH2F* geoRef_;
  
  TH1F* qie_adc_[4];
  TH1F* num_digi_[4];
  TH1F* qie_capid_[4];
  TH1F* qie_dverr_[4];

  TH2F* ProblemDigiCells;
  TH2F* ProblemDigiCells_DEPTH[4];
  
  TH2F* ProblemDigiCellsHB;
  TH2F* ProblemDigiCellsHE;
  TH2F* ProblemDigiCellsHO;
  TH2F* ProblemDigiCellsHF;
  
  TH2F* ProblemDigiCellsHB_DEPTH[4];
  TH2F* ProblemDigiCellsHE_DEPTH[4];
  TH2F* ProblemDigiCellsHO_DEPTH[4];
  TH2F* ProblemDigiCellsHF_DEPTH[4];
  
  double errorFrac_;

  /*
  TH2F* RawPedestalMean[4];
  TH2F* RawPedestalRMS[4];
  TH2F* SubPedestalMean[4]; 
  TH2F* SubPedestalRMS[4]; 
  */

  TH2F* HBRawPedestalMean[4];
  TH2F* HBRawPedestalRMS[4];
  TH2F* HBSubPedestalMean[4];
  TH2F* HBSubPedestalRMS[4];

  TH2F* HERawPedestalMean[4];
  TH2F* HERawPedestalRMS[4];
  TH2F* HESubPedestalMean[4];
  TH2F* HESubPedestalRMS[4];

  TH2F* HORawPedestalMean[4];
  TH2F* HORawPedestalRMS[4];
  TH2F* HOSubPedestalMean[4];
  TH2F* HOSubPedestalRMS[4];

  TH2F* HFRawPedestalMean[4];
  TH2F* HFRawPedestalRMS[4];
  TH2F* HFSubPedestalMean[4];
  TH2F* HFSubPedestalRMS[4];

  TH1F* HBRawPedestalMean_1D[4];
  TH1F* HBRawPedestalRMS_1D[4];
  TH1F* HBSubPedestalMean_1D[4];
  TH1F* HBSubPedestalRMS_1D[4];

  TH1F* HERawPedestalMean_1D[4];
  TH1F* HERawPedestalRMS_1D[4];
  TH1F* HESubPedestalMean_1D[4];
  TH1F* HESubPedestalRMS_1D[4];

  TH1F* HORawPedestalMean_1D[4];
  TH1F* HORawPedestalRMS_1D[4];
  TH1F* HOSubPedestalMean_1D[4];
  TH1F* HOSubPedestalRMS_1D[4];

  TH1F* HFRawPedestalMean_1D[4];
  TH1F* HFRawPedestalRMS_1D[4];
  TH1F* HFSubPedestalMean_1D[4];
  TH1F* HFSubPedestalRMS_1D[4];

};

#endif
