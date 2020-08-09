
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include "GMS_pmc_meter.h"

gms::system
::PMCMeter::PMCMeter(const std::size nruns,
		     const fptr ptfunc,
		     const char * __restrict func_name,
		     const char * __restrict file_name,
		     const char * __restrict event1_name,
		     const char * __restrict event2_name,
		     const char * __restrict event3_name,
		     const char * __restrict event4_name,
		     const bool  warmup,
		     const int32_t thnum)
:
m_nruns{nruns},
m_ptfunc{ptfunc},
m_func_name{func_name},
m_file_name{file_name},
m_event1_name{event1_name},
m_event2_name{event2_name},
m_event3_name{event3_name},
m_event4_name{event4_name},
m_warmup{warmup}{


 
   // Begin an allocation
   m_instr_issued = std::vector<PFC_CNT>(0LL,m_nruns);
   m_core_cycles  = std::vector<PFC_CNT>(0LL,m_nruns);
   m_ref_cycles   = std::vector<PFC_CNT>(0LL,m_nruns);
   m_hw_event1    = std::vector<PFC_CNT>(0LL,m_nruns);
   m_hw_event2    = std::vector<PFC_CNT>(0LL,m_nruns);
   m_hw_event3    = std::vector<PFC_CNT>(0LL,m_nruns);
   m_hw_event4    = std::vector<PFC_CNT>(0LL,m_nruns);
}
  




gms::system
::PMCMeter
::~PMCMeter() {
   if(m_ptfunc != nullptr) m_ptfunc = nullptr;
   if(m_event1_name != nullptr) m_event1_name = nullptr;
   if(m_event2_name != nullptr) m_event2_name = nullptr;
   if(m_event3_name != nullptr) m_event3_name = nullptr;
   if(m_event4_name != nullptr) m_event4_name = nullptr;
}


int32_t
gms::system
::PMCMeter::
run_benchmark() {

    // Attempt to initialize lbpfc PMC access driver
    int32_t err = 9999;
    err = pfcPinThread(thnum);
    if(err != 0) {
       std::cerr << __PRETTY_FUNCTION__ << ": ***Fatal*** -- could not set Core affinity: " << err << std::endl;
       return (err);
    }
    err = 9999;
    err = pfcInit();
    if(err != 0) {
        std::cerr << __PRETTY_FUNCTION__ << ": ***Fatal*** -- could not initialize the PMC driver: " << err << std::endl;
        return (err);
     }
     // Could not reach this block location in theory
     // Begin PMC initialization.
     static const PFC_CNT ZERO_CNT[7] = {0LL,0LL,0LL,0LL,0LL,0LL,0LL};
     PFC_CNT CNT[7] = {0LL,0LL,0LL,0LL,0LL,0LL,0LL};
     PFC_CFG CFG[7] = {2,2,2,0,0,0,0};
     int32_t err;
     CFG[3] = pfcParseCfg(m_event1_name);
     CFG[4] = pfcParseCfg(m_event2_name);
     CFG[5] = pfcParseCfg(m_event3_name);
     CFG[6] = pfcParseCfg(m_event4_name);
     err = -9;
     if(m_warmup) {
         m_ptfunc();
     }
     for(std::size_t i = 0; i != m_nruns; ++i) {
        //
        err = pfcWrCfgs(0,7,CFG);
        if(err != 0) {
              std::cerr << __PRETTY_FUNC__ << ": ***FATAL*** -- failed to configure a PMC: " << err << std::endl;
	      return (err);
         }
         pfcWrCnts(0,7,ZERO_CNT);
         memset(CNT,0,sizeof(CNT));
         PFCSTART(CNT);
	 m_ptfunc();
	 PFCEND(CNT);
	 pfcRemoveBias(CNT,1);
	 m_instr_issued[i] = CNT[0];
	 m_core_cycles[i]  = CNT[1];
	 m_ref_cycles[i]   = CNT[2];
	 m_hw_event1[i]    = CNT[3];
	 m_hw_event2[i]    = CNT[4];
	 m_hw_event3[i]    = CNT[5];
	 m_hw_event4[i]    = CNT[6];
     }
     pfcFini();
     return (0);
}
