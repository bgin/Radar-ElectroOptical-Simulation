
#ifndef __GMS_PMC_EVENTS_NAMES_H__
#define __GMS_PMC_EVENTS_NAMES_H__



//
// Global data structure representing currently
// programmed (by msr-tools) performance events
// This structure should be initialized at the program
// startup.
//


typedef struct PmcEventsNames {

   char * pmc_event0;
   char * pmc_event1;
   char * pmc_event2;
   char * pmc_event3;
   char * pmc_event4;
   char * pmc_event5;
   char * pmc_event6;
   char * pmc_event7;
  
} PmcEventsNames_t;



#endif /*__GMS_PMC_EVENTS_NAMES_H__*/
