


               uint64_t prog_counters_start[4] = {};                
               uint64_t prog_counters_end[4]   = {};                
               uint64_t tsc_start,tsc_stop;                         
               uint64_t act_cyc_start,act_cyc_stop;                 
               uint64_t ref_cyc_start,ref_cyc_stop;                 
               uint64_t [[maybe_unused]] dummy1;
               uint64_t [[maybe_unused]] dummy2;
               uint64_t [[maybe_unused]] dummy3;     
               int32_t core_counter_width;                          
               double utilization,nom_ghz,avg_ghz;                  
               dummy1 = rdtsc();                                    
               dummy2 = rdtsc();                                    
               dummy3 = rdpmc(0);                                   
               core_counter_width = get_core_counter_width();       
               for(int32_t i = 0; i != 4; ++i) 
               {                    
                   prog_counters_start[i] = rdpmc(i);               
               }                                                    
               act_cyc_start = rdpmc_actual_cycles();               
               ref_cyc_start = rdpmc_reference_cycles();            
               tsc_start = rdtsc();    