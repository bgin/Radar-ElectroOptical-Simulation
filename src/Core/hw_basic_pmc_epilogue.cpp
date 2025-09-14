

                for(int32_t i = 0; i != 4; ++i)                                                                                                                               
                {                                                                                                                                                             
                     prog_counters_stop[i] = rdpmc(i);                                                                                                                        
                }                                                                                                                                                                                                                                                                                                                                                                    \
                act_cyc_stop = rdpmc_actual_cycles();                                                                                                                         
                ref_cyc_stop = rdpmc_reference_cycles();                                                                                                                      
                tsc_stop = rdtsc();                                                                                                                                           
	            core_counter_width = get_core_counter_width();                                                                                                                            
	            nom_ghz = get_TSC_frequency()/1.0e9;                                                                                                                          
	            utilization = (double)(ref_cyc_stop-ref_cyc_start)/(double)(tsc_stop-tsc_start-18ull);                                                                 
	            avg_ghz = (double)(act_cyc_stop-act_cyc_start)/(double)(tsc_stop-tsc_start-18ull)*nom_ghz;     
                 syslog(LOG_INFO,"%-10s:\n", __PRETTY_FUNCTION__);                                                                                                             
	            syslog(LOG_INFO, "*************** Hardware Counters -- Dump Begin **************");                                                                           
	            syslog(LOG_INFO,"Core utilization                      : %f\n",utilization );                                                                                 
	            syslog(LOG_INFO,"Core average frequency                : %f\n",avg_ghz);                                                                                      
	            syslog(LOG_INFO,"Reference cycles                      : %20lld\n",ref_cyc_stop-ref_cyc_start);                                                                
	            syslog(LOG_INFO,"Actual cycles                         : %20lld\n",act_cyc_stop-act_cyc_start);                                                               
	            syslog(LOG_INFO,"%20lld\n",corrected_pmc_delta(prog_counters_stop[0],prog_counters_start[0],core_counter_width));      
	            syslog(LOG_INFO,"%20lld\n",corrected_pmc_delta(prog_counters_stop[1],prog_counters_start[1],core_counter_width));      
	            syslog(LOG_INFO,"%20lld\n",corrected_pmc_delta(prog_counters_stop[2],prog_counters_start[2],core_counter_width));      
	            syslog(LOG_INFO,"%20lld\n",corrected_pmc_delta(prog_counters_stop[3],prog_counters_start[3],core_counter_width));          
	            syslog(LOG_INFO, "*************** Hardware Counters -- Dump End   **************");                 
                                                                                                                                       
	             