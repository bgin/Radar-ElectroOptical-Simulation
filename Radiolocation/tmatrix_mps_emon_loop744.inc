

emon_loop744:  select case(pmc_types)
                   case(1)
		      result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                        -C "INST_RETIRED.ANY,CPU_CLK_UNHALTED.THREAD,CPU_CLK_UNHALTED.REF_TSC, &
                                         CPU_CLK_UNHALTED.THREAD_ANY" ")
		               if(result == .false.) then
			           ret = GETLASTERRORQQ()
				   print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
			       end if
	           case(2)
		      result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                        -C "LD_BLOCKS.STORE_FORWARD,LD_BLOCKS.NO_SR,MISALIGN_MEM_REF.LOADS, 
                                        MISALIGN_MEM_REF.STORES" ")
			       if(result == .false.) then
			            ret = GETLASTERRORQQ()
				    print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
			       end if
	           case(3)
		      result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                        -C LD_BLOCKS_PARTIAL.ADDRESS_ALIAS,DTLB_LOAD_MISSES.MISS_CAUSES_A_WALK &
                                         DTLB_LOAD_MISSES.WALK_COMPLETED_4K,DTLB_LOAD_MISSES.WALK_COMPLETED_2M_4M")
			       if(result == .false.) then
			            ret = GETLASTERRORQQ()
				    print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
			       end if
	           case(4)
		      result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                        -C DTLB_LOAD_MISSES.WALK_COMPLETED,DTLB_LOAD_MISSES.WALK_DURATION,
                                         DTLB_LOAD_MISSES.STLB_HIT_4K,DTLB_LOAD_MISSES.STLB_HIT")
			       if(result == .false.) then
			            ret = GETLASTERRORQQ()
				    print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
			       end if
		   case(5)
		      result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                        -C DTLB_LOAD_MISSES.PDE_CACHE_MISS,INT_MISC.RECOVERY_CYCLES,  &
                                         INT_MISC.RECOVERY_CYCLES_ANY,UOPS_ISSUED.ANY")
			       if(result == .false.) then
			            ret = GETLASTERRORQQ()
				    print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
			       end if
		   case(6)
		      result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                        -C UOPS_ISSUED.STALL_CYCLES,UOPS_ISSUED.CORE_STALL_CYCLES,    &
                                        UOPS_ISSUED.SLOW_LEA,UOPS_ISSUED.SINGLE_MUL")
			       if(result == .false.) then
			            ret = GETLASTERRORQQ()
				    print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
			       end if
	           case(7)
		       result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                         -C ARITH.DIVIDER_UOPS, L2_RQSTS.DEMAND_DATA_RD_MISS,         &
                                          L2_RQSTS.RFO_MISS, L2_RQSTS.CODE_RD_MISS, L2_RQSTS.ALL_DEMAND_MISS ")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(8)
		       result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                         -C L2_RQSTS.L2_PF_MISS, L2_RQSTS.MISS, L2_RQSTS.DEMAND_DATA_RD_HIT, &
                                          L2_RQSTS.RFO_HIT")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(9)
		       result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                         -C L2_RQSTS.CODE_RD_HIT, L2_RQSTS.L2_PF_HIT, L2_RQSTS.ALL_DEMAND_DATA_RD, &
                                          L2_RQSTS.ALL_RFO")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(10)
		       result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                         -C L2_RQSTS.ALL_CODE_RD, L2_RQSTS.ALL_DEMAND_REFERENCES, L2_RQSTS.ALL_PF, &
                                          L2_RQSTS.REFERENCES")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(11)
		       result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                         -C LONGEST_LAT_CACHE.MISS, LONGEST_LAT_CACHE.REFERENCE,            &
                                          CPU_CLK_UNHALTED.THREAD_P, CPU_CLK_UNHALTED.THREAD_P_ANY ")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(12)
		        result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                          -C CPU_CLK_THREAD_UNHALTED.REF_XCLK, CPU_CLK_THREAD_UNHALTED.REF_XCLK_ANY, &
					   CPU_CLK_UNHALTED.REF_XCLK, CPU_CLK_UNHALTED.REF_XCLK_ANY ")
                                if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(13)
		        result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                          -C CPU_CLK_THREAD_UNHALTED.ONE_THREAD_ACTIVE,  &
                                           CPU_CLK_UNHALTED.ONE_THREAD_ACTIVE ")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(14)
		        result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                          -C L1D_PEND_MISS.PENDING, L1D_PEND_MISS.PENDING_CYCLES,           &
                                           L1D_PEND_MISS.PENDING_CYCLES_ANY,L1D_PEND_MISS.REQUEST_FB_FULL ")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(15)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C L1D_PEND_MISS.FB_FULL,DTLB_STORE_MISSES.MISS_CAUSES_A_WALK, &
                                            DTLB_STORE_MISSES.WALK_COMPLETED_4K,DTLB_STORE_MISSES.WALK_COMPLETED ")
			       if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(16)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C LOAD_HIT_PRE.SW_PF,LOAD_HIT_PRE.HW_PF,EPT.WALK_CYCLES,  &
                                           L1D.REPLACEMENT ")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(17)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C MOVE_ELIMINATION.INT_ELIMINATED,MOVE_ELIMINATION.SIMD_ELIMINATED")
                                if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(18)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C RS_EVENTS.EMPTY_CYCLES,RS_EVENTS.EMPTY_END,IDQ.EMPTY,IDQ.MITE_UOPS ")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(19)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C IDQ.MITE_CYCLES,IDQ.DSB_UOPS,IDQ.DSB_CYCLES,IDQ.MS_DSB_UOPS")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(20)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C IDQ.MS_DSB_CYCLES,IDQ.MS_DSB_OCCUR,IDQ.ALL_DSB_CYCLES_4_UOPS, &
                                            IDQ.ALL_DSB_CYCLES_ANY_UOPS")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(21)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C IDQ.MS_MITE_UOPS,IDQ.ALL_MITE_CYCLES_4_UOPS,IDQ.ALL_MITE_CYCLES_ANY_UOPS,
                                            IDQ.MS_UOPS")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(22)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C IDQ.MS_CYCLES,IDQ.MS_SWITCHES,IDQ.MITE_ALL_UOPS")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(23)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C ICACHE.HIT,ICACHE.MISSES,ICACHE.IFETCH_STALL,ICACHE.IFDATA_STALL")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(24)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C ITLB_MISSES.MISS_CAUSES_A_WALK,ITLB_MISSES.WALK_COMPLETED_4K, &
                                            ITLB_MISSES.WALK_COMPLETED,ITLB_MISSES.WALK_DURATION")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(25)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C ITLB_MISSES.STLB_HIT,ILD_STALL.LCP,ILD_STALL.IQ_FULL")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(26)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C BR_INST_EXEC.NON_TAKEN_CONDITIONAL,BR_INST_EXEC.TAKEN_CONDITIONAL, &
                                            BR_INST_EXEC.TAKEN_DIRECT_JUMP,BR_INST_EXEC.TAKEN_INDIRECT_JUMP_NON_CALL_RET")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(27)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C BR_INST_EXEC.TAKEN_INDIRECT_NEAR_RETURN,BR_INST_EXEC.TAKEN_INDIRECT_NEAR_CALL,
                                            BR_INST_EXEC.TAKEN_DIRECT_NEAR_CALL,BR_INST_EXEC.ALL_CONDITIONAL")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(28)
		         result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                           -C BR_INST_EXEC.ALL_DIRECT_JMP,BR_INST_EXEC.ALL_INDIRECT_JUMP_NON_CALL_RET, &
                                            BR_INST_EXEC.ALL_INDIRECT_NEAR_RETURN,BR_INST_EXEC.ALL_DIRECT_NEAR_CALL")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(29)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C BR_INST_EXEC.ALL_BRANCHES,BR_MISP_EXEC.NONTAKEN_CONDITIONAL,  &
                                             BR_MISP_EXEC.TAKEN_CONDITIONAL,BR_MISP_EXEC.TAKEN_INDIRECT_JUMP_NON_CALL_RET")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(30)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C BR_MISP_EXEC.TAKEN_RETURN_NEAR,BR_MISP_EXEC.TAKEN_INDIRECT_NEAR_CALL, &
                                             BR_MISP_EXEC.ALL_CONDITIONAL,BR_MISP_EXEC.ALL_BRANCHES")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(31)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C IDQ_UOPS_NOT_DELIVERED.CORE,IDQ_UOPS_NOT_DELIVERED.CYCLES_0_UOPS_DELIV_CORE, &
                                             IDQ_UOPS_NOT_DELIVERED.CYCLES_LE_1_UOP_DELIV.CORE,IDQ_UOPS_NOT_DELIVERED.CYCLES_LE_2_UOP_DELIV.CORE")
                                if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
                   case(32)
                          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
			                    -C IDQ_UOPS_NOT_DELIVERED.CYCLES_LE_3_UOP_DELIV.CORE,  &
                                             IDQ_UOPS_NOT_DELIVERED.CYCLES_FE_WAS_OK")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(33)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_0,UOPS_EXECUTED_PORT.PORT_0_CORE,      &
                                             UOPS_DISPATCHED_PORT.PORT_0")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(34)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_1,UOPS_EXECUTED_PORT.PORT_1_CORE,  &
                                             UOPS_DISPATCHED_PORT.PORT_1")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(35)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_2,UOPS_EXECUTED_PORT.PORT_2_CORE, &
                                            UOPS_DISPATCHED_PORT.PORT_2 ")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(36)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_3,UOPS_EXECUTED_PORT.PORT_3_CORE, &
                                             UOPS_DISPATCHED_PORT.PORT_3")
			         if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				 end if
		   case(37)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_4,UOPS_EXECUTED_PORT.PORT_4_CORE, &
                                             UOPS_DISPATCHED_PORT.PORT_4")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(38)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_5,UOPS_EXECUTED_PORT.PORT_5_CORE, &
                                            UOPS_DISPATCHED_PORT.PORT_5")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(39)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_6,UOPS_EXECUTED_PORT.PORT_6_CORE, &
                                             UOPS_DISPATCHED_PORT.PORT_6")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(40)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED_PORT.PORT_7,UOPS_EXECUTED_PORT.PORT_7_CORE, &
                                             UOPS_DISPATCHED_PORT.PORT_7")
			        if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(41)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C RESOURCE_STALLS.ANY,RESOURCE_STALLS.RS,RESOURCE_STALLS.SB, &
                                             RESOURCE_STALLS.ROB")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if
		   case(42)
		          result = SYSTEMQQ("emon -v -M -m -d -f emon_tmatrix_mps_driver_loop744.txt -t0.0079 &
                                            -C UOPS_EXECUTED.STALL_CYCLES,UOPS_EXECUTED.CORE_CYCLES_NONE, &
                                             UOPS_EXECUTED.CORE")
				if(result == .false.) then
				     ret = GETLASTERRORQQ()
				     print*, "SYSTEMQQ: Failed to execute emon command -- reason: ", ret
				end if