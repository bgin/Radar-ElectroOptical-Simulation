
vtune_analysis:   select case(analysis_type)
                     case(1)
		         result = SYSTEMQQ("amplxe-cl -collect hotspots -knob sampling-mode=hw -knob sampling-interval=1                            &
                                            -knob enable-stack-collection=true -knob enable-characterization-insights=true                          &
                                             -no-allow-multiple-runs -call-stack-mode=user-only &
                                            data-limit=500 -d 1  -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
		     case(2)
		          result = SYSTEMQQ("amplxe-cl -collect memory-consumption -knob mem-object-size-min-thres=2048          &
                                            -no-allow-multiple-runs -call-stack-mode=user-only &
                                            data-limit=500 -d 1  -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
	             case(3)
		          result = SYSTEMQQ("amplxe-cl -collect-with runsa -cpu-mask=0-7 -no-allow-multiple-runs -call-stack-mode=user-only &
                                            data-limit=500 -d 1  -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
		     case(4)
		          result = SYSTEMQQ("amplxe-cl -collect uarch-exploration -knob sampling-interval=1 -knob collect-frontend-bound=true &
                                             -knob collect-bad-speculation=true -knob collect-memory-bound=true -knob collect-core-bound &
                                             -knob collect-retiring=true -knob collect-memory-bandwidth=true -knob pmu-collection-mode &
                                             -no-allow-multiple-runs -call-stack-mode=user-only &
                                              data-limit=500 -d 1 -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
		     case(5)
		          result = SYSTEMQQ("amplxe-cl -collect memory-access -knob sampling-interval=1 -knob analyze-mem-objects=true &
                                             -knob mem-object-size-min-thres=2048 -knob dram-bandwidth-limits=true 
                                             -no-allow-multiple-runs -call-stack-mode=user-only &
                                            data-limit=500 -d 1 -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
		     case(6)
		          result = SYSTEMQQ("amplxe-cl -collect threading -knob sampling-and-waits=hw  -no-allow-multiple-runs -call-stack-mode=user-only &
                                            data-limit=500 -d 1 -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
		    case(7)
		          result = SYSTEMQQ("amplxe-cl -collect hpc-performance -knob sampling-interval=1 -knob enable-stack-collection=true &
                                            -knob collect-memory-bandwidth=true -knob collect-affinity=false -knob analyze-openmp=true  &
                                            -no-allow-multiple-runs -call-stack-mode=user-only &
                                            data-limit=500 -d 1 -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
		    case(8)
		          result = SYSTEMQQ("amplxe-cl -io -no-allow-multiple-runs -call-stack-mode=user-only &
                                            data-limit=500 -d 1 -inline-mode=on -loop-mode=loop-and-function target-pid"//pid//" &
					    -v")
		          if(result == .false.) then
                            ret = GETLASTERRORQQ()
                            print*, "SYSTEMQQ: Failed to launch vtune profiler  -- reason: ", ret, "vtune_analysis"
                          end if
		    case default
		          print*, "vtune_analysis -- Invalid select case argument: ", analysis_type
	    end select vtune_analysis