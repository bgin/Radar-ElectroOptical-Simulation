

perf_loop1097:   select case(pmc_loop1097)
                   case(1)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set1.txt -er203 -er803 -er105 -er107 &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
		   case(2)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set2.txt -er100e -er0214 -er3f24 -er003c &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
		   case(3)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set3.txt -er0148 -er0248 -er0149 -er024c &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
		   case(4)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set4.txt -er0151 -er0279 -er0180 -er0280 &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
		   case(5)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set5.txt -er0480 -er0185 -er4188 -er019c &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
		   case(6)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set6.txt -er01a1 -er02a1 -er04a1 -er08a1 &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
		   case(7)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set7.txt -er10a1 -er20a1 -er80a1 -er01a2 &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
		   case(8)
		         result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1097_set8.txt -er11bc -er01c2  &
                                           -C7  -a sleep 0.1")
                         if(result == .false.) then
                              ret = GETLASTERRORQQ()
                              print*, "SYSTEMQQ: perf_loop1097 -- Failed to execute perf command -- reason: ", ret, "loop1097"
                         end if
	           case default
		         print*, "perf_loop1097: Invalid select case argument: ", pmc_loop1097
	           end select perf_loop1097