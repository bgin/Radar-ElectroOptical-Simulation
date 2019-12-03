

perf_loop1828:  select case(pmc_loop1828)
                   case(1)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set1.txt -er100 -er300 -er203 -er803 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(2)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set2.txt -er105 -er108 -er100e -er200e &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(3)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set3.txt -er400e -er0214 -er2424 -ere424 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(4)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set4.txt -er003c -er0148 -er248 -er0149 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(5)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set5.txt -er0151 -er015e -er0263 -er0279 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(6)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set6.txt -er0879 -er1879 -er0180 -er0280 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(7)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set7.txt -er0480 -er0185 -er0187 -er0487 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(8)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set8.txt -er4188 -er8188 -erff88 -er4189 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(9)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set9.txt -er019c -er01a1 -er02a1 -er04a1 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(10)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set10.txt -er08a1 -er10a1 -er20a1 -er40a1 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(11)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set11.txt -er01a2 -er08a2 -er10a2 -er02b1 &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case(12)
		       result = SYSTEMQQ("perf stat -o tmatrix_mps_driver_loop1828_set12.txt -er00c0 -er08c1 -er10c1 -er01c2 -er076c &
                                        -C 7  -a sleep 1")
                       if(result == .false.) then
                          ret = GETLASTERRORQQ()
                          print*, "SYSTEMQQ: perf_loop1828 -- Failed to execute perf command -- reason: ", ret, "loop1828"
                       end if
		   case default
		       print*, "perf_loop1828 -- Invalid select case argument: ", pmc_loop1828
		   end select  perf_loop1828