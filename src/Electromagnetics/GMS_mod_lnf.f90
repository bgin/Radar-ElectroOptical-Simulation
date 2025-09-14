

module mod_lnf

     use mod_kinds, only : int4,dp


     contains

       real(kind=dp) function lnf(z)
         real(kind=dp)  :: z
         ! Locals
         real(kind=dp), dimension(11), parameter :: c0 = [0.16427423239836267e+5_dp, -0.48589401600331902e+5_dp, &
                                                          0.55557391003815523e+5_dp, -0.30964901015912058e+5_dp, &
                                                          0.87287202992571788e+4_dp, -0.11714474574532352e+4_dp, &
                                                          0.63103078123601037e+2_dp, -0.93060589791758878e+0_dp, &
                                                          0.13919002438227877e-2_dp,-0.45006835613027859e-8_dp,&
                                                          0.13069587914063262e-9_dp]
         real(kind=dp) :: a,b,cp
         integer(kind=int4) :: i
         ! Exec code
         a = 1.0_dp
         cp = 2.5066282746310005_dp
         b = z+10.5_dp
         b = (z+0.5_dp)*log(b)-b
         do i=1, 11
            z = z+1.0_dp
            a = a+c0(i)/z
         end do
         lnf = b+log(cp*a)
         
    end function lnf


end module mod_lnf
