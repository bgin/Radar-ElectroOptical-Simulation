mkdir -p ./lib
mkdir -p ./objs
( cd ./src; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f ../lib/libmudpack.a ../objs/cud2.o ../objs/cud24.o ../objs/cud24cr.o ../objs/cud24sp.o ../objs/cud2cr.o ../objs/cud2sp.o ../objs/cud3.o ../objs/cud34.o ../objs/cud34sp.o ../objs/cud3cr.o ../objs/cud3ln.o ../objs/cud3pn.o ../objs/cud3sp.o ../objs/cudcom.o ../objs/cuh2.o ../objs/cuh24.o ../objs/cuh24cr.o ../objs/cuh2cr.o ../objs/cuh3.o ../objs/cuh34.o ../objs/mud2.o ../objs/mud24.o ../objs/mud24cr.o ../objs/mud24sp.o ../objs/mud2cr.o ../objs/mud2sa.o ../objs/mud2sp.o ../objs/mud3.o ../objs/mud34.o ../objs/mud34sp.o ../objs/mud3cr.o ../objs/mud3ln.o ../objs/mud3pn.o ../objs/mud3sa.o ../objs/mud3sp.o ../objs/mudcom.o ../objs/muh2.o ../objs/muh24.o ../objs/muh24cr.o ../objs/muh2cr.o ../objs/muh3.o ../objs/muh34.o ../objs/resc2.o ../objs/resc2cr.o ../objs/resc2sp.o ../objs/resc3.o ../objs/resc3sp.o ../objs/resm2.o ../objs/resm2cr.o ../objs/resm2sp.o ../objs/resm3.o ../objs/resm3cr.o ../objs/resm3sp.o 
gfortran  -c cud2.f -o ../objs/cud2.o
gfortran  -c cud24.f -o ../objs/cud24.o
gfortran  -c cud24cr.f -o ../objs/cud24cr.o
gfortran  -c cud24sp.f -o ../objs/cud24sp.o
gfortran  -c cud2cr.f -o ../objs/cud2cr.o
gfortran  -c cud2sp.f -o ../objs/cud2sp.o
gfortran  -c cud3.f -o ../objs/cud3.o
gfortran  -c cud34.f -o ../objs/cud34.o
gfortran  -c cud34sp.f -o ../objs/cud34sp.o
gfortran  -c cud3cr.f -o ../objs/cud3cr.o
gfortran  -c cud3ln.f -o ../objs/cud3ln.o
gfortran  -c cud3pn.f -o ../objs/cud3pn.o
gfortran  -c cud3sp.f -o ../objs/cud3sp.o
gfortran  -c cudcom.f -o ../objs/cudcom.o
gfortran  -c cuh2.f -o ../objs/cuh2.o
gfortran  -c cuh24.f -o ../objs/cuh24.o
gfortran  -c cuh24cr.f -o ../objs/cuh24cr.o
gfortran  -c cuh2cr.f -o ../objs/cuh2cr.o
gfortran  -c cuh3.f -o ../objs/cuh3.o
gfortran  -c cuh34.f -o ../objs/cuh34.o
gfortran  -c mud2.f -o ../objs/mud2.o
gfortran  -c mud24.f -o ../objs/mud24.o
gfortran  -c mud24cr.f -o ../objs/mud24cr.o
gfortran  -c mud24sp.f -o ../objs/mud24sp.o
gfortran  -c mud2cr.f -o ../objs/mud2cr.o
gfortran  -c mud2sa.f -o ../objs/mud2sa.o
gfortran  -c mud2sp.f -o ../objs/mud2sp.o
gfortran  -c mud3.f -o ../objs/mud3.o
gfortran  -c mud34.f -o ../objs/mud34.o
gfortran  -c mud34sp.f -o ../objs/mud34sp.o
gfortran  -c mud3cr.f -o ../objs/mud3cr.o
gfortran  -c mud3ln.f -o ../objs/mud3ln.o
gfortran  -c mud3pn.f -o ../objs/mud3pn.o
gfortran  -c mud3sa.f -o ../objs/mud3sa.o
gfortran  -c mud3sp.f -o ../objs/mud3sp.o
gfortran  -c mudcom.f -o ../objs/mudcom.o
gfortran  -c muh2.f -o ../objs/muh2.o
gfortran  -c muh24.f -o ../objs/muh24.o
gfortran  -c muh24cr.f -o ../objs/muh24cr.o
gfortran  -c muh2cr.f -o ../objs/muh2cr.o
gfortran  -c muh3.f -o ../objs/muh3.o
gfortran  -c muh34.f -o ../objs/muh34.o
gfortran  -c resc2.f -o ../objs/resc2.o
gfortran  -c resc2cr.f -o ../objs/resc2cr.o
gfortran  -c resc2sp.f -o ../objs/resc2sp.o
gfortran  -c resc3.f -o ../objs/resc3.o
gfortran  -c resc3sp.f -o ../objs/resc3sp.o
gfortran  -c resm2.f -o ../objs/resm2.o
gfortran  -c resm2cr.f -o ../objs/resm2cr.o
gfortran  -c resm2sp.f -o ../objs/resm2sp.o
gfortran  -c resm3.f -o ../objs/resm3.o
gfortran  -c resm3cr.f -o ../objs/resm3cr.o
gfortran  -c resm3sp.f -o ../objs/resm3sp.o
/usr/bin/ar -rv ../lib/libmudpack.a ../objs/cud2.o ../objs/cud24.o ../objs/cud24cr.o ../objs/cud24sp.o ../objs/cud2cr.o ../objs/cud2sp.o ../objs/cud3.o ../objs/cud34.o ../objs/cud34sp.o ../objs/cud3cr.o ../objs/cud3ln.o ../objs/cud3pn.o ../objs/cud3sp.o ../objs/cudcom.o ../objs/cuh2.o ../objs/cuh24.o ../objs/cuh24cr.o ../objs/cuh2cr.o ../objs/cuh3.o ../objs/cuh34.o ../objs/mud2.o ../objs/mud24.o ../objs/mud24cr.o ../objs/mud24sp.o ../objs/mud2cr.o ../objs/mud2sa.o ../objs/mud2sp.o ../objs/mud3.o ../objs/mud34.o ../objs/mud34sp.o ../objs/mud3cr.o ../objs/mud3ln.o ../objs/mud3pn.o ../objs/mud3sa.o ../objs/mud3sp.o ../objs/mudcom.o ../objs/muh2.o ../objs/muh24.o ../objs/muh24cr.o ../objs/muh2cr.o ../objs/muh3.o ../objs/muh34.o ../objs/resc2.o ../objs/resc2cr.o ../objs/resc2sp.o ../objs/resc3.o ../objs/resc3sp.o ../objs/resm2.o ../objs/resm2cr.o ../objs/resm2sp.o ../objs/resm3.o ../objs/resm3cr.o ../objs/resm3sp.o 
ar: creating archive ../lib/libmudpack.a
a - ../objs/cud2.o
a - ../objs/cud24.o
a - ../objs/cud24cr.o
a - ../objs/cud24sp.o
a - ../objs/cud2cr.o
a - ../objs/cud2sp.o
a - ../objs/cud3.o
a - ../objs/cud34.o
a - ../objs/cud34sp.o
a - ../objs/cud3cr.o
a - ../objs/cud3ln.o
a - ../objs/cud3pn.o
a - ../objs/cud3sp.o
a - ../objs/cudcom.o
a - ../objs/cuh2.o
a - ../objs/cuh24.o
a - ../objs/cuh24cr.o
a - ../objs/cuh2cr.o
a - ../objs/cuh3.o
a - ../objs/cuh34.o
a - ../objs/mud2.o
a - ../objs/mud24.o
a - ../objs/mud24cr.o
a - ../objs/mud24sp.o
a - ../objs/mud2cr.o
a - ../objs/mud2sa.o
a - ../objs/mud2sp.o
a - ../objs/mud3.o
a - ../objs/mud34.o
a - ../objs/mud34sp.o
a - ../objs/mud3cr.o
a - ../objs/mud3ln.o
a - ../objs/mud3pn.o
a - ../objs/mud3sa.o
a - ../objs/mud3sp.o
a - ../objs/mudcom.o
a - ../objs/muh2.o
a - ../objs/muh24.o
a - ../objs/muh24cr.o
a - ../objs/muh2cr.o
a - ../objs/muh3.o
a - ../objs/muh34.o
a - ../objs/resc2.o
a - ../objs/resc2cr.o
a - ../objs/resc2sp.o
a - ../objs/resc3.o
a - ../objs/resc3sp.o
a - ../objs/resm2.o
a - ../objs/resm2cr.o
a - ../objs/resm2sp.o
a - ../objs/resm3.o
a - ../objs/resm3cr.o
a - ../objs/resm3sp.o
( cd ./test; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f  tcud2.exe tcud24.exe tcud24cr.exe tcud24sp.exe tcud2cr.exe tcud2sp.exe tcud3.exe tcud34.exe tcud34sp.exe tcud3cr.exe tcud3sp.exe tcuh2.exe tcuh24.exe tcuh24cr.exe tcuh2cr.exe tcuh3.exe tcuh34.exe tmud2.exe tmud24.exe tmud24cr.exe tmud24sp.exe tmud2cr.exe tmud2sa.exe tmud2sp.exe tmud3.exe tmud34.exe tmud34sp.exe tmud3cr.exe tmud3sa.exe tmud3sp.exe tmuh2.exe tmuh24.exe tmuh24cr.exe tmuh2cr.exe tmuh3.exe tmuh34.exe
rm -f tcud2.exe
gfortran  tcud2.f -o tcud2.exe -L../lib -l mudpack
./tcud2.exe


 cud2 test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp =  3 jyq =  2 iex =  5 jey =  6
 nx =  49 ny =  65 iguess =  0 maxcy =  1
 method =  0 work space estimate =   42776

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cud2 intl =  0
 ierror =  0 minimum work space =   34192

 approximation call to cud2
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.615E-03
rm -f tcud24.exe
gfortran  tcud24.f -o tcud24.exe -L../lib -l mudpack
./tcud24.exe


 cud2 test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp =  3 jyq =  2 iex =  5 jey =  6
 nx =  49 ny =  65 iguess =  0 maxcy =  3
 method =  0 work space estimate =   34192

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cud2 intl =  0
 ierror =  0 minimum work space =   34192

 approximation call to cud2
 intl =  1 method =  0 iguess =  0 maxcy =  3
 ierror =  0
 maximum error  =   0.624E-03

 cud24 test  ierror =  0
 maximum error  =   0.179E-04
rm -f tcud24cr.exe
gfortran  tcud24cr.f -o tcud24cr.exe -L../lib -l mudpack
./tcud24cr.exe


 cud2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp =  3 jyq =  2 iex =  5 jey =  6
 nx =  49 ny =  65 iguess =  0 maxcy =  1
 method =  2 work space estimate =   64474

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cud2cr intl =  0
 ierror =  0 minimum work space =   64474

 approximation call to cud2cr
 intl =  1 method =  2 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.188E-03

 approximation call to cud2cr
 intl =  1 method =  2 iguess =  1 maxcy =  2
 ierror =  0
 maximum error  =   0.190E-03

 cud24cr test  ierror =  0
 maximum error  =   0.731E-05
rm -f tcud24sp.exe
gfortran  tcud24sp.f -o tcud24sp.exe -L../lib -l mudpack
./tcud24sp.exe


 cud2sp test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp =  2 jyq =  3 iex =  7 jey =  6
 nx = 129 ny =  97 iguess =  0 maxcy =  3
 method =  0 work space estimate =   48976

 multigrid option arguments 
 kcycle =  0
 iprer = -1
 ipost = -1
 intpol = **

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cud2sp intl =  0
 ierror =  0 minimum work space =   48976

 approximation call to cud2sp 
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.331E-03

 cud24sp test  ierror =  0
 maximum error  =   0.179E-04
rm -f tcud2cr.exe
gfortran  tcud2cr.f -o tcud2cr.exe -L../lib -l mudpack
./tcud2cr.exe


 cud2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp =  3 jyq =  2 iex =  5 jey =  6
 nx =  49 ny =  65 iguess =  0 maxcy =  1
 method =  2 work space estimate =   67426

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cud2cr intl =  0
 ierror =  0 minimum work space =   64474

 approximation call to cud2cr
 intl =  1 method =  2 iguess =  0
 ierror =  0
 maximum error  =   0.188E-03
rm -f tcud2sp.exe
gfortran  tcud2sp.f -o tcud2sp.exe -L../lib -l mudpack
./tcud2sp.exe


 cud2sp test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp =  2 jyq =  3 iex =  7 jey =  6
 nx = 129 ny =  97 iguess =  0 maxcy =  3
 method =  0 work space estimate =   64825

 multigrid option arguments 
 kcycle =  0
 iprer = -1
 ipost = -1
 intpol = **

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cud2sp intl =  0
 ierror =  0 minimum work space =   48976

 approximation call to cud2sp 
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.331E-03
rm -f tcud3.exe
gfortran  tcud3.f -o tcud3.exe -L../lib -l mudpack
./tcud3.exe


 cud3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  0 nzf =  0
 ixp =  2 jyq =  2 kzr =  3
 iex =  5 jey =  5 kez =  6
 nx =  33 ny =  33 nz =  97 iguess =  0
 maxcy =  1
 method =  3 work space length input = 1940400
 xa =  0.50 xb =  1.00
 yc =  0.50 yd =  1.00
 ze =  0.00 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to cud3 intl =  0
 ierror =  0 minimum work space = 1854498

 approximation call to cud3 
 intl =  1 method =  3 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.000E+00
rm -f tcud34.exe
gfortran  tcud34.f -o tcud34.exe -L../lib -l mudpack
./tcud34.exe


 cud3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  0 nzf =  0
 ixp =  2 jyq =  2 kzr =  3
 iex =  5 jey =  5 kez =  6
 nx =  33 ny =  33 nz =  97 iguess =  0
 maxcy =  1
 method =  3 work space length input = 1854498
 xa =  0.50 xb =  1.00
 yc =  0.50 yd =  1.00
 ze =  0.00 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to cud3 intl =  0
 ierror =  0 minimum work space = 1854498

 approximation call to cud3 
 intl =  1 method =  3 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.000E+00

 approximation call to cud3 
 intl =  1 method =  3 iguess =  1 maxcy =  3
 ierror =  0
 maximum error  =   0.000E+00

 cud24 test  ierror =  0
 maximum error  =   0.000E+00
rm -f tcud34sp.exe
gfortran  tcud34sp.f -o tcud34sp.exe -L../lib -l mudpack
./tcud34sp.exe


 cud3sp test 

 input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  2 nyd =  1
 nze =  2 nzf =  1
 ixp =  3 jyq =  2 kzr =  3
 iex =  5 jey =  5 kez =  5
 nx =  49 ny =  33 nz =  49 iguess =  0
 maxcy =  5
 method =  0 work space length input =  291605
 xa =  0.25 xb =  0.75
 yc =  0.33 yd =  0.67
 ze =  0.20 zf =  0.80
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to cud3sp intl =  0
 ierror =  0 minimum work space =  291605

 approximation call to cud3sp 
 intl =  1 method =  0 iguess =  0 maxcy =  5
 ierror =  0
 maximum error  =   0.325E-04

 cud34sp test  ierror =  0
 maximum error  =   0.169E-04
rm -f tcud3cr.exe
gfortran  tcud3cr.f -o tcud3cr.exe -L../lib -l mudpack
./tcud3cr.exe


 cud3cr test 

 input arguments 
 intl =  0
 nxa =  1 nxb =  2
 nyc =  0 nyd =  0
 nze =  1 nzf =  2
 ixp =  3 jyq =  2 kzr =  3
 iex =  3 jey =  6 kez =  3
 nx =  13 ny =  65 nz =  13
 iguess =  0 maxcy =  2
 method =  0 work space input =  173394
 xa =  0.00 xb =  1.00
 yc =  0.00 yd =  6.28
 ze =  0.00 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 mgopt(1)  =  0

 new cud3cr arguments 
 icrs 
    1    0    1
 tol = 0.0010
 maxit =  10

 initial call
 intl =  0
 iguess =  0
 ierror = -4
 minimum required work space length =   173394

 approximation call 
 intl =  1
 iguess =  0
 ierror =   0
 number of outer iterations executed =   3
 relative difference profile:
1.0000  1.0000  0.0000
 exact least squares error =        NaN
rm -f tcud3sp.exe
gfortran  tcud3sp.f -o tcud3sp.exe -L../lib -l mudpack
./tcud3sp.exe


 cud3sp test 

 input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  2 nyd =  1
 nze =  2 nzf =  1
 ixp =  3 jyq =  2 kzr =  3
 iex =  5 jey =  5 kez =  5
 nx =  49 ny =  33 nz =  49 iguess =  0
 maxcy =  3
 method =  0 work space length input =  318622
 xa =  0.25 xb =  0.75
 yc =  0.33 yd =  0.67
 ze =  0.20 zf =  0.80
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to cud3sp intl =  0
 ierror =  0 minimum work space =  291605

 approximation call to cud3sp 
 intl =  1 method =  0 iguess =  0 maxcy =  3
 ierror =  0
 maximum error  =   0.344E-04
rm -f tcuh2.exe
gfortran  tcuh2.f -o tcuh2.exe -L../lib -l mudpack
./tcuh2.exe


 cuh2 test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp = 11 jyq =  7 iex =  3 jey =  4
 nx =  45 ny =  57 iguess =  0 maxcy =  1
 method =  0 work space estimate =   34797

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cuh2 intl =  0
 ierror =  0 minimum work space =   30053

 approximation call to cuh2
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.805E-03
rm -f tcuh24.exe
gfortran  tcuh24.f -o tcuh24.exe -L../lib -l mudpack
./tcuh24.exe


 cuh2 test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp = 11 jyq =  7 iex =  3 jey =  4
 nx =  45 ny =  57 iguess =  0 maxcy =  3
 method =  0 work space estimate =   30053

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to cuh2 intl =  0
 ierror =  0 minimum work space =   30053

 approximation call to cuh2
 intl =  1 method =  0 iguess =  0 maxcy =  3
 ierror =  0
 maximum error  =   0.814E-03

 cuh24 test ierror =  0
 maximum error  =   0.104E-04
rm -f tcuh24cr.exe
gfortran  tcuh24cr.f -o tcuh24cr.exe -L../lib -l mudpack
./tcuh24cr.exe


 cuh2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp = 11 jyq =  9 iex =  3 jey =  4
 nx =  45 ny =  73 iguess =  0 maxcy =  2
 method =  2 work space estimate =   69661

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.250 xb =  1.000 yc =  0.000 yd =  0.500
 tolerance (error control) =    0.000E+00

 discretization call to cuh2cr intl =  0
 ierror =  0 minimum work space =   69661

 approximation call to cuh2cr
 intl =  1 method =  2 iguess =  0 maxcy =  2
 ierror =  0
 maximum error  =   0.112E-04

 cuh24cr test  ierror =  0
 maximum error  =   0.285E-06
rm -f tcuh2cr.exe
gfortran  tcuh2cr.f -o tcuh2cr.exe -L../lib -l mudpack
./tcuh2cr.exe


 cuh2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp = 11 jyq =  9 iex =  3 jey =  4
 nx =  45 ny =  73 iguess =  0 maxcy =  1
 method =  2 work space estimate =   72545

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.250 xb =  1.000 yc =  0.000 yd =  0.500
 tolerance (error control) =    0.000E+00

 discretization call to cuh2cr intl =  0
 ierror =  0 minimum work space =   69661

 approximation call to cuh2cr
 intl =  1 method =  2 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.110E-04
rm -f tcuh3.exe
gfortran  tcuh3.f -o tcuh3.exe -L../lib -l mudpack
./tcuh3.exe


 cuh3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  0 nzf =  0
 ixp =  5 jyq =  5 kzr =  7
 iex =  3 jey =  3 kez =  4
 nx =  21 ny =  21 nz =  57 iguess =  0
 maxcy =  1
 method =  3 work space length input =  499376
 xa =  0.50 xb =  1.00
 yc =  0.50 yd =  1.00
 ze =  0.00 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to cuh3 intl =  0
 ierror = -4 minimum work space =  488841

 approximation call to cuh3 
 intl =  1 method =  3 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.000E+00
rm -f tcuh34.exe
gfortran  tcuh34.f -o tcuh34.exe -L../lib -l mudpack
./tcuh34.exe


 cuh3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  0 nzf =  0
 ixp =  5 jyq =  5 kzr =  7
 iex =  3 jey =  3 kez =  4
 nx =  21 ny =  21 nz =  57 iguess =  0
 maxcy =  2
 method =  3 work space length input =  488841
 xa =  0.50 xb =  1.00
 yc =  0.50 yd =  1.00
 ze =  0.00 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to cuh3 intl =  0
 ierror = -4 minimum work space =  488841

 approximation call to cuh3 
 intl =  1 method =  3 iguess =  0 maxcy =  2
 ierror =  0
 maximum error  =   0.000E+00

 cuh34 test  ierror =  0
 maximum error  =   0.000E+00
rm -f tmud2.exe
gfortran  tmud2.f -o tmud2.exe -L../lib -l mudpack
./tmud2.exe


 mud2 test 

 integer input arguments 
 intl =  0 nxa =  2 nxb =  1 nyc =  1nyd =  2
 ixp =  3 jyq =  3 iex =  5 jey =  6
 nx =  49 ny =  97 iguess =  0 maxcy =  1
 method =  2 work space estimate =   83964

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to mud2 intl =  0
 ierror =  0 minimum work space =   70048

 approximation call to mud2
 intl =  1 method =  2 iguess =  0
 ierror =  0
 maximum error  =   0.331E-03
rm -f tmud24.exe
gfortran  tmud24.f -o tmud24.exe -L../lib -l mudpack
./tmud24.exe


 mud2 test 

 integer input arguments 
 intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp =  3 jyq =  3 iex =  5 jey =  6
 nx =  49 ny =  97 iguess =  0 maxcy =  2
 method =  2 work space estimate =   70048

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to mud2 intl =  0
 ierror =  0 minimum work space =   70048

 approximation call to mud2
 intl =  1 method =  2 iguess =  0
 ierror =  0
 maximum error  =   0.338E-03

 mud24 test  ierror =  0
 maximum error  =   0.513E-05
rm -f tmud24cr.exe
gfortran  tmud24cr.f -o tmud24cr.exe -L../lib -l mudpack
./tmud24cr.exe


 mud2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp =  3 jyq =  2 iex =  5 jey =  6
 nx =  49 ny =  65 iguess =  0 maxcy =  3
 method =  0 work space estimate =   51496

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to mud2cr intl =  0
 ierror =  0 minimum work space =   51496

 approximation call to mud2cr
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.624E-03

 mud24cr test  ierror =  0
 maximum error  =   0.389E-05
rm -f tmud24sp.exe
gfortran  tmud24sp.f -o tmud24sp.exe -L../lib -l mudpack
./tmud24sp.exe


 mud2sp test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp =  2 jyq =  3 iex =  6 jey =  5
 nx =  65 ny =  49 iguess =  0 maxcy =  3
 method =  0 work space estimate =   13264

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to mud2sp intl =  0
 ierror =  0 minimum work space =   13264

 approximation call to mud2sp
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.510E-04

 mud24sp test  ierror =  0
 maximum error  =   0.280E-05
rm -f tmud2cr.exe
gfortran  tmud2cr.f -o tmud2cr.exe -L../lib -l mudpack
./tmud2cr.exe


 mud2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp =  3 jyq =  2 iex =  5 jey =  6
 nx =  49 ny =  65 iguess =  0 maxcy =  1
 method =  0 work space estimate =   54686

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to mud2cr intl =  0
 ierror =  0 minimum work space =   51496

 approximation call to mud2cr
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.623E-03
rm -f tmud2sa.exe
gfortran  tmud2sa.f -o tmud2sa.exe -L../lib -l mudpack
./tmud2sa.exe


 mud2sa test

 integer input parameters 
intl =  0 nta =  1 ntb =  1 npc =  0 npd =  0
 itp =  5 jpq =  5 iet =  5 jep =  6
 nt =  81 np = 161 iguess =  0 maxcy =  1
 method =  2 work space estimate =  260820

 multigrid option parameters 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 ta =  0.000 tb =  3.142 pc =  0.000 pd =  6.283
 tolerance (error control) =   0.000E+00

 discretization call to mud2sa  intl =  0
 ierror =  0 minimum work space =  225200

 approximation call to mud2sa 
 intl =  1 method =  2 iguess =  0 maxcy =  1
 tolmax =  0.00
 ierror =  0
 maximum error  =  0.654E-04
rm -f tmud2sp.exe
gfortran  tmud2sp.f -o tmud2sp.exe -L../lib -l mudpack
./tmud2sp.exe


 mud2sp test 

 integer input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
 ixp =  2 jyq =  3 iex =  6 jey =  5
 nx =  65 ny =  49 iguess =  0 maxcy =  1
 method =  0 work space estimate =   17065

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to mud2sp intl =  0
 ierror =  0 minimum work space =   13264

 approximation call to mud2sp
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.520E-04
rm -f tmud3.exe
gfortran  tmud3.f -o tmud3.exe -L../lib -l mudpack
./tmud3.exe


 mud3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  2 nzf =  2
 ixp =  2 jyq =  3 kzr =  2
 iex =  5 jey =  4 kez =  6
 nx =  33 ny =  25 nz =  65 iguess =  0
 maxcy =  1
 method =  3 work space length input =  949725
 xa =  0.50 xb =  1.00
 yc =  1.00 yd =  2.00
 ze =  0.25 zf =  0.75
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to mud3 intl =  0
 ierror =  0 minimum work space =  824224

 approximation call to mud3 
 intl =  1 method =  3 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.117E-04
rm -f tmud34.exe
gfortran  tmud34.f -o tmud34.exe -L../lib -l mudpack
./tmud34.exe


 mud3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  2 nzf =  2
 ixp =  2 jyq =  3 kzr =  2
 iex =  5 jey =  4 kez =  6
 nx =  33 ny =  25 nz =  65 iguess =  0
 maxcy =  1
 method =  3 work space length input =  824224
 xa =  0.50 xb =  1.00
 yc =  1.00 yd =  2.00
 ze =  0.25 zf =  0.75
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to mud3 intl =  0
 ierror =  0 minimum work space =  824224

 approximation call to mud3 
 intl =  1 method =  3 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.117E-04

 approximation call to mud3 
 intl =  1 method =  3 iguess =  1 maxcy =  2
 ierror =  0
 maximum error  =   0.107E-04

 mud34 test  ierror =  0
 maximum error  =   0.739E-05
rm -f tmud34sp.exe
gfortran  tmud34sp.f -o tmud34sp.exe -L../lib -l mudpack
./tmud34sp.exe


 mud3sp test 

 input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  2 nyd =  1
 nze =  2 nzf =  1
 ixp =  2 jyq =  3 kzr =  5
 iex =  6 jey =  5 kez =  4
 nx =  65 ny =  49 nz =  41 iguess =  0
 maxcy =  3
 method =  0 work space length input =  514258
 xa =  0.50 xb =  1.00
 yc =  0.50 yd =  1.00
 ze =  0.50 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to mud3sp intl =  0
 ierror =  0 minimum work space =  472557

 approximation call to mud3sp 
 intl =  1 method =  0 iguess =  0 maxcy =  3
 ierror =  0
 maximum error  =   0.719E-05

 mud34sp test  ierror =  0
 maximum error  =   0.775E-06
rm -f tmud3cr.exe
gfortran  tmud3cr.f -o tmud3cr.exe -L../lib -l mudpack
./tmud3cr.exe


 mud3cr test 

 input arguments 
 intl =  0
 nxa =  1 nxb =  2
 nyc =  0 nyd =  0
 nze =  1 nzf =  2
 ixp =  3 jyq =  2 kzr =  3
 iex =  4 jey =  7 kez =  4
 nx =  25 ny = 129 nz =  25
 iguess =  0 maxcy =  2
 method =  0 work space input = 1176431
 xa =  0.00 xb =  1.00
 yc =  0.00 yd =  6.28
 ze =  0.00 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 mgopt(1)  =  0

 new mud3cr arguments 
 icrs 
    1    0    1
 tol = 0.0010
 maxit =  10

 initial call
 intl =  0
 iguess =  0
 ierror =  0
 minimum required work space length =  1176431

 approximation call 
 intl =  1
 iguess =  0
 ierror =   0
 number of outer iterations executed =   8
 relative difference profile:
0.2694  0.0853  0.0357  0.0135  0.0075  0.0031  0.0016  0.0008
 exact least squares error =  0.228E-03
rm -f tmud3sa.exe
gfortran  tmud3sa.f -o tmud3sa.exe -L../lib -l mudpack
./tmud3sa.exe


 mud3sa test

 input arguments
intl =  0 nra =  1 nrb =  1 ntc =  1 ntd =  1
 npe =  0 npf =  0
 irp =  3 jtq =  2 kpr =  2
 ier =  5 jet =  5 kep =  6
 nr =  49 nt =  33 np =  65 iguess =  0
 maxcy =  1
 method =  1 work space length input = 1598164
 ra =  0.50 rb =  1.00
 tc =  0.79 td =  2.36
 pe =  0.00 pf =  6.28
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to mud3sa intl =  0
 ierror =  0 minimum work space = 1598164

 approximation call to mud3sa 
 intl =  1 method =  1 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.369E-03
rm -f tmud3sp.exe
gfortran  tmud3sp.f -o tmud3sp.exe -L../lib -l mudpack
./tmud3sp.exe


 mud3sp test 

 input arguments 
intl =  0 nxa =  2 nxb =  1 nyc =  2 nyd =  1
 nze =  2 nzf =  1
 ixp =  2 jyq =  3 kzr =  5
 iex =  6 jey =  5 kez =  4
 nx =  65 ny =  49 nz =  41 iguess =  0
 maxcy =  3
 method =  0 work space length input =  514258
 xa =  0.50 xb =  1.00
 yc =  0.50 yd =  1.00
 ze =  0.50 zf =  1.00
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to mud3sp intl =  0
 ierror =  0 minimum work space =  472557

 approximation call to mud3sp 
 intl =  1 method =  0 iguess =  0 maxcy =  3
 ierror =  0
 maximum error  =   0.719E-05
rm -f tmuh2.exe
gfortran  tmuh2.f -o tmuh2.exe -L../lib -l mudpack
./tmuh2.exe


 muh2 test

 integer input parameters 
intl =  0 nta =  1 ntb =  1 npc =  0 npd =  0
 itp =  9 jpq =  9 iet =  4 jep =  5
 nt =  73 np = 145 iguess =  0 maxcy =  1
 method =  2 work space estimate =  214396

 multigrid option parameters 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 ta =  0.000 tb =  3.142 pc =  0.000 pd =  6.283
 tolerance (error control) =   0.000E+00

 discretization call to muh2  intl =  0
 ierror =  0 minimum work space =  186661

 approximation call to muh2 
 intl =  1 method =  2 iguess =  0 maxcy =  1
 tolmax =  0.00
 ierror =  0
 maximum error  =  0.795E-04
rm -f tmuh24.exe
gfortran  tmuh24.f -o tmuh24.exe -L../lib -l mudpack
./tmuh24.exe


 muh2 test

 integer input parameters 
intl =  0 nta =  1 ntb =  1 npc =  0 npd =  0
 itp =  9 jpq =  9 iet =  4 jep =  5
 nt =  73 np = 145 iguess =  0 maxcy =  1
 method =  2 work space estimate =  214396

 multigrid option parameters 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 ta =  0.000 tb =  3.142 pc =  0.000 pd =  6.283
 tolerance (error control) =   0.000E+00

 discretization call to muh2  intl =  0
 ierror =  0 minimum work space =  186661

 approximation call to muh2 
 intl =  1 method =  2 iguess =  0 maxcy =  1
 tolmax =  0.00
 ierror =  0
 maximum error  =  0.795E-04

 approximation call to muh2 
 intl =  1 method =  2 iguess =  1 maxcy =  2
 tolmax =  0.00
 ierror =  0
 maximum error  =  0.839E-04

 muh24 test  ierror =  0
 maximum error  =  0.212E-06
rm -f tmuh24cr.exe
gfortran  tmuh24cr.f -o tmuh24cr.exe -L../lib -l mudpack
./tmuh24cr.exe


 muh2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp = 15 jyq =  9 iex =  3 jey =  5
 nx =  61 ny = 145 iguess =  0 maxcy =  3
 method =  0 work space estimate =  149055

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to muh2cr intl =  0
 ierror =  0 minimum work space =  149055

 approximation call to muh2cr
 intl =  1 method =  0 iguess =  0 maxcy =  3
 ierror =  0
 maximum error  =   0.411E-03

 muh24cr test  ierror =  0
 maximum error  =   0.218E-05
rm -f tmuh2cr.exe
gfortran  tmuh2cr.f -o tmuh2cr.exe -L../lib -l mudpack
./tmuh2cr.exe


 muh2cr test 

 integer input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
 ixp = 15 jyq =  9 iex =  3 jey =  5
 nx =  61 ny = 145 iguess =  0 maxcy =  1
 method =  0 work space estimate =  151335

 multigrid option arguments 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 floating point input parameters 
 xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
 tolerance (error control) =    0.000E+00

 discretization call to muh2cr intl =  0
 ierror =  0 minimum work space =  149055

 approximation call to muh2cr
 intl =  1 method =  0 iguess =  0
 ierror =  0
 maximum error  =   0.411E-03
rm -f tmuh3.exe
gfortran  tmuh3.f -o tmuh3.exe -L../lib -l mudpack
./tmuh3.exe


 muh3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  2 nzf =  2
 ixp =  5 jyq =  5 kzr =  7
 iex =  4 jey =  3 kez =  4
 nx =  41 ny =  21 nz =  57 iguess =  0
 maxcy =  1
 method =  3 work space length input =  779587
 xa =  0.50 xb =  1.00
 yc =  1.00 yd =  2.00
 ze =  0.25 zf =  0.75
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to muh3 intl =  0
 ierror =  0 minimum work space =  776711

 approximation call to muh3 
 intl =  1 method =  3 iguess =  0 maxcy =  1
 ierror =  0
 maximum error  =   0.200E-04
rm -f tmuh34.exe
gfortran  tmuh34.f -o tmuh34.exe -L../lib -l mudpack
./tmuh34.exe


 muh3 test 

 input arguments 
intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  1
 nze =  2 nzf =  2
 ixp =  5 jyq =  5 kzr =  7
 iex =  4 jey =  3 kez =  4
 nx =  41 ny =  21 nz =  57 iguess =  0
 maxcy =  2
 method =  3 work space length input =  776771
 xa =  0.50 xb =  1.00
 yc =  1.00 yd =  2.00
 ze =  0.25 zf =  0.75
 tolmax =  0.000E+00

 multigrid options 
 kcycle =  2
 iprer =  2
 ipost =  1
 intpol =  3

 discretization call to muh3 intl =  0
 ierror =  0 minimum work space =  776711

 approximation call to muh3 
 intl =  1 method =  3 iguess =  0 maxcy =  2
 ierror =  0
 maximum error  =   0.212E-04

 muh34 test ierror =  0
 maximum error  =   0.727E-05
