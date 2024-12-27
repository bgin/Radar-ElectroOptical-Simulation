mkdir -p ./lib
mkdir -p ./objs
( cd ./src; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f ../lib/libfish90.a ../objs/fish.o ../objs/blktri.o ../objs/cblktri.o ../objs/cmgnbn.o ../objs/comf.o ../objs/fftpack.o ../objs/genbun.o ../objs/gnbnaux.o ../objs/hstcrt.o ../objs/hstcsp.o ../objs/hstcyl.o ../objs/hstplr.o ../objs/hstssp.o ../objs/hw3crt.o ../objs/hwscrt.o ../objs/hwscsp.o ../objs/hwscyl.o ../objs/hwsplr.o ../objs/hwsssp.o ../objs/pois3d.o ../objs/poistg.o ../objs/sepaux.o ../objs/sepeli.o ../objs/sepx4.o ../lib/fish.mod ../lib/FISH.mod fish.f
cp fish.F.orig fish.F
gfortran -I../lib -c fish.F -o ../objs/fish.o; mv FISH.mod ../lib
gfortran -I../lib -c blktri.f -o ../objs/blktri.o
gfortran -I../lib -c cblktri.f -o ../objs/cblktri.o
gfortran -I../lib -c cmgnbn.f -o ../objs/cmgnbn.o
gfortran -I../lib -c comf.f -o ../objs/comf.o
gfortran -I../lib -c fftpack.f -o ../objs/fftpack.o
gfortran -I../lib -c genbun.f -o ../objs/genbun.o
gfortran -I../lib -c gnbnaux.f -o ../objs/gnbnaux.o
gfortran -I../lib -c hstcrt.f -o ../objs/hstcrt.o
gfortran -I../lib -c hstcsp.f -o ../objs/hstcsp.o
gfortran -I../lib -c hstcyl.f -o ../objs/hstcyl.o
gfortran -I../lib -c hstplr.f -o ../objs/hstplr.o
gfortran -I../lib -c hstssp.f -o ../objs/hstssp.o
gfortran -I../lib -c hw3crt.f -o ../objs/hw3crt.o
gfortran -I../lib -c hwscrt.f -o ../objs/hwscrt.o
gfortran -I../lib -c hwscsp.f -o ../objs/hwscsp.o
gfortran -I../lib -c hwscyl.f -o ../objs/hwscyl.o
gfortran -I../lib -c hwsplr.f -o ../objs/hwsplr.o
gfortran -I../lib -c hwsssp.f -o ../objs/hwsssp.o
gfortran -I../lib -c pois3d.f -o ../objs/pois3d.o
gfortran -I../lib -c poistg.f -o ../objs/poistg.o
gfortran -I../lib -c sepaux.f -o ../objs/sepaux.o
gfortran -I../lib -c sepeli.f -o ../objs/sepeli.o
gfortran -I../lib -c sepx4.f -o ../objs/sepx4.o
/usr/bin/ar -rv ../lib/libfish90.a ../objs/fish.o ../objs/blktri.o ../objs/cblktri.o ../objs/cmgnbn.o ../objs/comf.o ../objs/fftpack.o ../objs/genbun.o ../objs/gnbnaux.o ../objs/hstcrt.o ../objs/hstcsp.o ../objs/hstcyl.o ../objs/hstplr.o ../objs/hstssp.o ../objs/hw3crt.o ../objs/hwscrt.o ../objs/hwscsp.o ../objs/hwscyl.o ../objs/hwsplr.o ../objs/hwsssp.o ../objs/pois3d.o ../objs/poistg.o ../objs/sepaux.o ../objs/sepeli.o ../objs/sepx4.o 
ar: creating archive ../lib/libfish90.a
a - ../objs/fish.o
a - ../objs/blktri.o
a - ../objs/cblktri.o
a - ../objs/cmgnbn.o
a - ../objs/comf.o
a - ../objs/fftpack.o
a - ../objs/genbun.o
a - ../objs/gnbnaux.o
a - ../objs/hstcrt.o
a - ../objs/hstcsp.o
a - ../objs/hstcyl.o
a - ../objs/hstplr.o
a - ../objs/hstssp.o
a - ../objs/hw3crt.o
a - ../objs/hwscrt.o
a - ../objs/hwscsp.o
a - ../objs/hwscyl.o
a - ../objs/hwsplr.o
a - ../objs/hwsssp.o
a - ../objs/pois3d.o
a - ../objs/poistg.o
a - ../objs/sepaux.o
a - ../objs/sepeli.o
a - ../objs/sepx4.o
( cd ./test; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f  tblktri.exe tcblktri.exe tcmgnbn.exe tgenbun.exe thstcrt.exe thstcsp.exe thstcyl.exe thstplr.exe thstssp.exe thw3crt.exe thwscrt.exe thwscsp.exe thwscyl.exe thwsplr.exe thwsssp.exe tpois3d.exe tpoistg.exe tsepeli.exe tsepx4.exe
rm -f tblktri.exe
gfortran -I../lib -I../lib tblktri.f -o tblktri.exe -L../lib -l fish90
./tblktri.exe
     BLKTRI TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.6478E-05
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.2737E-02
     The output from your computer is: 
     IERROR =           0  Discretization Error =   1.27372220E-02
rm -f tcblktri.exe
gfortran -I../lib -I../lib tcblktri.f -o tcblktri.exe -L../lib -l fish90
./tcblktri.exe
     CBLKTRI TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.6457E-05
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.2737E-02
     The output from your computer is: 
     IERROR =           0  Discretization Error =   1.27370311E-02
rm -f tcmgnbn.exe
gfortran -I../lib -I../lib tcmgnbn.f -o tcmgnbn.exe -L../lib -l fish90
./tcmgnbn.exe
     CMGNBN TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 9.1620E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 9.1801E-3
     The output from your computer is: 
     IERROR =           0  Discretization Error =   9.18554235E-03
rm -f tgenbun.exe
gfortran -I../lib -I../lib tgenbun.f -o tgenbun.exe -L../lib -l fish90
./tgenbun.exe
     GENBUN TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 9.6406E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 9.6556E-3
     The output from your computer is: 
     IERROR =           0  Discretization Error =   9.65690613E-03
rm -f thstcrt.exe
gfortran -I../lib -I../lib thstcrt.f -o thstcrt.exe -L../lib -l fish90
./thstcrt.exe
     HSTCRT TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.2600E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.2586E-3
     The output from your computer is: 
     IERROR =           0  Discretization Error =   1.25920773E-03
rm -f thstcsp.exe
gfortran -I../lib -I../lib thstcsp.f -o thstcsp.exe -L../lib -l fish90
./thstcsp.exe
     HSTCSP TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 5.5843E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 5.5845E-3
     The output from your computer is: 
     IERROR =           0  Discretization Error =   5.58453798E-03
rm -f thstcyl.exe
gfortran -I../lib -I../lib thstcyl.f -o thstcyl.exe -L../lib -l fish90
./thstcyl.exe
     HSTCYL TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  PERTRB =-4.4311E-4
     Discretization Error = 7.5280E-5 
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  PERTRB =-4.4321E-4
     Discretization Error = 7.3557E-5
     The output from your computer is: 
     IERROR =           0  PERTRB =  -4.43209137E-04
     Discretization Error =   7.31476321E-05
rm -f thstplr.exe
gfortran -I../lib -I../lib thstplr.f -o thstplr.exe -L../lib -l fish90
./thstplr.exe
     HSTPLR TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.1303E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 1.1300E-3
     The output from your computer is: 
     IERROR =           0  Discretization Error =   1.12986565E-03
rm -f thstssp.exe
gfortran -I../lib -I../lib thstssp.f -o thstssp.exe -L../lib -l fish90
./thstssp.exe
     HSTSSP TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  PERTRB = 6.35830E-4
     discretization error = 3.37523E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  PERTRB = 6.35919E-4
     discretization error = 3.38144E-3
     The output from your computer is: 
     IERROR =           0  PERTRB =   6.35743549E-04
     discretization error =   3.37797590E-03
rm -f thw3crt.exe
gfortran -I../lib -I../lib thw3crt.f -o thw3crt.exe -L../lib -l fish90
./thw3crt.exe
     HW3CRT TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 9.6480E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 9.6480E-3
     The output from your computer is: 
     IERROR =           0  Discretization Error =   9.64805484E-03
rm -f thwscrt.exe
gfortran -I../lib -I../lib thwscrt.f -o thwscrt.exe -L../lib -l fish90
./thwscrt.exe
     HWSCRT TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 5.36508-4
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 4.9305E-4
     The output from your computer is: 
     IERROR =           0  Discretization Error =   5.11169434E-04
rm -f thwscsp.exe
gfortran -I../lib -I../lib thwscsp.f -o thwscsp.exe -L../lib -l fish90
./thwscsp.exe
  HWSCSP TEST RUN, EXAMPLE 1 *** 
  Previous 64 bit floating point arithmetic result 
  ierror = 0
  discretization error = 7,9984E-4 
  Previous 32 bit floating point arithmetic result 
  ierror = 0
  discretization error = 7.9907E-4 
  The output from your computer is: 
  IERROR =           0
  Discretization Error =  7.99238682E-04
  ********** 
  ********** 
  HWSCSP TEST RUN, EXAMPLE 2 *** 
  Previous 64 bit floating point arithmetic result 
  ierror = 0
  discretization error = 5.8682E-5 
  Previous 32 bit floating point arithmetic result 
  ierror = 0
  discretization error = 5.9962E-5 
  The output from your computer is: 
  IERROR =           0
  Discretization Error =  6.08861446E-05
rm -f thwscyl.exe
gfortran -I../lib -I../lib thwscyl.f -o thwscyl.exe -L../lib -l fish90
./thwscyl.exe
     HWSCYL TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  PERTRB = 2.2674E-4
     Discretization Error = 3.7367E-4 
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  PERTRB = 2.26976-4
     Discretization Error = 3.5554E-4
     The output from your computer is: 
     IERROR =           0  PERTRB =   2.26975826E-04
     Discretization Error =   3.62455845E-04
rm -f thwsplr.exe
gfortran -I../lib -I../lib thwsplr.f -o thwsplr.exe -L../lib -l fish90
./thwsplr.exe
     HWSPLR TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 6.19134E-4
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 6.20723E-4
     The output from your computer is: 
     IERROR =           0  Discretization Error =   6.22600317E-04
rm -f thwsssp.exe
gfortran -I../lib -I../lib thwsssp.f -o thwsssp.exe -L../lib -l fish90
./thwsssp.exe
     HWSSSP TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 3.38107E-3
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 3.3650E-3
     The output from your computer is: 
     IERROR =           0  Discretization Error =    1.0304666    
rm -f tpois3d.exe
gfortran -I../lib -I../lib tpois3d.f -o tpois3d.exe -L../lib -l fish90
./tpois3d.exe
     POIS3D TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 2.93277E-2
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 2.93390E-2
     The output from your computer is: 
     IERROR =           0  Discretization Error =   2.93388367E-02
rm -f tpoistg.exe
gfortran -I../lib -I../lib tpoistg.f -o tpoistg.exe -L../lib -l fish90
./tpoistg.exe
     POISTG TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 5.6417E-4
     Previous 32 bit floating point arithmetic result 
     IERROR = 0,  Discretization Error = 5.6183E-4
     The output from your computer is: 
     IERROR =           0  Discretization Error =   5.63159585E-04
rm -f tsepeli.exe
gfortran -I../lib -I../lib tsepeli.f -o tsepeli.exe -L../lib -l fish90
./tsepeli.exe
     SEPELI TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0
     Second Order Discretization Error = 9.7891E-5
     Fourth Order Discretization Error = 1.4735E-6
     Previous 32 bit floating point arithmetic result 
     IERROR = 0
     Second Order Discretization Error = 1.2708E-4
     Fourth Order Discretization Error = 3.1948E-5
     The output from your computer is: 
     IERROR =           0
     Second Order Discretization Error =  1.24454498E-04
     Fourth Order Discretization Error =  2.92062759E-05
rm -f tsepx4.exe
gfortran -I../lib -I../lib tsepx4.f -o tsepx4.exe -L../lib -l fish90
./tsepx4.exe
     SEPEX4 TEST RUN *** 
     Previous 64 bit floating point arithmetic result 
     IERROR = 0
     Second Order Discretization Error = 1.5985E-4
     Fourth Order Discretization Error = 1.8575E-6
     Previous 32 bit floating point arithmetic result 
     IERROR = 0
     Second Order Discretization Error = 1.5044E-4
     Fourth Order Discretization Error = 1.5736E-5
     The output from your computer is: 
     IERROR =           0
     Second Order Discretization Error =  1.55925751E-04
     Fourth Order Discretization Error =  1.62124634E-05
