
include ../make.inc

SRC=fish.f blktri.f cblktri.f cmgnbn.f comf.f fftpack.f \
    genbun.f gnbnaux.f hstcrt.f hstcsp.f hstcyl.f hstplr.f \
    hstssp.f hw3crt.f hwscrt.f hwscsp.f hwscyl.f hwsplr.f \
    hwsssp.f pois3d.f poistg.f sepaux.f sepeli.f sepx4.f

OBJ=$(subst .f,.o,$(SRC))
OBJS=$(addprefix ../objs/,$(OBJ))


$(LIB) : $(OBJS)
	$(AR) -rv $@ $? 

ifeq ($(UNAMES),AIX)

fish.f : fish.F
	$(CPP) $< ; mv F$@ $@

endif

ifeq ($(UNAMES),SunOS)

fish.f : fish.F
	$(CPP) $<

endif

ifeq ($(UNAMES),Linux)

fish.f : fish.F
	$(CPP) $< > $@

endif

ifeq ($(UNAMES),IRIX64)

fish.f : fish.F
	$(CPP) $< > $@

../objs/fish.o : fish.f
	$(F90) -c $< -o $@; mv FISH.mod ../lib

../objs/%.o : %.f
	$(F90) -c $< -o $@

endif

ifeq ($(UNAMES),Darwin)

fish.F : fish.F.orig
	cp fish.F.orig fish.F

../objs/fish.o : fish.F
	$(F90) -c $< -o $@; mv FISH.mod ../lib

../objs/%.o : %.f
	$(F90) -c $< -o $@

endif

ifneq ($(UNAMES),IRIX64)

../objs/%.o : %.f
	$(F90) -c $< -o $@

endif

clean:
	rm -f $(LIB) $(OBJS) ../lib/fish.mod ../lib/FISH.mod fish.f
