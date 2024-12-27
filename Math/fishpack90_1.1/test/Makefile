
include ../make.inc

SRC=tblktri.f tcblktri.f tcmgnbn.f tgenbun.f \
    thstcrt.f thstcsp.f thstcyl.f thstplr.f thstssp.f \
    thw3crt.f thwscrt.f thwscsp.f thwscyl.f thwsplr.f \
    thwsssp.f tpois3d.f tpoistg.f tsepeli.f tsepx4.f

EXES=$(subst .f,.exe, $(SRC))

all : $(EXES)

$(EXES) : $(SRC)

%.exe : %.f
	rm -f $@
	$(F90) -I../lib $< -o $@ -L../lib -l fish90
	./$@


clean :
	rm -f $(EXES)
