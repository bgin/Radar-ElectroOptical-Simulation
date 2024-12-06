
include ../make.inc

SRC=tcud2.f tcud24.f tcud24cr.f tcud24sp.f tcud2cr.f tcud2sp.f tcud3.f \
    tcud34.f tcud34sp.f tcud3cr.f tcud3sp.f tcuh2.f tcuh24.f tcuh24cr.f \
    tcuh2cr.f tcuh3.f tcuh34.f tmud2.f tmud24.f tmud24cr.f tmud24sp.f \
    tmud2cr.f tmud2sa.f tmud2sp.f tmud3.f tmud34.f tmud34sp.f tmud3cr.f \
    tmud3sa.f tmud3sp.f tmuh2.f tmuh24.f tmuh24cr.f tmuh2cr.f tmuh3.f \
    tmuh34.f

EXES=$(subst .f,.exe, $(SRC))

all : $(EXES)

$(EXES) : $(SRC)

%.exe : %.f
	rm -f $@
	$(F90) $< -o $@ -L../lib -l mudpack
	./$@


clean :
	rm -f $(EXES)
