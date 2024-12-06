
include ../make.inc

SRC=cud2.f cud24.f cud24cr.f cud24sp.f cud2cr.f cud2sp.f cud3.f cud34.f \
    cud34sp.f cud3cr.f cud3ln.f cud3pn.f cud3sp.f cudcom.f cuh2.f cuh24.f \
    cuh24cr.f cuh2cr.f cuh3.f cuh34.f mud2.f mud24.f mud24cr.f mud24sp.f \
    mud2cr.f mud2sa.f mud2sp.f mud3.f mud34.f mud34sp.f mud3cr.f mud3ln.f \
    mud3pn.f mud3sa.f mud3sp.f mudcom.f muh2.f muh24.f muh24cr.f muh2cr.f \
    muh3.f muh34.f resc2.f resc2cr.f resc2sp.f resc3.f resc3sp.f resm2.f \
    resm2cr.f resm2sp.f resm3.f resm3cr.f resm3sp.f

OBJ=$(subst .f,.o,$(SRC))
OBJS=$(addprefix ../objs/,$(OBJ))

$(LIB) : $(OBJS)
	$(AR) -rv $@ $? 

../objs/%.o : %.f
	$(F90) -c $< -o $@

clean:
	rm -f $(LIB) $(OBJS) 
