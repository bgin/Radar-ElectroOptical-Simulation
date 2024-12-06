
include make.inc

all: libfish testfish

libfish:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) )

testfish:
	( cd ./test; $(MAKE) clean; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )
