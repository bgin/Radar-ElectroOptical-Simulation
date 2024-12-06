
include make.inc

all: libmudpack testmudpack

libmudpack:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) )

testmudpack:
	( cd ./test; $(MAKE) clean; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )
