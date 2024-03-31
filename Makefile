CC= cc
FC= gfortran
CLIBS= -lm

OBJ = nitfread.o rwfits.o psf2d.o fitlm.o linsolve.o sortc.o

# Uncomment these lines to build without CFITSIO support.
#CFLAGS= -O -g -Wall 
#LIBS = -lm

# Uncomment these lines to include CFITSIO support.
CFLAGS= -O -g -Wall -DUSECFITSIO -I$(ATLAS_HOME)/include

# Uncomment these lines to include CFITSIO support and OpenMP multithreading
# CFLAGS= -O -g -Wall -DUSECFITSIO -DUSEMP -fopenmp

OBJ += cfitsrw.o
#LIBS = -lcfitsio -lm
LIBS = -lm

all: tphot tphot_mp

tphot_mp.o: tphot.c psf2d.h
	$(CC) $(CFLAGS) -DUSEMP -fopenmp -c tphot.c -o tphot_mp.o

.c.o:
	$(CC) $(CFLAGS) -c $<

tphot: $(OBJ) tphot.o psf2d.h
	$(CC) $(CFLAGS) tphot.o $(OBJ) $(LIBS) -o tphot $(ATLAS_HOME)/lib/libcfitsio.a 

tphot_mp: $(OBJ) tphot_mp.o psf2d.h
	$(CC) $(CFLAGS) -DUSEMP -fopenmp tphot_mp.o $(OBJ) $(LIBS) -o tphot_mp $(ATLAS_HOME)/lib/libcfitsio.a 

tphot35: tphot35.o nitfread.o rwfits.o psf2d.o psf35mm.o fitlm.o linsolve.o sortc.o
	$(CC) $(CFLAGS) tphot35.o nitfread.o rwfits.o \
	 psf2d.o psf35mm.o fitlm.o linsolve.o sortc.o -lm -o tphot35

install:
	install -m 0775 tphot tphot_mp $(ATLAS_HOME)/bin/
	install -m 0644 tphot.man $(ATLAS_HOME)/man/man1/tphot.1

clean:
	rm -f *.o tphot tphot_mp tphot35
