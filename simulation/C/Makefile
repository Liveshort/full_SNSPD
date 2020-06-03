CC     = gcc
FC     = gfortran
CFLAGS = -g -I -Wall -Wextra -fopenmp -std=c11 -fPIC -Ofast
FFLAGS = -g -fopenmp

DEPDIR = include
OBJDIR = obj
LIBDIR = lib

CBLSLB = $(LIBDIR)/libcblas.a
BLSLB  = $(LIBDIR)/librefblas.a
LPKLB  = $(LIBDIR)/liblapack.a
LPKELB = $(LIBDIR)/liblapacke.a

OPNBLS = $(LIBDIR)/libopenblas.a

ifneq ("$(wildcard $(OPNBLS))","")
BLSDIR = include_openblas
LIBS   = $(OPNBLS) -lm -lpthread
else
BLSDIR = include_netlib
LIBS   = -lm $(LPKELB) $(LPKLB) $(BLSLB) $(CBLSLB)
endif

_DEPS  = snspd.h thermal.h electrical.h helper.h transmission.h types.h linalg.h
DEPS   = $(patsubst %,$(DEPDIR)/%,$(_DEPS))

_OBJ   = types.o linalg.o helper.o transmission.o thermal.o electrical.o snspd.o simulation.o
OBJ    = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

$(OBJDIR)/%.o: %.c $(DEPS)
	$(CC) -I $(DEPDIR) -I $(BLSDIR) $(CFLAGS) -c -o $@ $<

all: simulation

simulation: $(OBJ) $(LIBS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o bld/$@ $^

run:
	./bld/simulation ${args}

.PHONY: clean all run

clean:
	rm -f $(OBJDIR)/*.o bld/snspd bld/snspd.so bld/simulation
