CC     = gcc
COPT   = -O3 -I$(HOME)/include -L$(HOME)/lib
EXE    = protrusion
OFILES = protrusion.o
LIBS   = -lbiop -lgen -lm -lxml2
LFILES = 

$(EXE) : $(OFILES) $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f $(OFILES) $(LFILES)

distclean : clean
	\rm -f $(EXE)

