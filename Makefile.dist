EXE    = pseudotorsion
OFILES = pseudotorsion.o
LFILES = bioplib/ReadPDB.o \
         bioplib/SelectCaPDB.o \
         bioplib/chindex.o \
         bioplib/padterm.o \
         bioplib/FindNextResidue.o \
         bioplib/fsscanf.o \
         bioplib/phi.o
LIBS   = -lm
CC     = gcc
LOPTS  =

$(EXE) : $(OFILES) $(LFILES)
	$(CC) $(LOPTS) -o $(EXE) $(OFILES) $(LFILES) $(LIBS)

.c.o   :
	$(CC) -o $@ -c $<

