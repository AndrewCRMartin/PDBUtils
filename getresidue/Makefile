LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include

CC = cc
OFILES = getresidue.o
LIBS   = -lbiop -lgen -lm -lxml2
#CFLAGS = -g -ansi -Wall -DDEBUG=1
#CFLAGS = -g -ansi -Wall
CFLAGS = -O3 -ansi -Wall

getresidue : $(OFILES)
	$(CC) $(CFLAGS) -o $@ $(OFILES) -L $(LIBDIR) $(LIBS)

.c.o :
	$(CC) $(CFLAGS) -c $< -I $(INCDIR)

clean :
	rm $(OFILES)

install :
	cp getresidue $(HOME)/bin
