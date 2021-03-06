# ===========================================================================
#  Makefile for MUSIC                                 Chun Shen Dec. 12, 2014
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++-mp-5
CFLAGS= -Wall -O3 -fopenmp
#CFLAGS= -Wall -g -fopenmp

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	mpihydro
endif

SRC		=	main.cpp music.cpp grid.cpp eos.cpp init.cpp \
            util.cpp evolve.cpp advance.cpp u_derivative.cpp \
			reconst.cpp dissipative.cpp minmod.cpp grid_info.cpp \
			cornelius.cpp read_in_parameters.cpp

INC		= 	music.h grid.h eos.h init.h util.h data.h \
			evolve.h advance.h u_derivative.h reconst.h dissipative.h \
			minmod.h grid_info.h cornelius.h read_in_parameters.h emoji.h \
			field.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
