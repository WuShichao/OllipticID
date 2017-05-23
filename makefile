#
# Olliptic makefile
#

# ---------------------------------------------------------------
# Source, header and object file names
# ---------------------------------------------------------------

SRCS  = main.cc \
	interface.cc \
	domain.cc \
	ffunction.cc \
	interpolation.cc\
	scalarfields.cc \
	multigrid.cc \
	elliptic.cc \
	operator.cc \
	output.cc \
	boundary.cc \
	diff_2nd.cc \
	diff_4th.cc \
	diff_6th.cc \
	diff_8th.cc


# header files

HDRS  = domain.h \
	ffunction.h\
	elliptic.h\
	scalarfields.h \
	multigrid.h \
	interface.h

# object files

OBJS  = main.o \
	interface.o \
	domain.o  \
	ffunction.o \
	interpolation.o\
	scalarfields.o \
	multigrid.o\
	elliptic.o \
	operator.o \
	output.o \
	boundary.o \
	diff_2nd.o \
	diff_4th.o \
	diff_6th.o \
	diff_8th.o

# ---------------------------------------------------------------
# Set here exe name
# ---------------------------------------------------------------

EXE = Olliptic

# ---------------------------------------------------------------
# Set here compiler options
# ---------------------------------------------------------------

#CPP      = mpic++
#CPP      = icpc
CPP	= g++

CPPFLAGS = -O3 -ferror-limit=1  ###-pg 
debug: CPPFLAGS = -O0 -g -pg -Wall


#MPIFLAG = -DOLLIN_MPI 

#LDFLAGS = -L/opt/intel/Compiler/11.1/064/lib/intel64 -lifcore\
#	-L/opt/intel/Compiler/11.1/064/mkl/lib/em64t\
#	-L/opt/intel/Compiler/11.1/064/mkl/lib/em64t -lmkl -lguide -lpthread

LDFLAGS = -lm 

# ---------------------------------------------------------------
# Dir
# ---------------------------------------------------------------

BASE = $(shell /bin/pwd)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src
EXED = $(BASE)/exe
HDRSD = $(BASE)/include

vpath %.o $(OBJD)
vpath %.cc $(SRCD)
vpath %.h $(HDRSD)



# ---------------------------------------------------------------
# target 1 : run exe
# ---------------------------------------------------------------

all: intro setdirs run


intro:
	@echo ''
	@echo '===> Making Olliptic ...'	
	@echo ''	

setdirs:
	@mkdir -p $(EXED) $(OBJD) 

run: $(OBJS) $(HDRS)
	@echo ''
	@echo '===> Building the executable'	
	@echo ''	
	@cd $(OBJD); $(CPP) $(CPPFLAGS) $(LDFLAGS) -o $(EXED)/$(EXE).x $(OBJS)
	@echo ''
	@echo '============ All Done! ============'
	@echo ' '

# ---------------------------------------------------------------
# Target 2 : debug
# ---------------------------------------------------------------

debug: intro clean $(OBJS) $(HDRS)
	@echo ''
	@echo '===> Building the executable'	
	@echo ''
	@cd $(OBJD); $(CPP) $(CPPFLAGS) $(LDFLAGS) -o $(EXED)/$(EXE).x  $(OBJS)      
	@echo ''
	@echo '============ All Done! ============'
	@echo ''
	@echo 'DEBUG:'
	@echo ' gdb [exe] [corefile]'
	@echo ''
	@echo 'Remember the settings for the core output !!!'
	@echo '  (sh/bash)    ulimit -c unlimited'
	@echo '  (tcsh/csh)   limit coredumpsize unlimited'
	@echo ' '
	@echo 'PROFILE:'
	@echo ' time [exe]'
	@echo ' gprof -b [exe] gmon.out'
	@echo ' '

# ---------------------------------------------------------------
# Objects
# ---------------------------------------------------------------

%.o $(OBJD)/%.o: %.cc
	@echo ''
	@echo '===> Building' $*.o
	@echo ''
	@cd $(OBJD); $(CPP) $(CPPFLAGS) -c $(MPIFLAG) -I $(HDRSD) $(CCFLAGS) $(SRCD)/$(*).cc

# ---------------------------------------------------------------
# Target 3 : clean tmp files
# ---------------------------------------------------------------

clean: intro
	@rm -f $(BASE)/obj/*.o
	@rm -f $(BASE)/*~
	@rm -f $(SRCD)/*~
	@rm -f $(HDRSD)/*~
	@rm -f $(BASE)/*.dat
	@rm -f $(BASE)/*.xg
	@echo " === clean done ! === "

# ---------------------------------------------------------------
# Target 4 : delete
# ---------------------------------------------------------------

delete: intro
	@rm -f $(OBJD)/*.o
	@rm -f $(SRCD)/*~
	@rm -f $(HDRSD)/*~
	@rm -f $(BASE)/*~
	@rm -f $(BASE)/*\#
	@rm -f $(BASE)/*.x
	@rm -f $(BASE)/*.out
	@rm -f $(BASE)/*.dat
	@rm -f $(BASE)/*.xg
	@echo " === delete done ! === "

# ---------------------------------------------------------------
# target 5 : help
# ---------------------------------------------------------------

help: 
	@echo ' '
	@echo ' ---------------------------'
	@echo ' Makefile - HELP'
	@echo ' ---------------------------'
	@echo '  '
	@echo ' make        : compile with debug options'
	@echo ' '
	@echo ' make run    : compile with optimizations'
	@echo ' '
	@echo ' make clean  : clean object, ~ and data files'
	@echo ' '
	@echo ' make delete : delete object, ~, data and exe file'
	@echo ' '

