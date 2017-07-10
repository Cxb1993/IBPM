#Makefile for ibfs on mac

FFILES = parameters grid myfft user variables operators operators_pressure main
MFILES = getdata makeplot

IBFILES = ${FFILES:=.o}

COMPILE = gfortran-mp-6
FLAGS1 = -m64 -I/usr/local/fftw3/include -C
FLAGS2 = -L/usr/local/fftw3/lib -lfftw3 -lblas -llapack -o

BASE = /Users/andresgoza/Documents/CFPG/IB_2d_code/PC_code/gilead/test2
SRC = src
OBJ = obj
BIN = bin
OUT = output
IN  = input

DIR = ${SRC} ${OBJ} ${BIN} ${OUT} ${IN}
DIRS = ${DIR:=.dir}

IBFILES_SRC := $(foreach p, $(IBFILES), $(SRC)/$(p))
IBFILES_OBJ := $(foreach p, $(IBFILES), $(OBJ)/$(p))

ib :  ${IBFILES_OBJ}
	$(COMPILE) $(FLAGS1) $(IBFILES_OBJ) $(FLAGS2) $(BIN)/ib -J$(OBJ) 

setup: ${DIRS}
	@for i in ${DIR} ; do \
               mkdir -p $${i} ;\
	done
	@for i in ${FFILES} ; do \
               echo "cp  -p ${BASE}/${SRC}/$${i}.f90 ${SRC}/$${i}.f90 " ;\
               cp  -p ${BASE}/${SRC}/$${i}.f90 ${SRC}/$${i}.f90 ;\
	done
	@for i in ${MFILES} ; do \
               echo "cp  -p ${BASE}/${OUT}/$${i}.m ${OUT}/$${i}.m " ;\
               cp  -p ${BASE}/${OUT}/$${i}.m ${OUT}/$${i}.m ;\
	done

realclean: ${DIRS} # 
	@for i in ${DIR} ; do \
               rm -rf $${i} ;\
	done

clean: 
	rm -f $(BIN)/ib ${OBJ}/*.o ${OBJ}/*.mod *.mod

$(OBJ)/%.o: $(SRC)/%.f90
	$(COMPILE) $(FLAGS1) -o $@ -c $<

%.dir : 	
	$< 
