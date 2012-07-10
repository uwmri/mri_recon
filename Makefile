VER = 22
GE = 0
STATIC=0
DICOM =0

ifeq ($(GE),1)
	CC = g++32
else
	CC = g++
endif


ifeq ($(STATIC),1)
	SLINK = -static
else
	SLINK = 
endif

ifeq ($(GE),1)
	DCMTK_DIR = /usr/local/olinux32/dcmtk
else
	DCMTK_DIR = /usr/local/dcmtk
endif
	
#include $(DCMTK_DIR)/Makefile.def
#include /export/home/kvigen/linux/src/dicom/Makefile.def


LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
  CFLAGS_LOCAL  = -c -g -DLINUX -m64 -DRECON_64BIT -Wall -D_FILE_OFFSET_BITS=64 -fopenmp -O2
else
	ifeq ($(GE),1)
   		CFLAGS_LOCAL  = -c -g -DLINUX -ansi-pedantic -Wall -I$(DCMTK_DIR)/include -I/usr/local/olinux32/include -L$(DCMTK_DIR)/lib -L/usr/local/olinux32/lib
	else
		CFLAGS_LOCAL  = -c -g -DLINUX -Wall 
	endif
endif

ifeq ($(DICOM),1)
CFLAGS_DCMTK  = -D_REENTRANT -D_XOPEN_SOURCE_EXTENDED -D_XOPEN_SOURCE=500            \
                -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_POSIX_C_SOURCE=199506L  \
                -DHAVE_CONFIG_H -D UW_DICOM
else
CFLAGS_DCMTK = -D_REENTRANT -D_XOPEN_SOURCE_EXTENDED -D_XOPEN_SOURCE=500            \
                -D_BSD_SOURCE -D_BSD_COMPAT -D_OSF_SOURCE -D_POSIX_C_SOURCE=199506L  \
                -DHAVE_CONFIG_H 
endif


CFLAGS = $(CFLAGS_LOCAL) $(CFLAGS_DCMTK) -Wno-deprecated 

INC_DIRS = -I$(DCMTK_DIR)/include -IGEHEADER 

ifeq ($(LBITS),64)
  LIB_DIRS = -L$(DCMTK_DIR)/lib64 -L/usr/lib64 
else
  	ifeq ($(GE),1)
  		LIB_DIRS = -L$(DCMTK_DIR)/lib -L/usr/local/olinux32/lib 
	else
		LIB_DIRS = -L$(DCMTK_DIR)/lib -L/usr/lib 
	endif
endif

FFTW3_LIBS = -lfftw3f_threads -lfftw3f 
LOCAL_LIBS =  $(FFTW3_LIBS) -lpthread -lm -lz -fopenmp 
STATIC_LIBS = 

ifeq ($(DICOM),1)
	DICOM_LIBS = -luw_lnx$(VER) -ldcmdata -lofstd
else
	DICOM_LIBS = 
endif

RECON_VERSION = ESE$(VER)_RECON


LIBS = $(DICOM_LIBS) $(LOCAL_LIBS)

# ML libraries
LOECHHOME=/export/home/loecher
LIB_DIRS += -L$(LOECHHOME)/linux/arma322/usr/lib64/
INC_DIRS += -I$(LOECHHOME)/linux/arma322/usr/include/
LIB_DIRS += -L$(LOECHHOME)/linux/acml440/gfortran64/lib/ -L$(LOECHHOME)/linux/acml440/gfortran64_mp/lib/
INC_DIRS += -I$(LOECHHOME)/linux/acml440/gfortran64/include/ -I$(LOECHHOME)/linux/acml440/gfortran64_mp/include/
LIBS += -lgfortran -lacml_mp

# For Main COmpiles
#RECON_OBJECTS = gridFFT.o wavelet3D.o polynomial_fitting.o gating_lib.o tornado_lib.o trajectory_lib.o io_lib.o matrix_lib.o master_lib.o iterative_lib.o csi_lib.o pcvipr_gradwarp.o master_recon.o 
RECON_OBJECTS =  gridFFT.o wavelet3D.o recon.o recon_lib.o softthreshold.o ge_pfile_lib.o spirit.o


# ditto
recon = recon$(VER)

all:  $(recon)

# recon-library-links compiles
recon$(VER): $(RECON_OBJECTS)
	$(CC) $(STATIC_LIBS) -o recon$(VER) $(RECON_OBJECTS) $(LIB_DIRS) $(LIBS) -D $(RECON_VERSION) $(SLINK) 

# recon compiles
recon.o: recon.cxx
	$(CC) $(INC_DIRS) $(CFLAGS)  -D $(RECON_VERSION) $(SLINK) recon.cxx 


# library compile
master_lib.o: master_lib.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) master_lib.c
pcvipr_gradwarp.o: pcvipr_gradwarp.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) pcvipr_gradwarp.c
matrix_lib.o: matrix_lib.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) matrix_lib.c
csi_lib.o: csi_lib.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) csi_lib.c
iterative_lib.o: iterative_lib.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) iterative_lib.c
io_lib.o: io_lib.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) io_lib.c
gating_lib.o: gating_lib.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) gating_lib.cpp
wavelet3D.o: wavelet3D.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) wavelet3D.cpp
recon_lib.o: recon_lib.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) recon_lib.cpp
softthreshold.o: softthreshold.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) softthreshold.cpp
gridFFT.o: gridFFT.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) gridFFT.cpp
ge_pfile_lib.o: ge_pfile_lib.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) ge_pfile_lib.cpp
tornado_lib.o: tornado_lib.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) tornado_lib.c
trajectory_lib.o: trajectory_lib.c
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) trajectory_lib.c
polynomial_fitting.o: polynomial_fitting.c
	$(CC) $(DEB_HOST_GNU_CPU) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) polynomial_fitting.c
spirit.o: spirit.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) -D $(RECON_VERSION) $(SLINK) spirit.cpp


clean:
	rm -f $(RECON_OBJECTS) $(recon) $(filter) $(FILTER_OBJECTS) $(dixon) $(DIXON_OBJECTS)
	
clean_o:
	rm *.o -f  



