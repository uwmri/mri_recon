STATIC=0
DICOM =0
CC = g++

ifeq ($(STATIC),1)
	SLINK = -static
else
	SLINK = 
endif

LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
  	CFLAGS_LOCAL  = -c -g -DLINUX -m64 -DRECON_64BIT -Wall -D_FILE_OFFSET_BITS=64 -fopenmp -O3  
else
	CFLAGS_LOCAL  = -c -g -DLINUX -Wall 
endif

CFLAGS = $(CFLAGS_LOCAL) $(CFLAGS_DCMTK) -Wno-deprecated 

INC_DIRS = 

ifeq ($(LBITS),64)
  LIB_DIRS = -L$(DCMTK_DIR)/lib64 -L/usr/lib64 
else
  	ifeq ($(GE),1)
  		LIB_DIRS = -L$(DCMTK_DIR)/lib -L/usr/local/olinux32/lib 
	else
		LIB_DIRS = -L$(DCMTK_DIR)/lib -L/usr/lib 
	endif
endif

FFTW3_LIBS = -lfftw3f_omp -lfftw3f 
LOCAL_LIBS =  $(FFTW3_LIBS) -lpthread -lm -fopenmp 
STATIC_LIBS = 
LIBS = $(DICOM_LIBS) $(LOCAL_LIBS)

# Armadillo/ACML libraries (coil compression and spirit stuff)
MRFLOWHOME=/export/home/mrflow
LIB_DIRS += -L$(MRFLOWHOME)/linux/arma322/usr/lib64/
INC_DIRS += -I$(MRFLOWHOME)/linux/arma322/usr/include/
LIB_DIRS += -L$(MRFLOWHOME)/linux/acml440/gfortran64/lib/ -L$(MRFLOWHOME)/linux/acml440/gfortran64_mp/lib/
INC_DIRS += -I$(MRFLOWHOME)/linux/acml440/gfortran64/include/ -I$(MRFLOWHOME)/linux/acml440/gfortran64_mp/include/
LIBS += -lacml_mp -lacml_mv

# Blitz++
LIB_DIRS += -L$ /export/home/kmjohnso/linux/lib/
INC_DIRS += -I$ /export/home/kmjohnso/linux/include/

RUNPATH=$(MRFLOWHOME)/linux/arma322/usr/lib64/:$(MRFLOWHOME)/linux/acml440/gfortran64/lib/:$(MRFLOWHOME)/linux/acml440/gfortran64_mp/lib/

# For Main COmpiles
RECON_OBJECTS =  mri_data.o gridFFT.o spirit.o temporal_diff.o wavelet3D.o recon.o recon_lib.o softthreshold.o #ge_pfile_lib.o

# ditto
all:  recon_binary

# recon-library-links compiles
recon_binary: $(RECON_OBJECTS)
	$(CC) $(STATIC_LIBS) -o recon_binary $(RECON_OBJECTS) $(LIB_DIRS) -Wl,-rpath,$(RUNPATH) $(LIBS) $(SLINK) 

# recon compiles
recon.o: recon.cxx
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) recon.cxx 
wavelet3D.o: wavelet3D.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) wavelet3D.cpp
spirit.o: spirit.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) spirit.cpp
temporal_diff.o: temporal_diff.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) temporal_diff.cpp
recon_lib.o: recon_lib.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) recon_lib.cpp
mri_data.o: mri_data.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) mri_data.cpp
softthreshold.o: softthreshold.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) softthreshold.cpp
gridFFT.o: gridFFT.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) gridFFT.cpp
ge_pfile_lib.o: ge_pfile_lib.cpp
	$(CC) $(INC_DIRS) $(CFLAGS) $(SLINK) ge_pfile_lib.cpp

clean:
	rm -f $(RECON_OBJECTS) $(recon) 
	
clean_o:
	rm *.o -f  



