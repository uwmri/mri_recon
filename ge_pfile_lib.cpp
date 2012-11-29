/************************************************
Libraries for GE P-files

Initial Author: Kevin M. Johnson
Description: This code contains files to read p-file data

*************************************************/
#include "ge_pfile_lib.h"

void PFILE::read_header(char Ifilename[]){
  
  strcpy(filename,Ifilename);

  FILE *fphead;
  char header[RDB_HEADER_SIZE_BYTES];
  
  printf("Read P-file:: Byte size %d\n",(int)RDB_HEADER_SIZE_BYTES);
  fflush(stdout);
  if ((fphead = fopen64(filename,"r")) != NULL) {
    fseeko64(fphead, (long long) 0, (int) 0);
    if( fread(header,sizeof(char),RDB_HEADER_SIZE_BYTES, fphead) !=RDB_HEADER_SIZE_BYTES){
    	printf("Failed to read in raw P-file header.\n"); 
    	exit(1);
    }
    fclose(fphead);
  } else {
    printf("Failed to open in raw P-file header.\n"); 
    exit(1);
  }
  
  // printf("Copy Header\n");
  fflush(stdout);
  memcpy(&rdbhead, header+RDB_HDR_OFF, RDB_HDR_SIZE);
  memcpy(&acq_tab, header+RDB_DATAACQ_OFF, RDB_DATAACQ_SIZE);
  memcpy(&examhead, header+RDB_EXAMDATATYPE_OFF, RDB_EXAMDATATYPE_SIZE);
  memcpy(&serieshead, header+RDB_SERIESDATATYPE_OFF, RDB_SERIESDATATYPE_SIZE);
  memcpy(&imagehead, header+RDB_MRIMAGEDATATYPE_OFF, RDB_MRIMAGEDATATYPE_SIZE);

  // Raw Data 
  cout << "Read P-file:: Slices " << rdbhead.rdb_hdr_nslices << " Views = " << rdbhead.rdb_hdr_da_yres - 1 << " Xres= " <<rdbhead.rdb_hdr_da_xres <<endl;
  RawData.alloc(rdbhead.rdb_hdr_nslices,rdbhead.rdb_hdr_da_yres - 1,rdbhead.rdb_hdr_da_xres);
}

//  Blindly reads all data for one coil - No Sorting for Encodig,echo,etc
void PFILE::read_data(int coil){
	 
  // Open File
  long long total_frames = rdbhead.rdb_hdr_da_yres*rdbhead.rdb_hdr_nslices;
  FILE *fd=NULL;
  if (!(fd=fopen64(filename,"r"))) {
    fprintf(stderr,"Can't open raw data file %s for read!\n",filename);
    exit(1);
  }
  
  // Allocate Memory for Xres
  short int *buf=NULL;
  int *buf_int=NULL;
  if(rdbhead.rdb_hdr_point_size==2){
  	buf = new short int[2*rdbhead.rdb_hdr_da_xres];
  }else{
  	buf_int = new int[2*rdbhead.rdb_hdr_da_xres];
  }
  
  // Skip to Coil
  fseeko64(fd, (long long)RDB_HEADER_SIZE_BYTES, SEEK_SET);
  long long skip = 2*(long long)rdbhead.rdb_hdr_point_size*(long long)rdbhead.rdb_hdr_da_xres*total_frames*(long long)coil;
    
    	
  // Read in all the data			
  int pts;
  for(int k= 0; k< rdbhead.rdb_hdr_nslices; k++) {
	 // Skip baselines 
	 skip = 2*(long long)rdbhead.rdb_hdr_point_size*(long long)rdbhead.rdb_hdr_da_xres;
     fseeko64(fd, skip, SEEK_CUR);
	 
	 for(int j=0; j< rdbhead.rdb_hdr_da_yres - 1; j++){
	 	 // Now read all slices
	 	if(rdbhead.rdb_hdr_point_size==2){
	 		if ((pts=fread(buf,2*sizeof(short),rdbhead.rdb_hdr_da_xres,fd))!=rdbhead.rdb_hdr_da_xres){
       			cout << "Slice " << k << " view = " << j << " only " << pts << " points " << endl;
       			exit(1);
     		}
			
			for(int i=0; i< rdbhead.rdb_hdr_da_xres;i++){
				RawData[k][j][i]=complex<int>( (int)buf[2*i],(int)buf[2*i+1]); 
			}
	 	}else{
	 		if ((pts=fread(buf_int,2*sizeof(int),rdbhead.rdb_hdr_da_xres,fd))!=rdbhead.rdb_hdr_da_xres){
       			cout << "Slice " << k << " view = " << j << " only " << pts << " points " << endl;
       			exit(1);
     		}
			
			for(int i=0; i< rdbhead.rdb_hdr_da_xres;i++){
				RawData[k][j][i]=complex<int>( buf_int[2*i], buf_int[2*i+1]); 
			}
	 	}
	 }// Yres
  }// Slice
  fclose(fd);
}

