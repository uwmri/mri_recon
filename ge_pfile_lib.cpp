#include "ge_pfile_lib.hpp"

using namespace NDarray;

void PFILE::read_header(char Ifilename[]){
  
  strcpy(filename,Ifilename);

  FILE *fphead;
  char header[RDB_HEADER_SIZE_BYTES];
  
  printf("Read P-file:: Byte size %d\n",(int)RDB_HEADER_SIZE_BYTES);
  fflush(stdout);
  if ((fphead = fopen64(filename,"r")) != NULL) {
    fseeko64(fphead, (long long) 0, (int) 0);
    if( fread(header,sizeof(char),RDB_HEADER_SIZE_BYTES, fphead) !=RDB_HEADER_SIZE_BYTES){
    	throw std::runtime_error("Couldn't read Pfile header");
	}
    fclose(fphead);
  } else {
    throw std::runtime_error("Failed to open P-file");
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
}


void PFILE::read_data( complex<float> *dat,int coil, int offset, int number, int stride){
	read_data(dat,coil,offset,number, stride,0);
}


void PFILE::read_data( complex<float> *dat,int coil, int offset, int number, int stride, int start){
  
 // cout << "Read" << coil << "," << offset << "," << number << "," << stride << endl;
  int i, j;
  short int *buf=0;
  int *buf_int=0;
  int total_frames;
  long long skip;
  FILE *fd;
  
  int xres =  rdbhead.rdb_hdr_da_xres - start;
   
  /**Seek to coil*/
  total_frames = rdbhead.rdb_hdr_da_yres* rdbhead.rdb_hdr_nslices;
  if (!(fd=fopen64(filename,"r"))) {
    fprintf(stderr,"Can't open raw data file %s for read!\n",filename);
    exit(1);
  }
  
  if(rdbhead.rdb_hdr_point_size==2){
  	buf = (short int *)malloc(rdbhead.rdb_hdr_da_xres*2*sizeof(short));
  }else{
  	buf_int = (int *)malloc(rdbhead.rdb_hdr_da_xres*2*sizeof(int));
  }
  
  /*Seek to Projection*/
  long long seek_proj = (long long)offset;
  long long seek_base = (long long)(ceil( (float)seek_proj / (float)( rdbhead.rdb_hdr_da_yres-1)) ) ; 
  
  /*Skip Header*/
  fseeko64(fd, (long long)RDB_HEADER_SIZE_BYTES, 1);
  
  /*Skip to Coil*/
  skip = 2*rdbhead.rdb_hdr_point_size*(long long)rdbhead.rdb_hdr_da_xres*(long long)total_frames*(long long)coil;
  fseeko64(fd, skip, 1);
  
  /*Skip to Encoding*/
  skip =  2*(long long)rdbhead.rdb_hdr_point_size*(long long)rdbhead.rdb_hdr_da_xres * (seek_proj + seek_base); 
  fseeko64(fd, skip, 1);
  
  /*Now Read All the data*/
  for (i= 0; i< number; i++) {
	 
	 if ((seek_proj%( rdbhead.rdb_hdr_da_yres-1) ) == 0){
      skip = 2*(long long)rdbhead.rdb_hdr_point_size*(long long)rdbhead.rdb_hdr_da_xres; 
      fseeko64(fd, skip, 1);
	 }
	 	 
	 if(rdbhead.rdb_hdr_point_size==2){
	 	if ((j=fread(buf,2*sizeof(short),rdbhead.rdb_hdr_da_xres,fd))!=rdbhead.rdb_hdr_da_xres){
       		fprintf(stderr,"Not enough data, projection %d has %d samples\n",i, j);
       		exit(1);
     	}
	 	
	 	for (j=start; j < rdbhead.rdb_hdr_da_xres; j++){
			dat[i*xres + (j-start)] += complex<float>( (float)buf[j*2],(float)buf[j*2+1]);
	 	}
	 }else{
	 	if ((j=fread(buf_int,2*sizeof(int),rdbhead.rdb_hdr_da_xres,fd))!=rdbhead.rdb_hdr_da_xres){
       		fprintf(stderr,"Not enough data, projection %d has %d samples\n",i, j);
       		exit(1);
     	}
	 
	 	for (j=start; j < rdbhead.rdb_hdr_da_xres; j++){
			dat[i*xres + (j-start)] += complex<float>( (float)buf_int[j*2],(float)buf_int[j*2+1]);
      	}
	 }
	 seek_proj++;	 
	 
	 /*Account for Stride*/
	 skip  = (stride-1)*2*(long long)rdbhead.rdb_hdr_point_size*(long long)rdbhead.rdb_hdr_da_xres; 
	 skip += (long long)( (ceil( (float)seek_proj / (float) ( rdbhead.rdb_hdr_da_yres-1)))- (ceil( (float)(seek_proj+stride-1) / (float)( rdbhead.rdb_hdr_da_yres-1))) );
	 seek_proj+=(stride-1);
	 fseeko64(fd, skip, 1);
	 
  }
  
  /* Cleanup*/
  fclose(fd);
  if(rdbhead.rdb_hdr_point_size==2){
    free(buf);
  }else{
    free(buf_int);
  }   
 return;

}

void PFILE::write_dicom( NDarray::Array< complex<float>,3> &X, int lx,const string series_description, float max_all,float zoom_x,float zoom_y, float zoom_z){

		
	printf("Initialize Dicom Header\n");  
	DATA_INFO DI;

	/*TEMP*/
	// strcpy(serieshead.prtcl,"");

	/* we will zero image to get it square if x != y */
	imagehead.imatrix_X=  X.length(firstDim);
	imagehead.imatrix_Y = X.length(secondDim);
	imagehead.vasflags &= ~2;
	imagehead.slthick = zoom_z * imagehead.dfov/ (float)X.length(thirdDim);
	imagehead.dfov *= zoom_x;
	
	cout << "Slthick is " << imagehead.slthick << endl;

	/*Ge Dicom Routine*/
	initGEImage( &DI, &examhead, &serieshead, &imagehead, &rdbhead );
	DI.image_xres = X.length(firstDim);
	DI.image_yres = X.length(secondDim);
	DI.numslices =  X.length(thirdDim);
	DI.series_numimages = X.length(thirdDim);
	DI.slthick = imagehead.slthick;

	imagehead.slquant = X.length(thirdDim);
	serieshead.se_numimages = X.length(thirdDim);
	
	float *temp = new float[ X.length(firstDim)*X.length(secondDim)];
	float scale = 4095/max_all;
	short *signa_image;
	
	char *a=new char[series_description.size()+1];
	a[series_description.size()]=0;
	memcpy(a,series_description.c_str(),series_description.size());
	
	sprintf(DI.series_description,"%s",a );

	imagehead.im_seno    =lx;
	DI.new_series_number = lx;
	DI.base_series_number= lx;	

	for (int k = 0; k < X.length(thirdDim); k++){
			DI.image_number = k+1;
			strcpy(DI.outfilename, "");

			int pos = 0;
			for(int j=0; j< X.length(secondDim); j++){
				for(int i=0; i< X.length(firstDim); i++){
					temp[pos] = scale*abs( X(i,j,k));
					pos++;
			}}

			if (procGEImage(temp, &signa_image, &DI, &rdbhead,&examhead,&serieshead,&imagehead) != 0) {
				printf("Error processing magnitude image data.\n");
				exit(1);
			}

			if(writeDCM(DI.outfilename, rdbhead, examhead, serieshead, imagehead, signa_image, 0, 0, 0) != 0) {
				printf("Error writing magnitude image file.\n");
				exit (1);
			}
	}
	delete [] temp;
	
	imagehead.dfov /= zoom_x;


}

void PFILE::write_dicom( NDarray::Array< float,3> &X, int lx,const string series_description, float max_all,float zoom_x,float zoom_y, float zoom_z){

		
	printf("Initialize Dicom Header\n");  
	DATA_INFO DI;

	/*TEMP*/
	// strcpy(serieshead.prtcl,"");

	/* we will zero image to get it square if x != y */
	imagehead.imatrix_X=  X.length(firstDim);
	imagehead.imatrix_Y = X.length(secondDim);
	imagehead.vasflags &= ~2;
	imagehead.slthick = imagehead.dfov/ (float)X.length(thirdDim);
	cout << "Slthick is " << imagehead.slthick << endl;
   

	/*Ge Dicom Routine*/
	initGEImage( &DI, &examhead, &serieshead, &imagehead, &rdbhead );
	
	// Size of Image
	DI.image_xres = X.length(firstDim);
	DI.image_yres = X.length(secondDim);
	DI.numslices =  X.length(thirdDim);
	
	DI.series_numimages = X.length(thirdDim);
	DI.slthick = imagehead.slthick;

	imagehead.slquant = X.length(thirdDim);
	serieshead.se_numimages = X.length(thirdDim);
	
	float *temp = new float[ X.length(firstDim)*X.length(secondDim)];
	float scale = 4095/max_all;
	short *signa_image;
	
	char *a=new char[series_description.size()+1];
	a[series_description.size()]=0;
	memcpy(a,series_description.c_str(),series_description.size());
	
	sprintf(DI.series_description,"%s",a );

	imagehead.im_seno    =lx;
	DI.new_series_number = lx;
	DI.base_series_number= lx;	

	for (int k = 0; k < X.length(thirdDim); k++){
			DI.image_number = k+1;
			strcpy(DI.outfilename, "");

			int pos = 0;
			for(int j=0; j< X.length(secondDim); j++){
				for(int i=0; i< X.length(firstDim); i++){
					temp[pos] = scale*abs( X(i,j,k));
					pos++;
			}}

			if (procGEImage(temp, &signa_image, &DI, &rdbhead,&examhead,&serieshead,&imagehead) != 0) {
				printf("Error processing magnitude image data.\n");
				exit(1);
			}

			if(writeDCM(DI.outfilename, rdbhead, examhead, serieshead, imagehead, signa_image, 0, 0, 0) != 0) {
				printf("Error writing magnitude image file.\n");
				exit (1);
			}
	}
	delete [] temp;

}

	

/*********************************************
get_unit_vector_physical

Input: 
	Imagehead		image header from scan 
	scn_tbl 		table storing positions
	
Output:NA
		
Description:
	Gets the physical coordinate system to image translation matrix.
	This routine just copies info from image header. It mostly
	here to make things easier to understand. IE what corner to 
	start from, direction,etc. This is unlikely to work for planes
	that haven't been sorted out. 
	
	System was designed for 2D so the vectors are slightly confusing
*********************************************/

void PFILE::calc_unit_vectors(int rcxres,int rcyres,int rczres,float zoom_x,float zoom_y,float zoom_z){
	
	float delx, dely, delz;
	float actual_fovx,actual_fovy,actual_fovz;
	float prescribed_fovx,prescribed_fovy,prescribed_fovz;
		
	
	/*Get FOV (temp)*/
	prescribed_fovx = imagehead.dfov;
	prescribed_fovy = imagehead.dfov;
	prescribed_fovz = imagehead.dfov;
	 
	actual_fovx = prescribed_fovx * zoom_x; 
	actual_fovy = prescribed_fovy * zoom_y; 
	actual_fovz = prescribed_fovz * zoom_z; 
		
	delx= actual_fovx / (float)rcxres;
	dely= actual_fovy / (float)rcyres;
	delz= actual_fovz / (float)rczres;
	
	/*Corner Position*/
	physical_tbl.sx =  ( acq_tab[0].gw_point2[0]); /*Start Location No Shifts or Zoom*/
	physical_tbl.sy =  ( acq_tab[0].gw_point2[1]); /*Start Location No Shifts or Zoom*/
	physical_tbl.sz =  ( acq_tab[0].gw_point2[2]); /*Start Location No Shifts or Zoom*/
	// printf("CORNER = %f %f %f\n",physical_tbl.sx,physical_tbl.sy,physical_tbl.sz);
	
	/*Get Unit Vectors*/
	physical_tbl.ix =  ( acq_tab[0].gw_point1[0] - acq_tab[0].gw_point2[0] ) / imagehead.dfov * delx;
	physical_tbl.iy =  ( acq_tab[0].gw_point1[1] - acq_tab[0].gw_point2[1] ) / imagehead.dfov * delx;
	physical_tbl.iz =  ( acq_tab[0].gw_point1[2] - acq_tab[0].gw_point2[2] ) / imagehead.dfov * delx;
		
	physical_tbl.jx =  ( acq_tab[0].gw_point3[0] - acq_tab[0].gw_point1[0] ) / imagehead.dfov *dely;
	physical_tbl.jy =  ( acq_tab[0].gw_point3[1] - acq_tab[0].gw_point1[1] ) / imagehead.dfov *dely;
	physical_tbl.jz =  ( acq_tab[0].gw_point3[2] - acq_tab[0].gw_point1[2] ) / imagehead.dfov *dely;
	
	physical_tbl.kx =  -( physical_tbl.iy*physical_tbl.jz - physical_tbl.iz*physical_tbl.jy)/ dely/ delx * delz; 
	physical_tbl.ky =  -( physical_tbl.iz*physical_tbl.jx - physical_tbl.ix*physical_tbl.jz)/ dely/ delx * delz; 
	physical_tbl.kz =  -( physical_tbl.ix*physical_tbl.jy - physical_tbl.iy*physical_tbl.jx)/ dely/ delx * delz; 
	
	/*Need to fix offsets for Zooming */	
	physical_tbl.sx -=  physical_tbl.ix/delx*(actual_fovx - prescribed_fovx)/2.0;
	physical_tbl.sy -=  physical_tbl.iy/delx*(actual_fovx - prescribed_fovx)/2.0;
	physical_tbl.sz -=  physical_tbl.iz/delx*(actual_fovx - prescribed_fovx)/2.0;

	physical_tbl.sx -=  physical_tbl.jx/dely*(actual_fovy - prescribed_fovy)/2.0;
	physical_tbl.sy -=  physical_tbl.jy/dely*(actual_fovy - prescribed_fovy)/2.0;
	physical_tbl.sz -=  physical_tbl.jz/dely*(actual_fovy - prescribed_fovy)/2.0;
	
	physical_tbl.sx -=  physical_tbl.kx/delz*(actual_fovz - prescribed_fovz)/2.0;
	physical_tbl.sy -=  physical_tbl.ky/delz*(actual_fovz - prescribed_fovz)/2.0;
	physical_tbl.sz -=  physical_tbl.kz/delz*(actual_fovz - prescribed_fovz)/2.0;
	
	/****Make Sure SZ is Zero (Assumes Table Was Shifted)**/
	physical_tbl.sz = -0.5* ( (float)rczres*physical_tbl.kz + (float)rcyres*physical_tbl.jz  + (float)rcxres*physical_tbl.iz );


		/*Now Convert to Physical Units --- 
	 Currently using the rotation:
	 	[-1  0  0;
		  0 -1  0
		  0  0  1]  *Head First*
	
	 	[ 1  0  0;
		  0 -1  0
		  0  0 -1]  *Feet First*
		  
	For physical to logical. These are based on empirical derivation.
	*/	
	
	
	if(rdbhead.rdb_hdr_entry==1){
		/*Head First*/
		logical_tbl.ix *= -1; 
		logical_tbl.iy *= -1;
		logical_tbl.iz *=  1;
		
		logical_tbl.jx *= -1; 
		logical_tbl.jy *= -1;
		logical_tbl.jz *=  1;
		
		logical_tbl.kx *= -1; 
		logical_tbl.ky *= -1;
		logical_tbl.kz *=  1;
				
		logical_tbl.sx *= -1; 
		logical_tbl.sy *= -1;
		logical_tbl.sz *=  1;
		
		logical_tbl.sz -=  rdbhead.rdb_hdr_scancent;
	}else{
		/*Feet First*/
		logical_tbl.ix *=  1; 
		logical_tbl.iy *= -1;
		logical_tbl.iz *= -1;
		
		logical_tbl.jx *=  1; 
		logical_tbl.jy *= -1;
		logical_tbl.jz *= -1;
		
		logical_tbl.kx *=  1; 
		logical_tbl.ky *= -1;
		logical_tbl.kz *= -1;
				
		logical_tbl.sx *=  1; 
		logical_tbl.sy *= -1;
		logical_tbl.sz *= -1;
		
		logical_tbl.sz +=  rdbhead.rdb_hdr_scancent;
	}
			
	/*Read in Fov Info */ 
	printf("Center %f %f %f \n", imagehead.ctr_R,  imagehead.ctr_A,  imagehead.ctr_S);
	printf("Normal %f %f %f \n", imagehead.norm_R, imagehead.norm_A, imagehead.norm_S);
	printf("Top L  %f %f %f \n", imagehead.tlhc_R, imagehead.tlhc_A, imagehead.tlhc_S);
	printf("Top R  %f %f %f \n", imagehead.trhc_R, imagehead.trhc_A, imagehead.trhc_S);
	printf("Bot R  %f %f %f \n", imagehead.brhc_R, imagehead.brhc_A, imagehead.brhc_S);	
	printf("Pass Number %d\n",acq_tab[0].pass_number);
	printf("Swap PF %d\n",imagehead.swappf);
	printf("Patient Position %d\n",serieshead.position);
	printf("Slice Thickness %f\n",imagehead.slthick);
	printf("Scan Center %f\n",rdbhead.rdb_hdr_scancent); 
	printf("Patient Position %d\n",rdbhead.rdb_hdr_position);
	printf("Patient Entry %d\n",rdbhead.rdb_hdr_entry);
	printf("Landmark %f\n",rdbhead.rdb_hdr_lmhor);
	
	/***INFO FROM ACQ ***/
	printf("GW_PT1 %f %f %f\n",acq_tab[0].gw_point1[0],acq_tab[0].gw_point1[1],acq_tab[0].gw_point1[2]);
	printf("GW_PT2 %f %f %f\n",acq_tab[0].gw_point2[0],acq_tab[0].gw_point2[1],acq_tab[0].gw_point2[2]);
	printf("GW_PT3 %f %f %f\n",acq_tab[0].gw_point3[0],acq_tab[0].gw_point3[1],acq_tab[0].gw_point3[2]);
	printf("Pass Number %d\n",acq_tab[0].pass_number);
	printf("Transpose %d\n",acq_tab[0].transpose);
	printf("Rotate %d\n",acq_tab[0].rotate);
	
	printf("Physical Matrix=\n");
	printf(" | %5f %5f %5f |[i] %5f\n",physical_tbl.ix,physical_tbl.jx,physical_tbl.kx,physical_tbl.sx);
	printf(" | %5f %5f %5f |[j]+%5f\n",physical_tbl.iy,physical_tbl.jy,physical_tbl.ky,physical_tbl.sy);
	printf(" | %5f %5f %5f |[k] %5f\n",physical_tbl.iz,physical_tbl.jz,physical_tbl.kz,physical_tbl.sz);
	printf("\n");
	
	

}


