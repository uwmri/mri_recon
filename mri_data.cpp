
void MRI_DATA::read_kmap( ){

	FILE *fid;
	float *map_exportx = (float *)malloc(sizeof(float)*flag->xres*flag->nproj);		
	float *map_exporty = (float *)malloc(sizeof(float)*flag->xres*flag->nproj);		
	float *map_exportz = (float *)malloc(sizeof(float)*flag->xres*flag->nproj);		
	float *map_exportw = (float *)malloc(sizeof(float)*flag->xres*flag->nproj);		
			
	char fname[80];
	sprintf(fname,"%sKMAPX_VD_%d.dat",flag->data_folder,vdt);
	fid=fopen(fname,"r");
	fread(map_exportx,sizeof(float),flag->xres*flag->nproj,fid);
	fclose(fid);
	
	sprintf(fname,"%sKMAPY_VD_%d.dat",flag->data_folder,vdt);
	fid=fopen(fname,"r");
	fread(map_exporty,sizeof(float),flag->xres*flag->nproj,fid);
	fclose(fid);
		
	if(flag->ss_2d==1){
 	memset(map_exportz,0,(size_t)( 	flag->xres*flag->nproj*sizeof(float)));
	}else{
	sprintf(fname,"%sKMAPZ_VD_%d.dat",flag->data_folder,vdt);
	fid=fopen(fname,"r");
	fread(map_exportz,sizeof(float),flag->xres*flag->nproj,fid);
	fclose(fid);
	}
	
	sprintf(fname,"%sKWEIGHT.dat",flag->data_folder,vdt);
	fid=fopen(fname,"r");
	fread(map_exportw,sizeof(float),flag->xres*flag->nproj,fid);
	fclose(fid);
			
	for (int j=0; j < flag->nproj; j++){
     for (int i=0; i < flag->xres; i++){
    	flag->kmap->m[j][i].i = map_exportx[flag->xres*j + i];
		flag->kmap->m[j][i].j = map_exporty[flag->xres*j + i];
		flag->kmap->m[j][i].k = map_exportz[flag->xres*j + i];
		flag->kmap->m[j][i].wt= map_exportw[flag->xres*j + i];
	}}
		
	free(map_exportx);
	free(map_exporty);
	free(map_exportz);
	free(map_exportw);
}