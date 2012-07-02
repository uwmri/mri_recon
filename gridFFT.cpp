/************************************************
Gridding and FFT Libraries for K-space to Image Domain Tranformation

Initial Author: 
	Kevin M. Johnson

Usage Example: 
	// Do Once
	gridFFT gridding;
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(64,64,64,3,KAISER);
	
	// Use for transform
	gridding.forward()


*************************************************/

#include "gridFFT.h"

//----------------------------------------
// Deconstructor - Free's Memory
//----------------------------------------

gridFFT::~gridFFT(){
  delete [] winx;
  delete [] winy;
  delete [] winz;
  delete [] grid_filterX;
  delete [] grid_filterY;
  delete [] grid_filterZ;
}

//----------------------------------------
// Constructor - Sets Default Vals
//----------------------------------------
gridFFT::gridFFT(){
	grid_x=1.375;
  	grid_y=1.375;
  	grid_z=1.375;
	overgrid = 1.375;
	kernel_type = KAISER;
  	betaX=12;
	betaY=12;
	betaZ=12;
	dwinX=-1;
	dwinY=-1;
	dwinZ=-1;
	winx=NULL;
	winy=NULL;
	winz=NULL;
	grid_filterX=NULL;
	grid_filterY=NULL;
	grid_filterZ=NULL;
	fft_plan=NULL;
	ifft_plan=NULL;
	k_rad=9999.0;
}

//----------------------------------------
// Allocate Memory for Gridding
//----------------------------------------

void gridFFT::alloc_grid(){
	k3d_grid.alloc(Sz,Sy,Sx);
	k3dref.alloc(Nz,Ny,Nx);
}

//----------------------------------------
// Return Array
//----------------------------------------

array3D< complex<float> >gridFFT::return_array( void){
	return( k3dref );		 
}

//----------------------------------------
// FFT Planning - Based on FFTW Library
//----------------------------------------

void gridFFT::plan_fft( void ){
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(THREADS);
    	
	//- Load old Plan if Possible
	FILE *fid;
	FILE *fid2;
	char fname[300];
	char hname[300];
	char com[1300];
	gethostname( hname,299);
#ifdef SCANNER
	sprintf(fname,"/usr/g/research/pcvipr/fft_wisdom_host_%s_x%d_y%d_z%d.dat",hname,Sx,Sy,Sz);
#else
	sprintf(fname,"/export/home/kmjohnso/FFT_PLANS/fft_wisdom_host_%s_x%d_y%d_z%d.dat",hname,Sx,Sy,Sz);
#endif
	printf("The FFT File will be %s\n",fname);
	
				
	if(  (fid=fopen(fname,"r")) != NULL){
		fftwf_import_wisdom_from_file(fid);
		fclose(fid);
	}	
	
	cout << " Planning FFT " << endl;
	fft_plan = fftwf_plan_dft_3d(Sz, Sy, Sx,(fftwf_complex *)k3d_grid[0][0],(fftwf_complex*)k3d_grid[0][0],FFTW_FORWARD, FFTW_MEASURE);
	
	cout << " Planning Inverse FFT" << endl;
	ifft_plan = fftwf_plan_dft_3d(Sz, Sy, Sx,(fftwf_complex *)k3d_grid[0][0],(fftwf_complex*)k3d_grid[0][0],FFTW_BACKWARD, FFTW_MEASURE);
		
	/*In case New Knowledge Was Gained*/	
	if( (fid2 = fopen(fname, "w")) == NULL){
		printf("Could Not Export FFT Wisdom\n");
	}else{
		fftwf_export_wisdom_to_file(fid2);
		fclose(fid2);
		sprintf(com,"chmod 777 %s",fname);
		system(com); 
	}
	
	return;
}


//----------------------------------------
// Parse Command Line Args
//----------------------------------------

void gridFFT::read_commandline(int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if (strcmp("-h", pstring[pos] ) == 0) {
	  	printf("\n*********************************************\n");
	  	printf("Gridding Control:\n");
	  	printf("*********************************************\n");
	  	
		printf("-kaiser                :use kaiser bessel kernel\n");
	  	printf("-triangle              :use triangle kernel\n");
	  	printf("-gauss                 :use gaussian kernel\n");
	  	printf("-gridx #               :overgrid factor in x\n");
	  	printf("-gridy #               :overgrid factor in y\n");
	  	printf("-gridz #               :overgrid factor in z\n");
	  	printf("-fast_grid             :no overgriding with narrow kernel\n");
	  	printf("-moderate_grid         :moderate kernel and overgrid\n");
		printf("-large_grid            :use large oversampling with smaller kernel\n");
		
		float_flag("-grid_x",grid_x);
		float_flag("-grid_y",grid_y);
		float_flag("-grid_z",grid_z);
		trig_flag(KAISER,"-kaiser",kernel_type);
		trig_flag(TRIANGLE,"-triangle",kernel_type);
		
	// Special Copies
	}else if(strcmp("-fast_grid", pstring[pos]) == 0) {
	  	overgrid = 1.0;
		kernel_type = TRIANGLE;
		dwinX = 1.0;
		dwinY = 1.0;
		dwinZ = 1.0;
		
	}else if(strcmp("-moderate_grid", pstring[pos]) == 0) {
	  	overgrid = 1.25;
		dwinX = 1.5;
		dwinY = 1.5;
		dwinZ = 1.5;
				
	}else if(strcmp("-large_grid", pstring[pos]) == 0) {
	  	overgrid = 2.0;
		dwinX = 4;
		dwinY = 4;
		dwinZ = 4;
	}
  }
}    


//----------------------------------------
//    Setup for Gridding 
//----------------------------------------

void gridFFT::precalc_gridding(int NzT,int NyT,int NxT, int directions){
  
  // How many pts in kernel per delta k
  grid_modX = 600;
  grid_modY = 600;
  grid_modZ = 600;
  
  // Copy to Class
  Nx = NxT;
  Ny = NyT;
  Nz = NzT;
    
  // Get rounded Gridding ratio*
  grid_x =  16.0*ceil( ( overgrid * (float)Nx )/16.0	) / (float)Nx;
  grid_y =  16.0*ceil( ( overgrid * (float)Ny )/16.0	) / (float)Ny;
  grid_z =  16.0*ceil( ( overgrid * (float)Nz )/16.0	) / (float)Nz;
  
  // Compute Grid Size 
  Sz = (int)(grid_z *Nz);
  Sy = (int)(grid_y *Ny);
  Sx = (int)(grid_x *Nx);
  
  // ------------------------------------------------
  //    Kernel Calculations
  // ------------------------------------------------
  
  switch(kernel_type){
  		case(TRIANGLE):{
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 1.0) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 1.0) : ( dwinY );
		dwinZ   = (dwinZ == -1 ) ? ( 1.0) : ( dwinZ );
		
		// Grid Length for precomputed kernel
		int grid_lengthX = (int)( (float)dwinX*(float)grid_modX);
		int grid_lengthY = (int)( (float)dwinY*(float)grid_modY);
		int grid_lengthZ = (int)( (float)dwinZ*(float)grid_modZ);
		
		// Alloc Structs for Gridding
		grid_filterX= new float[ grid_lengthX+10];
		grid_filterY= new float[ grid_lengthY+10];
		grid_filterZ= new float[ grid_lengthZ+10];
		memset(grid_filterX, 0, (size_t)( (int)(grid_lengthX+10)*sizeof(float)));
  		memset(grid_filterY, 0, (size_t)( (int)(grid_lengthY+10)*sizeof(float)));
  		memset(grid_filterZ, 0, (size_t)( (int)(grid_lengthZ+10)*sizeof(float)));
    
		// Compute Seperable Kernel
		for(int i=0; i<(grid_lengthX+1); i++){
			float grid_pos = (float)i / (float)grid_lengthX;
			grid_filterX[i]  = 1.0 - grid_pos;
	 	}
	
		for(int i=0; i<(grid_lengthY+1); i++){
			float grid_pos = (float)i / (float)grid_lengthY;
			grid_filterY[i]  = 1.0 - grid_pos;
	 	}
	
		for(int i=0; i<(grid_lengthZ+1); i++){
			float grid_pos = (float)i / (float)grid_lengthZ;
			grid_filterZ[i]  = 1.0 - grid_pos;
	 	}
	}break;
	
	case(KAISER):{
		
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 2.5) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 2.5) : ( dwinY );
		dwinZ   = (dwinZ == -1 ) ? ( 2.5) : ( dwinZ );
	
		// Grid Length for precomputed kernel
		int grid_lengthX = (int)( (float)dwinX*(float)grid_modX);
		int grid_lengthY = (int)( (float)dwinY*(float)grid_modY);
		int grid_lengthZ = (int)( (float)dwinZ*(float)grid_modZ);
		
		// Alloc Structs for Gridding
		grid_filterX= new float[ grid_lengthX+10];
		grid_filterY= new float[ grid_lengthY+10];
		grid_filterZ= new float[ grid_lengthZ+10];
		memset(grid_filterX, 0, (size_t)( (int)(grid_lengthX+10)*sizeof(float)));
  		memset(grid_filterY, 0, (size_t)( (int)(grid_lengthY+10)*sizeof(float)));
  		memset(grid_filterZ, 0, (size_t)( (int)(grid_lengthZ+10)*sizeof(float)));
    
		// Get optimal Beta per Beatty et al
	 	betaX = PI*sqrtf(  (dwinX*dwinX)/(grid_x*grid_x)*(grid_x -0.5)*(grid_x-0.5) - 0.8);
	 	betaY = PI*sqrtf(  (dwinY*dwinY)/(grid_y*grid_y)*(grid_y -0.5)*(grid_y-0.5) - 0.8);
	 	betaZ = PI*sqrtf(  (dwinZ*dwinZ)/(grid_z*grid_z)*(grid_z -0.5)*(grid_z-0.5) - 0.8);
		
		float beta_minX = sqrt(pow(PI*2*dwinX/grid_x,2.0) -  9.6752);
		float beta_minY = sqrt(pow(PI*2*dwinY/grid_y,2.0) -  9.6752);
		float beta_minZ = sqrt(pow(PI*2*dwinZ/grid_z,2.0) -  9.6752);
		betaX = ( beta_minX > betaX) ? ( beta_minX ) : ( betaX);
		betaY = ( beta_minY > betaY) ? ( beta_minY ) : ( betaY);
		betaZ = ( beta_minZ > betaZ) ? ( beta_minZ ) : ( betaZ);
	
		// Compute Seperable Kernels
		for(int i=0; i<(grid_lengthX+1); i++){
			float grid_pos=  ( (float)i )/( (float)grid_lengthX);
			float grid_arg = sqrtf( 1.0 - grid_pos*grid_pos );
			grid_filterX[i] = bessi0(betaX*grid_arg)/bessi0(betaX);
		}
	
		for(int i=0; i<(grid_lengthY+1); i++){
			float grid_pos=  ( (float)i )/( (float)grid_lengthY);
			float grid_arg = sqrtf( 1.0 - grid_pos*grid_pos );
			grid_filterY[i] = bessi0(betaY*grid_arg)/bessi0(betaY);
		}
		
		for(int i=0; i<(grid_lengthZ+1); i++){
			float grid_pos=  ( (float)i )/( (float)grid_lengthZ);
			float grid_arg = sqrtf( 1.0 - grid_pos*grid_pos );
			grid_filterZ[i] = bessi0(betaZ*grid_arg)/bessi0(betaZ);
		}
		
	}break;
  }
  
  // ------------------------------------------------
  //    Image Domain Calcs (crop + deapp)
  // ------------------------------------------------
  
  
  // Calculations to Determine Crop Positions
  og_sx =  (int)( (float)Nx*(grid_x - 1)/ 2.0);
  og_sy =  (int)( (float)Ny*(grid_y - 1)/ 2.0);
  og_sz =  (int)( (float)Nz*(grid_z - 1)/ 2.0);
  og_ex =  og_sx + Nx;
  og_ey =  og_sy + Ny;
  og_ez =  og_sz + Nz;
  
  printf("\n\nGridding Kernel Info\n");
  printf("Dwin 		%f %f %f\n",dwinX,dwinY,dwinZ); 
  printf("Mod 		%f %f %f\n",grid_modX,grid_modY,grid_modZ); 
  printf("Og %d-%d x %d-%d x %d-%d\n",og_sx,og_ex,og_sy,og_ey,og_sz,og_ez);  
  
  // Deapp Windows
  winx = new float[Sx];
  winy = new float[Sy];
  winz = new float[Sz];

  for( int i = 0; i < Sx;i++){
  	winx[i] = 0.0;
  	float ipos = i - (float)Sx/2.0;
	for(int grid_pos = 0; grid_pos < dwinX*grid_modX; grid_pos++){ 
		winx[i] += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modX / (float)Sx)*grid_filterX[grid_pos];
	}
	winx[i] = (float)grid_modX/winx[i];
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winx[i]*=fact / Sx;
  }

  for( int i = 0; i < Sy;i++){
  	winy[i] = 0.0;
  	float ipos = i - (float)Sy/2.0;
	for(int grid_pos = 0; grid_pos < dwinY*grid_modY; grid_pos++){ 
		winy[i] += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modY / (float)Sy)*grid_filterY[grid_pos];
	}
	winy[i] = (float)grid_modY/winy[i];
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winy[i]*=fact / Sy;
  }
   
  for( int i = 0; i < Sz;i++){
  	winz[i] = 0.0;
  	float ipos = i - (float)Sz/2.0;
	for(int grid_pos = 0; grid_pos < dwinZ*grid_modZ; grid_pos++){ 
		winz[i] += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modZ / (float)Sz)*grid_filterZ[grid_pos];
	}
	winz[i]  = (float)grid_modZ/winz[i];
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winz[i]*=fact / Sz;
  }
  
  // Allocate Memory
  cout << "Alloc Grid" << endl;
  alloc_grid();
  
  // Setup FFT	
  plan_fft();
  
}

#include "sys/time.h"
double gettime(void){
   	timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec % 1000000000) + tv.tv_usec/1000000.0;
};
  
//----------------------------------------
//    Forward Transform
//----------------------------------------

void gridFFT::forward(complex<float> *data, float *kx, float *ky, float *kz, float *kw,int Npts){
	k3d_grid.zero();
	grid_forward(data,kx,ky,kz,kw,Npts);
	chop();
	fftwf_execute(fft_plan);
	deapp_chop_crop();
}	

void gridFFT::backward(complex<float> *data, float *kx, float *ky, float *kz, float *kw,int Npts){
	k3d_grid.zero();
	icrop_deapp_chop();
	fftwf_execute(ifft_plan);
	chop();
	grid_backward(data,kx,ky,kz,kw,Npts);
}

void gridFFT::chop(void){
	#pragma omp parallel for 
	for(int k=0; k < Sz; k++){
      for (int j=0; j< Sy; j++){
    	for (int i=0; i < Sx; i++){
			float fact =  ((float)( 2*(( i + j + k )%2) - 1));
			k3d_grid[k][j][i] *= fact;
	}}}
}

//----------------------------------------
//    Crop from Gridding Matrix to Image
//----------------------------------------
void gridFFT::deapp_chop_crop(){
	
	#pragma omp parallel for 
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz[k+og_sz];
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy[j+og_sy];
		for(int i=0; i<Nx; i++) {
			k3dref[k][j][i] = ( k3d_grid[k+og_sz][j+og_sy][i+og_sx]*(wty*winx[i+og_sx]));
	}}}
}

//----------------------------------------
//    Crop from Gridding Matrix to Image
//----------------------------------------
void gridFFT::icrop_deapp_chop(){
	#pragma omp parallel for 
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz[k+og_sz];
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy[j+og_sy];
		for(int i=0; i<Nx; i++) {
			k3d_grid[k+og_sz][j+og_sy][i+og_sx] = ( k3dref[k][j][i]*(wty*winx[i+og_sx]));
	}}}
}

//-----------------------------------------------------
//  Multiple by Precalculated Deappodization Window
//-----------------------------------------------------
void gridFFT::deapp_chop(){
 	#pragma omp parallel for 
	for(int k=0; k<Sz; k++){ 
	 float wtz = winz[k];
     for(int j=0; j<Sy; j++){ 
   	  float wty = wtz*winy[j];
	  for(int i=0; i< Sx; i++) {
    	k3d_grid[k][j][i] *= winx[i]*wty;
	}}}
}

// -------------------------------------------------------
//  This is the main function for Gridding.  Assumes data
//  is already density compensated, etc.
// -------------------------------------------------------

void gridFFT::grid_forward( complex<float> *data, float *kx, float *ky, float *kz, float *kw,int Npts){

	float cx = Sx/2 -1;
	float cy = Sy/2 -1;
	float cz = Sz/2 -1;
	
	#pragma omp parallel for schedule(dynamic,1024) 
	for (int i=0; i < Npts; i++) {
      	
		complex<float>temp =data[i];
				
		// Do not grid zeros
     	if( temp==complex<float>(0,0)) continue;
		
		// Density Comp
		temp *= kw[i];
				
	    // Calculate the exact kspace sample point in 
	    // dimension flag->grid* kspace that this point (i,j)
	    // is contributing too.
	   		
		// Compute Coordinates + Check
		float dkx = kx[i]*grid_x + cx;
		int sx = (int)ceil( dkx - dwinX);
		if(sx <0)   continue;
		int ex = (int)floor(dkx + dwinX);
		if(ex >= Sx) continue;  
	    
		// Compute Coordinates + Check
		float dky = ky[i]*grid_y + cy;
		int sy = (int)ceil( dky - dwinY);
		if(sy <0)   continue;
		int ey = (int)floor(dky + dwinY);
		if(ey >= Sy) continue;  
	    
		// Compute Coordinates + Check
		float dkz = kz[i]*grid_z + cz;
		int sz = (int)ceil( dkz - dwinZ);
		if(sz <0)   continue;
		int ez = (int)floor(dkz + dwinZ);
		if(ez >= Sz) continue;  
		
		float kr = kx[i]*kx[i] + ky[i]*ky[i] + kz[i]*kz[i];
		temp *= exp( -kr / (2.0*k_rad*k_rad) );

		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    		float delz = fabs(grid_modZ*(dkz -(float)lz));
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ[(int)delz]*( 1.0-dz) + grid_filterZ[(int)delz +1]*dz;
			
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly));
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY[(int)dely]*( 1.0-dy) + grid_filterY[(int)dely +1]*dy );
			 	 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx));
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX[(int)delx]*( 1.0-dx) + grid_filterX[(int)delx +1]*dx );
			 		
					complex<float>temp2 = wtx*temp;
					float RD = real(temp2);
					float ID = imag(temp2);
					float *I = reinterpret_cast<float *>(&k3d_grid.vals[lz][ly][lx]);
					float *R = I++;
					
					// Prevent Race conditions in multi-threaded
					#pragma omp atomic
					*R+=RD;
					
					#pragma omp atomic
					*I+=ID;
																								
					/*This Memory Access is the Bottleneck - Also not thread safe!*/	 			 
			 		// k3d_grid.vals[lz][ly][lx]+=temp2;
			}/* end lz loop */
	  	  }/* end ly */
		 }/* end lx */
	}/* end data loop */
	return;
}
	
void gridFFT::grid_backward( complex<float> *data, float *kx, float *ky, float *kz,float *kw,int Npts){

	float cx = Sx/2 -1;
	float cy = Sy/2 -1;
	float cz = Sz/2 -1;
	
	#pragma omp parallel for 
	for (int i=0; i < Npts; i++) {
      	
		// Calculate the exact kspace sample point in 
	    // dimension flag->grid* kspace that this point (i,j)
	    // is contributing too.
	   	data[i] = complex<float>(0,0);
				
		// Compute Coordinates + Check
		float dkx = kx[i]*grid_x + cx;
		int sx = (int)ceil( dkx - dwinX);
		if(sx <0)   continue;
		int ex = (int)floor(dkx + dwinX);
		if(ex >= Sx) continue;  
	    
		// Compute Coordinates + Check
		float dky = ky[i]*grid_y + cy;
		int sy = (int)ceil( dky - dwinY);
		if(sy <0)   continue;
		int ey = (int)floor(dky + dwinY);
		if(ey >= Sy) continue;  
	    
		// Compute Coordinates + Check
		float dkz = kz[i]*grid_z + cz;
		int sz = (int)ceil( dkz - dwinZ);
		if(sz <0)   continue;
		int ez = (int)floor(dkz + dwinZ);
		if(ez >= Sz) continue;  
		
		complex<float>temp(0,0);
					
		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    		float delz = fabs(grid_modZ*(dkz -(float)lz));
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ[(int)delz]*( 1.0-dz) + grid_filterZ[(int)delz +1]*dz;
			
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly));
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY[(int)dely]*( 1.0-dy) + grid_filterY[(int)dely +1]*dy );
			 	 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx));
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX[(int)delx]*( 1.0-dx) + grid_filterX[(int)delx +1]*dx );
			 											
					/*This Memory Access is the Bottleneck*/	 			 
			 		temp += wtx*k3d_grid[lz][ly][lx];
			  	 
	    	}/* end lz loop */
	  	  }/* end ly */
		 }/* end lx */
		 
		 //float kr = kx[i]*kx[i] + ky[i]*ky[i] + kz[i]*kz[i];
		 //temp *= exp( -kr / (2.0*k_rad*k_rad) );
		 //temp *= kw[i];
		 
		 data[i] = temp;
	
	}/* end data loop */
	return;
}	
	


// For Kaiser Bessel Window
float gridFFT::bessi0(float x)
{
  float ax, ans;
  double y;


  if ((ax=fabs(x)) < 3.75) {
    y = x/3.75;
    y*=y;
    ans = 1.0 + y*(3.5156229+y*(3.0899424+y*(1.2067492 +
	  y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    ans=(exp((double)ax)/sqrt((double)ax))*(0.39894228+y*(0.1328592e-1 +
	 y*(0.225319e-2+y*(-0.157565e-2 + y*(0.916281e-2 +
	 y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1 +
         y*0.392377e-2))))))));
  }
  return(ans);
}







