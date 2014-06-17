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
#include "io_templates.hpp"
#include "tictoc.hpp"
using namespace NDarray;

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
	kernel_type = KAISER_KERNEL;
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
	fft_in_z=1;
  	fft_in_y=1;
  	fft_in_z=1;
  	grid_in_x=1;
  	grid_in_y=1;
  	grid_in_z=1;
	k_rad=9999.0;
	time_grid =0 ;
	double_grid =0;
}

//----------------------------------------
// Allocate Memory for Gridding
//----------------------------------------

void gridFFT::alloc_grid(){
	k3d_grid.setStorage( ColumnMajorArray<3>());
	k3d_grid.resize(Sx,Sy,Sz);
	k3d_grid = 0;

	Array< complex<float>,3>image2 = k3d_grid( Range(og_sx,og_ex-1),Range(og_sy,og_ey-1),Range(og_sz,og_ez-1));
	image.reference(image2);
}


//----------------------------------------
// FFT Planning - Based on FFTW Library
//----------------------------------------

void gridFFT::plan_fft( void ){
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    cout << "FFT Planning with " << omp_get_max_threads() << " threads"<< endl;
		
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
	
	cout << "Test" << endl;
	fftwf_complex *ptr = reinterpret_cast<fftwf_complex*>(k3d_grid.data());
		
	cout << " Planning FFT " << endl << flush; 
	fft_plan = fftwf_plan_dft_3d(Sz,Sy,Sx,ptr,ptr,FFTW_FORWARD, FFTW_MEASURE);
	
	cout << " Planning Inverse FFT" << endl << flush;
	ifft_plan = fftwf_plan_dft_3d(Sz,Sy,Sx,ptr,ptr,FFTW_BACKWARD, FFTW_MEASURE);
		
	/*In case New Knowledge Was Gained*/	
	if( (fid2 = fopen(fname, "w")) == NULL){
		printf("Could Not Export FFT Wisdom\n");
	}else{
		fftwf_export_wisdom_to_file(fid2);
		fclose(fid2);
		sprintf(com,"chmod 777 %s",fname);
		if( system(com) != 1 ){
			cout << "Failed to Change FFT Plan Permissions" << endl;
		}
	}
	
	return;
}


// ----------------------
// Help Message
// ----------------------
void gridFFT::help_message(void){
	cout << "----------------------------------------------" << endl;
	cout << "   Gridding Control " << endl;
	cout << "----------------------------------------------" << endl;
	
	cout<<"Control" << endl;
	help_flag("-kaiser","use kaiser bessel kernel");
 	help_flag("-triangle","use triangle kernel");
 	help_flag("-sinc","use sinc kernel");
 	help_flag("-overgrid []","overgrid by factor []");
 	help_flag("-fast_grid","no overgrid,traingle kernel");
	help_flag("-time_grid","output times for gridding");
			
	cout<<"Fine Control" << endl;
	help_flag("-grid_in_x []","0=use nearest neighbor interpolation in x");
	help_flag("-grid_in_y []","0=use nearest neighbor interpolation in y");
	help_flag("-grid_in_z []","0=use nearest neighbor interpolation in z");
	help_flag("-fft_in_x []","0=no fft in x");
	help_flag("-fft_in_y []","0=no fft in y");
	help_flag("-fft_in_z []","0=no fft in z");
	help_flag("-dwinX []","Size of kernel in x");
	help_flag("-dwinY []","Size of kernel in y");
	help_flag("-dwinZ []","Size of kernel in z");
	help_flag("-double_grid","Use complex<double> for forward gridding (slow)");
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
	 
		float_flag("-overgrid",overgrid);
		float_flag("-dwinX",dwinX);
		float_flag("-dwinY",dwinY);
		float_flag("-dwinZ",dwinZ);
		
		int_flag("-grid_in_x",grid_in_x);
		int_flag("-grid_in_y",grid_in_y);
		int_flag("-grid_in_z",grid_in_z);
		int_flag("-fft_in_x",fft_in_x);
		int_flag("-fft_in_y",fft_in_y);
		int_flag("-fft_in_z",fft_in_z);
				
		trig_flag(KAISER_KERNEL,"-kaiser",kernel_type);
		trig_flag(TRIANGLE_KERNEL,"-triangle",kernel_type);
		trig_flag(SINC_KERNEL,"-sinc",kernel_type);
		trig_flag(1,"-time_grid",time_grid);
		trig_flag(1,"-double_grid",double_grid);
	// Special Copies
	}else if(strcmp("-fast_grid", pstring[pos]) == 0) {
	  	overgrid = 1.0;
		kernel_type = TRIANGLE_KERNEL;
		dwinX = 1.0;
		dwinY = 1.0;
		dwinZ = 1.0;
	}
  }
}    


//----------------------------------------
//    Setup for Gridding 
//----------------------------------------

void gridFFT::precalc_kernel(int NzT,int NyT,int NxT, int directions){
  
    // Copy to Class
   Nx = NxT;
   Ny = NyT;
   Nz = NzT;
  
  // How many pts in kernel per delta k for kernel lookup table
  grid_modX = 600;
  grid_modY = 600;
  grid_modZ = 600;
 
  // Get rounded Gridding ratio*
  if(grid_in_x ==1){
  	grid_x =  16.0*ceil( ( overgrid * (float)Nx )/16.0	) / (float)Nx;
  }else{
  	grid_x = 1;
  }
  
  if(grid_in_y ==1){
  	grid_y =  16.0*ceil( ( overgrid * (float)Ny )/16.0	) / (float)Ny;
  }else{
  	grid_y = 1;
  }
  
  if(grid_in_z ==1){
  	grid_z =  16.0*ceil( ( overgrid * (float)Nz )/16.0	) / (float)Nz;
  }else{
  	grid_z = 1;
  }
  
  // Compute Grid Size 
  Sz = (int)(grid_z *Nz);
  Sy = (int)(grid_y *Ny);
  Sx = (int)(grid_x *Nx);
  
  // ------------------------------------------------
  //    Kernel Calculations
  // ------------------------------------------------
  
  switch(kernel_type){
  		case(TRIANGLE_KERNEL):{
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 1.0) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 1.0) : ( dwinY );
		dwinZ   = (dwinZ == -1 ) ? ( 1.0) : ( dwinZ );
		
		// Grid Length for precomputed kernel
		int grid_lengthX = (int)( (float)dwinX*(float)grid_modX);
		int grid_lengthY = (int)( (float)dwinY*(float)grid_modY);
		int grid_lengthZ = (int)( (float)dwinZ*(float)grid_modZ);
		
		// Alloc Lookup Table Structs for Gridding
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
	
	case(KAISER_KERNEL):{
		
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
	
	case(SINC_KERNEL):{
		
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 15.0) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 15.0) : ( dwinY );
		dwinZ   = (dwinZ == -1 ) ? ( 15.0) : ( dwinZ );
	
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
    		
		
		// Compute Seperable Kernels
		for(int i=0; i<(grid_lengthX+1); i++){
			double grid_pos=  ( (double)i )/( (double)grid_lengthX);
			if(grid_pos < 0.01){
				grid_filterX[i] = 1.0;
			}else{
				grid_filterX[i]=(float)( sin(grid_pos*PI)/(grid_pos*PI));
			}
		}
		
		for(int i=0; i<(grid_lengthY+1); i++){
			double grid_pos=  ( (double)i )/( (double)grid_lengthY);
			if(grid_pos < 0.01){
				grid_filterY[i] = 1.0;
			}else{
				grid_filterY[i]=(float)( sin(grid_pos*PI)/(grid_pos*PI));
			}
		}
	
		for(int i=0; i<(grid_lengthZ+1); i++){
			double grid_pos=  ( (double)i )/( (double)grid_lengthZ);
			if(grid_pos < 0.01){
				grid_filterZ[i] = 1.0;
			}else{
				grid_filterZ[i]=(float)( sin(grid_pos*PI)/(grid_pos*PI));
			}
		}
		
	}break;
	
  }
  
}
 
void gridFFT::precalc_gridding(int NzT,int NyT,int NxT, int directions){
  
    
  precalc_kernel(NzT,NyT,NxT,directions);
  
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
  for( int i = 0; i < Sx;i++){
  	winx[i] = 0.0;
  	float ipos = i - (float)Sx/2.0;
	for(int grid_pos = 0; grid_pos < dwinX*grid_modX; grid_pos++){ 
		// Fourier Transform of Kernel
		winx[i] += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modX / (float)Sx)*grid_filterX[grid_pos];
	}
	winx[i] = (float)grid_modX/winx[i];
	
	// Put chopping + zeroing in window to save time
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winx[i]*=fact / Sx;
	winx[i]*=( i < og_sx) ? ( 0.0 ) : ( 1.0);
	winx[i]*=( i > og_ex-1) ? ( 0.0 ) : ( 1.0);
	
	if(grid_in_x==0){
		winx[i]=1.0;
	}
  }

  winy = new float[Sy];
  for( int i = 0; i < Sy;i++){
  	winy[i] = 0.0;
  	float ipos = i - (float)Sy/2.0;
	for(int grid_pos = 0; grid_pos < dwinY*grid_modY; grid_pos++){ 
		winy[i] += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modY / (float)Sy)*grid_filterY[grid_pos];
	}
	winy[i] = (float)grid_modY/winy[i];
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winy[i]*=fact / Sy;
  	winy[i]*=( i < og_sy) ? ( 0.0 ) : ( 1.0);
	winy[i]*=( i > og_ey-1) ? ( 0.0 ) : ( 1.0);
	
	if(grid_in_y==0){
		winy[i]=1.0;
	}
  }
   
  winz = new float[Sz];
  for( int i = 0; i < Sz;i++){
  	winz[i] = 0.0;
  	float ipos = i - (float)Sz/2.0;
	for(int grid_pos = 0; grid_pos < dwinZ*grid_modZ; grid_pos++){ 
		winz[i] += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modZ / (float)Sz)*grid_filterZ[grid_pos];
	}
	winz[i]  = (float)grid_modZ/winz[i];
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winz[i]*=fact / Sz;
	winz[i]*=( i < og_sz) ? ( 0.0 ) : ( 1.0);
	winz[i]*=( i > og_ez-1) ? ( 0.0 ) : ( 1.0);
	
	if(grid_in_z==0){
		winz[i]=1.0;
	}
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
//    Gridding Calls
//----------------------------------------

void gridFFT::forward( Array<complex<float>,3>&X,\
					   const Array<complex<float>,3>&smap,\
					   const Array<complex<float>,3>&data,\
					   const Array<float,3>&kx,\
					   const Array<float,3>&ky,\
					   const Array<float,3>&kz,\
					   const Array<float,3>&kw){
	tictoc T;
	if(time_grid) T.tic(); 
	k3d_grid=0; // Zero K-Space
	if(time_grid) cout << "Forward::zero:"<< T << endl << flush;
	
	if(time_grid) T.tic(); 
	chop_grid_forward(data,kx,ky,kz,kw); // Grid data to K-Space
	if(time_grid)cout << ",grid:"<< T << endl << flush;
	
	if(time_grid) T.tic(); 
	fftwf_execute(fft_plan); //FFT 
	if(time_grid)cout << ",fft:"<< T << endl << flush;
	
	if(time_grid) T.tic(); 
	forward_image_copy(X,smap); // Copy result to X
	if(time_grid)cout << ",copy:"<< T << endl<< flush;
}

// No Sensitivity map
void gridFFT::forward( Array<complex<float>,3>&X,\
					   const Array<complex<float>,3>&data,\
					   const Array<float,3>&kx,\
					   const Array<float,3>&ky,\
					   const Array<float,3>&kz,\
					   const Array<float,3>&kw){
	tictoc T;
	if(time_grid) T.tic(); 
	k3d_grid=0; // Zero K-Space
	if(time_grid) cout << "Forward::zero:"<< T << flush;
	
	if(time_grid) T.tic(); 
	chop_grid_forward(data,kx,ky,kz,kw); // Grid data to K-Space
	if(time_grid)cout << ",grid:"<< T<< flush;
	
	if(time_grid) T.tic(); 
	fftwf_execute(fft_plan); //FFT 
	if(time_grid)cout << ",fft:"<< T<< flush;
	
	if(time_grid) T.tic(); 
	forward_image_copy(X); // Copy result to X
	if(time_grid)cout << ",copy:"<< T << endl<< flush;
}

void gridFFT::backward(const Array<complex<float>,3>&X,\
					   const Array<complex<float>,3>&smap,\
					   Array<complex<float>,3>&data,\
					   const Array<float,3>&kx,\
					   const Array<float,3>&ky,\
					   const Array<float,3>&kz,\
					   const Array<float,3>&kw){
	
	tictoc T;
	if(time_grid) T.tic(); 
	backward_image_copy(X,smap); // Copy image to gridding 
	if(time_grid) cout << "Backward::copy:"<< T;
	
	if(time_grid) T.tic(); 
	fftwf_execute(ifft_plan); // Do FFT 
	if(time_grid)cout << ",ifft:"<< T;
	
	if(time_grid) T.tic(); 
	chop_grid_backward(data,kx,ky,kz,kw); // Inverse gridding
	if(time_grid)cout << ",igrid:"<< T << endl;
	
}

// Same but without smap
void gridFFT::backward(const Array<complex<float>,3>&X,\
					   Array<complex<float>,3>&data,\
					   const Array<float,3>&kx,\
					   const Array<float,3>&ky,\
					   const Array<float,3>&kz,\
					   const Array<float,3>&kw){
	
	tictoc T;
	if(time_grid) T.tic(); 
	backward_image_copy(X); // Copy image to gridding 
	if(time_grid) cout << "Backward::copy:"<< T;
	
	if(time_grid) T.tic(); 
	fftwf_execute(ifft_plan); // Do FFT 
	if(time_grid)cout << ",ifft:"<< T;
	
	if(time_grid) T.tic();
	chop_grid_backward(data,kx,ky,kz,kw); // Inverse gridding
	if(time_grid)cout << ",igrid:"<< T << endl;
}

//----------------------------------------
//    Crop from Gridding Matrix to Image
//----------------------------------------
void gridFFT::forward_image_copy(Array<complex<float>,3>&X){
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz[k+og_sz];
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy[j+og_sy];
		for(int i=0; i<Nx; i++) {
			X(i,j,k) = ( image(i,j,k)*wty*winx[i+og_sx]); /* Don't just sum coils when no sense map is given*/
	}}}
}

void gridFFT::forward_image_copy(Array<complex<float>,3>&X,const Array<complex<float>,3>&smap){
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz[k+og_sz];
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy[j+og_sy];
		for(int i=0; i<Nx; i++) {
			X(i,j,k) += ( image(i,j,k)*wty*winx[i+og_sx]*conj(smap(i,j,k)));
	}}}
}

//----------------------------------------
//    Copy Image to gridding Matrix, Multiply by Smaps, and Zero 
//----------------------------------------

void gridFFT::backward_image_copy(const Array<complex<float>,3>&X){
	
	if(overgrid > 1.0){
		k3d_grid = 0;
	}
	
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz[k+og_sz];
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy[j+og_sy];
		for(int i=0; i<Nx; i++) {
			image(i,j,k) = X(i,j,k)*wty*winx[i+og_sx];
	}}}
}

void gridFFT::backward_image_copy(const Array<complex<float>,3>&X,const Array<complex<float>,3>&smap){
	
	if(overgrid > 1.0){
		k3d_grid = 0;
	}
	
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz[k+og_sz];
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy[j+og_sy];
		for(int i=0; i<Nx; i++) {
			image(i,j,k) =  X(i,j,k)*wty*winx[i+og_sx]*smap(i,j,k);
	}}}
}

// -------------------------------------------------------
//  This is the main function for Gridding.  Assumes data
//  is already density compensated, etc.
// -------------------------------------------------------

void gridFFT::chop_grid_forward( const Array<complex<float>,3>&dataA, const Array<float,3>&kxA,const Array<float,3>&kyA,const Array<float,3>&kzA,const Array<float,3>&kwA){

	float cx = Sx/2;
	float cy = Sy/2;
	float cz = Sz/2;
	
	if( !dataA.isStorageContiguous()){
		cout << "Non-contiguous storage doesn't work yet" << endl;
		exit(1);
	}
	
	int Npts = dataA.numElements();
	const complex<float> *data = dataA.data();
	const float *kx = kxA.data();
	const float *ky = kyA.data();
	const float *kz = kzA.data();
	const float *kw = kwA.data();
	
	// For testing
	Array< complex<double>,3>tempD;
	if(double_grid==1){
		cout << "WARNING::gridding using double" << endl; 
		tempD.setStorage( ColumnMajorArray<3>());
		tempD.resize(Sx,Sy,Sz);
		tempD = 0;
	}
					
	#pragma omp parallel for 
	for (int i=0; i < Npts; i++) {
      	
		complex<float>temp =data[i];
				
		// Density Comp
		temp *= kw[i];
		
		// Do not grid zeros
     	if( temp==complex<float>(0,0)) continue;
		
				
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
		int sz,ez;
		if(grid_in_z==1){
			sz = (int)ceil( dkz - dwinZ);
			ez = (int)floor(dkz + dwinZ);
		}else{
			sz = (int)( dkz);
			ez = sz;
		}
		if(sz <0)   continue;
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
			 		
					wtx *=  ((float)( 2*(( lx + ly + lz )%2) - 1)); // Chop in  gridding now
					
					if(double_grid){
						complex<float>temp2 = wtx*temp;
						double RD = (double)real(temp2);
						double ID = (double)imag(temp2);
						double *I = reinterpret_cast<double *>(&tempD(lx,ly,lz));
						double *R = I++;
					
						// Prevent Race conditions in multi-threaded
						#pragma omp atomic
						*R+=RD;
					
						#pragma omp atomic
						*I+=ID;
					
					}else{
						complex<float>temp2 = wtx*temp;
						float RD = real(temp2);
						float ID = imag(temp2);
						float *I = reinterpret_cast<float *>(&k3d_grid(lx,ly,lz));
						float *R = I++;
					
						// Prevent Race conditions in multi-threaded
						#pragma omp atomic
						*R+=RD;
					
						#pragma omp atomic
						*I+=ID;
					}																	
					/*This Memory Access is the Bottleneck - Also not thread safe!*/	 			 
			 		// k3d_grid.vals[lz][ly][lx]+=temp2;
			}/* end lz loop */
	  	  }/* end ly */
		 }/* end lx */
	}/* end data loop */
	
	// Copy back if using double
	if(double_grid){
	for(int k=0; k< Sz; k++){ 
	  for(int j=0; j<Sy; j++){ 
	   for(int i=0; i<Sx; i++) {
			k3d_grid(i,j,k) =  tempD(i,j,k);
	}}}
	}
	
	
	return;
}


// -------------------------------------------------------
//  This is the main function for Gridding.  Assumes data
//  is already density compensated, etc.
// -------------------------------------------------------

void gridFFT::grid_backward( const Array<float,3>&imX, Array<float,3>&dataA, const Array<float,3>&kxA,const Array<float,3>&kyA,const Array<float,3>&kzA){

	float cx = imX.length(firstDim)/2;
	float cy = imX.length(secondDim)/2;
	float cz = imX.length(thirdDim)/2;
		
	if( !dataA.isStorageContiguous()){
		cout << "Non-contiguous storage doesn't work yet" << endl;
		exit(1);
	}
	
	int Npts = dataA.numElements();
	float *data = dataA.data();
	const float *kx = kxA.data();
	const float *ky = kyA.data();
	const float *kz = kzA.data();
	
	#pragma omp parallel for 
	for (int i=0; i < Npts; i++) {
      	
		float temp =0.0;
					
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
		int sz,ez;
		if(grid_in_z==1){
			sz = (int)ceil( dkz - dwinZ);
			ez = (int)floor(dkz + dwinZ);
		}else{
			sz = (int)( dkz);
			ez = sz;
		}
		if(sz <0)   continue;
		if(ez >= Sz) continue;  
		
		
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
			 		
					temp+= wtx*imX(lx,ly,lz);
					
					
			}/* end lz loop */
	  	  }/* end ly */
		 }/* end lx */
		 data[i] = temp;
		 	 
		 
	}/* end data loop */
	
	return;
}

void gridFFT::grid_forward( Array<float,3>&imX, const Array<float,3>&dataA, const Array<float,3>&kxA,const Array<float,3>&kyA,const Array<float,3>&kzA){

	float cx = imX.length(firstDim)/2;
	float cy = imX.length(secondDim)/2;
	float cz = imX.length(thirdDim)/2;
	
	if( !dataA.isStorageContiguous()){
		cout << "Non-contiguous storage doesn't work yet" << endl;
		exit(1);
	}
	
	int Npts = dataA.numElements();
	const float *data = dataA.data();
	const float *kx = kxA.data();
	const float *ky = kyA.data();
	const float *kz = kzA.data();
	
	#pragma omp parallel for 
	for (int i=0; i < Npts; i++) {
      	
		float temp =data[i];
					
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
		int sz,ez;
		if(grid_in_z==1){
			sz = (int)ceil( dkz - dwinZ);
			ez = (int)floor(dkz + dwinZ);
		}else{
			sz = (int)( dkz);
			ez = sz;
		}
		if(sz <0)   continue;
		if(ez >= Sz) continue;  
		
		
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
			 		
					float *image_temp = &imX(lx,ly,lz);					
					float temp2 = wtx*temp;
					#pragma omp atomic 
					*image_temp+= temp2;
					
					
			}/* end lz loop */
	  	  }/* end ly */
		 }/* end lx */
	}/* end data loop */
	
	return;
}




	
void gridFFT::chop_grid_backward(Array<complex<float>,3>&dataA, const Array<float,3>&kxA,const Array<float,3>&kyA,const Array<float,3>&kzA,const Array<float,3>&kwA){

	float cx = Sx/2;
	float cy = Sy/2;
	float cz = Sz/2;
	
	if( !dataA.isStorageContiguous()){
		cout << "Non-contiguous storage doesn't work yet" << endl;
		exit(1);
	}
	
	int Npts = dataA.numElements();
	complex<float> *data = dataA.data();
	const float *kx = kxA.data();
	const float *ky = kyA.data();
	const float *kz = kzA.data();
	const float *kw = kwA.data();
	
	#pragma omp parallel for 
	for (int i=0; i < Npts; i++) {
      	
		
		// Do not grid zeros
     		if( kw[i]==0.0) continue;
		
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
		int sz,ez;
		if(grid_in_z==1){
			sz = (int)ceil( dkz - dwinZ);
			ez = (int)floor(dkz + dwinZ);
		}else{
			sz = (int)( dkz);
			ez = sz;
		}
		if(sz <0)   continue;
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
			 		wtx *=  ((float)( 2*(( lx + ly + lz )%2) - 1)); // Chop in inverse gridding now
					temp += wtx*k3d_grid(lx,ly,lz);
			  	 
	    	}/* end lz loop */
	  	  }/* end ly */
		 }/* end lx */
		 data[i] += temp;
	
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







