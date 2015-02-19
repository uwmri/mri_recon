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



}

//----------------------------------------
// Constructor - Sets Default Vals
//----------------------------------------
gridFFT::gridFFT(){
	grid_x= -1;
  	grid_y= -1;
  	grid_z= -1;
	overgrid = 1.375;
	kernel_type = KAISER_KERNEL;
  	betaX=12;
	betaY=12;
	betaZ=12;
	dwinX=-1;
	dwinY=-1;
	dwinZ=-1;
	fft_plan=NULL;
	ifft_plan=NULL;
	fft_in_x=1;
  	fft_in_y=1;
  	fft_in_z=1;
  	grid_in_x=1;
  	grid_in_y=1;
  	grid_in_z=-1;
	grid_scale_x = 1.0;
	grid_scale_y = 1.0;
	grid_scale_z = 1.0;
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
	
	
	// Get FFTW Plan Folder
	char const* tmp = getenv("FFT_PLAN_PATH");
	string fft_wisdom_folder;
	if(tmp == NULL){
		cout << "Warning no FFTW Plan folder " << endl;
		cout << "  -- Set FFT_PLAN_PATH or plans may be slow" << endl;
	}else{
		fft_wisdom_folder = string(tmp);
	}
		
	// Create Name
	ostringstream stringStream;
	stringStream << fft_wisdom_folder << "/fft_wisdom_host_" << hname << "_x" << Sx << "_y" << Sy << "_z" << Sz << ".dat";
	string fft_name =  stringStream.str();
	

	cout << "The FFT File will be" << fft_name << endl;
				
	if(  (fid=fopen(fft_name.c_str(),"r")) != NULL){
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
		float_flag("-grid_x",grid_x);
		float_flag("-grid_y",grid_y);
		float_flag("-grid_z",grid_z);
		
		float_flag("-grid_scale_x",grid_scale_x);
		float_flag("-grid_scale_y",grid_scale_y);
		float_flag("-grid_scale_z",grid_scale_z);
		
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

void gridFFT::precalc_kernel(void){
   
  
  // How many pts in kernel per delta k for kernel lookup table
  grid_modX = 600;
  grid_modY = 600;
  grid_modZ = 600;
 
  
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
		grid_filterX.resize( grid_lengthX+10);
		grid_filterY.resize( grid_lengthY+10);
		grid_filterZ.resize( grid_lengthZ+10);
    	grid_filterX = 0.0;
		grid_filterY = 0.0;
		grid_filterZ = 0.0;
		
    
		// Compute Seperable Kernel
		for(int i=0; i<(grid_lengthX+1); i++){
			float grid_pos = (float)i / (float)grid_lengthX;
			grid_filterX(i)  = 1.0 - grid_pos;
	 	}
	
		for(int i=0; i<(grid_lengthY+1); i++){
			float grid_pos = (float)i / (float)grid_lengthY;
			grid_filterY(i)  = 1.0 - grid_pos;
	 	}
	
		for(int i=0; i<(grid_lengthZ+1); i++){
			float grid_pos = (float)i / (float)grid_lengthZ;
			grid_filterZ(i)  = 1.0 - grid_pos;
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
		
		// Alloc Lookup Table Structs for Gridding
		grid_filterX.resize( grid_lengthX+10);
		grid_filterY.resize( grid_lengthY+10);
		grid_filterZ.resize( grid_lengthZ+10);
    	grid_filterX = 0.0;
		grid_filterY = 0.0;
		grid_filterZ = 0.0;

		// Get optimal Beta per Beatty et al
		float act_grid_x = grid_x*grid_scale_x;
		float act_grid_y = grid_y*grid_scale_y;
		float act_grid_z = grid_z*grid_scale_z;
				
	 	betaX = PI*sqrtf(  (dwinX*dwinX)/(act_grid_x*act_grid_x)*(act_grid_x -0.5)*(act_grid_x-0.5) - 0.8);
	 	betaY = PI*sqrtf(  (dwinY*dwinY)/(act_grid_y*act_grid_y)*(act_grid_y -0.5)*(act_grid_y-0.5) - 0.8);
	 	betaZ = PI*sqrtf(  (dwinZ*dwinZ)/(act_grid_z*act_grid_z)*(act_grid_z -0.5)*(act_grid_z-0.5) - 0.8);
		
		
		float beta_minX = sqrt(pow(PI*2*dwinX/act_grid_x,2.0) -  9.6752);
		float beta_minY = sqrt(pow(PI*2*dwinY/act_grid_y,2.0) -  9.6752);
		float beta_minZ = sqrt(pow(PI*2*dwinZ/act_grid_z,2.0) -  9.6752);
		betaX = ( beta_minX > betaX) ? ( beta_minX ) : ( betaX);
		betaY = ( beta_minY > betaY) ? ( beta_minY ) : ( betaY);
		betaZ = ( beta_minZ > betaZ) ? ( beta_minZ ) : ( betaZ);
		
		
		// Compute Seperable Kernels
		for(int i=0; i<(grid_lengthX+1); i++){
			float grid_pos=  ( (float)i )/( (float)grid_lengthX);
			float grid_arg = sqrtf( 1.0 - grid_pos*grid_pos );
			grid_filterX(i) = bessi0(betaX*grid_arg)/bessi0(betaX);
		}
	
		for(int i=0; i<(grid_lengthY+1); i++){
			float grid_pos=  ( (float)i )/( (float)grid_lengthY);
			float grid_arg = sqrtf( 1.0 - grid_pos*grid_pos );
			grid_filterY(i) = bessi0(betaY*grid_arg)/bessi0(betaY);
		}
		
		for(int i=0; i<(grid_lengthZ+1); i++){
			float grid_pos=  ( (float)i )/( (float)grid_lengthZ);
			float grid_arg = sqrtf( 1.0 - grid_pos*grid_pos );
			grid_filterZ(i) = bessi0(betaZ*grid_arg)/bessi0(betaZ);
		}
		
	}break;
	
	case(SINC_KERNEL):{
		
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 2.5) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 2.5) : ( dwinY );
		dwinZ   = (dwinZ == -1 ) ? ( 2.5) : ( dwinZ );
	
		// Grid Length for precomputed kernel
		int grid_lengthX = (int)( (float)dwinX*(float)grid_modX);
		int grid_lengthY = (int)( (float)dwinY*(float)grid_modY);
		int grid_lengthZ = (int)( (float)dwinZ*(float)grid_modZ);
		
		// Alloc Lookup Table Structs for Gridding
		grid_filterX.resize( grid_lengthX+10);
		grid_filterY.resize( grid_lengthY+10);
		grid_filterZ.resize( grid_lengthZ+10);
    	grid_filterX = 0.0;
		grid_filterY = 0.0;
		grid_filterZ = 0.0;
				
		// Compute Seperable Kernels
		for(int i=0; i<(grid_lengthX+1); i++){
			float grid_pos=  ( (float)i * (float)dwinX / (float)grid_scale_x)/( (float)grid_lengthX);
			if(grid_pos < 0.01){
				grid_filterX(i) = 1.0;
			}else{
				grid_filterX(i)=(float)( sin(grid_pos*PI)/(grid_pos*PI));
			}
			grid_filterX(i) *=  ( 0.54  + 0.46*cos( PI*(float)i / ( (float)grid_lengthZ))); 
		}
		
		for(int i=0; i<(grid_lengthY+1); i++){
			float grid_pos=  ( (float)i * (float)dwinY/(float)grid_scale_y)/( (float)grid_lengthY);
			if(grid_pos < 0.01){
				grid_filterY(i) = 1.0;
			}else{
				grid_filterY(i)=(float)( sin(grid_pos*PI)/(grid_pos*PI));
			}
			grid_filterY(i) *=   ( 0.54  + 0.46*cos( PI*(float)i / ( (float)grid_lengthY)));
		}
	
		for(int i=0; i<(grid_lengthZ+1); i++){
			float grid_pos=  ( (float)i * (float)dwinZ / (float)grid_scale_z)/( (float)grid_lengthZ);
			if(grid_pos < 0.01){
				grid_filterZ(i) = 1.0;
			}else{
				grid_filterZ(i)=(float)( sin(grid_pos*PI)/(grid_pos*PI));
			}
			grid_filterZ(i) *= ( 0.54  + 0.46*cos( PI*(float)i / ( (float)grid_lengthZ))); 
		}
		
	}break;
	
  }
  
  // Normalize
  grid_filterX *= 0.5*grid_modX /  sum(grid_filterX);
  grid_filterY *= 0.5*grid_modY /  sum(grid_filterY);
  grid_filterZ *= 0.5*grid_modZ /  sum(grid_filterZ);

}


void gridFFT::precalc_gridding(int NzT,int NyT,int NxT, TrajDim trajectory_dims, TrajType trajectory_type ){
  
  Nx = NxT;
  Ny = NyT;
  Nz = NzT;
  
  
  // ---------------------------------------------
  // Determine what needs to be grid
  // ---------------------------------------------
  
  if( (trajectory_dims==TWOD) || (trajectory_type!= THREEDNONCARTESIAN) ){
  	if(grid_in_z ==-1){
		grid_in_z = 0;	
	}
  }
  if(grid_in_z == -1){
  	grid_in_z =1;
  }
  
  if(trajectory_type==CARTESIAN){
  	grid_in_y = 0;
  }
  
   // Get rounded Gridding ratio*
  if(grid_in_x ==1){
  	if(grid_x==-1){
		grid_x =  16.0*ceil( ( overgrid * (float)Nx )/16.0	) / (float)Nx;
  	}
  }else{
  	grid_x = 1;
  }
  
  if(grid_in_y ==1){
  	if(grid_y==-1){
		grid_y =  16.0*ceil( ( overgrid * (float)Ny )/16.0	) / (float)Ny;
	}
  }else{
  	grid_y = 1;
  }
  
  if(grid_in_z ==1){
  	if(grid_z==-1){
		grid_z =  16.0*ceil( ( overgrid * (float)Nz )/16.0	) / (float)Nz;
	}
  }else{
  	grid_z = 1;
  }
   // Compute Grid Size 
  Sz = (int)(grid_z *Nz);
  Sy = (int)(grid_y *Ny);
  Sx = (int)(grid_x *Nx);  
      
  precalc_kernel();
  
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
  printf("Grid in x=%d, y=%d, z=%d\n",grid_in_x,grid_in_y,grid_in_z);
  
  // Deapp Windows
  winx.resize(Sx);
  if( (fft_in_x==1) && (grid_in_x==1) ){
  winx = 0.0;
  for( int i = 0; i < Sx;i++){
  	float ipos = i - (float)Sx/2.0;
	for(int grid_pos = 0; grid_pos < dwinX*grid_modX; grid_pos++){ 
		// Fourier Transform of Kernel
		winx(i) += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modX / (float)Sx)*grid_filterX(grid_pos);
	}
	winx(i) = (float)grid_modX/winx(i);
	
	// Put chopping + zeroing in window to save time
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winx(i)*=fact / Sx;
	winx(i)*=( i < og_sx) ? ( 0.0 ) : ( 1.0);
	winx(i)*=( i > og_ex-1) ? ( 0.0 ) : ( 1.0);
	
	if(grid_in_x==0){
		winx(i)=1.0;
	}
  }
  winx /= max(winx);
  }else{
  	winx = 1.0;
  }
  
  winy.resize(Sy);
  if( (fft_in_y==1) && (grid_in_y==1) ){
  
  winy = 0.0;
  for( int i = 0; i < Sy;i++){
  	float ipos = i - (float)Sy/2.0;
	for(int grid_pos = 0; grid_pos < dwinY*grid_modY; grid_pos++){ 
		winy(i) += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modY / (float)Sy)*grid_filterY(grid_pos);
	}
	
	winy(i) = (float)grid_modY/winy(i);
	float fact =  ((float)( 2*(( i  )%2) - 1));
	winy(i)*=fact / Sy;
  	winy(i)*=( i < og_sy) ? ( 0.0 ) : ( 1.0);
	winy(i)*=( i > og_ey-1) ? ( 0.0 ) : ( 1.0);
	
	if(grid_in_y==0){
		winy(i)=1.0;
	}
  }
  winy /= max(winy);
  }else{
  	winy=1.0;
  }
   
  winz.resize(Sz);
  if( (fft_in_z==1) && (grid_in_z==1) ){
  	winz = 0.0;
  	for( int i = 0; i < Sz;i++){
  		float ipos = i - (float)Sz/2.0;
		for(int grid_pos = 0; grid_pos < dwinZ*grid_modZ; grid_pos++){ 
			winz(i) += 2*cos( 2*PI*ipos* grid_pos / (float)grid_modZ / (float)Sz)*grid_filterZ(grid_pos);
		}
		winz(i)  = (float)grid_modZ/winz(i);
		float fact =  ((float)( 2*(( i  )%2) - 1));
		winz(i)*=fact / Sz;
		winz(i)*=( i < og_sz) ? ( 0.0 ) : ( 1.0);
		winz(i)*=( i > og_ez-1) ? ( 0.0 ) : ( 1.0);
	
		if(grid_in_z==0){
		winz(i)=1.0;
		}	
  	}
  	winz /= max(winz);  
  }else{
  	winz = 1.0;
  }
  
  // Allocate Memory
  cout << "Alloc Grid" << endl;
  alloc_grid();
  
  // Setup 3D FFT
  if( (fft_in_z==1) && (fft_in_y==1) && (fft_in_x==1) ){	
  	plan_fft();
  }
  
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
	do_fft();
	if(time_grid)cout << ",fft:"<< T << endl << flush;
	
	if(time_grid) T.tic(); 
	forward_image_copy(X,smap); // Copy result to X
	if(time_grid)cout << ",copy:"<< T << endl<< flush;
}

void gridFFT::do_fft( void ){

	if( (fft_in_z==1) && (fft_in_y==1) && (fft_in_x==1) ){
		fftwf_execute(fft_plan); //FFT 
	}else{
		if(fft_in_x){
			fft3( k3d_grid,0,FFTW_FORWARD,0);
		}
		
		if(fft_in_y){
			fft3( k3d_grid,1,FFTW_FORWARD,0);
		}
		
		if(fft_in_z){
			fft3( k3d_grid,2,FFTW_FORWARD,0);
		}
	}
}

void gridFFT::do_ifft( void ){

	if( (fft_in_z==1) && (fft_in_y==1) && (fft_in_x==1) ){
		fftwf_execute(ifft_plan); //FFT 
	}else{
		if(fft_in_x){
			fft3( k3d_grid, 0,FFTW_BACKWARD,0);
		}
		
		if(fft_in_y){
			fft3( k3d_grid, 1,FFTW_BACKWARD,0);
		}
		
		if(fft_in_z){
			fft3( k3d_grid, 2,FFTW_BACKWARD,0);
		}
	}
}


// Float Sensitivty Map (intensity correction)
void gridFFT::forward( Array<complex<float>,3>&X,\
					   const Array<float,3>&smap,\
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
	do_fft();
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
	do_fft();
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
	do_ifft();
	if(time_grid)cout << ",ifft:"<< T;
	
	if(time_grid) T.tic(); 
	chop_grid_backward(data,kx,ky,kz,kw); // Inverse gridding
	if(time_grid)cout << ",igrid:"<< T << endl;
	
}

// Float sensitivty map (intensity correction)
void gridFFT::backward(const Array<complex<float>,3>&X,\
					   const Array<float,3>&smap,\
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
	do_ifft();
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
	do_ifft();
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
	  float wtz = winz(k+og_sz);
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy(j+og_sy);
		for(int i=0; i<Nx; i++) {
			X(i,j,k) = ( image(i,j,k)*wty*winx(i+og_sx)); /* Don't just sum coils when no sense map is given*/
	}}}
}

void gridFFT::forward_image_copy(Array<complex<float>,3>&X,const Array<complex<float>,3>&smap){
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz(k+og_sz);
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy(j+og_sy);
		for(int i=0; i<Nx; i++) {
			X(i,j,k) += ( image(i,j,k)*wty*winx(i+og_sx)*conj(smap(i,j,k)));
	}}}
}


void gridFFT::forward_image_copy(Array<complex<float>,3>&X,const Array< float,3>&smap){
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz(k+og_sz);
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy(j+og_sy);
		for(int i=0; i<Nx; i++) {
			X(i,j,k) += ( image(i,j,k)*wty*winx(i+og_sx)*smap(i,j,k));
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
	  float wtz = winz(k+og_sz);
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy(j+og_sy);
		for(int i=0; i<Nx; i++) {
			image(i,j,k) = X(i,j,k)*wty*winx(i+og_sx);
	}}}
}

void gridFFT::backward_image_copy(const Array<complex<float>,3>&X,const Array<complex<float>,3>&smap){
	
	if(overgrid > 1.0){
		k3d_grid = 0;
	}
	
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz(k+og_sz);
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy(j+og_sy);
		for(int i=0; i<Nx; i++) {
			image(i,j,k) =  X(i,j,k)*wty*winx(i+og_sx)*smap(i,j,k);
	}}}
}


void gridFFT::backward_image_copy(const Array<complex<float>,3>&X,const Array<float,3>&smap){
	
	if(overgrid > 1.0){
		k3d_grid = 0;
	}
	
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){ 
	  float wtz = winz(k+og_sz);
	  for(int j=0; j<Ny; j++){ 
	    float wty = wtz*winy(j+og_sy);
		for(int i=0; i<Nx; i++) {
			image(i,j,k) =  X(i,j,k)*wty*winx(i+og_sx)*smap(i,j,k);
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
	
	#pragma omp parallel for schedule(dynamic,1024)
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
		int sx;
		int ex;
		if(grid_in_x){
			sx = max( (int)ceil( dkx - dwinX),0);
			ex = min( (int)floor(dkx + dwinX),Sx-1);
		}else{
			sx = (int)( dkx);
			ex = sx;
		}
		if(sx >= Sx) continue;
		if(ex < 0) continue;  
				
		float dky = ky[i]*grid_y + cy;
		int sy;
		int ey;
		if(grid_in_y){
			sy = max( (int)ceil( dky - dwinY),0);
			ey = min( (int)floor(dky + dwinY),Sy-1);
		}else{
			sy = (int)(dky);
			ey = sy;
		}
		if(sy >= Sy) continue;
		if(ey < 0) continue;  
			
		
		float dkz = kz[i]*grid_z + cz;
		int sz;
		int ez;
		if(grid_in_z){
			sz = max( (int)ceil( dkz - dwinZ),0);
			ez = min( (int)floor(dkz + dwinZ),Sz-1);
		}else{
			sz = (int)(dkz);
			ez = sz;
		}
		if(sz >= Sz) continue;
		if(ez < 0) continue;  
			
		
		float kr=0.0;
		if( fft_in_x){
			kr+=(kx[i]*kx[i]);
		}
		
		if( fft_in_y){
			kr+=(ky[i]*ky[i]);
		}
		
		if( fft_in_z){
			kr+=(kz[i]*kz[i]);
		}
		temp *= exp( -kr / (2.0*k_rad*k_rad) );

		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    			float delz = fabs(grid_modZ*(dkz -(float)lz));
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
			if(!grid_in_z){
				wtz =1.0;
			}
			
			if(fft_in_z){
				wtz *=( (float)(2*(lz%2) -1 ));
			}
				
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly));
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy );
			 	if(fft_in_y){
					wty *=( (float)(2*(ly%2) -1 ));
				}
				 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx));
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX( (int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx );
			 		
					if(fft_in_x){
						wtx *=( (float)(2*(lx%2) -1 ));
					}
					
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
	
	int SizeX = imX.length(firstDim);
	int SizeY = imX.length(secondDim);
	int SizeZ = imX.length(thirdDim);
		
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
		if(ex >= SizeX) continue;  
	    
		// Compute Coordinates + Check
		float dky = ky[i]*grid_y + cy;
		int sy;
		int ey;
		if(grid_in_y==1){
			sy = (int)ceil( dky - dwinY);
			ey = (int)floor(dky + dwinY);
		}else{
			sy = (int)( dky);
			ey = sy;
		}
		if(sy <0)   continue;
		if(ey >= SizeY) continue;  
		
		
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
		if(ez >= SizeZ) continue;  
		
		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    		float delz = fabs(grid_modZ*(dkz -(float)lz));
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
			
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly));
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy );
			 	 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx));
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX((int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx );
			 		
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
	
	int SizeX = imX.length(firstDim);
	int SizeY = imX.length(secondDim);
	int SizeZ = imX.length(thirdDim);
		
	if( !dataA.isStorageContiguous()){
		cout << "Non-contiguous storage doesn't work yet" << endl;
		exit(1);
	}
	
	int Npts = dataA.numElements();
	const float *data = dataA.data();
	const float *kx   = kxA.data();
	const float *ky   = kyA.data();
	const float *kz   = kzA.data();
	
	/*
	cout << "Grid Size = " << SizeX << " x " << SizeY  << " x " << SizeZ << endl;
	cout << "Center X = " << cx << " x " << cy  << " x " << cz << endl;
	cout << "Dwin = " << dwinX << " x " << dwinY  << " x " << dwinZ << endl;
	cout << "Grid Length = " << grid_filterX.length(firstDim) << " x " << grid_filterY.length(firstDim)  << " x " << grid_filterZ.length(firstDim) << endl;
	cout << "Grid Mod  = " << grid_modX << " x " << grid_modY  << " x " << grid_modZ << endl;
	cout << "Grid In   = " << grid_in_x << " x " << grid_in_y  << " x " << grid_in_z << endl;
	*/
	
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
		if(ex >= SizeX) continue;  
	    
		// Compute Coordinates + Check
		float dky = ky[i]*grid_y + cy;
		int sy;
		int ey;
		if(grid_in_y==1){
			sy = (int)ceil( dky - dwinY);
			ey = (int)floor(dky + dwinY);
		}else{
			sy = (int)( dky);
			ey = sy;
		}
		if(sy <0)   continue;
		if(ey >= SizeY) continue;  
		
		
		
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
		if(ez >= SizeZ) continue;  
		
		
		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    		float delz = fabs(grid_modZ*(dkz -(float)lz));
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
			
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly));
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy );
			 	 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx));
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX((int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx );
			 		
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
	
	#pragma omp parallel for schedule(dynamic,1024)
	for (int i=0; i < Npts; i++) {
      	
		
		// Do not grid zeros
     	if( kw[i]==0.0) continue;
					
	    // Calculate the exact kspace sample point in 
	    // dimension flag->grid* kspace that this point (i,j)
	    // is contributing too.
	   		
		// Compute Coordinates + Check
		float dkx = kx[i]*grid_x + cx;
		int sx;
		int ex;
		if(grid_in_x){
			sx = max( (int)ceil( dkx - dwinX),0);
			ex = min( (int)floor(dkx + dwinX),Sx-1);
		}else{
			sx = (int)( dkx);
			ex = sx;
		}
		if(sx >= Sx) continue;
		if(ex < 0) continue;  
		
		
		float dky = ky[i]*grid_y + cy;
		int sy;
		int ey;
		if(grid_in_y){
			sy = max( (int)ceil( dky - dwinY),0);
			ey = min( (int)floor(dky + dwinY),Sy-1);
		}else{
			sy = (int)(dky);
			ey = sy;
		}
		if(sy >= Sy) continue;
		if(ey < 0) continue;  
			
		
		float dkz = kz[i]*grid_z + cz;
		int sz;
		int ez;
		if(grid_in_z){
			sz = max( (int)ceil( dkz - dwinZ),0);
			ez = min( (int)floor(dkz + dwinZ),Sz-1);
		}else{
			sz = (int)(dkz);
			ez = sz;
		}
		if(sz >= Sz) continue;
		if(ez < 0) continue;  

		complex<float>temp(0,0);
					
		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    		float delz = fabs(grid_modZ*(dkz -(float)lz));
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
			if(!grid_in_z){
				wtz =1.0;
			}
			
			if(fft_in_z){
				wtz *=( (float)(2*(lz%2) -1 ));
			}
			
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly));
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy );
			 	if(fft_in_y){
					wty *=( (float)(2*(ly%2) -1 ));
				}
								 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx));
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX((int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx );
			 		
					if(fft_in_x){
						wtx *=( (float)(2*(lx%2) -1 ));
					}
														
					/*This Memory Access is the Bottleneck*/	 			 
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







