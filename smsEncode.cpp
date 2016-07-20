/************************************************
Gridding and FFT Libraries for K-space to Image Domain Tranformation

Initial Author: 
	Kevin M. Johnson

Usage Example: 	

	This is a speciallized gridding routine especially for simultaneous multislice with sliding slices

	Creates block structured images

*************************************************/

#include "smsEncode.h"
#include "io_templates.hpp"
#include "tictoc.hpp"
using namespace NDarray;


//----------------------------------------
// Constructor - Sets Default Vals
//----------------------------------------
smsEncode::smsEncode(){
	grid_x= -1;
  	grid_y= -1;
  	overgrid = 1.5;
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
  	grid_in_x=1;
  	grid_in_y=1;
  	grid_in_z=1;
	grid_scale_x = 1.0;
	grid_scale_y = 1.0;
	grid_scale_z = 1.0;
	time_grid =0 ;
	double_grid =0;
	
	sms_factor =2;
	sms_slice_phase = 0.0;
	sms_fwhm = 2.0;
	
}

//----------------------------------------
// Allocate Memory for Gridding
//----------------------------------------

void smsEncode::alloc_grid(){
		
	// Allocate a 5D container to grid to
	cout << "Allocating Gridding Array: " << Sx << " x " << Sy << " x " << Sz << " :: " << (Nt+sms_factor-1) << " x " << Ne << endl;
	Array< Array< complex<float>, 3>,2> temp =  Alloc5DContainer< complex<float> >(Sx, Sy, Sz, Nt+sms_factor-1,Ne);
	k3d_grid.reference(temp);
}

// ----------------------
// Help Message
// ----------------------
void smsEncode::help_message(void){
	cout << "----------------------------------------------" << endl;
	cout << "   SMS Encode Control " << endl;
	cout << "----------------------------------------------" << endl;
	
	cout<<"Control" << endl;
	help_flag("-kaiser","use kaiser bessel kernel");
 	help_flag("-triangle","use triangle kernel");
 	help_flag("-overgrid []","overgrid by factor []");
 	help_flag("-fast_grid","no overgrid,traingle kernel");
	help_flag("-time_grid","output times for gridding");
			
	cout<<"Fine Control" << endl;
	help_flag("-grid_in_x []","0=use nearest neighbor interpolation in x");
	help_flag("-grid_in_y []","0=use nearest neighbor interpolation in y");
	help_flag("-grid_in_z []","0=use nearest neighbor interpolation in z");
	help_flag("-fft_in_x []","0=no fft in x");
	help_flag("-fft_in_y []","0=no fft in y");
	help_flag("-dwinX []","Size of kernel in x");
	help_flag("-dwinY []","Size of kernel in y");
	help_flag("-dwinZ []","Size of kernel in z");
	
}
	
 
//----------------------------------------
// Parse Command Line Args
//----------------------------------------


void smsEncode::read_commandline(int numarg, char **pstring){

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
		float_flag("-sms_fwhm",sms_fwhm);
		
		float_flag("-grid_scale_x",grid_scale_x);
		float_flag("-grid_scale_y",grid_scale_y);
		float_flag("-grid_scale_z",grid_scale_z);
		
		int_flag("-grid_in_x",grid_in_x);
		int_flag("-grid_in_y",grid_in_y);
		int_flag("-grid_in_z",grid_in_z);
		int_flag("-fft_in_x",fft_in_x);
		int_flag("-fft_in_y",fft_in_y);
		
		float_flag("-sms_slice_phase",sms_slice_phase);		
				
		trig_flag(KAISER_KERNEL,"-kaiser",kernel_type);
		trig_flag(TRIANGLE_KERNEL,"-triangle",kernel_type);
		trig_flag(SINC_KERNEL,"-sinc",kernel_type);
		trig_flag(POLY_KERNEL,"-poly_kernel",kernel_type);
		
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

void smsEncode::precalc_kernel(void){
   
  
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
		
		// Grid Length for precomputed kernel
		int grid_lengthX = (int)( (float)dwinX*(float)grid_modX);
		int grid_lengthY = (int)( (float)dwinY*(float)grid_modY);
		
		// Alloc Lookup Table Structs for Gridding
		grid_filterX.resize( grid_lengthX+10);
		grid_filterY.resize( grid_lengthY+10);
		grid_filterX = 0.0;
		grid_filterY = 0.0;
		
    
		// Compute Seperable Kernel
		for(int i=0; i<(grid_lengthX+1); i++){
			float grid_pos = (float)i / (float)grid_lengthX;
			grid_filterX(i)  = 1.0 - grid_pos;
	 	}
	
		for(int i=0; i<(grid_lengthY+1); i++){
			float grid_pos = (float)i / (float)grid_lengthY;
			grid_filterY(i)  = 1.0 - grid_pos;
	 	}
	
	}break;
	
	case(KAISER_KERNEL):{
		
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 2.5) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 2.5) : ( dwinY );
	
		// Grid Length for precomputed kernel
		int grid_lengthX = (int)( (float)dwinX*(float)grid_modX);
		int grid_lengthY = (int)( (float)dwinY*(float)grid_modY);

		// Alloc Lookup Table Structs for Gridding
		grid_filterX.resize( grid_lengthX+10);
		grid_filterY.resize( grid_lengthY+10);
    	grid_filterX = 0.0;
		grid_filterY = 0.0;

		// Get optimal Beta per Beatty et al
		float act_grid_x = grid_x*grid_scale_x;
		float act_grid_y = grid_y*grid_scale_y;
				
	 	betaX = PI*sqrtf(  (dwinX*dwinX)/(act_grid_x*act_grid_x)*(act_grid_x -0.5)*(act_grid_x-0.5) - 0.8);
	 	betaY = PI*sqrtf(  (dwinY*dwinY)/(act_grid_y*act_grid_y)*(act_grid_y -0.5)*(act_grid_y-0.5) - 0.8);
		
		float beta_minX = sqrt(pow(PI*2*dwinX/act_grid_x,2.0) -  9.6752);
		float beta_minY = sqrt(pow(PI*2*dwinY/act_grid_y,2.0) -  9.6752);
		betaX = ( beta_minX > betaX) ? ( beta_minX ) : ( betaX);
		betaY = ( beta_minY > betaY) ? ( beta_minY ) : ( betaY);
		
		
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
		
	}break;
	
	
	case(POLY_KERNEL):{
		
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 2) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 2) : ( dwinY );
			
		// Grid Length for precomputed kernel
		int grid_lengthX = (int)( (float)dwinX*(float)grid_modX);
		int grid_lengthY = (int)( (float)dwinY*(float)grid_modY);
				
		// Alloc Lookup Table Structs for Gridding
		grid_filterX.resize( grid_lengthX+10);
		grid_filterY.resize( grid_lengthY+10);
		
    	grid_filterX = 0.0;
		grid_filterY = 0.0;
		
		loadKernelTable( grid_filterX);
		loadKernelTable( grid_filterY);
		
	}break;
	
	
  }
  
  
  // Z has a distinct kernel
  dwinZ   = sms_fwhm;
  cout << "SMS FWHM = " << sms_fwhm << endl;
  float sigma = (float)sms_fwhm / 2.354820045030949;
  int grid_lengthZ = (int)( (float)dwinZ*(float)grid_modZ);
  grid_filterZ.resize( grid_lengthZ+10);
  grid_filterZ = 0.0;
  for(int i=0; i<(grid_lengthZ+10); i++){
  	float temp = (float)i/ grid_modZ;
	grid_filterZ(i) = exp(-temp*temp/ (2.*sigma*sigma));
  }
  
  // Normalize
  grid_filterX *= 0.5 /  ( sum(grid_filterX) / (float)grid_modX );
  grid_filterY *= 0.5 /  ( sum(grid_filterY) / (float)grid_modY );
  grid_filterZ *= 0.5 /  ( sum(grid_filterZ) / (float)grid_modZ );

  ArrayWrite(grid_filterZ,"GridFilterZ.dat");

}


void smsEncode::precalc_gridding(int NzT,int NyT,int NxT, int NtT, int NeT, float sms_factorT,TrajDim trajectory_dims, TrajType trajectory_type ){
  
  Nx = NxT;
  Ny = NyT;
  Nz = NzT;
  Ne = NeT;
  Nt = NtT;
  sms_factor = sms_factorT;
  
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
  	if(grid_x == -1){
		grid_x =   overgrid;
  	}
  }else{
  	grid_x = 1;
  }
  
  if(grid_in_y ==1){
  	if(grid_y==-1){
		grid_y =  overgrid;
	}
  }else{
  	grid_y = 1;
  }
  
   // Compute Grid Size 
  Sz = Nz;
  Sy = (int)(grid_y *Ny);
  Sx = (int)(grid_x *Nx);  
      
  precalc_kernel();
  
  // ------------------------------------------------
  //    Image Domain Calcs (crop + deapp)
  // ------------------------------------------------
  
  // Calculations to Determine Crop Positions
  og_sx =  (int)( (float)Nx*(grid_x - 1)/ 2.0);
  og_sy =  (int)( (float)Ny*(grid_y - 1)/ 2.0);
  og_ex =  og_sx + Nx;
  og_ey =  og_sy + Ny;
  
  printf("\n\nGridding Kernel Info\n");
  printf("Dwin 		%f %f %f\n",dwinX,dwinY,dwinZ); 
  printf("Mod 		%f %f %f\n",grid_modX,grid_modY,grid_modZ); 
  printf("Og %d-%d x %d-%d \n",og_sx,og_ex,og_sy,og_ey);  
  printf("Grid in x=%d, y=%d, z=%d \n",grid_in_x,grid_in_y,grid_in_z);
  printf("Size =%d x %d x %d x %d \n",Sx,Sy,Sz,Ne);
  
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
  	}
  	winy /= max(winy);
  }else{
  	winy=1.0;
  }

  // Allocate Memory
  cout << "Alloc Grid" << endl;
  alloc_grid();
  
}


void smsEncode::do_fft( void ){
	for(int t =0; t < k3d_grid.length(firstDim); t++){
		for(int e =0; e < k3d_grid.length(secondDim); e++){
	
			Array< complex<float>,3> temp = k3d_grid(t,e);
			if(fft_in_x){
				fft3(temp, 0, FFTW_FORWARD, 0);
			}
		
			if(fft_in_y){
				fft3(temp, 1, FFTW_FORWARD, 0);
			}
		}
	}
}

void smsEncode::do_ifft( void ){
	for(int t =0; t < k3d_grid.length(firstDim); t++){
		for(int e =0; e < k3d_grid.length(secondDim); e++){
			Array< complex<float>,3> temp = k3d_grid(t,e);
			if(fft_in_y){
				fft3(temp, 1, FFTW_BACKWARD, 0);
			}
			if(fft_in_x){
				fft3(temp, 0, FFTW_BACKWARD, 0);
			}
		
		}
	}
}


/**
 * Forward gridding 
 * @param X image to be accumulated
 * @param sensitivity map
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */					   
void smsEncode::forward( Array<Array<complex<float>,3>,2>&xdf,\
					    Array<complex<float>,3>&smap,\
					    Array< Array<complex<float>,3>,2>&data,\
					    Array< Array<float,3>,2>&kx,\
					    Array< Array<float,3>,2>&ky,\
					    Array< Array<float,3>,2>&kz,\
					    Array< Array<float,3>,2>&kw,\
						Array< Array<float,3>,3>&z){
	
	
	tictoc T;
	if(time_grid) T.tic(); 
	for(int t =0; t < k3d_grid.length(firstDim); t++){
		for(int e =0; e < k3d_grid.length(secondDim); e++){
			k3d_grid(t,e)=complex<float>(0.0,0.0); // Zero K-Space
	}}
	if(time_grid) cout << "Forward::zero:" << T << flush;
	
	if(time_grid) T.tic(); 
	chop_grid_forward(data,kx,ky,kz,kw,z); // Grid data to K-Space
	if(time_grid)cout << ",grid:"<< T << flush;
	
	if(time_grid) T.tic(); 
	do_fft();
	if(time_grid)cout << ",fft:"<< T <<  flush;
	
	if(time_grid) T.tic(); 
	accumulate(xdf,smap); // Deapp,multiply by sensitivity map, and copy
	if(time_grid)cout << ",accumulate:"<< T << endl<< flush;
	
}

/**
 * Backward gridding with sensitivity map
 * @param X image to be accumulated
 * @param sensitivity map
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */					   
void smsEncode::backward( Array<Array<complex<float>,3>,2>&X,\
					   Array<complex<float>,3>&smap,\
					   Array< Array<complex<float>,3>,2>&data,\
					   Array< Array<float,3>,2>&kx,\
					   Array< Array<float,3>,2>&ky,\
					   Array< Array<float,3>,2>&kz,\
					   Array< Array<float,3>,2>&kw,\
					   Array< Array<float,3>,3>&z){
	
	tictoc T;
	
	if(time_grid) T.tic(); 
	for(int t =0; t < k3d_grid.length(firstDim); t++){
		for(int e =0; e < k3d_grid.length(secondDim); e++){
			k3d_grid(e)=0; // Zero K-Space
	}}
	if(time_grid) cout << "Backward::zero:"<< T << flush;
		
	if(time_grid) T.tic(); 
	set_image(X,smap); // Copy image to gridding 
	if(time_grid) cout << "Backward::copy:"<< T << flush;
		
	if(time_grid) T.tic(); 
	do_ifft();
	if(time_grid)cout << ",ifft:"<< T << flush;
	
	if(time_grid) T.tic(); 
	chop_grid_backward(data,kx,ky,kz,kw,z);
	if(time_grid)cout << ",igrid:"<< T << endl << flush;
	
}


//----------------------------------------
//    Crop from Gridding Matrix to Image
//----------------------------------------
		
void smsEncode::accumulate( Array< Array< complex<float>, 3>,2>&X,  Array< complex<float>, 3>&smap ){
	
	for(int t=0; t < X.length(firstDim); t++){
		for(int e=0; e < X.length(secondDim); e++){
	
		#pragma omp parallel for
		for(int k=0; k< Nz; k++){ 
	  		for(int j=0; j<Ny; j++){ 
	  			float wty = winy(j+og_sy);
	  			for(int i=0; i<Nx; i++){
					float wt = wty*winx(i+og_sx);
					X(t,e)(i,j,k) += wt*k3d_grid(t,e)(i+og_sx,j+og_sy,k)*conj(smap(i,j,k));
		}}}
	}}
}


void smsEncode::set_image(  Array< Array< complex<float>, 3>,2>&X,  Array< complex<float>, 3>&smap){

	for(int t=0; t < X.length(firstDim); t++){
		for(int e=0; e < X.length(secondDim); e++){
	
		#pragma omp parallel for
		for(int k=0; k< Nz; k++){ 
			for(int j=0; j<Ny; j++){ 
	  			float wty = winy(j+og_sy);
				for(int i=0; i<Nx; i++){
					float wt = wty*winx(i+og_sx);
					k3d_grid(t,e)(i+og_sx,j+og_sy,k)= wt*X(t,e)(i,j,k)*smap(i,j,k);
		}}}
	}}
}



// -------------------------------------------------------
//  This is the main function for Gridding.  Assumes data
//  is already density compensated, etc.
// -------------------------------------------------------
void nested_workaround( long index, int *N,int *idx, int total){
	long tempi = index;
	for( int pos=0; pos < total; pos++){
		idx[pos] = tempi % N[pos];
		tempi = (tempi-idx[pos]) / N[pos];		
	}
}

void smsEncode::chop_grid_forward( Array< Array<complex<float>,3>,2>&data,\
					   Array< Array<float,3>,2>&kxA,\
					   Array< Array<float,3>,2>&kyA,\
					   Array< Array<float,3>,2>&kzA,\
					   Array< Array<float,3>,2>&kwA,\
					   Array< Array<float,3>,3>&zA){

	
	float cx = Sx/2;
	float cy = Sy/2;
	
	if(time_grid){
		cout << "Range Kx = " << min(kxA(0)) << " to " << max(kxA(0)) << endl;
		cout << "Range Ky = " << min(kyA(0)) << " to " << max(kyA(0)) << endl;
		cout << "Range Kz = " << min(kzA(0)) << " to " << max(kzA(0)) << endl;
		cout << "Range Z =  " << min(zA(0)) << " to " << max(zA(0)) << endl;
		cout << "Range Kw = " << min(kwA(0)) << " to " << max(kwA(0)) << endl;
		cout << "Max Kdata = " << max(abs(data(0))) << endl;
		cout << "Encodes = " << Ne << endl;
		cout << "Nt =      " << Nt << endl;
	}

	
	long total_images = Nt*Ne*sms_factor;
	int *N = new int[3];
	N[0] = Nt;
	N[1] = Ne;
	N[2] = sms_factor;
	
	#pragma omp parallel for
	for( int pos=0; pos < total_images; pos++){
		
		// Get the actual position
		int *I = new int[3];
		nested_workaround(pos,N,I,3);
		int t = I[0];
		int e = I[1];
		int sms_pos = I[2];
		delete [] I;
		
		float cz = Sz*0.5;
		
		for(int kk =0; kk < data(t,e).length(thirdDim); kk++){
		for(int jj =0; jj < data(t,e).length(secondDim); jj++){
		for(int ii =0; ii < data(t,e).length(firstDim); ii++){
		
			float kx = kxA(t,e)(ii,jj,kk);
			float ky = kyA(t,e)(ii,jj,kk);
			float kz = kzA(t,e)(ii,jj,kk);
			float zz = zA(t,e,sms_pos)(ii,jj,kk);
			float kw = kwA(t,e)(ii,jj,kk);
			complex<float>temp =data(t,e)(ii,jj,kk);
				
			// Density Comp
			temp *= kw;
			temp *= polar<float>(1.0,sms_slice_phase*sms_pos); 
			
			// Do not grid zeros
     		if( temp==complex<float>(0,0)) continue;
		
				
	 		// Calculate the exact kspace sample point in 
			// dimension flag->grid* kspace that this point (i,j)
			// is contributing too.
	   		
			// X position
			float dkx = kx*grid_x + cx;
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
				
			// Y position	
			float dky = ky*grid_y + cy;
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
			
			// Z position
			float dkz = zz + cz;
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
			
			
			// Now loop to put in the data
			for(int lz =sz; lz<=ez; lz++){
    			float delz = fabs(grid_modZ*(dkz -(float)lz));
				float dz = delz - (float)((int)delz);
				float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
				if(!grid_in_z){
					wtz =1.0;
				}
			
				// Combined Kz - z phase
				complex<float> Cwtz = polar<float>(wtz,-2.0*PI*zz*kz);
							
				for(int ly =sy; ly<=ey; ly++){
        		
					float dely = fabs(grid_modY*(dky -(float)ly));
					float dy = dely - (float)((int)dely);
					complex<float> wty = Cwtz*((float)( grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy ));
			 		if(fft_in_y){
						wty *=( (float)(2*(ly%2) -1 ));
					}
				 
					for(int lx =sx; lx<=ex; lx++){
			 			float delx = fabs(grid_modX*(dkx -(float)lx));
			 			float dx = delx - (float)((int)delx);
						complex<float> wtx = wty*((float)(  grid_filterX( (int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx ));
			 		
						if(fft_in_x){
							wtx *=( (float)(2*(lx%2) -1 ));
						}
					
						complex<float>temp2 = wtx*temp;
						float RD = real(temp2);
						float ID = imag(temp2);
						float *I = reinterpret_cast<float *>(&k3d_grid(t+sms_pos,e)(lx,ly,lz));
						float *R = I++;
					
						// Prevent Race conditions in multi-threaded
						#pragma omp atomic
						*R+=RD;
					
						#pragma omp atomic
						*I+=ID;
																	
					}/* end lz loop */
	  	 		}/* end ly */
		 	}/* end lx */
	
		}}} // Data Loop
	
	
	}/* Time , Encode, Sms Factor */
	
	delete [] N;
	return;
}




	
void smsEncode::chop_grid_backward( Array< Array<complex<float>,3>,2>&data,\
					   Array< Array<float,3>,2>&kxA,\
					   Array< Array<float,3>,2>&kyA,\
					   Array< Array<float,3>,2>&kzA,\
					   Array< Array<float,3>,2>&kwA,\
					   Array< Array<float,3>,3>&zA){


	float cx = Sx/2;
	float cy = Sy/2;
	
	long total_images = Nt*Ne*sms_factor;
	int *N = new int[3];
	N[0] = Nt;
	N[1] = Ne;
	N[2] = sms_factor;
	
	#pragma omp parallel for
	for( int pos=0; pos < total_images; pos++){
		
		// Get the actual position
		int *I = new int[3];
		nested_workaround(pos,N,I,3);
		int t = I[0];
		int e = I[1];
		int sms_pos = I[2];
		delete [] I;
		
		float cz = Sz*0.5;
				
		Array< complex<float>,3>dataA = data(t,e);
		
		for(int kk =0; kk < data(t,e).length(thirdDim); kk++){
		for(int jj =0; jj < data(t,e).length(secondDim); jj++){
		for(int ii =0; ii < data(t,e).length(firstDim); ii++){
		
			float kx = kxA(t,e)(ii,jj,kk);
			float ky = kyA(t,e)(ii,jj,kk);
			float kz = kzA(t,e)(ii,jj,kk);
			float zz = zA(t,e,sms_pos)(ii,jj,kk);
			float kw = kwA(t,e)(ii,jj,kk);
			
			// Do not grid zeros
     		if( kw==0) continue;
		
	 		// Calculate the exact kspace sample point in 
			// dimension flag->grid* kspace that this point (i,j)
			// is contributing too.
	   		
			// Compute Coordinates + Check
			float dkx = kx*grid_x + cx;
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
				
			float dky = ky*grid_y + cy;
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
			
		
			float dkz = zz + cz;
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
			
			complex<float>temp(0.0,0.0);
			/*This is the main loop - most time is spent here*/
			for(int lz =sz; lz<=ez; lz++){
    			
				float delz = fabs(grid_modZ*(dkz -(float)lz));
				float dz = delz - (float)((int)delz);
				float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
				if(!grid_in_z){
					wtz =1.0;
				}
				
				complex<float> Cwtz = polar<float>(wtz,2*PI*kz*zz);
			
				for(int ly =sy; ly<=ey; ly++){
        		
					float dely = fabs(grid_modY*(dky -(float)ly));
					float dy = dely - (float)((int)dely);
					complex<float> wty = Cwtz*((float)(  grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy ));
			 		if(fft_in_y){
						wty *=( (float)(2*(ly%2) -1 ));
					}
				 	
					for(int lx =sx; lx<=ex; lx++){
			 			float delx = fabs(grid_modX*(dkx -(float)lx));
			 			float dx = delx - (float)((int)delx);
						complex<float> wtx =wty*((float)(  grid_filterX( (int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx ));
			 		
						if(fft_in_x){
							wtx *=( (float)(2*(lx%2) -1 ));
						}
										
						temp += wtx*k3d_grid(t+sms_pos,e)(lx,ly,lz);
					
										
					}/* end lx loop */
	  	 		}/* end ly */
		 	}/* end lx */
		 	temp *= polar<float>(1.0,-sms_slice_phase*sms_pos); 
		 	dataA(ii,jj,kk) += temp;
		}}}/* end data loop */

		
	}// Nested loop

	delete [] N;
	
	return;
}	
	


// For Kaiser Bessel Window
float smsEncode::bessi0(float x)
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


/* The kernel's radius FOV product is the length,
 * in pixels, to the first truncation point.
 */

void  smsEncode::loadKernelTable(Array<float,1> & out)
{
	int len = out.numElements();
	double c0 = 1.;
  	double c1 = 0.04522831;
  	double c2 = -3.36020304;
  	double c3 = 1.12417012;
  	double c4 = 2.82448025;
  	double c5 = -1.63447764;

  	for (int i = 0; i < len; i++) {
    	double x = double(i)/double(len);
    	double x2 = x*x;
    	double x4 = x2*x2;
   	 	double x3 = x*x2;
    	double x5 = x*x4;
    	out(i) = (float)( c0 + c1*x + c2*x2 + c3*x3 + c4*x4 + c5*x5);
    }
	
    return;
} 





