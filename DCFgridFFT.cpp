/************************************************
Gridding and FFT Libraries for K-space to Image Domain Tranformation

Initial Author: 
	Kevin M. Johnson

Usage Example: 
	// Do Once
	DCFgridFFT gridding;
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(64,64,64,3,KAISER);
	
	// Use for transform
	gridding.forward()


*************************************************/

#include "DCFgridFFT.h"
#include "io_templates.hpp"
#include "tictoc.hpp"
using namespace NDarray;


//----------------------------------------
// Constructor - Sets Default Vals
//----------------------------------------
DCFgridFFT::DCFgridFFT(){
	grid_x= -1;
  	grid_y= -1;
  	grid_z= -1;
	overgrid = 1.375;
	kernel_type = POLY_KERNEL;
  
	dwinX=-1;
	dwinY=-1;
	dwinZ=-1;
	
  	grid_in_x=1;
  	grid_in_y=1;
  	grid_in_z=-1;
	grid_scale_x = 1.0;
	grid_scale_y = 1.0;
	grid_scale_z = 1.0;
	time_grid =0 ;
	
	acc = 4;

}

// ----------------------
// Help Message
// ----------------------
void DCFgridFFT::help_message(void){
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
	help_flag("-dwinX []","Size of kernel in x");
	help_flag("-dwinY []","Size of kernel in y");
	help_flag("-dwinZ []","Size of kernel in z");
}
	
 
//----------------------------------------
// Parse Command Line Args
//----------------------------------------

void DCFgridFFT::read_commandline(int numarg, char **pstring){

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

		trig_flag(POLY_KERNEL,"-poly_kernel",kernel_type);
		
		trig_flag(1,"-time_grid",time_grid);

	// Special Copies
	}
  }
}    



/* The kernel's radius FOV product is the length,
 * in pixels, to the first truncation point.
 */

void  DCFgridFFT::loadKernelTable(Array<float,1> & out)
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

//----------------------------------------
//    Setup for Gridding 
//----------------------------------------

void DCFgridFFT::precalc_kernel(void){
   
  
  // How many pts in kernel per delta k for kernel lookup table
  grid_modX = 600;
  grid_modY = 600;
  grid_modZ = 600;
 
  
  // ------------------------------------------------
  //    Kernel Calculations
  // ------------------------------------------------
  
  switch(kernel_type){
  		
	default:
	case(POLY_KERNEL):{
		
		// Kernel Half Size 
		dwinX   = (dwinX == -1 ) ? ( 2) : ( dwinX );
		dwinY   = (dwinY == -1 ) ? ( 2) : ( dwinY );
		dwinZ   = (dwinZ == -1 ) ? ( 2) : ( dwinZ );
	
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
		
		loadKernelTable( grid_filterX);
		loadKernelTable( grid_filterY);
		loadKernelTable( grid_filterZ);
		
	}break;
	
	
  }
  
  // Normalize
  grid_filterX *= 0.5*grid_modX /  sum(grid_filterX);
  grid_filterY *= 0.5*grid_modY /  sum(grid_filterY);
  grid_filterZ *= 0.5*grid_modZ /  sum(grid_filterZ);

}


// -------------------------------------------------------
//  This is the main function for Gridding.  Assumes data
//  is already density compensated, etc.
// -------------------------------------------------------

void DCFgridFFT::forward( Array< float,3> & image, const Array<float,3>&dataA, const Array<float,3>&kxA,const Array<float,3>&kyA,const Array<float,3>&kzA){
	int Sx = image.length(firstDim);
	int Sy = image.length(secondDim);
	int Sz = image.length(thirdDim);
	
	float cx = Sx/2;
	float cy = Sy/2;
	float cz = Sz/2;
	
	long Npts = dataA.numElements();
		
	//long stride_x = 1;
	long stride_y = dataA.length(firstDim);
	long stride_z = dataA.length(secondDim);
			
	#pragma omp parallel for
	for (long index=0; index < Npts; index++) {
      	
		
		// Nested parallelism workaround
		int ii =   index % stride_y;
		long tempi = (index-ii) / stride_y;
		int jj =   tempi % stride_z;
		int kk =  ( tempi - jj) / stride_z;
		
		float kx = kxA(ii,jj,kk);
		float ky = kyA(ii,jj,kk);
		float kz = kzA(ii,jj,kk);
		float temp =dataA(ii,jj,kk);
		
		float kr = sqrt(kx*kx + ky*ky + kz*kz);
		float scale = kr/128.0*(acc-1) + 1;
		
		
		// Do not grid zeros
     	if( temp==0.0) continue;
						
	 	// Calculate the exact kspace sample point in 
		// dimension flag->grid* kspace that this point (i,j)
		// is contributing too.
	   		
		// Compute Coordinates + Check
		float dkx = kx*grid_x + cx;
		int sx;
		int ex;
		if(grid_in_x){
			sx = max( (int)ceil( dkx - scale*dwinX),0);
			ex = min( (int)floor(dkx + scale*dwinX),Sx-1);
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
			sy = max( (int)ceil( dky - scale*dwinY),0);
			ey = min( (int)floor(dky + scale*dwinY),Sy-1);
		}else{
			sy = (int)(dky);
			ey = sy;
		}
		if(sy >= Sy) continue;
		if(ey < 0) continue;  
			
		
		float dkz = kz*grid_z + cz;
		int sz;
		int ez;
		if(grid_in_z){
			sz = max( (int)ceil( dkz - scale*dwinZ),0);
			ez = min( (int)floor(dkz + scale*dwinZ),Sz-1);
		}else{
			sz = (int)(dkz);
			ez = sz;
		}
		if(sz >= Sz) continue;
		if(ez < 0) continue;  
			
		
		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    		float delz = fabs(grid_modZ*(dkz -(float)lz))/scale;
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
			if(!grid_in_z){
				wtz =1.0;
			}
			
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly))/scale;
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy );
							 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx))/scale;
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX( (int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx );
			 		
					float temp2 = wtx*temp;
					
					#pragma omp atomic 
					image(lx,ly,lz) += temp2;
			
				}/* end lz loop */
	  	 	}/* end ly */
		 }/* end lx */
	}/* end data loop */
	

	return;
}



void DCFgridFFT::scale_kw(Array< float,3>&dataA, const Array<float,3>&kxA,const Array<float,3>&kyA,const Array<float,3>&kzA){

	// Scale Kw by the kernel size
	long Npts = dataA.numElements();
		
	//long stride_x = 1;
	long stride_y = dataA.length(firstDim);
	long stride_z = dataA.length(secondDim);
			
	for (long index=0; index < Npts; index++) {
      	
		// Nested parallelism workaround
		int ii =   index % stride_y;
		long tempi = (index-ii) / stride_y;
		int jj =   tempi % stride_z;
		int kk =  ( tempi - jj) / stride_z;
		
		float kx = kxA(ii,jj,kk);
		float ky = kyA(ii,jj,kk);
		float kz = kzA(ii,jj,kk);
		
		float kr = sqrt(kx*kx + ky*ky + kz*kz);
		float scale = kr/128.0*(acc-1) + 1;
		
		dataA(ii,jj,kk) *= pow(scale,3);
	}
		
}
	
void DCFgridFFT::backward(Array< float,3> & image,Array< float,3>&dataA, const Array<float,3>&kxA,const Array<float,3>&kyA,const Array<float,3>&kzA){
	int Sx = image.length(firstDim);
	int Sy = image.length(secondDim);
	int Sz = image.length(thirdDim);
	
	
	float cx = Sx/2;
	float cy = Sy/2;
	float cz = Sz/2;
	
	long Npts = dataA.numElements();
	//long stride_x = 1;
	long stride_y = dataA.length(firstDim);
	long stride_z = dataA.length(secondDim);
	
	#pragma omp parallel for 
	for (int index=0; index < Npts; index++) {
      	
		// Nested parallelism workaround
		int ii =  index % stride_y;
		long tempi = (index-ii) / stride_y;
		int jj =  tempi % stride_z;
		int kk =  ( tempi - jj) / stride_z;
				
		float kx = kxA(ii,jj,kk);
		float ky = kyA(ii,jj,kk);
		float kz = kzA(ii,jj,kk);
		
		float kr = sqrt(kx*kx + ky*ky + kz*kz);
		float scale = kr/128*(acc-1) + 1;
				
		
		dataA(ii,jj,kk) = 0.0;
		
					
	    // Calculate the exact kspace sample point in 
	    // dimension flag->grid* kspace that this point (i,j)
	    // is contributing too.
	   		
		// Compute Coordinates + Check
		float dkx = kx*grid_x + cx;
		int sx;
		int ex;
		if(grid_in_x){
			sx = max( (int)ceil( dkx - scale*dwinX),0);
			ex = min( (int)floor(dkx + scale*dwinX),Sx-1);
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
			sy = max( (int)ceil( dky - scale*dwinY),0);
			ey = min( (int)floor(dky + scale*dwinY),Sy-1);
		}else{
			sy = (int)(dky);
			ey = sy;
		}
		if(sy >= Sy) continue;
		if(ey < 0) continue;  
			
		
		float dkz = kz*grid_z + cz;
		int sz;
		int ez;
		if(grid_in_z){
			sz = max( (int)ceil( dkz - scale*dwinZ),0);
			ez = min( (int)floor(dkz + scale*dwinZ),Sz-1);
		}else{
			sz = (int)(dkz);
			ez = sz;
		}
		if(sz >= Sz) continue;
		if(ez < 0) continue;  

		float temp = 0.0;
					
		/*This is the main loop - most time is spent here*/
		for(int lz =sz; lz<=ez; lz++){
    		float delz = fabs(grid_modZ*(dkz -(float)lz))/scale;
			float dz = delz - (float)((int)delz);
			float wtz = grid_filterZ((int)delz)*( 1.0-dz) + grid_filterZ((int)delz +1)*dz;
			if(!grid_in_z){
				wtz =1.0;
			}
			
			for(int ly =sy; ly<=ey; ly++){
        		
				float dely = fabs(grid_modY*(dky -(float)ly))/scale;
				float dy = dely - (float)((int)dely);
				float wty =wtz*(  grid_filterY((int)dely)*( 1.0-dy) + grid_filterY((int)dely +1)*dy );
								 
				for(int lx =sx; lx<=ex; lx++){
			 		float delx = fabs(grid_modX*(dkx -(float)lx))/scale;
			 		float dx = delx - (float)((int)delx);
					float wtx =wty*(  grid_filterX((int)delx)*( 1.0-dx) + grid_filterX((int)delx +1)*dx );
			 		
					/*This Memory Access is the Bottleneck*/	 			 
			 		temp += wtx*image(lx,ly,lz);
			  	 
	    	}/* end lz loop */
	  	  }/* end ly */
		 }/* end lx */
		
		 dataA(ii,jj,kk) = temp;
	
	}/* end data loop */
	return;
}	
	







