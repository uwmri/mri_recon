#include "spirit.h"
#include "io_templates.cpp"

using arma::cx_fmat;
using arma::fvec;
using arma::uvec;
using namespace std;

SPIRIT::SPIRIT() 
{    
    // Kernel Size Defaults 
    kr_f = 0.0;
	krx_f = 3.0;  
	kry_f = 3.0;  
	krz_f = 3.0;  
	
	// Calibration Size Defaults 
	cr = 0;
	crx = 12;
	cry = 12;
	crz = 12;
	
	// Shape of Calibration Region
	shape = SP_CIRCLE;
	
	// Shrinkage operations for memory saving
	mapshrink = 2;
	mapthresh = 0.0;
	
	// Phase type
	phase_type = SP_COIL_PHASE;

	// Write out lots of images
	debug =0;
}


void SPIRIT::init(int xres, int yres, int zres, int nc){

    // Size of Final Images
    rcxres = xres;
	rcyres = yres;
	rczres = zres;
	ncoils = nc;
	
	// For Radial Kernel  
	if (kr_f > 0.0) {
	    krx_f = kr_f;  
		kry_f = kr_f;  
		krz_f = kr_f;  
	}
    krx = (int)ceil(krx_f);
	kry = (int)ceil(kry_f);
	krz = (int)ceil(krz_f);
	krz = (zres ==1) ? ( 0 ) : ( krz); // for 2D Kernel
	
	// For Calibration Size
	if (cr > 0) {
	    crx = cr;  
		cry = cr;  
		crz = cr;  
	}
	crz = (zres == 1) ? ( 0 ) : ( crz); // for 2D Kernel
	
    
}

// ----------------------
// Help Message
// ----------------------
void SPIRIT::help_message(void){
    cout << "----------------------------------------------" << endl;
	cout << "  Control for ESPIRIT " << endl;
	cout << "----------------------------------------------" << endl;
	
	help_flag("-sp_krx_f []","kernel size in x");
	help_flag("-sp_kry_f []","kernel size in y");
	help_flag("-sp_krz_f []","kernel size in z");
	help_flag("-sp_kr_f []","kernel size in r");
	
	help_flag("-sp_crx []","auto cal size in x");
	help_flag("-sp_cry []","auto cal size in y");
	help_flag("-sp_crz []","auto cal size in z");
	help_flag("-sp_cr []","auto cal size in r");
	
	help_flag("-sp_square","square cal shape");
	help_flag("-sp_circle","circle cal shape");
	
	help_flag("-sp_mapshrink []","shrink cal res for eigen decomp");
	help_flag("-sp_mapthresh []","threshold of eigen vals");
	help_flag("-sp_debug","write out intermediate images");
	help_flag("-sp_coil_phase","use phase from one coil");

}

void SPIRIT::read_commandline(int numarg, char **pstring){
    
#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
    
	for(int pos=0; pos < numarg; pos++){
	    
		if(strcmp("-h", pstring[pos] ) == 0) {
		    
			// Kernel Size
			float_flag("-sp_kr_f",kr_f);
			float_flag("-sp_krx_f",krx_f);
			float_flag("-sp_kry_f",kry_f);
			float_flag("-sp_krz_f",krz_f);

			// ACS Size
			int_flag("-sp_cr",cr);
			int_flag("-sp_crx",crx);
			int_flag("-sp_cry",cry);
			int_flag("-sp_crz",crz);
			
			// ACS Shape
			trig_flag(SP_SQUARE,"-sp_square",shape);
			trig_flag(SP_CIRCLE,"-sp_circle",shape);
			
			trig_flag(1,"-sp_debug",debug);		
			trig_flag(SP_COIL_PHASE,"-sp_coil_phase",phase_type);
		
		    int_flag("-sp_mapshrink",mapshrink);
		    float_flag("-sp_mapthresh",mapthresh);

		}
	}
}    

//-----------------------------------------------------
//  E-Spirit Call 
//-----------------------------------------------------
void SPIRIT::generateEigenCoils( Array< complex<float>,4 > &smaps){

    if(debug){
		ArrayWriteMag(smaps,"InputSmaps.dat");
    }

    // FFT Smaps Back to K-Space
    cout<< "FFT back to K-space" << endl << flush;
    for(int coil=0; coil< smaps.length(fourthDim); coil++){
		Array< complex<float>,3>SmapRef = smaps(Range::all(),Range::all(),Range::all(),coil);
		ifft(SmapRef); // In Place FFT
    }

    if(debug){
		ArrayWriteMag(smaps,"SmapKspace.dat");
    }

    // Array Reference for Low with Blitz		  
    cout<< "Calibrate" << endl << flush;
    calibrate_ellipsoid(smaps);

    // Code to get Sense maps from kernel		  
    prep(); // FFTs		  
    getcoils(smaps);

    // Phase Correction for Sense Maps
    phase_correct(smaps);

    // Interpolate back to high resolution size
    if(debug==1){
		ArrayWriteMag(smaps,"FinalSmaps.dat");
    }

}


//-----------------------------------------------------
//  E-Spirit Call 
//   Input:
//  	kdata -- k-space data in Cartesian Coordinates (Nx,Ny,Nz,Coils)
//   Description:
//		Does 
//-----------------------------------------------------
void SPIRIT::calibrate_ellipsoid(Array< complex<float>,4 > &kdata)
{
    cout << "Calibrating..." << endl << flush;
    
	// Centers of k-space
    int cx, cy, cz;
    cx = kdata.length(firstDim)/2;
    cy = kdata.length(secondDim)/2;
    cz = kdata.length(thirdDim)/2;
	cout << "Calibration Centered on " << cx << "," << cy << "," << cz << endl;
	
    // Subtract kernel borders from ACS size
    crx -= krx;
    cry -= kry;
    crz -= krz;

    int num_ACS = 0;
    int num_kernel = 0;

    if (shape==SP_CIRCLE) { // ellipse mode check
		// Figure out how many acs there are in the ellipse
		for (int z = -crz; z <= crz; z++) {
	    for (int y = -cry; y <= cry; y++) {
		for (int x = -crx; x <= crx; x++) {
		    
			float rad; 
			if(crz ==0){
				rad = (float)(x*x)/(crx*crx) + (float)(y*y)/(cry*cry);
		    }else{
				rad = (float)(x*x)/(crx*crx) + (float)(y*y)/(cry*cry) + (float)(z*z)/(crz*crz);
		    }
			if (rad <= 1) {num_ACS++;}      
		}}}

		// Figure out how many points are in the kernel
		for (int z = -krz; z <= krz; z++) {
	    for (int y = -kry; y <= kry; y++) {
		for (int x = -krx; x <= krx; x++) {
		    float rad; 
			if(krz_f ==0){
				rad = (float)(x*x)/(krx_f*krx_f) + (float)(y*y)/(kry_f*kry_f);
		    }else{
				rad = (float)(x*x)/(krx_f*krx_f) + (float)(y*y)/(kry_f*kry_f) + (float)(z*z)/(krz_f*krz_f);
		    }
			if (rad <= 1) {num_kernel++;}      
		}}}
		num_kernel *= ncoils;
	} else {
		num_ACS = (2*crx+1)*(2*cry+1)*(2*crz+1);
		num_kernel = (2*krx+1)*(2*kry+1)*(2*krz+1)*ncoils;
    }
	cout << "num_ACS: " << num_ACS << " num_kernel: " << num_kernel << endl;

    // Kernel Matrix
    cx_fmat A(num_ACS,num_kernel);
    	
    cout << "Populating Kernel: " << endl << flush;
    int Arow = 0;
	for (int iz = -crz; iz <= crz; iz++) { 
	for (int iy = -cry; iy <= cry; iy++) {
	for (int ix = -crx; ix <= crx; ix++) {

		    float rad;
			if( shape == SP_SQUARE){
				rad =0.0;
			}else if(crz==0){
				rad = (float)(ix*ix)/(crx*crx) + (float)(iy*iy)/(cry*cry);
			}else{
				rad = (float)(ix*ix)/(crx*crx) + (float)(iy*iy)/(cry*cry) + (float)(iz*iz)/(crz*crz);
			}
		    if (rad > 1)
			continue;

		    
		    // Get Neighborhood Values
		    int Acol = 0;
		    for (int jc = 0; jc < ncoils; jc++) {
			for (int jz = -krz; jz <= krz; jz++) {
			for (int jy = -kry; jy <= kry; jy++) {
			for (int jx = -krx; jx <= krx; jx++) {

				    float krad;
					if( shape == SP_SQUARE){
						krad = 0.0;
					}else if(krz_f == 0 ){
						krad =  (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f);
					}else{
						krad =  (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f) + (float)(jz*jz)/(krz_f*krz_f);
					}
				    if (krad <= 1) {
			
						// Get New Coord
						int kern_xx = jx+ix+cx;
						int kern_yy = jy+iy+cy;
						int kern_zz = jz+iz+cz;

						A(Arow,Acol) = (complex<float>)kdata(kern_xx,kern_yy,kern_zz,jc);
						Acol++;
					
				    }else{
						continue;
					}

			}}}} // Block of cal matrix
		    Arow++;

	} // ix
	} // iy
	} // iz

	// -------------------------------
	// Do SVD Decomposition
	// -------------------------------  
	
	cout << "A is " << A.n_rows << " x " << A.n_cols << endl;
		
	cout << "SVD on Kernel: "  << endl << flush;
  	cx_fmat U;
	fvec s;
	cx_fmat V;
	arma::svd_econ(U,s,V,A,'r');
	
	cout << "U is " << U.n_rows << " x " << U.n_cols << endl;
	cout << "S is " << s.n_rows << " x " << s.n_cols << endl;
	cout << "V is " << V.n_rows << " x " << V.n_cols << endl;
	
	V.save("V.txt",arma::raw_binary);
	s.save("S.txt",arma::raw_binary);
		
	cout << "Finding Number of Singular Values "  << endl << flush;
	float thresh = s(0)*sqrt(0.001);
	nV=1; 
	while( s(nV) > thresh ){
		nV++;	
	}
	nV--;
	cout << "Using " << nV << " singular values of kernel" << endl;
	
	cout << "Alloc Kernel: "  << endl << flush;
	
	// Blitz Alocation of Kernel (ColumnMajor Now!)
	k.setStorage( ColumnMajorArray<5>());
	k.resize(krx*2+1,kry*2+1,krz*2+1,ncoils,nV);
	k = complex<float>(0.0,0.0);
			
	// -------------------------------
	// Put the results into the right place in the kernel matrix
	// -------------------------------  

	for(int v=0; v<nV; v++){
	
		int xrow = 0;
		for (int jc = 0; jc < ncoils; jc++) {
	    for (int jz = -krz; jz <= krz; jz++) {
		for (int jy = -kry; jy <= kry; jy++) {
		for (int jx = -krx; jx <= krx; jx++) {
			float krad;
			if( shape == SP_SQUARE){
				krad = 0.0;
			}else if(krz_f == 0 ){
				krad =  (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f);
			}else{
				krad =  (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f) + (float)(jz*jz)/(krz_f*krz_f);
			}
			
			
			if (krad <= 1) {
			
			  // Get New Coord
			  int kern_xx = jx + krx;
			  int kern_yy = jy + kry;
			  int kern_zz = jz + krz;
			  
			  k(kern_xx,kern_yy,kern_zz,jc,v) = V(xrow,v);
			  xrow++;
			}else{
				continue;
			}
		}}}} // Block 	        
	} // nV 
    
	if(debug){
	    ArrayWriteMag(k,"KernelKspace.dat");
	}
	return;    
}

/**
 * FFT the k-space kernel into image space.
 */
void SPIRIT::prep() {
    
	cout << "FFT to SPIRiT Hybrid Space" << endl;
	
	// Hybrid Storage to save memory
	int sx = k.length(firstDim);
	int sy = k.length(secondDim);
	int sz = rczres;
	int nc = ncoils;
	
	// Blitz Alocation of Kernel (ColumnMajor Now!)
	im.setStorage( ColumnMajorArray<5>());
	im.resize(sx,sy,sz,nc,nV);
	im = 0;  
	
	cout << "Hybrid Kernel Size = " << sx << " x " << sy << " x " << sz << " x " << nc << " x " << nV << endl;
	
	/////////////
	// Copy k into the right place into im
	int z_off = (im.length(thirdDim)/2);
	
	Range all = Range::all();
	for (int v = 0; v < nV; v++) {
	    for (int coil = 0; coil < nc; coil++) {
		    Array< complex<float>,3>KernelCoil = k(all,all,all,coil,v);
		    Array< complex<float>,3>ImCoil = im(all,all,all,coil,v);
		    ImCoil(all,all,Range(-krz+z_off,krz+z_off)) = KernelCoil;
		    ifft(ImCoil,2); // Only FFT in third Dimension
	    }
	}
    
	if(debug){
	    ArrayWriteMag(im,"HybridKernelImage.dat");
	}
}


void SPIRIT::gs_orthogonalization(cx_fmat &m){
  
  for(int i=0; i< (int)m.n_rows; i++){
 	
	// Actual GS
	cx_fmat v = m.col(i);
	for(int j=0; j< i; j++){
		complex<float>scale(0,0);
		for(int k=0; k< (int)m.n_cols; k++){
			scale+= v(k)*conj(m(k,j));
		}
		v = v - scale*m.col(j);	
	}
	
	// Normalize
	float vsum =0;
	for(int k=0; k< (int)m.n_cols; k++){
		vsum += norm(v(k));
	}
	vsum = sqrt(vsum);
	
	// Put back
	for(int k=0; k< (int)m.n_cols; k++){
		m(k,i)=v(k)/vsum;
	}
  }
}

// Inplace orthognal iteration
void SPIRIT::orthogonal_iteration(cx_fmat &m){
	cx_fmat Q = arma::eye<cx_fmat>(m.n_rows,m.n_rows);
	cx_fmat A = Q;
	for(int iter=0; iter<30; iter++){
		A = Q;
		Q = m*A;	
		gs_orthogonalization(Q);
	}	
	m = Q;
}


/**
 * Generate coils from SPIRiT kernel.
 */
void SPIRIT::getcoils(Array< complex<float>,4 > &LR)
{
    int Nx = LR.length(firstDim);
    int Ny = LR.length(secondDim);
    int Nz = LR.length(thirdDim);
	
	Range all = Range::all();
    
	// Kernel Size
    int kNx = im.length(firstDim);
//    int kNy = im.length(secondDim);
//    int kNz = im.length(thirdDim);
    
	int cx = (Nx/2);
	int cy = (Ny/2);
	
	tictoc T;
	T.tic();
    cout << "Starting eig" << endl;

//    int count = 0;
	
	// Storage for Hybrid kx-y-coil-vector Kernel
	Array< complex<float>,4>temp;
	temp.setStorage( ColumnMajorArray<4>() );
	temp.resize(kNx,Ny,ncoils,nV);
	
	// Storage for x-coil-vector Kernel
	Array< complex<float>,3>x_line_full;
	x_line_full.setStorage( ColumnMajorArray<3>() );
	x_line_full.resize(Nx,ncoils,nV);
    
	for(int iz = 0; iz < Nz; iz ++) {

		cout << iz << " of " << Nz << endl << flush;
			
		// Get one Slice
		temp=complex<float>(0.0,0.0);
		Array< complex<float>,4>KernelTemp= im(all,all,iz,all,all);
		temp( all,Range(-kry+cy,kry+cy),all,all) = KernelTemp;
	
		// FFT to Full resolution in kx - y - z space
		for(int v=0; v<nV; v++){
			Array< complex<float>,3>fft_temp = temp(all,all,all,v);
			ifft(fft_temp,1);
		}
			
		for(int iy = 0; iy < Ny; iy ++) {
     	
		// Grab one row at a time
		Array< complex<float>,3>x_line = temp(all,iy,all,all); // fixed y and z
		x_line_full=complex<float>(0.0,0.0);
		x_line_full(Range(-krx+cx,krx+cx),all,all) = x_line;
		ifft(x_line_full,0);
		
		#pragma omp parallel for
	 	for(int ix = 0; ix < Nx; ix ++) {
		    
	 	   // Copy to A 
		   cx_fmat A( ncoils ,nV);
		   for(int v = 0; v < nV; v++) {
			for(int coil = 0; coil < ncoils; coil++) {
				A(coil,v) = x_line_full(ix,coil,v);
		   }}
		   
		   // Compute AtA
		   //cx_fmat AtA = A*A.t();
		   
		   /*
		   fvec eigval;
		   cx_fmat eigvec;
		   eig_sym(eigval, eigvec, AtA); 
		   
		     Copy Back
           for(int jj = 0; jj < ncoils; jj++) {
                  LR(ix,iy,iz,jj)=conj( eigvec(jj,ncoils-1));
           }
           LR(ix,iy,iz,ncoils)=eigval(ncoils-1);
		   */
		   
		   // SVD_econ is much faster than eig_sym and orthogonal iteration (but still slow). 
		   cx_fmat U;
		   fvec s;
		   cx_fmat V;
		   svd_econ(U, s, V, A, 'l');
		   for(int jj = 0; jj < ncoils; jj++) {
                  LR(ix,iy,iz,jj)=conj( U(jj,0));
           }
           //LR(ix,iy,iz,ncoils)=s(0);
		
			/*
		   cx_fmat AtA = A*A.t();
		   orthogonal_iteration(AtA);
		   for(int jj = 0; jj < ncoils; jj++) {
             LR(ix,iy,iz,jj)= conj( AtA(jj,0));
           }*/
		   
		}}
	  
	 }// Slice
	 cout<< "Done with eig: took " << T << endl;
	    
	 if(debug){
		ArrayWriteMag(LR,"CoilsEig.dat");
	 }
}


/** 
 * Rotates the spirit coils to match up with low res gridded data.  
 * The rotation is weighted by the norm of the data's magnitude.  This still ends 
 * up leaving some small discontinuities, but not many.
 */
void SPIRIT::phase_correct(Array< complex<float>,4 > &maps)
{
    
	if(phase_type==SP_SMOOTH){
	    
		// Try to find a map that provides smooth combination 
		
		Array< complex<float>,4> Sblur;
		Sblur.setStorage( ColumnMajorArray<4>());
		Sblur.resize(maps.length(firstDim),maps.length(secondDim),maps.length(thirdDim),maps.length(fourthDim));
		Sblur = complex<float>(1.0f,0.0);
		
		int Nc = maps.length(fourthDim)-1;
		
		if(debug==1){
		    FILE *fid = fopen("SpiritCoilPhase.dat","a");
			for(int coil=0; coil < Nc; coil++){
			    for(int j= 0; j<maps.extent(1); j++){
				for(int i=0; i<maps.extent(0);i++){
				    float b = arg(maps(i,j,maps.length(thirdDim)/2)*maps(i,j,maps.length(thirdDim)/2,coil));
					fwrite(&b,1,sizeof(float),fid);
				}}}
		    fclose(fid);
		}	
	    
		for(int iter=0; iter<1000; iter++){
		    
			float cx = maps.length(firstDim)/2;
			float cy = maps.length(secondDim)/2;
			float cz = maps.length(thirdDim)/2;
			float ksx = 2; //0.1.*(krx_f*krx_f);
		    float ksy = 2; //0.1.*(kry_f*kry_f);
		    float ksz = 2; //1.*(krz_f*krz_f);
		    
			// Get a blurred sensitivity map
			cout << "Blurring" << endl << flush;
			Sblur = maps;
			
			for( int coil=0; coil<maps.length(fourthDim); coil++){
			    Array<complex<float>,3> Sblur_ref = Sblur(Range::all(),Range::all(),Range::all(),coil);
				fft(Sblur_ref);
			}
		    
			for( int coil=0; coil<Sblur.length(fourthDim); coil++){
			    for (int k = 0; k < Sblur.length(thirdDim); k++) {
				for (int j = 0; j < Sblur.length(secondDim); j++) {
				    for (int i = 0; i < Sblur.length(firstDim); i++) {
					float x = (float)i - cx;
					    float y = (float)j - cy;
					    float z = (float)k - cz;
					    Sblur(i,j,k,coil) *= exp( -( x*x/ksx + y*y/ksy + z*z/ksz));
				    }}}				
			}
		    
			for( int coil=0; coil<maps.length(fourthDim); coil++){
			    Array<complex<float>,3> Sblur_ref = Sblur(Range::all(),Range::all(),Range::all(),coil);
				ifft(Sblur_ref);
			}		
		    
			// Now Phase 
			cout << "Aligning Iteration = " << iter << endl << flush;
			for (int k = 0; k < maps.length(thirdDim); k++) {
			    for (int j = 0; j < maps.length(secondDim); j++) {
				for (int i = 0; i < maps.length(firstDim); i++) {		
				    complex<float>p(0,0);
					for( int coil=0; coil<maps.length(fourthDim); coil++){
					    p+= maps(i,j,k,coil)*conj(Sblur(i,j,k,coil));
					}
				    p = polar(1.0f,arg(p));
					
					for( int coil=0; coil<maps.length(fourthDim); coil++){
					    maps(i,j,k,coil)*= conj(p);
					}
				}}}
		    
			
			if(debug==1){					
			    FILE *fid = fopen("SpiritCoilPhaseBlur.dat","a");
				for(int coil=0; coil < Nc; coil++){
				    for(int j= 0; j<maps.extent(1); j++){
					for(int i=0; i<maps.extent(0);i++){
					    float b = arg(Sblur(i,j,maps.length(thirdDim)/2,coil));
						fwrite(&b,1,sizeof(float),fid);
					}}}
			    fclose(fid);			// Export


			    fid = fopen("SpiritCoilPhase.dat","a");
			    for(int coil=0; coil < Nc; coil++){
				for(int j= 0; j<maps.extent(1); j++){
				    for(int i=0; i<maps.extent(0);i++){
					float b = arg(maps(i,j,maps.length(thirdDim)/2,coil));
					fwrite(&b,1,sizeof(float),fid);
				    }}}
			    fclose(fid);
			}
		}				
	}else if(phase_type == SP_COIL_PHASE){	
	 	cout << "Phase Correct with one coil" << endl;
	    #pragma omp parallel for
		for (int z = 0; z < maps.length(thirdDim); z++) {
		    for (int y = 0; y < maps.length(secondDim); y++) {
			for (int x = 0; x < maps.length(firstDim); x++) {
				
					// Pick one coil
				    complex<float>CI = polar(1.0f, -arg( maps(x,y,z,0)));
					
					// Get Normalization
					float mag=0.0;
					for (int coil = 0; coil < maps.length(fourthDim); coil++) {
					    mag += norm(maps(x,y,z,coil));
					}
				    CI /= sqrt(mag);
					
					// Apply
					for (int coil = 0; coil < maps.length(fourthDim); coil++) {
					    maps(x,y,z,coil)*=CI;
					}
				
			}}}
	    
	}
    
}

