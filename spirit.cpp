#include "spirit.h"
#include "io_templates.cpp"

using arma::cx_mat;
using arma::vec;
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
  crx = 14;
  cry = 14;
  crz = 14;

  // Shape of Calibration Region
  shape = SP_CIRCLE;
  calib_type = SP_TIK;
  calib_lam = .01;
  
  // Shrinkage operations for speed
  mapshrink = 2;
  mapthresh = 0.0;
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
  
  // For Calibration Size
  if (cr > 0) {
    crx = cr;  
    cry = cr;  
    crz = cr;  
  } 
  
  // Blitz Alocation of Kernel (ColumnMajor Now!)
  k.setStorage( ColumnMajorArray<5>());
  k.resize(krx*2+1,kry*2+1,krz*2+1,nc,nc);

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
	
	help_flag("-sp_tik","regularized kernel cal");
	help_flag("-sp_tsvd","truncated svd kernel cal");
	
	help_flag("-sp_mapshrink []","shrink cal res for eigen decomp");
	help_flag("-sp_mapthresh []","threshold of eigen vals");
}

void SPIRIT::read_commandline(int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if(strcmp("-h", pstring[pos] ) == 0) {
	  	
		float_flag("-sp_kr_f",kr_f);
		float_flag("-sp_krx_f",krx_f);
		float_flag("-sp_kry_f",kry_f);
		float_flag("-sp_krz_f",krz_f);
	  	
		int_flag("-sp_cr",cr);
		int_flag("-sp_crx",crx);
		int_flag("-sp_cry",cry);
		int_flag("-sp_crz",crz);
		
		trig_flag(SP_SQUARE,"-sp_square",shape);
		trig_flag(SP_CIRCLE,"-sp_circle",shape);
		
		trig_flag(SP_TIK,"-sp_tik",calib_type);
		trig_flag(SP_TSVD,"-sp_tsvd",calib_type);
		
		int_flag("-sp_mapshrink",mapshrink);
		float_flag("-sp_mapthresh",mapthresh);
		
	}
  }
}    

//-----------------------------------------------------
//  E-Spirit Call 
//-----------------------------------------------------
void SPIRIT::generateEigenCoils( Array< complex<float>,4 > &smaps){
		  
		   // FFT Smaps Back to K-Space
		  for(int coil=0; coil< smaps.length(fourthDim); coil++){
		  	Array< complex<float>,3>SmapRef = smaps(Range::all(),Range::all(),Range::all(),coil);
			ifft(SmapRef); // In Place FFT
		  }
		   
		  // Array Reference for Low with Blitz		  
		  calibrate_ellipsoid(smaps);

		  // Code to get Sense maps from kernel		  
		  prep();		  
		  Array< complex<float>,4>LRmaps(rcxres/mapshrink,rcyres/mapshrink,rczres/mapshrink,smaps.length(fourthDim)+1);
		  getcoils(LRmaps);

		  // Phase Correction for Sense Maps
		  for(int coil=0; coil< smaps.length(fourthDim); coil++){
		  	Array< complex<float>,3>SmapRef = smaps(Range::all(),Range::all(),Range::all(),coil);
			ifft(SmapRef); // In Place FFT
		  }
		  Array< complex<float>,4>PhaseRef = smaps(Range(fromStart,toEnd,mapshrink),Range(fromStart,toEnd,mapshrink),Range(fromStart,toEnd,mapshrink),Range::all());
		  rotateCoils(LRmaps,PhaseRef);
		  
		  // Interpolate back to high resolution size
		  interpMaps(LRmaps,smaps);
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
  cout << "Calibrating..." << endl;
  // Centers of k-space
  int cx, cy, cz;
  cx = kdata.length(firstDim)/2;
  cy = kdata.length(secondDim)/2;
  cz = kdata.length(thirdDim)/2;
  
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
          float rad = (float)(x*x)/(crx*crx) + (float)(y*y)/(cry*cry) + (float)(z*z)/(crz*crz);
          if (rad <= 1) {num_ACS++;}      
    }}}
    
    // Figure out how many points are in the kernel
    for (int z = -krz; z <= krz; z++) {
      for (int y = -kry; y <= kry; y++) {
        for (int x = -krx; x <= krx; x++) {
          float rad = (float)(x*x)/(krx_f*krx_f) + (float)(y*y)/(kry_f*kry_f) + (float)(z*z)/(krz_f*krz_f);
          if (rad <= 1) {num_kernel++;}      
    }}}
    num_kernel *= ncoils;
    num_kernel -= 1; // Identity point
  } else {
    num_ACS = (2*crx+1)*(2*cry+1)*(2*crz+1);
    num_kernel = (2*krx+1)*(2*kry+1)*(2*krz+1)*ncoils-1;
  }
  
  cout << "num_ACS: " << num_ACS << " num_kernel: " << num_kernel << endl;
  
  // Allocate a bunch of matrices for openMP threads
  //int nthreads = omp_get_max_threads();
  cx_mat* A = new cx_mat[ncoils];
  cx_mat* b = new cx_mat[ncoils];
  cx_mat* Ap = new cx_mat[ncoils];
  cx_mat* ApA = new cx_mat[ncoils];
  cx_mat* S = new cx_mat[ncoils];
  cx_mat* x = new cx_mat[ncoils];
  
  cx_mat* U = new cx_mat[ncoils];
  vec* s = new vec[ncoils];
  cx_mat* V = new cx_mat[ncoils];
  
  #pragma omp parallel for
  for (int ic = 0; ic < ncoils; ic++) {
    int Arow,Acol;
    //int tid = omp_get_thread_num();
    int tid = ic;
    A[tid].zeros(num_ACS,num_kernel);
    b[tid].zeros(num_ACS,1);
    
    Arow = 0;
    for (int iz = -crz; iz <= crz; iz++) { 
      for (int iy = -cry; iy <= cry; iy++) {
        for (int ix = -crx; ix <= crx; ix++) {
          
          float rad = (shape == SP_SQUARE) ? 0 : (float)(ix*ix)/(crx*crx) + (float)(iy*iy)/(cry*cry) + (float)(iz*iz)/(crz*crz);
          if (rad <= 1) {
          // Old Index :: ki = (cx+ix) + kdata.dim[0]*(cy+iy) + kdata.dim[0]*kdata.dim[1]*(cz+iz) + kdata.dim[0]*kdata.dim[1]*kdata.dim[2]*ic; // Index in kdata
          b[tid](Arow,0) = (complex<double>)kdata(ix+cx,iy+cy,iz+cz,ic);
		  
  		  // Get Neighborhood
          Acol = 0;
          for (int jc = 0; jc < ncoils; jc++) {
            for (int jz = -krz; jz <= krz; jz++) {
              for (int jy = -kry; jy <= kry; jy++) {
                for (int jx = -krx; jx <= krx; jx++) {
                  
                  float krad = (shape == SP_SQUARE) ? 0 : (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f) + (float)(jz*jz)/(krz_f*krz_f);
                  if (krad <= 1) {
                  // Get New Coord
				  
				  int kern_xx = cx+jx+ix;
				  int kern_yy = cy+jy+iy;
				  int kern_zz = cz+jz+iz;
				  
				  if( (jc==ic) && (jz==0) && (jy==0) && (jx==0)){
				  }else{				  				  
                    A[tid](Arow,Acol) = (complex<double>)kdata(kern_xx,kern_yy,kern_zz,jc);
                    Acol++;
                  }
                  } else {continue;}
                  
                }
              }
            }
          } // jc
          Arow++;
          }else{continue;}
          
        } // ix
      }// iy
    } // iz
    
	// -------------------------------
	// Do the spirit calibration
    // -------------------------------  
	if (calib_type==SP_TIK) {
      Ap[tid] = arma::trans(A[tid]);
      ApA[tid] = Ap[tid] * A[tid];
  
      float beta = arma::norm(ApA[tid],"fro")*(calib_lam);
      S[tid].eye(ApA[tid].n_rows,ApA[tid].n_cols);
      S[tid] *= beta;
      
      ApA[tid] += S[tid];
      
      x[tid] = arma::solve(ApA[tid],Ap[tid]) * b[tid];
    } else if (calib_type==SP_TSVD) {
      arma::svd_econ(U[tid],s[tid],V[tid],A[tid]);
      uint tsvd_i = 0;
      double tsvd_max = s[tid].max();
      while( s[tid](tsvd_i)>calib_lam*tsvd_max && tsvd_i < A[tid].n_cols-1) { tsvd_i++; }
      cout << "Last SV index = " << tsvd_i << " out of " << s[tid].n_elem << endl;
      for(uint ti = 0; ti < A[tid].n_cols; ti++) {
        if (ti<tsvd_i) {s[tid](ti) = 1.0/s[tid](ti);}
        else {s[tid](ti) = 0;}
      }
      x[tid] = V[tid]*arma::diagmat(s[tid])*U[tid].t()*b[tid];
    }
    
    
	// -------------------------------
	// DPut the results of x into the right place in the kernel matrix
    // -------------------------------  
	
    // Kernel centers
    int ckx = k.length(firstDim)/2;
    int cky = k.length(secondDim)/2;
    int ckz = k.length(thirdDim)/2;
	
	// Copy to Kernel Matrix
    int xrow = 0;
    for (int jc = 0; jc < ncoils; jc++) {
      for (int jz = -krz; jz <= krz; jz++) {
        for (int jy = -kry; jy <= kry; jy++) {
          for (int jx = -krx; jx <= krx; jx++) {
            
            float grad = (shape == SP_SQUARE) ? 0 : (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f) + (float)(jz*jz)/(krz_f*krz_f);
            if (grad <= 1) {
            
			// Reverse x and y and z so it works right
            int kern_xx = ckx - jx;
			int kern_yy = cky - jy;
			int kern_zz = ckz - jz;
						
			if( (jc==ic) && (jz==0) && (jy==0) && (jx==0)){
			}else{				  				  
               k(kern_xx,kern_yy,kern_zz,jc,ic) = x[tid](xrow,0);
               xrow++;
            }        
            } else {continue;}
            
          }
        }
      }
    } // jc
          
  } // ic
  
    delete[] A;
    delete[] b;
    delete[] Ap;
    delete[] ApA;
    delete[] S;
    delete[] x;
    delete[] U;
    delete[] s;
    delete[] V;
  
  return;    
}

/**
* FFT the k-space kernel into image space.
*/
void SPIRIT::prep() {

  cout << "FFT to SPIRiT Image space" << endl;
  
  int sx = rcxres/mapshrink;
  int sy = rcyres/mapshrink;
  int sz = rczres/mapshrink;
  int nc = ncoils;
  
  // Blitz Alocation of Kernel (ColumnMajor Now!)
  im.setStorage( ColumnMajorArray<5>());
  im.resize(sx,sy,sz,nc,nc);
    
  
  /////////////
  // Copy k into the right place into im
  int x_off = (im.length(firstDim)/2);
  int y_off = (im.length(secondDim)/2);
  int z_off = (im.length(thirdDim)/2);
  
  int ksx = ((k.length(firstDim)-1)/2);
  int ksy = ((k.length(secondDim)-1)/2);
  int ksz = ((k.length(secondDim)-1)/2);
  
  int neg_mod = 0 - (im.length(firstDim)/2+im.length(secondDim)/2+im.length(thirdDim)/2)%2;
  if (neg_mod==0) {neg_mod = 1;}
  
  Range all = Range::all();
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nc; i++) {
      
	  Array< complex<float>,3>KernelCoil = k(all,all,all,i,j);
	  Array< complex<float>,3>ImCoil = im(all,all,all,i,j);
	  
	  for (int iz = -ksz; iz <= ksz; iz++) {
        for (int iy = -ksy; iy <= ksy; iy++) { 
          for (int ix = -ksx; ix <= ksx; ix++) {   
            ImCoil(ix+x_off,iy+y_off,iz+z_off) = (float)neg_mod*KernelCoil(ix+ksx,iy+ksy,iz+ksz);
          }
        }
      }
      
	  fft(ImCoil);
	  
    }
  }
}


// using Eigen::MatrixXcd;

/**
* Generate coils from SPIRiT kernel.
*/
void SPIRIT::getcoils(Array< complex<float>,4 > &LR)
{
  int ncoils = im.length(fourthDim);
  
  tictoc T;
  T.tic();
  cout << "Starting eig" << endl;
  
  #pragma omp parallel for
  for(int iz = 0; iz < im.extent(thirdDim); iz ++) {
    
	// Arrays for storage
  	cx_mat A;
	A.zeros(ncoils,ncoils);
  	cx_mat AtA;
	AtA.zeros(ncoils,ncoils);

	for(int iy = 0; iy < im.extent(secondDim); iy ++) {
      for(int ix = 0; ix < im.extent(firstDim); ix ++) {
  		
		// Copy to A		
		for(int jj = 0; jj < ncoils; jj++) {
          for(int ii = 0; ii < ncoils; ii++) {
            A(jj,ii) = (complex<double>)im(ix,iy,iz,ii,jj);
			AtA(ii,jj) = A(jj,ii);
          }
        }
		
		// Compute AtA
		A *= AtA;
  		
		// Eigen Vector
		vec eigval;
		cx_mat eigvec;
		eig_sym(eigval,eigvec,A);
		
		// Copy Back
		for(int jj = 0; jj < ncoils; jj++) {
			LR(ix,iy,iz,jj)=eigvec(jj,ncoils-1);
		}
		LR(ix,iy,iz,ncoils)=eigval(ncoils-1);
  }}}
  cout<< "Done with eig: took " << T << endl;
  // ArrayWriteMag(LR,"CoilsEig.dat");

}

void SPIRIT::interpMaps(Array< complex<float>,4 > &LR, Array< complex<float>,4 > &out) {
 
  int ncoils = out.length(fourthDim);
  
  Array< complex<float>,3>EigVal=LR(Range::all(),Range::all(),Range::all(),ncoils);
  mapthresh *= max(abs(EigVal));
  cout << "Interpolating espirit maps" << endl;
  
  // Linear interpolation up to the full resolution
  for(int i = 0; i < ncoils; i++) {
    Array<complex<float>,3>CoilRef =LR(Range::all(),Range::all(),Range::all(),i); 
		  
	#pragma omp parallel for
    for (int z = 0; z < out.length(thirdDim); z++) {
      int zz = z/mapshrink;
	  int zp = (zz + 1) % EigVal.length(thirdDim);
	  float dz = (float)(z%mapshrink)/(float)mapshrink;
		  
	  for (int y = 0; y < out.length(secondDim); y++) {
        int yy = y/mapshrink;
		int yp = (yy + 1) % EigVal.length(secondDim);
		float dy = (float)(y%mapshrink)/(float)mapshrink;
		  
		for (int x = 0; x < out.length(thirdDim); x++) {
          int xx = x/mapshrink;
		  int xp = (xx + 1) % EigVal.length(thirdDim);
		  float dx = (float)(x%mapshrink)/(float)mapshrink;
		    
		  
		  float svd_interp = abs(\
		   					     (1-dz)*(1-dy)*(1-dx)*EigVal(xx,yy,zz)+ \
                            	 (1-dz)*(1-dy)*(dx)*EigVal(xp,yy,zz)+ \
								 (1-dz)*(dy)*(1-dx)*EigVal(xx,yp,zz)+ \
                            	 (1-dz)*(dy)*(dx)*EigVal(xp,yp,zz)+ \
								 (dz)*(1-dy)*(1-dx)*EigVal(xx,yy,zp)+ \
                            	 (dz)*(1-dy)*(dx)*EigVal(xp,yy,zp)+ \
								 (dz)*(dy)*(1-dx)*EigVal(xx,yp,zp)+ \
                            	 (dz)*(dy)*(dx)*EigVal(xp,yp,zp));
								           
          if (svd_interp < mapthresh) {
           out(x,y,z,i) = 0;
          } else {
          			        
          out(x,y,z,i)  =  	     (1-dz)*(1-dy)*(1-dx)*CoilRef(xx,yy,zz)+ \
                            	 (1-dz)*(1-dy)*(dx)*CoilRef(xp,yy,zz)+ \
								 (1-dz)*(dy)*(1-dx)*CoilRef(xx,yp,zz)+ \
                            	 (1-dz)*(dy)*(dx)*CoilRef(xp,yp,zz)+ \
								 (dz)*(1-dy)*(1-dx)*CoilRef(xx,yy,zp)+ \
                            	 (dz)*(1-dy)*(dx)*CoilRef(xp,yy,zp)+ \
								 (dz)*(dy)*(1-dx)*CoilRef(xx,yp,zp)+ \
                            	 (dz)*(dy)*(dx)*CoilRef(xp,yp,zp);
           }
        }
      }
    }
    
    
  }
}

/** 
* Rotates the spirit coils to match up with low res gridded data.  
* The rotation is weighted by the norm of the data's magnitude.  This still ends 
* up leaving some small discontinuities, but not many.
*/
void SPIRIT::rotateCoils(Array< complex<float>,4 > &maps, Array< complex<float>,4 > &ref)
{

  #pragma omp parallel for
  for (int z = 0; z < ref.length(thirdDim); z++) {
    for (int y = 0; y < ref.length(secondDim); y++) {
      for (int x = 0; x < ref.length(firstDim); x++) {
        
		// Solve for scale factor + Normalization
		complex<float>CI(0,0);
		float mag=0.0;
		for (int coil = 0; coil < ref.length(fourthDim); coil++) {
			CI += ref(x,y,z,coil)*conj(maps(x,y,z,coil));
			mag += norm(maps(x,y,z,coil));
		}
		CI/=abs(CI);
		CI/=mag;
		
		// Scale
		for (int coil = 0; coil < ref.length(fourthDim); coil++) {
			maps(x,y,z,coil)*=CI;
		}
		
		  
		
		/*
		float diff = 0.0;
        float rot = 0.0;
        float scale = 0.0;
        for (int coil = 0; coil < ref.length(fourthDim); coil++) {
          diff = arg(ref(x,y,z,coil)*conj(maps(x,y,z,coil)));  // Switched to conj * (vs subtraction)
          if (coil > 0) { // Unwrap based on current guess of phase
           diff -= 2.0*PI*round((diff-rot/scale)/2.0/PI); 
          }
          rot += diff*norm(ref(x,y,z,coil));
          scale += norm(ref(x,y,z,coil));
        }
        rot /= scale;
        for (int coil = 0; coil < ref.length(fourthDim); coil++) {
          maps(x,y,z,coil) *= complex<float>(cos(rot),sin(rot)); 
        } 
		*/    
  }}}

}

