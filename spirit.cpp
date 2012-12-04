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
  
  // Phase type
  phase_type = SP_COIL_PHASE;
  
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
  
  // For Calibration Size
  if (cr > 0) {
    crx = cr;  
    cry = cr;  
    crz = cr;  
  } 
  
  // Blitz Alocation of Kernel (ColumnMajor Now!)
  k.setStorage( ColumnMajorArray<5>());
  k.resize(krx*2+1,kry*2+1,krz*2+1,nc,nc);
  k = 0;
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
	help_flag("-sp_debug","write out intermediate images");
	help_flag("-sp_coil_phase","use phase from one coil");
	
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

		trig_flag(1,"-sp_debug",debug);		
		trig_flag(SP_TIK,"-sp_tik",calib_type);
		trig_flag(SP_TSVD,"-sp_tsvd",calib_type);
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
		  prep();		  
		  Array< complex<float>,4>LRmaps(rcxres/mapshrink,rcyres/mapshrink,rczres/mapshrink,smaps.length(fourthDim)+1);
		  getcoils(LRmaps);

		  // Phase Correction for Sense Maps
		  for(int coil=0; coil< smaps.length(fourthDim); coil++){
		  	Array< complex<float>,3>SmapRef = smaps(Range::all(),Range::all(),Range::all(),coil);
			fft(SmapRef); // In Place FFT
		  }
		  Array< complex<float>,4>PhaseRef = smaps(Range(fromStart,toEnd,mapshrink),Range(fromStart,toEnd,mapshrink),Range(fromStart,toEnd,mapshrink),Range::all());
		  rotateCoils(LRmaps,PhaseRef);
		  
		  // Interpolate back to high resolution size
		  smaps=0;
		  interpMaps(LRmaps,smaps);
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
  
  // Kernel Matrix
  cx_mat A(num_ACS,num_kernel);
  cx_mat b(num_ACS,1);
  cx_mat x;
  
  cout << "Solving SPIRIT: ";
  for (int ic = 0; ic < ncoils; ic++) {
    
	cout << ic << " " << flush; 
	int Arow = 0;
    for (int iz = -crz; iz <= crz; iz++) { 
      for (int iy = -cry; iy <= cry; iy++) {
        for (int ix = -crx; ix <= crx; ix++) {
          
          float rad = (shape == SP_SQUARE) ? 0 : (float)(ix*ix)/(crx*crx) + (float)(iy*iy)/(cry*cry) + (float)(iz*iz)/(crz*crz);
          if (rad > 1)
		  	continue;
		   
		  // Data Point
		  b(Arow,0) = (complex<double>)kdata(ix+cx,iy+cy,iz+cz,ic);
		  
  		  // Get Neighborhood Values
          int Acol = 0;
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
                    A(Arow,Acol) = (complex<double>)kdata(kern_xx,kern_yy,kern_zz,jc);
                    Acol++;
                  }
                  } else {continue;}
                  
                }
              }
            }
          } // jc
          Arow++;
                   
        } // ix
      }// iy
    } // iz
    
	// -------------------------------
	// Do the spirit calibration
    // -------------------------------  
	
	if (calib_type==SP_TIK) {
      cx_mat AhA = arma::trans(A)*A;
      cx_mat AhB = arma::trans(A)*b;
	    
      float beta = arma::norm(AhA,"fro")*(calib_lam);
      for(unsigned int i=0; i< AhA.n_rows;i++){
	  	 AhA(i,i) += beta;
      }
      x = arma::solve(AhA,AhB);
    }else if (calib_type==SP_TSVD) {
	  cx_mat U;
	  vec s;
	  cx_mat V;
	  
      arma::svd_econ(U,s,V,A);
      uint tsvd_i = 0;
      double tsvd_max = s.max();
      while( s(tsvd_i)>calib_lam*tsvd_max && tsvd_i < A.n_cols-1) { tsvd_i++; }
      cout << "Last SV index = " << tsvd_i << " out of " << s.n_elem << endl;
      for(uint ti = 0; ti < A.n_cols; ti++) {
        if (ti<tsvd_i) {s(ti) = 1.0/s(ti);}
        else {s(ti) = 0;}
      }
      x = V*arma::diagmat(s)*U.t()*b;
    }
    
    
	// -------------------------------
	// Put the results into the right place in the kernel matrix
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
            if (grad > 1)
            	continue;
			
			// Reverse x and y and z so it works right
            int kern_xx = ckx - jx;
			int kern_yy = cky - jy;
			int kern_zz = ckz - jz;
						
			if( (jc==ic) && (jz==0) && (jy==0) && (jx==0)){
			}else{				  				  
               k(kern_xx,kern_yy,kern_zz,jc,ic) = x(xrow,0);
               xrow++;
            }        
            
          }
        }
      }
    } // jc
          
  } // ic

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
  im = 0;  
  
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
	  ImCoil(Range(-ksx+x_off,ksx+x_off),Range(-ksy+y_off,ksy+y_off),Range(-ksz+z_off,ksz+z_off)) = KernelCoil;
	  fft(ImCoil);
	  
    }
  }
  
  if(debug){
   ArrayWriteMag(im,"KernelImage.dat");
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
  	cx_mat At;
	At.zeros(ncoils,ncoils);

	for(int iy = 0; iy < im.extent(secondDim); iy ++) {
      for(int ix = 0; ix < im.extent(firstDim); ix ++) {
  		
		// Copy to A		
		for(int jj = 0; jj < ncoils; jj++) {
          for(int ii = 0; ii < ncoils; ii++) {
            A(ii,jj) = (complex<double>)im(ix,iy,iz,jj,ii);
          }
        }
		
		// Compute AtA
		At = A.t();
		AtA = At*A;
  		
		// Eigen Vector
		vec eigval;
		cx_mat eigvec;
		eig_sym(eigval,eigvec,AtA);
		
		// Copy Back
		for(int jj = 0; jj < ncoils; jj++) {
			LR(ix,iy,iz,jj)=eigvec(jj,ncoils-1);
		}
		LR(ix,iy,iz,ncoils)=eigval(ncoils-1);
  }}}
  
  cout<< "Done with eig: took " << T << endl;
  
  if(debug){
  	ArrayWriteMag(LR,"CoilsEig.dat");
  }
}

void SPIRIT::interpMaps(Array< complex<float>,4 > &LR, Array< complex<float>,4 > &out) {
 
  int ncoils = out.length(fourthDim);
  
  Array< complex<float>,3>EigVal=LR(Range::all(),Range::all(),Range::all(),ncoils);
  mapthresh *= max(abs(EigVal));
  cout << "Interpolating espirit maps" << endl;
  out =0;
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
		  
		for (int x = 0; x < out.length(firstDim); x++) {
          int xx = x/mapshrink;
		  int xp = (xx + 1) % EigVal.length(firstDim);
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

  if(1==1){ // phase_type==SP_SMOOTH){
  
  	// Try to find a map that provides smooth combination 
  	Array< complex<float>,3> GRAD;
	GRAD.setStorage( ColumnMajorArray<3>());
  	GRAD.resize(ref.length(firstDim),ref.length(secondDim),ref.length(thirdDim));
  	
	Array< complex<float>,3> P;
	P.setStorage( ColumnMajorArray<3>());
  	P.resize(ref.length(firstDim),ref.length(secondDim),ref.length(thirdDim));
  	P = complex<float>(1.0,0.0);
	
	int Nc = maps.length(fourthDim)-1;
		
	if(debug==1){
		FILE *fid = fopen("SpiritCoilPhase.dat","a");
		for(int coil=0; coil < Nc; coil++){
		for(int j= 0; j<maps.extent(1); j++){
	  	for(int i=0; i<maps.extent(0);i++){
			float b = arg( P(i,j,ref.length(thirdDim)/2) * maps(i,j,ref.length(thirdDim)/2,coil) );
			fwrite(&b,1,sizeof(float),fid);
		}}}
		fclose(fid);
	}
	
	#pragma omp parallel for
  	for (int z = 0; z < ref.length(thirdDim); z++) {
      for (int y = 0; y < ref.length(secondDim); y++) {
      for (int x = 0; x < ref.length(firstDim); x++) {
    	complex<float>CI = conj( maps(x,y,z,0));
		CI /= abs(CI);
		for (int coil = 0; coil < Nc; coil++) {
			maps(x,y,z,coil)*=CI;
		}
		
     }}}
	
	if(debug==1){
		FILE *fid = fopen("SpiritCoilPhase.dat","a");
		for(int coil=0; coil < Nc; coil++){
		for(int j= 0; j<maps.extent(1); j++){
	  	for(int i=0; i<maps.extent(0);i++){
			float b = arg(P(i,j,ref.length(thirdDim)/2)*maps(i,j,ref.length(thirdDim)/2,coil));
			fwrite(&b,1,sizeof(float),fid);
		}}}
		fclose(fid);
	}	
	
	for(int iter=0; iter<100; iter++){
		// Get Gradient
		GRAD=(complex<float>(0.0,0));
		int Nx = ref.length(firstDim);
		int Ny = ref.length(secondDim);
		int Nz = ref.length(thirdDim);
		for(int coil=0; coil < Nc; coil++){
				
		#pragma omp parallel for
		for (int z = 0; z < ref.length(thirdDim); z++) {
     	for (int y = 0; y < ref.length(secondDim); y++) {
      	for (int x = 0; x < ref.length(firstDim); x++) {
  				int xm = (x-1 + Nx) % Nx;
				int xp = (x+1 + Nx) % Nx;
				int ym = (y-1 + Ny) % Ny;
				int yp = (y+1 + Ny) % Ny;
				int zm = (z-1 + Nz) % Nz;
				int zp = (z+1 + Nz) % Nz;
				GRAD(x,y,z) += conj(maps(x,y,z,coil))*( (complex<float>(6.0,0.0))*maps(x,y,z,coil)*P(x,y,z) - maps(xm,y,z,coil)*P(xm,y,z) - maps(xp,y,z,coil)*P(xp,y,z) - maps(x,ym,z,coil)*P(x,ym,z) - maps(x,yp,z,coil)*P(x,yp,z) - maps(x,y,zm,coil)*P(x,y,zm) -maps(x,y,zp,coil)*P(x,y,zp));
				
		}}}
		}
    		
		// Step
		GRAD *= 0.15;
		P -=GRAD;
		float scale = ((float)(Nx*Ny*Nz))/sum(abs(P));
		P *= scale; 
		
		cout<< "Iter = " << iter << ": Error " << sum(abs(GRAD)) << endl;
		
		if(debug==1){					
			FILE *fid = fopen("SpiritPhase.dat","a");
			for(int j= 0; j<maps.extent(1); j++){
			for(int i=0; i<maps.extent(0);i++){
					float b = arg(P(i,j,Nz/2));
					fwrite(&b,1,sizeof(float),fid);
			}}
			fclose(fid);
		
			// Export
			fid = fopen("SpiritCoilPhase.dat","a");
			for(int coil=0; coil < Nc; coil++){
			for(int j= 0; j<maps.extent(1); j++){
	  		for(int i=0; i<maps.extent(0);i++){
				float b = arg(P(i,j,ref.length(thirdDim)/2)*maps(i,j,ref.length(thirdDim)/2,coil));
				fwrite(&b,1,sizeof(float),fid);
			}}}
			fclose(fid);
		}
	}				
	
	// Apply Phase
	for (int z = 0; z < ref.length(thirdDim); z++) {
    for (int y = 0; y < ref.length(secondDim); y++) {
    for (int x = 0; x < ref.length(firstDim); x++) {
  		float mag =0.0;
		for(int coil=0; coil < Nc; coil++){			
			mag+=norm(maps(x,y,z,coil));
		}
		
		float phase = arg(P(x,y,z));
		complex<float>temp(cos(phase),sin(phase));
		temp = temp/sqrtf(mag);
		for(int coil=0; coil < Nc; coil++){			
			maps(x,y,z,coil)*=temp;
		}
	}}}

  }else{	
  	
  #pragma omp parallel for
  for (int z = 0; z < ref.length(thirdDim); z++) {
    for (int y = 0; y < ref.length(secondDim); y++) {
      for (int x = 0; x < ref.length(firstDim); x++) {
    	
		if(phase_type == SP_LOW_RES_PHASE){
			// Solve for scale factor + Normalization
			complex<float>CI(0,0);
			float mag=0.0;
			for (int coil = 0; coil < ref.length(fourthDim); coil++) {
				CI += ref(x,y,z,coil)*conj(maps(x,y,z,coil));
				mag += norm(maps(x,y,z,coil));
			}
			CI/=abs(CI);
			CI/=sqrt(mag);
		
			// Scale
			for (int coil = 0; coil < ref.length(fourthDim); coil++) {
				maps(x,y,z,coil)*=CI;
			}
		}else{
			complex<float>CI = conj( maps(x,y,z,0));
			
			float mag=0.0;
			for (int coil = 0; coil < ref.length(fourthDim); coil++) {
				mag += norm(maps(x,y,z,coil));
			}
			
			CI /= abs(CI);
			CI /= sqrt(mag);
			for (int coil = 0; coil < ref.length(fourthDim); coil++) {
				maps(x,y,z,coil)*=CI;
			}
		}
  }}}

  }

}

