#include "spirit.h"
#include <Eigen/Dense>

using arma::cx_mat;
using arma::vec;
using arma::uvec;
using namespace std;

SPIRIT::SPIRIT() 
{    
    
  kr_f = 0.0;
  krx_f = 3.0;  
  kry_f = 3.0;  
  krz_f = 3.0;  
  
  cr = 0;
  crx = 14;
  cry = 14;
  crz = 14;

  shape = SP_CIRCLE;
  calib_type = SP_TIK;
  calib_lam = .01;
  
  mapshrink = 2;
  mapthresh = 0.0;
  
  /*
  krx_f = kx;
  kry_f = ky;
  krz_f = kz;
  
  krx = (int)ceil(kx);
  kry = (int)ceil(ky);
  krz = (int)ceil(kz);
  
  cout << "Allocating spirit: " << nc << " " << krz << " " << kry << " " << krx << endl;
  
  k.alloc(nc,nc,krz*2+1,kry*2+1,krx*2+1);
  
  ncoils = nc;
  
  shape = SP_CIRCLE;
  calib_type = SP_TIK;   
  */
}

void SPIRIT::init(int xres, int yres, int zres, int nc){
  rcxres = xres;
  rcyres = yres;
  rczres = zres;
  ncoils = nc;
  
  if (kr_f > 0.0) {
    krx_f = kr_f;  
    kry_f = kr_f;  
    krz_f = kr_f;  
  }
  
  if (cr > 0) {
    crx = cr;  
    cry = cr;  
    crz = cr;  
  }
  
  krx = (int)ceil(krx_f);
  kry = (int)ceil(kry_f);
  krz = (int)ceil(krz_f);
  
  k.alloc(nc,nc,krz*2+1,kry*2+1,krx*2+1);
  
}

void SPIRIT::read_commandline(int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if (strcmp("-h", pstring[pos] ) == 0) {
	  	printf("\n*********************************************\n");
	  	printf(" SPIRiT Control:\n");
	  	printf("*********************************************\n");
	  
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

void SPIRIT::generateEigenCoils(array4D< complex<float> > &smaps, gridFFT &shrinkgrid, MRI_DATA &data){
		  arrayND< complex<float> >kdata = smaps; // switch to arrayND (indexing is more compatible with my old code)
		  kdata.fft3(-1);
		  
		  array4D< complex<float> >LRmaps;
		  LRmaps.alloc(ncoils+1,rczres/mapshrink+1,rcyres/mapshrink+1,rcxres/mapshrink+1);
		  
		  calibrate_ellipsoid(kdata);
		  prep();		  
		  getcoils(LRmaps);
		  
		  array4D< complex<float> > phaseref;
		  phaseref.alloc(ncoils,rczres/mapshrink,rcyres/mapshrink,rcxres/mapshrink);
		  
	    shrinkgrid.precalc_gridding(rczres/mapshrink,rcyres/mapshrink,rcxres/mapshrink,3);
	    
		  shrinkgrid.k_rad = 16;
		  for(int coil=0; coil< ncoils; coil++){
        shrinkgrid.forward( data.kdata[coil][0][0],data.kx[0][0],data.ky[0][0],data.kz[0][0],data.kw[0][0],data.Num_Pts*data.Num_Readouts);
        phaseref[coil]= ( shrinkgrid.return_array() );
      }
		  
		  rotateCoils(LRmaps, phaseref);
		  
		  interpMaps(LRmaps,smaps);
		  
		  im.freeArray();
		  LRmaps.freeArray();
		  phaseref.freeArray();
		  kdata.freeArray();
}

/**
* Calibrate SPIRiT kernel.
*/
void SPIRIT::calibrate_ellipsoid(arrayND< complex<float> > &kdata)
{
  cout << "Calibrating..." << endl;
  // Centers of k-space
  int cx, cy, cz;
  cx = kdata.dim[0]/2;
  cy = kdata.dim[1]/2;
  cz = kdata.dim[2]/2;
  
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
    int Arow, ki, Acol, kj;
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
          ki = (cx+ix) + kdata.dim[0]*(cy+iy) + kdata.dim[0]*kdata.dim[1]*(cz+iz) + kdata.dim[0]*kdata.dim[1]*kdata.dim[2]*ic; // Index in kdata
          b[tid](Arow,0) = (complex<double>)kdata.vals[ki];
  
          Acol = 0;
          for (int jc = 0; jc < ncoils; jc++) {
            for (int jz = -krz; jz <= krz; jz++) {
              for (int jy = -kry; jy <= kry; jy++) {
                for (int jx = -krx; jx <= krx; jx++) {
                  
                  float krad = (shape == SP_SQUARE) ? 0 : (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f) + (float)(jz*jz)/(krz_f*krz_f);
                  if (krad <= 1) {
                  kj = (cx+ix+jx) + kdata.dim[0]*(cy+iy+jy) + kdata.dim[0]*kdata.dim[1]*(cz+iz+jz) + kdata.dim[0]*kdata.dim[1]*kdata.dim[2]*jc;
                  
                  if (kj != ki) {
                    A[tid](Arow,Acol) = (complex<double>)kdata.vals[kj];
                    Acol++;
                  }
                  } else {continue;}
                  
                }
              }
            }
          } // jc
          Arow++;
          } else {continue;}
          
        }
      }
    } // iz
    
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
    
    
    // Put the results of x into the right place in the kernel matrix
    // Kernel centers
    int ckx = k.dim[0]/2;
    int cky = k.dim[1]/2;
    int ckz = k.dim[2]/2;
    
    int cindex =  ckx + k.dim[0]*cky + k.dim[0]*k.dim[1]*ckz
                  + k.dim[0]*k.dim[1]*k.dim[2]*ic
                  + k.dim[0]*k.dim[1]*k.dim[2]*k.dim[3]*ic;
    
    int xrow = 0;
    int Gindex;
    for (int jc = 0; jc < ncoils; jc++) {
      for (int jz = -krz; jz <= krz; jz++) {
        for (int jy = -kry; jy <= kry; jy++) {
          for (int jx = -krx; jx <= krx; jx++) {
            
            float grad = (shape == SP_SQUARE) ? 0 : (float)(jx*jx)/(krx_f*krx_f) + (float)(jy*jy)/(kry_f*kry_f) + (float)(jz*jz)/(krz_f*krz_f);
            if (grad <= 1) {
            // Reverse x and y and z so it works right
            Gindex = (ckx-jx) + k.dim[0]*(cky-jy)
                    + k.dim[0]*k.dim[1]*(ckz-jz)
                    + k.dim[0]*k.dim[1]*k.dim[2]*jc
                    + k.dim[0]*k.dim[1]*k.dim[2]*k.dim[3]*ic;
            
                    
            if (Gindex != cindex) {
              k.vals[Gindex] = x[tid](xrow,0);
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
  
  im.alloc(nc,nc,sz,sy,sx);
  
  /////////////
  // Copy k into the right place into im
  int im_i, k_i;
  k_i = 0;
  int x_off = (im.dim[0]/2);
  int y_off = (im.dim[1]/2);
  int z_off = (im.dim[2]/2);
  
  int ksx = ((k.dim[0]-1)/2);
  int ksy = ((k.dim[1]-1)/2);
  int ksz = ((k.dim[2]-1)/2);
  
  int neg_mod = 0 - (im.dim[0]/2+im.dim[1]/2+im.dim[2]/2)%2;
  if (neg_mod==0) {neg_mod = 1;}
  
  for (int j = 0; j < k.dim[3]; j++) {
    for (int i = 0; i < k.dim[3]; i++) {
      for (int iz = -ksz; iz <= ksz; iz++) {
        for (int iy = -ksy; iy <= ksy; iy++) { 
          for (int ix = -ksx; ix <= ksx; ix++) {   
            im_i = (ix + x_off) + im.dim[0]*(iy + y_off) +
                   im.dim[0]*im.dim[1]*(iz + z_off) +
                   im.dim[0]*im.dim[1]*im.dim[2]*i +
                   im.dim[0]*im.dim[1]*im.dim[2]*im.dim[3]*j;
            im.vals[im_i] = k.vals[k_i]*(float)neg_mod;
            k_i++;
          }
        }
      }
    }
  }
  
  im.fft3(1);
  
}

using Eigen::MatrixXcd;

/**
* Generate coils from SPIRiT kernel.
*/
void SPIRIT::getcoils(array4D< complex<float> > &LR)
{
  tictoc T;
  T.tic();
  cout << "Starting eig" << endl;
   
  int imdim[] = {im.dim[0],im.dim[1],im.dim[2]};
  int ncoils = im.dim[3];
  
  int nthreads = omp_get_max_threads();

  MatrixXcd* A = new MatrixXcd[nthreads];
  MatrixXcd* AtA = new MatrixXcd[nthreads];
  Eigen::SelfAdjointEigenSolver<MatrixXcd>* EIG = new Eigen::SelfAdjointEigenSolver<MatrixXcd>[nthreads];
  MatrixXcd* V = new MatrixXcd[nthreads];
  Eigen::VectorXd* SV = new Eigen::VectorXd[nthreads];
  
  for (int i = 0; i < nthreads; i++) {
    A[i].resize(ncoils,ncoils);
    AtA[i].resize(ncoils,ncoils);
    V[i].resize(ncoils,ncoils);
    SV[i].resize(ncoils);
    EIG[i] = Eigen::SelfAdjointEigenSolver<MatrixXcd>(ncoils);
  }
  
  #pragma omp parallel for
  for(int iz = 0; iz < imdim[2]; iz ++) {
    for(int iy = 0; iy < imdim[1]; iy ++) {
      for(int ix = 0; ix < imdim[0]; ix ++) {
        
        int tid = omp_get_thread_num();
        
        for(int jj = 0; jj < ncoils; jj++) {
          for(int ii = 0; ii < ncoils; ii++) {
            A[tid](jj,ii) = (complex<double>)im(jj,ii,iz,iy,ix);
          }
        }
        
        AtA[tid] = A[tid].adjoint() * A[tid];
        EIG[tid].compute(AtA[tid]);
        
        V[tid] = EIG[tid].eigenvectors();
        SV[tid] = EIG[tid].eigenvalues();
        
        for(int recv_cnt = 0; recv_cnt < ncoils; recv_cnt++) {
          LR[recv_cnt][iz][iy][ix] = V[tid](recv_cnt,ncoils-1);
        }
        LR[ncoils][iz][iy][ix] = sqrt(SV[tid](ncoils-1)); //eigenvalue 
  }}}
  T.toc("Done with eig");
}

void SPIRIT::interpMaps(array4D< complex<float> > &LR, array4D< complex<float> > &out) {
  
  int ncoils = out.Nt;
  
  mapthresh *= abs(LR[ncoils].max()); 
  
  cout << "Interpolating espirit maps" << endl;
  
  // Linear interpolation up to the full resolution
  for(int i = 0; i < ncoils; i++) {
    #pragma omp parallel for
    for (int z = 0; z < out.Nz; z++) {
      int zz = z/mapshrink;
      float dz = (float)(z%mapshrink)/mapshrink;
      for (int y = 0; y < out.Ny; y++) {
        int yy = y/mapshrink;
        float dy = (float)(y%mapshrink)/mapshrink;
        for (int x = 0; x < out.Nx; x++) {
          int xx = x/mapshrink;
          float dx = (float)(x%mapshrink)/mapshrink;
          
          float svd_interp = abs((1-dz)*(1-dy)*(1-dx)*LR[ncoils][zz][yy][xx] + \
                            (1-dz)*(1-dy)*(dx)*LR[ncoils][zz][yy][xx+1] + \
                            (1-dz)*(dy)*(1-dx)*LR[ncoils][zz][yy+1][xx] + \
                            (1-dz)*(dy)*(dx)*LR[ncoils][zz][yy+1][xx+1] + \
                            (dz)*(1-dy)*(1-dx)*LR[ncoils][zz+1][yy][xx] + \
                            (dz)*(1-dy)*(dx)*LR[ncoils][zz+1][yy][xx+1] + \
                            (dz)*(dy)*(1-dx)*LR[ncoils][zz+1][yy+1][xx] + \
                            (dz)*(dy)*(dx)*LR[ncoils][zz+1][yy+1][xx+1]);
          
          if (svd_interp < mapthresh) {
           out[i][z][y][x] = 0;
          } else {
                            
          out[i][z][y][x] = (1-dz)*(1-dy)*(1-dx)*LR[i][zz][yy][xx] + \
                            (1-dz)*(1-dy)*(dx)*LR[i][zz][yy][xx+1] + \
                            (1-dz)*(dy)*(1-dx)*LR[i][zz][yy+1][xx] + \
                            (1-dz)*(dy)*(dx)*LR[i][zz][yy+1][xx+1] + \
                            (dz)*(1-dy)*(1-dx)*LR[i][zz+1][yy][xx] + \
                            (dz)*(1-dy)*(dx)*LR[i][zz+1][yy][xx+1] + \
                            (dz)*(dy)*(1-dx)*LR[i][zz+1][yy+1][xx] + \
                            (dz)*(dy)*(dx)*LR[i][zz+1][yy+1][xx+1];
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
void SPIRIT::rotateCoils(array4D< complex<float> > &maps, array4D< complex<float> > &ref)
{
  #pragma omp parallel for
  for (int z = 0; z < ref.Nz; z++) {
    for (int y = 0; y < ref.Ny; y++) {
      for (int x = 0; x < ref.Nx; x++) {
        float diff = 0.0;
        float rot = 0.0;
        float scale = 0.0;
        for (int coil = 0; coil < ref.Nt; coil++) {
          diff = arg(ref[coil][z][y][x]) - arg(maps[coil][z][y][x]);
          if (coil > 0) { // Unwrap based on current guess of phase
           diff -= 2.0*PI*round((diff-rot/scale)/2.0/PI); 
          }
          rot += diff*norm(ref[coil][z][y][x]);
          scale += norm(ref[coil][z][y][x]);
        }
        rot /= scale;
        for (int coil = 0; coil < ref.Nt; coil++) {
          maps[coil][z][y][x] *= complex<float>(cos(rot),sin(rot)); 
        }     
  }}}
}

/** 
* Perform a DFT to get a single point in image space. Note that u, v, and w are 
* supposed to come in scaled by the image size.  (After trying this is way too 
* slow to be useful.)
*/
complex<float> SPIRIT::getImPointDFT(float u, float v, float w, int c1, int c2) {
  complex<float> out = 0.0;
  complex<float> i1 = complex<float>(0.0,2.0*PI);
  
  for(int z = 0; z < k.dim[2]; z++) {
    for(int y = 0; y < k.dim[1]; y++) {
      for(int x = 0; x < k.dim[0]; x++) {
         int mod = ((x+y+z)%2) == 0 ? 1 : -1;
         out += (float)mod*k(c2,c1,z,y,x)*exp(i1*(u*x + v*y + w*z));
  }}}
  
  return out;
  
}

