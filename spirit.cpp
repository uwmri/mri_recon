#include "spirit.h"

using namespace arma;
using namespace std;

SPIRIT::SPIRIT(int kx, int ky, int kz, int nc) 
{
    k.alloc(nc,nc,kz,ky,kx);
}

void SPIRIT::calibrate(int cxs, int cys, int czs, 
                         arrayND< complex<float> > &kdata, double calib_lam)
{
  // Data has most likely come in image space, so do a FFT
  cout << "ifft for calibration" << endl;
  ifft3(kdata);
  
  cout << "Calibrating..." << endl;
  // Centers of k-space
  int cx, cy, cz;
  cx = kdata.dim[0]/2;
  cy = kdata.dim[1]/2;
  cz = kdata.dim[2]/2;
  
  // Size of cropped source data
  int wx = cxs-(k.dim[0]-1);
  int wy = cys-(k.dim[1]-1);
  int wz = czs-(k.dim[2]-1);
  
  // Number of ACS points we have (omits border points)
  int num_ACS = (wx*wy*wz);
  
  if (num_ACS < (k.dim[3]*k.dim[0]*k.dim[1]*k.dim[2]-1)) {
    printf("WARNING: Not enough calibration points, system is underdetermined\n");
  } 
  
  // Set up A and b
  int ncoils = k.dim[3];
  
  cx_mat* nA = new cx_mat[ncoils];
  cx_mat* nb = new cx_mat[ncoils];
  cx_mat* Ap = new cx_mat[ncoils];
  cx_mat* ApA = new cx_mat[ncoils];
  cx_mat* S = new cx_mat[ncoils];
  cx_mat* x = new cx_mat[ncoils];
  
  // Spirit kernel half width
  int Gwx, Gwy, Gwz;
  Gwx = ((k.dim[0]-1)/2);
  Gwy = ((k.dim[1]-1)/2);
  Gwz = ((k.dim[2]-1)/2);
  
  // ACS region in the source data
  int startx, starty, startz;
  int stopx, stopy, stopz;
  startx = cx - (cxs/2) + Gwx;
  stopx  = cx + (cxs/2) - Gwx;
  starty = cy - (cys/2) + Gwy;
  stopy  = cy + (cys/2) - Gwy;
  startz = cz - (czs/2) + Gwz;
  stopz  = cz + (czs/2) - Gwz;  
  
  if ((czs == 1) && (Gwz == 0)) {
    stopz++;
  }
  
  // Fill in A and b
  #pragma omp parallel for
  for (int ic = 0; ic < k.dim[3]; ic++) {
    int Abi, ki, Acol, kj;
    int tid = ic;
    nA[tid].zeros(num_ACS,k.dim[0]*k.dim[1]*k.dim[2]*k.dim[3]-1);
    nb[tid].zeros(num_ACS,1);
    //printf("Calibrating coil: %d\n",ic);
    for (int iz = startz; iz<stopz; iz++) { 
      for (int iy = starty; iy<stopy; iy++) {
        for (int ix = startx; ix<stopx; ix++) {
          
          Abi = (ix-startx) + wx*(iy-starty) + wx*wy*(iz-startz); // Row in A and b 
          ki = ix + kdata.dim[0]*iy + kdata.dim[0]*kdata.dim[1]*iz + kdata.dim[0]*kdata.dim[1]*kdata.dim[2]*ic; // Index in kdata
          
          // ACS points (b vector)
          nb[tid](Abi,0) = (complex<double>)kdata.vals[ki];

          // Grab surrounding source points to go in A
          Acol = 0; // Easier way to index which column we are in for the A matrix
          for (int jc = 0; jc < k.dim[3]; jc++) {
            for (int jz = iz-Gwz; jz <=  iz+Gwz; jz++) {
              for (int jy = iy-Gwy; jy <=  iy+Gwy; jy++) {
                for (int jx = ix-Gwx; jx <=  ix+Gwx; jx++) {
                  
                  kj = jx + kdata.dim[0]*jy + kdata.dim[0]*kdata.dim[1]*jz + kdata.dim[0]*kdata.dim[1]*kdata.dim[2]*jc;
                  
                  if (kj != ki) {
                    nA[tid](Abi,Acol) = (complex<double>)kdata.vals[kj];
                    Acol++;
                  }
                  
                }
              }
            }
          } // jc
          
        }
      }
    } // iz
    
    /*
    cx_mat U;
    U.zeros(A.n_rows,A.n_cols);
    vec s;
    s.zeros(A.n_cols);
    cx_mat V;
    V.zeros(A.n_cols,A.n_cols);
    svd_econ(U,s,V,A);
    cout << "SVD:" << endl;
    cout << s.max() << " " << s.min() << endl;
    
    // TSVD stuff
    int tsvd_i = 0;
    double tsvd_max = s.max();
    while( s(tsvd_i)>0.05*tsvd_max ) { tsvd_i++; }
    cout << "Last SV index = " << tsvd_i << " " << s(tsvd_i) << endl;
    for(int ti = 0; ti < A.n_cols; ti++) {
      if (ti<tsvd_i) {s(ti) = 1.0/s(ti);}
      else {s(ti) = 0;}
    }
    cx_mat x = V*diagmat(s)*U.t()*b;
    */
    
    Ap[tid] = trans(nA[tid]);
    ApA[tid] = Ap[tid] * nA[tid];

    float beta = norm(ApA[tid],"fro")*(calib_lam);
    S[tid].eye(ApA[tid].n_rows,ApA[tid].n_cols);
    S[tid] *= beta;
    
    ApA[tid] += S[tid];
    
    x[tid] = solve(ApA[tid],Ap[tid]) * nb[tid];
    
    
    // Now put the results of x into the right place in kdata
    int cindex =  Gwx + k.dim[0]*Gwy + k.dim[0]*k.dim[1]*Gwz
                  + k.dim[0]*k.dim[1]*k.dim[2]*ic
                  + k.dim[0]*k.dim[1]*k.dim[2]*k.dim[3]*ic;
    
    
    int xrow = 0;
    int Gindex;
    for (int kc = 0; kc < k.dim[3]; kc++) {
      for (int kz = 0; kz <  k.dim[2]; kz++) {
        for (int ky = 0; ky <  k.dim[1]; ky++) {
          for (int kx = 0; kx <  k.dim[0]; kx++) {
            
            // Reverse x and y and z so it works right
            Gindex = (k.dim[0]-(kx+1)) + k.dim[0]*(k.dim[1]-(ky+1))
                    + k.dim[0]*k.dim[1]*(k.dim[2]-(kz+1))
                    + k.dim[0]*k.dim[1]*k.dim[2]*kc
                    + k.dim[0]*k.dim[1]*k.dim[2]*k.dim[3]*ic;
            
            if (Gindex != cindex) {
              k.vals[Gindex] = x[tid](xrow,0);
              xrow++;
            }
            
          }
        }
      }
    } // kc
    
          
  } // ic
  return;    
}

void SPIRIT::prep(int sx, int sy, int sz, int nc) {
  
  im.alloc(nc,nc,sz,sy,sx);
  
  int imsize3 = sx*sy*sz;
  int fftdim[] = {im.dim[2], im.dim[1], im.dim[0]}; 
  fftwf_init_threads();
  fftwf_plan_with_nthreads(16);
  
  // FFTW plan
  fftwf_plan fftplan = fftwf_plan_many_dft(3, fftdim, im.dim[3]*im.dim[3],
                                  (fftwf_complex *) im.vals, fftdim,
                                  1, imsize3,
                                  (fftwf_complex *) im.vals, fftdim,
                                  1, imsize3,
                                  FFTW_FORWARD, FFTW_MEASURE);
  
  /////////////
  // Copy k into the right place into im
  int im_i, k_i;
  k_i = 0;
  int x_off = (im.dim[0]/2) - ((k.dim[0]-1)/2);
  int y_off = (im.dim[1]/2) - ((k.dim[1]-1)/2);
  int z_off = (im.dim[2]/2) - ((k.dim[2]-1)/2);
  
  int neg_mod = 0 - (im.dim[0]/2+im.dim[1]/2+im.dim[2]/2)%2;
  if (neg_mod==0) {neg_mod = 1;}
  
  for (int j = 0; j < k.dim[3]; j++) {
    for (int i = 0; i < k.dim[3]; i++) {
      for (int iz = 0; iz < k.dim[2];iz++) {
        for (int iy = 0; iy < k.dim[1];iy++) { 
          for (int ix = 0; ix < k.dim[0];ix++) {   
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
  
  /////////
  // FFT shift
  for (int i = 0; i < im.dim[3]*im.dim[3]; i++) {
    int offset = i * im.dim[0]*im.dim[1]*im.dim[2];
    fftshift3(im.dim[0],im.dim[1],im.dim[2],(im.vals+offset));
  }
  
  //////////////
  // FFT all the matrices    
  fftwf_execute(fftplan);
  
  ////////
  // FFT shift
  for (int i = 0; i < im.dim[3]*im.dim[3]; i++) {
    int offset = i * im.dim[0]*im.dim[1]*im.dim[2];
    fftshift3(im.dim[0],im.dim[1],im.dim[2],(im.vals+offset));
  }
  
}

void SPIRIT::getcoils(array4D< complex<float> > &out, int shrink, float thresh)
{
  cout << "Starting eig" << endl;
  
  int imdim[] = {im.dim[0],im.dim[1],im.dim[2]};
  int ncoils = im.dim[3];
  
  array4D< complex<float> >LR; // Low res coils, +1 in dimension for interp
  LR.alloc(im.dim[3]+1,im.dim[2]+1,im.dim[1]+1,im.dim[0]+1); //+1 in coils for svd val
  LR.zero();
  
  int nthreads = omp_get_max_threads();
  
  cx_mat* nA = new cx_mat[nthreads];
  cx_mat* nAtA = new cx_mat[nthreads];
  cx_mat* nU = new cx_mat[nthreads];
  vec* ns = new vec[nthreads];
  cx_mat* nV = new cx_mat[nthreads];
  
  // Initialize
  for (int i = 0; i < nthreads; i++) {
    nA[i].zeros(ncoils,ncoils);
    nAtA[i].zeros(ncoils,ncoils);
    nU[i].zeros(ncoils,ncoils);
    ns[i].zeros(ncoils,1);
    nV[i].zeros(ncoils,ncoils);
  }
  
  /*
  cx_mat A;
  A.zeros(ncoils*nthreads, ncoils);
  
  A = trans(nA[0]);
  
  cx_mat AtA;
  AtA.zeros(ncoils*nthreads, ncoils);
  
  cx_mat U;
  U.zeros(ncoils*nthreads, ncoils);
  vec s;
  s.zeros(ncoils*nthreads, 1);
  cx_mat V;
  V.zeros(ncoils*nthreads, ncoils);
  */
  
  
  #pragma omp parallel for
  for(int iz = 0; iz < imdim[2]; iz ++) {
    for(int iy = 0; iy < imdim[1]; iy ++) {
      for(int ix = 0; ix < imdim[0]; ix ++) {
        
        int tid = omp_get_thread_num();
        
        for(int jj = 0; jj < ncoils; jj++) {
          for(int ii = 0; ii < ncoils; ii++) {
            int ki = ix + imdim[0]*iy + imdim[0]*imdim[1]*iz + imdim[0]*imdim[1]*imdim[2]*ii + imdim[0]*imdim[1]*imdim[2]*ncoils*jj;
            nA[tid](jj,ii) = complex<double>(im.vals[ki]);
          }
        }
        
        nAtA[tid] = trans(nA[tid]) * nA[tid];
        svd_econ(nU[tid],ns[tid],nV[tid],nAtA[tid]);
        
        //uword index;
        //realval = real(eigval);
        //realval.max(index);
        
        for(int recv_cnt = 0; recv_cnt < ncoils; recv_cnt++) {
          LR[recv_cnt][iz][iy][ix] = nV[tid](recv_cnt,0);
        }
        LR[ncoils][iz][iy][ix] = ns[tid](0);
        //out[ncoils][iz][iy][ix] = s(0);
  }}}
  
  //thresh *= abs(LR[ncoils].max()); 
  
  // Linear interpolation up to the full resolution
  for(int i = 0; i < ncoils; i++) {
    #pragma omp parallel for
    for (int z = 0; z < out.Nz; z++) {
      int zz = z/shrink;
      float dz = (float)(z%shrink)/shrink;
      for (int y = 0; y < out.Ny; y++) {
        int yy = y/shrink;
        float dy = (float)(y%shrink)/shrink;
        for (int x = 0; x < out.Nx; x++) {
          int xx = x/shrink;
          float dx = (float)(x%shrink)/shrink;
          
          float svd_interp = abs((1-dz)*(1-dy)*(1-dx)*LR[ncoils][zz][yy][xx] + \
                            (1-dz)*(1-dy)*(dx)*LR[ncoils][zz][yy][xx+1] + \
                            (1-dz)*(dy)*(1-dx)*LR[ncoils][zz][yy+1][xx] + \
                            (1-dz)*(dy)*(dx)*LR[ncoils][zz][yy+1][xx+1] + \
                            (dz)*(1-dy)*(1-dx)*LR[ncoils][zz+1][yy][xx] + \
                            (dz)*(1-dy)*(dx)*LR[ncoils][zz+1][yy][xx+1] + \
                            (dz)*(dy)*(1-dx)*LR[ncoils][zz+1][yy+1][xx] + \
                            (dz)*(dy)*(dx)*LR[ncoils][zz+1][yy+1][xx+1]);
          
          if (svd_interp < thresh) {
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
  cout << "Done with eig" << endl;
}

//ifft (not just for ndims=3, but a 3D transform)
void ifft3(arrayND< complex<float> > &in){
  fftwf_init_threads();
  fftwf_plan_with_nthreads(omp_get_max_threads());
                 
  int imsize3 = in.dim[2]*in.dim[1]*in.dim[0];
  int fftdim[] = {in.dim[2], in.dim[1], in.dim[0]}; 
  
  int ntransforms = 1;
  for (int i = 3; i < in.ndims; i++) {ntransforms*=in.dim[i];}
  
  // FFTW plan
  fftwf_plan fftplan = fftwf_plan_many_dft(3, fftdim, ntransforms,
                        (fftwf_complex *) in.vals, fftdim,
                        1, imsize3,
                        (fftwf_complex *) in.vals, fftdim,
                        1, imsize3,
                        FFTW_BACKWARD, FFTW_ESTIMATE);
  
  for (int i = 0; i < in.dim[3]; i++) {
    int offset = i * imsize3;
    fftshift3(in.dim[0],in.dim[1],in.dim[2],(in.vals+offset));
  }
  
  fftwf_execute(fftplan);
  
  for (int i = 0; i < in.dim[3]; i++) {
    int offset = i * imsize3;
    fftshift3(in.dim[0],in.dim[1],in.dim[2],(in.vals+offset));
  }
}
		 


void fftshift3(int xs, int ys, int zs, complex<float>* data)
{
  int index, mod;
  #pragma omp parallel for default(none) shared(data,xs,ys,zs) private(index,mod)
  for (int k = 0; k < zs; k++) {
    for (int j = 0; j < ys; j++) {
      for (int i = 0; i < xs; i++) {
        index = i + j*xs + k*xs*ys;
        mod = ((i+j+k)%2) == 0 ? 1 : -1;
        data[index] = (float)mod*data[index];
      }
    }
  }
  return;
}
