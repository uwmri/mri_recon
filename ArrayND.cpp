// ND array class

#define MAXDIMS 16
#define OMPSPLIT 1024

#include <fftw3.h>

template < class C >
class arrayND{ 
	public:	
		C *vals; 
		
		int ndims;  // Number of active dimensions
		int imdims; // Unused so far, potentially to dilineate image dimensions (2D vs 3D) from other dims (for ffts and such)
		int dim[MAXDIMS]; // Fastest changing dimension first (i.e. x)
		size_t numel; // Total number of elements
		
		int MemExists;
		
		/*---------------------------------------------
		 *    Constructor / Decontructor / etc
		 *---------------------------------------------*/
		
		// Constructor
		arrayND< C >(){
			MemExists = 0;
			memset(dim,0,(size_t)(sizeof(int)*MAXDIMS));
		}
		
		arrayND< C >(const array4D<C> &copy){
			MemExists = 0;
			memset(dim,0,(size_t)(sizeof(int)*MAXDIMS));
			
			alloc(copy.Nt, copy.Nz, copy.Ny, copy.Nx);
			
			int i = 0;
			
			for (int c = 0; c < copy.Nt; c++) {
      for (int z = 0; z < copy.Nz; z++) {
      for (int y = 0; y < copy.Ny; y++) {
      for (int x = 0; x < copy.Nx; x++) {
        vals[i] = copy.vals[c][z][y][x];
        i++;
			}}}}
			
		}
		
		/*
		void copy_pointers( array4D< C > &temp ){
			Nz = temp->Nz;
			Ny = temp->Ny;
			Nx = temp->Nx;
			vals = temp->vals[0][0][0];
		}
		*/
		
		// Allocation with Size
		void alloc(int z,int y,int x){
			ndims = 3;
			dim[0] = x;
			dim[1] = y;
			dim[2] = z;
			
			numel = (size_t)x*y*z;
			
			vals = alloc_vals();
			MemExists = 1;
		}
		
		// Allocation with Size
		void alloc(int d4, int z,int y,int x){
			ndims = 4;
			dim[0] = x;
			dim[1] = y;
			dim[2] = z;
			dim[3] = d4;
			
			numel = (size_t)x*y*z*d4;
			
			vals = alloc_vals();
			MemExists = 1;
		}
		
		
		void alloc(int d5, int d4, int z,int y,int x){
			ndims = 5;
			dim[0] = x;
			dim[1] = y;
			dim[2] = z;
			dim[3] = d4;
			dim[4] = d5;
			
			numel = (size_t)x*y*z*d4*d5;
			
			vals = alloc_vals();
			MemExists = 1;
		}
		
		
		// Base Code to Allocate a Contigous Memory Block	
		C *alloc_vals( void ){
			C* temp;
			temp = new C[numel];
			memset(temp,0,(size_t)(sizeof(C)*numel));
			return(temp);		
		}	
			
		// Cleans Memory	
		void freeArray(){
      if(MemExists){
        free(vals);
      // cout << "Destroying" << endl;
      }
    }
		
    inline C& operator()(int z, int y, int x) {
     return vals[x+y*dim[0]+z*dim[0]*dim[1]]; 
    }
    
    inline C& operator()(int d4, int z, int y, int x) {
     return vals[x+y*dim[0]+z*dim[0]*dim[1]+d4*dim[0]*dim[1]*dim[2]]; 
    }
    
    inline C& operator()(int d5, int d4, int z, int y, int x) {
     return vals[x+y*dim[0]+z*dim[0]*dim[1]+d4*dim[0]*dim[1]*dim[2]+d5*dim[0]*dim[1]*dim[2]*dim[3]]; 
    }
    
    	/*
		void point_to_3D( array3D< C > *temp ){
			Nz = temp->Nz;
			Ny = temp->Ny;
			Nx = temp->Nx;
			vals = temp->vals;
		}
		*/
		
		/*---------------------------------------------
		 *    Utilities - Operators
		 *---------------------------------------------*/
		 
		 // Allow Indexing with arrayND M[k][j][i] (suppose to be bad but useful for coding)
		 /*
		 inline C** & operator[](const int i){
		 	return( vals[i]);		 
		 }
		 
		 inline C & operator ()(const int i){
		 	return( vals[0][0][i]);		 
		 }
		 */
		 
		 // Equals Sets Value		 
		 void operator = (C temp){
		 	#pragma omp parallel for 
			for(int i=0; i<numel; i++){
				vals[i]=temp;					
        }
       }
		 
		 // Scaling overloaded	
		 void operator *= (float X){
		 	#pragma omp parallel for 		 
			for(int i=0; i<numel; i++){
				vals[i] *= X;			
			}
		 }
		 
		 // Scaling overloaded	
		 void operator *= (complex<float> XX){
		 	#pragma omp parallel for 		 
			for(int i=0; i<numel; i++){
				vals[i] *= XX;			
			}	 
		 }
		 
		 // Point Wise Multiply of Same Type
		 template < class X>		 		 
		 void operator *= (X temp){
		 	check_array_dims(temp);
			
			#pragma omp parallel for 		 
			for(int i=0; i<numel; i++){
				vals[i] *= temp.vals[i];					
			}	 
		 }
		 
		 // Set Array to Zero Fast
		 void zero( void ){
		 	memset(vals,0,sizeof(C)*numel);
		 }
		 
		 // Pointwise equals 
		 void operator = (arrayND<C>temp){
		 	check_array_dims(temp);
			
			#pragma omp parallel for 
			for(int i=0; i<numel; i++){
				vals[i]=temp.vals[i];					
			} 
		 }
		 
		 // Sum of X^2
		 double energy (void){
		  
		 	double *temp = new double[OMPSPLIT];
		 	memset(temp,0,(size_t)(sizeof(double)*OMPSPLIT));
		 	size_t splitsize = (int)(numel/OMPSPLIT)+1;
		 	
		 	#pragma omp parallel for 		 
			for(int i=0; i< OMPSPLIT; i++){
			  double sum = 0.0;
			  size_t start = i*splitsize;
			  size_t stop = (i+1)*splitsize;
			  if (stop > numel) {stop=numel;}
			  
			  for(unsigned int j = start; j < stop; j++) {
			    sum += (double)norm(vals[j]);
			  }
			  
			  temp[i] = sum;	
			}
			
			
			double X=0.0;
			for(int i=0; i<OMPSPLIT; i++){
				X += temp[i];
			}	
			X *= 1.0/(float)numel;
			
			delete []temp;
		 	return(X);
		 } 
		 
		  // Sum of X
		 C max( void ){
			
		  C *temp = new C[OMPSPLIT];
		 	memset(temp,0,(size_t)(sizeof(C)*OMPSPLIT));
		 	size_t splitsize = (int)(numel/OMPSPLIT)+1;
		 	
		 	#pragma omp parallel for 		 
			for(int i=0; i< OMPSPLIT; i++){
			  C M = 0.0;
			  size_t start = i*splitsize;
			  size_t stop = (i+1)*splitsize;
			  if (stop > numel) {stop=numel;}
			  
			  for(unsigned int j = start; j < stop; j++) {
			    M = ( abs(M) > abs(vals[j])) ? ( M ) : ( vals[j]);
			  }
			  
			  temp[i] = M;	
			}
		 	
			
			C X = 0.0;
			for(int i=0; i<OMPSPLIT; i++){
				X = ( abs(temp[i]) > abs(X) ) ? ( temp[i] ) : ( X );
			}		
			
			delete []temp;
		 	return(X);		
		 } 
		 
		 // Sum of X
		 C sum( void ){
		   
      C *temp = new C[OMPSPLIT];
		 	memset(temp,0,(size_t)(sizeof(C)*OMPSPLIT));
		 	size_t splitsize = (int)(numel/OMPSPLIT)+1;
			
		 	#pragma omp parallel for 		 
			for(int i=0; i< OMPSPLIT; i++){
			  C sum = 0.0;
			  size_t start = i*splitsize;
			  size_t stop = (i+1)*splitsize;
			  if (stop > numel) {stop=numel;}
			  
			  for(unsigned int j = start; j < stop; j++) {
			    sum += vals[j];
			  }
			  
			  temp[i] = sum;	
			}
		 	
		 	C X=0.0;
			for(int i=0; i<OMPSPLIT; i++){
				X += temp[i];
			}	
			X *= 1.0/(float)numel;
			
			delete []temp;
		 	return(X);
		 } 
		 	
		 // In place Square		 
		 void sqr( void ){
		 	#pragma omp parallel for 
			for(int i=0; i<numel; i++){
				vals[i] *=vals[i];					
			}	 
		 }
		 		 		 
		 // Doesn't work/make sense for complex
		 void sqrt( void ){
		 	#pragma omp parallel for 
			for(int i=0; i<numel; i++){
				vals[i]=sqrtf(vals[i]);					
			}	 
		 }
		 
		 // Only works for complex 
		 void complex_conj( void ){
		 	#pragma omp parallel for 
			for(int i=0; i<numel; i++){
				vals[i]=conj(vals[i]);					
        }			 
       }
		 
		 void check_array_dims( arrayND<C>temp ){
		 	
		   if( temp.ndims != ndims){
		     cout << "Checking array dimensions failed on: # of dimensions (" << ndims << " != " << temp.ndims << ")" << endl;
		   }else{
		     for (int i = 0; i < ndims; i++) {
		       if( temp.dim[i] != dim[i]){
             cout << "Checking array dimensions failed on: dim[" << i << "] (" << dim[i] << " != " << temp.dim[i] << ")" << endl;
           }
		     }
		   }
		 }
		 
		 void conjugate_multiply(arrayND<C>temp){ 
		 	
			check_array_dims(temp);
			
			#pragma omp parallel for 		 
			for(int i=0; i<numel; i++){
				vals[i] *= conj(temp.vals[i]);					
			}	
		 
		 }
		 
		 // For Add calcs 
		 template < class X>		 		 
		 void operator += (X temp){
		 	check_array_dims(temp);
			
			#pragma omp parallel for 		 
			for(int i=0; i<numel; i++){
				vals[i] += temp.vals[i];					
			}	 
		 }
		 
		 template < class X>		 		 
		 void operator -= (X temp){
		 
		 	check_array_dims(temp);
			
		 	#pragma omp parallel for 		 
			for(int i=0; i<numel; i++){
				vals[i] -=  temp.vals[i];					
			}		 
		 }
		 		 		 	
		 void fft3(int dir) {
		   fftwf_init_threads();
       fftwf_plan_with_nthreads(omp_get_max_threads());
          
        //- Load old Plan if Possible
        FILE *fid;
        FILE *fid2;
        char fname[300];
        char hname[300];
        char com[1300];
        gethostname( hname,299);
      #ifdef SCANNER
        sprintf(fname,"/usr/g/research/pcvipr/fft_wisdom_host_%s_x%d_y%d_z%d.dat",hname,dim[0],dim[1],dim[2]);
      #else
        sprintf(fname,"/export/home/kmjohnso/FFT_PLANS/fft_wisdom_host_%s_x%d_y%d_z%d.dat",hname,dim[0],dim[1],dim[2]);
      #endif
        printf("The FFT File will be %s\n",fname);
       
        fftwf_plan fft_plan;
        fftwf_plan ifft_plan;
        
        fftwf_plan fft_plan_temp;
        fftwf_plan ifft_plan_temp;
        
      if(  (fid=fopen(fname,"r")) != NULL){
        fftwf_import_wisdom_from_file(fid);
        fclose(fid);
      }	
      
      int imsize = dim[0]*dim[1]*dim[2];
      
      // Plan the transforms in temp arrays so that data isnt overwritten
      complex<float> *temp = new complex<float>[imsize];
      fft_plan_temp =  fftwf_plan_dft_3d(dim[2], dim[1], dim[0],(fftwf_complex *)temp,(fftwf_complex*)temp,FFTW_FORWARD, FFTW_MEASURE);
      ifft_plan_temp = fftwf_plan_dft_3d(dim[2], dim[1], dim[0],(fftwf_complex *)temp,(fftwf_complex*)temp,FFTW_BACKWARD, FFTW_MEASURE);
      delete[] temp;
      
      int ntransforms = 1;
      for (int i = 3; i < ndims; i++) {ntransforms*=dim[i];}
      
      fftshift3();
      
      for (int i = 0; i < ntransforms; i++) {
        fft_plan =  fftwf_plan_dft_3d(dim[2], dim[1], dim[0],(fftwf_complex *)(vals+i*imsize),(fftwf_complex*)(vals+i*imsize),FFTW_FORWARD, FFTW_MEASURE);
        ifft_plan = fftwf_plan_dft_3d(dim[2], dim[1], dim[0],(fftwf_complex *)(vals+i*imsize),(fftwf_complex*)(vals+i*imsize),FFTW_BACKWARD, FFTW_MEASURE);
        if(dir==1) {fftwf_execute(fft_plan);}
        else if(dir==-1) {fftwf_execute(ifft_plan);}
      }
      
      fftshift3();
      
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
     }
		 
     void fftshift3()
    {
      int ntransforms = 1;
      for (int i = 3; i < ndims; i++) {ntransforms*=dim[i];}
      
      #pragma omp parallel for
      for (int kk = 0; kk < ntransforms; kk++) {
        for (int k = 0; k < dim[2]; k++) {
          for (int j = 0; j < dim[1]; j++) {
            for (int i = 0; i < dim[0]; i++) {
              int mod = ((i+j+k)%2) == 0 ? 1 : -1;
              (*this)(kk,k,j,i) = (float)mod*(*this)(kk,k,j,i);
      }}}}
      return;
    }
     
		 // Write binary file		 
		 void write(char *name){
		 	FILE *fp=fopen(name,"w");
     		fwrite(vals, numel, sizeof(C), fp);
     		fclose(fp); 
		 }
		
	private:
	
};

