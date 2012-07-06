// Array Class for MRI Reconstruction
//   Used to allow very specific control over memory
//   and parallelization of key operations 
// Main Classes: 
//		Array2D,Array3D	 - With contigous memory access
//		Array4D,Array5D	 - Non-contigous memory access (3D are contigous)

#ifndef hARRAY
#define hARRAY

#define MAXDIMS 16

using namespace std;

template < class C >
class arrayND{ 
	public:	
		C *vals; 
		
		int ndims;
		int dim[MAXDIMS];
		size_t numel;
		
		int MemExists;
		
		/*---------------------------------------------
		 *    Constructor / Decontructor / etc
		 *---------------------------------------------*/
		
		// Constructor
		arrayND< C >(){
			MemExists = 0;
			memset(dim,0,(size_t)(sizeof(int)*MAXDIMS));
		}
		
		// Allocation with Size
		void alloc(int z,int y,int x){
			ndims = 3;
			dim[0] = x;
			dim[1] = y;
			dim[2] = z;
			
			numel = x*y*z;
			
			vals = alloc_vals();
			MemExists =1;
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
		 
		 // Sum of X^2
		 double energy (void){
		  
		  int nthreads = omp_get_num_threads();
		 	double *temp = new double[nthreads];
		 	memset(temp,0,(size_t)(sizeof(double)*nthreads));
		 	
		 	#pragma omp parallel for 		 
			for(int i=0; i<numel; i++){
			  int tid = omp_get_thread_num();
			  temp[tid] += (double)norm(vals[i]);		
			}
			
			double X=0.0;
			for(int i=0; i<nthreads; i++){
				X += temp[i];
			}	
			X *= 1.0/numel;
			
			delete []temp;
		 	return(X);
			
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
		 void operator = (array3D<C>temp){
		 	check_array_dims(temp);
			
			#pragma omp parallel for 
			for(int i=0; i<numel; i++){
				vals[i]=temp.vals[i];					
			} 
		 }
		 
		  // Sum of X
		 C max( void ){
			C *temp = new C[Nz];
			
			#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
				// Y Energy
				C jM =0;
				for(int j=0; j<Ny; j++){
					for(int i=0; i<Nx; i++){
						jM = ( abs(jM) > abs(vals[k][j][i])) ? ( jM ) : ( vals[k][j][i]);			
					}	
				}
				temp[k] = jM; // Add energy for whole slice
			}
			
			C X=0.0;
			for(int k=0; k<Nz; k++){
				X= ( abs(temp[k]) > abs(X) ) ? ( temp[k] ) : ( X );
			}		 
		 	return(X);		
		 } 
		 
		 // Sum of X
		 C sum( void ){
			C *temp = new C[Nz];
			
			#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
				float scale = 1.0/( (float)(Nx*Ny*Nz));
				// Y Energy
				C jE =0;
				for(int j=0; j<Ny; j++){
					// X Energy
					C iE =0;
					for(int i=0; i<Nx; i++){
						iE += vals[k][j][i];			
					}	
					jE+=iE*scale; // Add energy for whole row
				}
				temp[k] = jE; // Add energy for whole slice
			}
			
			C X=0.0;
			for(int k=0; k<Nz; k++){
				X+= temp[k];
			}		 
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
		 
		 void check_array_dims( array3D<C>temp ){
		 	if( temp.Nx != Nx){
				cout << "Conjugate Multi:Conflicting Dims " << Nx << " != " << temp.Nx << endl;
			}
		 	
			if( temp.Ny != Ny){
				cout << "Conjugate Multi:Conflicting Dims " << Ny << " != " << temp.Ny << endl;
			}
		 	
			if( temp.Nz != Nz){
				cout << "Conjugate Multi:Conflicting Dims " << Nz << " != " << temp.Nz << endl;
			}
		 }
		 
		 void conjugate_multiply(array3D<C>temp){ 
		 	
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
		 		 		 		 
		 // Write binary file		 
		 void write(char *name){
		 	FILE *fp=fopen(name,"w");
     		fwrite(vals, numel, sizeof(C), fp);
     		fclose(fp); 
		 }
		
		 
		 /* For Maximum calcs (prevent overloaded template)
		float Xmax( float A, float B){	return( (A>B) ? (A) : (B)); }
		double Xmax( double A, double B){	return( (A>B) ? (A) : (B)); }
		short Xmax( short A, short B){	return( (A>B) ? (A) : (B)); }
		int Xmax( int A, int B){	return( (A>B) ? (A) : (B)); }
		complex<float> Xmax( complex<float> A, complex<float> B){ return( (abs(A)>abs(B)) ? (A) : (B)); }
		complex<double> Xmax( complex<double> A, complex<double> B){ return( (abs(A)>abs(B)) ? (A) : (B)); }
		complex<short> Xmax( complex<short> A, complex<short> B){ return( (abs(A)>abs(B)) ? (A) : (B)); }
		complex<int> Xmax( complex<int> A, complex<int> B){ return( (abs(A)>abs(B)) ? (A) : (B)); }
		
		  Maximum 
		 C max(void){
		 	C m=vals[0][0][0];
			for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				m = Xmax(m,vals[k][j][i]);
			}}}
		 	return(m);
		 }*/
		
	private:
	
};



#endif


