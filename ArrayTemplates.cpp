// Array Class for MRI Reconstruction
//   Used to allow very specific control over memory
//   and parallelization of key operations 
// Main Classes: 
//		Array2D,Array3D	 - With contigous memory access
//		Array4D,Array5D	 - Non-contigous memory access (3D are contigous)

#ifndef hARRAY
#define hARRAY

using namespace std;

template < class C >
class array3D{ 
	public:	
		C ***vals; 
		int Nx;
		int Ny;
		int Nz;
		int Numel;
		int MemExists;
		
		/*---------------------------------------------
		 *    Constructor / Decontructor / etc
		 *---------------------------------------------*/
		
		// Constructor
		array3D< C >(){
			MemExists = 0;
		}
		
		// Allocation with Size
		void alloc(int z,int y,int x){
			Nz = z;
			Ny = y;
			Nx = x;
			Numel = x*y*z;
			
			// Allocate 3D Data arbitrary type 
			// cout << "Alloc " << Sz << " x " << Sy << " x " << Sx << endl;
			vals = alloc3D();
			MemExists =1;
		}
		
		// Base Code to Allocate a Contigous Memory Block	
		C ***alloc3D( void ){
			C*** temp;
			
			temp = new C**[Nz];
			temp[0] = new C*[Nz*Ny];
			temp[0][0] = new C[Nz*Ny*Nx];
			for(int j=1; j<Ny;j++){
				temp[0][j]= temp[0][j-1] + Nx;
			}
			for(int k=1; k<Nz;k++){
				temp[k]= temp[k-1] + Ny;
				temp[k][0]= temp[k-1][0] + Ny*Nx;
				for(int j=1;j<Ny;j++){
					temp[k][j] = temp[k][j-1] + Nx;
				}
			}
			memset(temp[0][0],0,(size_t)(sizeof(C)*Nz*Ny*Nx));
			return(temp);		 
		}	
			
		// Cleans Memory	
		void freeArray(){
		    if(MemExists){
				free(vals[0][0]);
    			free(vals[0]);
    			free(vals);
				// cout << "Destroying" << endl;
			}
    	}
		
		void point_to_3D( array3D< C > *temp ){
			Nz = temp->Nz;
			Ny = temp->Ny;
			Nx = temp->Nx;
			vals = temp->vals;
		}
		
		/*---------------------------------------------
		 *    Utilities - Operators
		 *---------------------------------------------*/
		 
		 // Allow Indexing with arrayND M[k][j][i] (suppose to be bad but useful for coding)
		 inline C** & operator[](const int i){
		 	return( vals[i]);		 
		 }
		 
		 inline C & operator ()(const int i){
		 	return( vals[0][0][i]);		 
		 }
		 
		 
		 // Equals Sets Value		 
		 void operator = (C temp){
		 	#pragma omp parallel for 
			for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i]=temp;					
			}}}		 
		 }
		 
		 // Scaling overloaded	
		 void operator *= (float X){
		 	#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i] *= X;			
			}}}		 
		 }
		 
		 // Scaling overloaded	
		 void operator *= (complex<float> XX){
		 	#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i] *= XX;			
			}}}		 
		 }
		 
		 // Sum of X^2
		 double energy (void){
		 	double *temp = new double[Nz];
		 	
			#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
				
				float scale = 1.0/( (float)(Nx*Ny*Nz));
				
				// Y Energy
				double jE =0;
				for(int j=0; j<Ny; j++){
					// X Energy
					double iE =0;
					for(int i=0; i<Nx; i++){
						iE += (double)norm(vals[k][j][i]);			
					}	
					jE+=iE*scale;
				}
				temp[k] = jE; // Add energy for whole slice [race condition]
			}
			
			double X=0.0;
			for(int k=0; k<Nz; k++){
				X+= temp[k];
			}		 
			delete []temp;
		 	return(X);
		 } 
		 
		 // Point Wise Multiply of Same Type
		 template < class X>		 		 
		 void operator *= (X temp){
		 	check_array_dims(temp);
			
			#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i] *= temp.vals[k][j][i];					
			}}}		 
		 }
		 
		 // Set Array to Zero Fast
		 void zero( void ){
		 	memset(vals[0][0],0,sizeof(C)*Nx*Ny*Nz);
		 }
		 
		 // Pointwise equals 
		 void operator = (array3D<C>temp){
		 	check_array_dims(temp);
			
			#pragma omp parallel for 
			for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i]=temp.vals[k][j][i];					
			}}}		 
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
		 
		  // min of X
		 C min( void ){
			C *temp = new C[Nz];
			
			#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
				// Y Energy
				C jM =0;
				for(int j=0; j<Ny; j++){
					for(int i=0; i<Nx; i++){
						jM = ( abs(jM) < abs(vals[k][j][i])) ? ( jM ) : ( vals[k][j][i]);			
					}	
				}
				temp[k] = jM; // Add energy for whole slice
			}
			
			C X=0.0;
			for(int k=0; k<Nz; k++){
				X= ( abs(temp[k]) < abs(X) ) ? ( temp[k] ) : ( X );
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
			for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i] *=vals[k][j][i];					
			}}}			 
		 }
		 		 		 
		 // Doesn't work/make sense for complex
		 void sqrt( void ){
		 	#pragma omp parallel for 
			for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i]=sqrtf(vals[k][j][i]);					
			}}}			 
		 }
		 
		 // Only works for complex 
		 void complex_conj( void ){
		 	#pragma omp parallel for 
			for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i]=conj(vals[k][j][i]);					
			}}}			 
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
		 	for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i] *= conj(temp.vals[k][j][i]);					
			}}}		
		 
		 }
		 
		 // For Add calcs 
		 template < class X>		 		 
		 void operator += (X temp){
		 	check_array_dims(temp);
			
			#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i] += temp.vals[k][j][i];					
			}}}		 
		 }
		 
		 template < class X>		 		 
		 void operator -= (X temp){
		 
		 	check_array_dims(temp);
			
		 	#pragma omp parallel for 		 
		 	for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			for(int i=0; i<Nx; i++){
				vals[k][j][i] -=  temp.vals[k][j][i];					
			}}}		 
		 }
		 		 		 		 		 
		 // Write binary file		 
		 void write(char *name){
		 	FILE *fp=fopen(name,"w");
     		fwrite( vals[0][0],Nx*Ny*Nz, sizeof(C), fp);
     		fclose(fp); 
		 }
		 
		 // Write binary file (magnitude)		 
		 void write_mag(char *name){
		 	float *temp = new float[Nx];
			
			FILE *fp=fopen(name,"w");
     		for(int k=0; k<Nz; k++){
			for(int j=0; j<Ny; j++){
			
				for(int i=0; i<Nx; i++){
					temp[i] = abs(vals[k][j][i]);	
				}
			
				fwrite( temp,Nx, sizeof(float), fp);
     		}}
			fclose(fp); 
		 }
		
		 // Write binary file (magnitude, slice)		 
		 void write_mag(char *name,int slice,const char *type){
		 	float *temp = new float[Nx];
			
			FILE *fp=fopen(name,type);
     		for(int j=0; j<Ny; j++){
				for(int i=0; i<Nx; i++){
					temp[i] = abs(vals[slice][j][i]);	
				}
				fwrite( temp,Nx, sizeof(float), fp);
     		}
			fclose(fp);
			delete []temp; 
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

// 4D is really a container of 3D arrays
template < class C >
class array4D{ 
	public:	
		array3D<C> *vals; 
		int Nx;
		int Ny;
		int Nz;
		int Nt;
		int Numel;
		int MemExists;
		
		/*---------------------------------------------
		 *    Constructor / Decontructor / etc
		 *---------------------------------------------*/
		
		// Constructor
		array4D< C >(){
			MemExists = 0;
		}
				
		// Allocation with Size
		void alloc(int t,int z,int y,int x){
			Nz = z;
			Ny = y;
			Nx = x;
			Nt = t;
			Numel=t*z*y*x;
			vals = new array3D<C>[Nt];
			for(t=0; t<Nt; t++){
				vals[t].alloc(Nz,Ny,Nx);			
			}
		}
				
		// Deconstructor		
		void freeArray(){
		    if(MemExists){
				for(int i=0; i<Nt;i++){
					vals[i].freeArray();
				}
				delete [] vals;
			}
    	}
		
		void point_to_3D( array3D<C> *temp){
			Nt=1;
			Nz=temp->Nz;
			Ny=temp->Ny;
			Nx=temp->Nx;
			vals = new array3D<C>[Nt];
			for(int t=0; t<Nt; t++){
				vals[t].point_to_3D( temp);
			}
		}
		
		void point_to_4D( array4D<C> *temp){
			Nt=temp->Nt;
			Nz=temp->Nz;
			Ny=temp->Ny;
			Nx=temp->Nx;
			vals = temp->vals;
		}
		
		void samesize( array4D<C> *temp){
			Nt=temp->Nt;
			Nz=temp->Nz;
			Ny=temp->Ny;
			Nx=temp->Nx;
			Nt=temp->Nt;
			vals = new array3D<C>[Nt];
			for(int t=0; t<Nt; t++){
				vals[t].alloc(Nz,Ny,Nx);			
			}
		}
		
		
		/*---------------------------------------------
		 *    Utilities - Operators
		 *---------------------------------------------*/
		 
		 // Allow Indexing with arrayND M[k][j][i] (suppose to be bad but useful for FFT coding)
		 inline array3D<C> & operator[](const int i){
		 	return(vals[i]);		 
		 }
		 
		 inline C & operator ()(const int i){
		 	int tt = (int)((float)i / (float)vals[0].Numel);
			int ii = i % vals[0].Numel;
			return( vals[tt](ii) );		 
		 }
		 
		 void operator = (array4D<C>temp){
		 	if( temp.Nt != Nt){
				cout << "Add:Conflicting Dims " << Nt << " != " << temp.Nt << endl;
			}
			for(int t=0; t<Nt; t++){
				vals[t]=temp.vals[t];					
			}		 
		 }
		 
		 void zero(void){
		 	for(int t=0; t<Nt; t++){
				vals[t].zero();					
			}		 
		 }
		
		 
		void conjugate_multiply(array4D<C>temp){ 
			if( temp.Nt != Nt){
				cout << "Add:Conflicting Dims " << Nt << " != " << temp.Nt << endl;
			}
			for(int t=0; t<Nt; t++){
				vals[t].conjugate_multiply(temp.vals[t]);					
			}		 
		 }
		
		
 		void operator -= (array4D<C>temp){
		 	if( temp.Nt != Nt){
				cout << "Add:Conflicting Dims " << Nt << " != " << temp.Nt << endl;
			}
			for(int t=0; t<Nt; t++){
				vals[t]-=temp.vals[t];					
			}		 
		 }
		 
		 double energy( void){
		 	double temp=0;
			for(int t=0; t<Nt; t++){
				temp+= vals[t].energy();  					
			}
			return(temp);		 
		 }
		 
		 C sum( void){
		 	C temp=0;
			for(int t=0; t<Nt; t++){
				temp+= vals[t].sum();  					
			}
			return(temp);		 
		 }
		 
		 C max( void){
		 	C m=0;
			for(int t=0; t<Nt; t++){
				C temp = vals[t].max();
				m = ( abs(m) > abs(temp) ) ? ( m ) : ( temp );					
			}
			return(m);		 
		 }
		
		C min( void){
		 	C m=0;
			for(int t=0; t<Nt; t++){
				C temp = vals[t].min();
				m = ( abs(m) < abs(temp) ) ? ( m ) : ( temp );					
			}
			return(m);		 
		 }
		 
		 
		void operator += (array4D<C>temp){
		 	if( temp.Nt != Nt){
				cout << "Add:Conflicting Dims " << Nt << " != " << temp.Nt << endl;
			}
			for(int t=0; t<Nt; t++){
				vals[t]+=temp.vals[t];					
			}		 
		 }
		 
		 void operator *= (array4D<C>temp){
		 	if( temp.Nt != Nt){
				cout << "Add:Conflicting Dims " << Nt << " != " << temp.Nt << endl;
			}
			for(int t=0; t<Nt; t++){
				vals[t]*=temp.vals[t];					
			}		 
		 }
		 
		 void operator *= (C temp){
		 	for(int t=0; t<Nt; t++){
				vals[t]*=temp;					
			}		 
		 }
				  
	private:
	
};

// 5D is really a container of 4D arrays
template < class C >
class array5D{ 
	public:	
		array4D<C> *vals; 
		int MemExists;
		
		int Nx;
		int Ny;
		int Nz;
		int Nt;
		int Ne;
		int Numel;
		/*---------------------------------------------
		 *    Constructor / Decontructor / etc
		 *---------------------------------------------*/
		
		// Constructor
		array5D< C >(){
			MemExists = 0;
		}
		
		// Allocation with Size
		void alloc(int e,int t,int z,int y,int x){
			Ne=e;
			Nt=t;
			Nz=z;
			Ny=y;
			Nx=x;
			Numel=e*t*z*y*x;
			vals = new array4D<C>[Ne];
			for(t=0; t<Ne; t++){
				vals[t].alloc(Nt,Nz,Ny,Nx);
			}
		}
				
		// Deconstructor		
		void freeArray(){
		    if(MemExists){
				for(int i=0; i<Ne;i++){
					vals[i].freeArray();
				}
				delete [] vals;
			}
    	}
		
		void point_to_3D( array3D<C> *temp){
			Ne=1;
			Nt=1;
			Nz=temp->Nz;
			Ny=temp->Ny;
			Nx=temp->Nx;
			vals = new array4D<C>[Ne];
			for(int t=0; t<Ne; t++){
				vals[t].point_to_3D( temp);
			}
					
		}
		
		void point_to_4D( array4D<C> *temp){
			Ne=1;
			Nt=temp->Nt;
			Nz=temp->Nz;
			Ny=temp->Ny;
			Nx=temp->Nx;
			vals = new array4D<C>[Ne];
			for(int t=0; t<Ne; t++){
				vals[t].point_to_4D( temp );
			}
		}
		
		
		/*---------------------------------------------
		 *    Utilities - Operators
		 *---------------------------------------------*/
		 
		 // Allow Indexing with arrayND M[k][j][i] (suppose to be bad but useful for FFT coding)
		 inline array4D<C> & operator[](const int i){
		 	return(vals[i]);		 
		 }
		 
		 inline C & operator ()(const int i){
		 	int tt = (int)((float)i / (float)vals[0].Numel);
			int ii = i % vals[0].Numel;
			return( vals[tt](ii) );		 
		 }
		 
		 void operator = (array5D<C>temp){
		 	if( temp.Ne !=  Ne){
				cout << "Add:Conflicting Dims " << Ne << " != " << temp.Ne << endl;
			}
			for(int t=0; t< Ne; t++){
				vals[t]=temp.vals[t];					
			}		 
		 }
		 
		 void zero(void){
		 	for(int t=0; t<Ne; t++){
				vals[t].zero();					
			}		 
		 }
		
		 
		void conjugate_multiply(array5D<C>temp){ 
			if( temp.Ne != Ne){
				cout << "Add:Conflicting Dims " << Ne << " != " << temp.Ne << endl;
			}
			for(int t=0; t<Ne; t++){
				vals[t].conjugate_multiply(temp.vals[t]);					
			}		 
		 }
		
		
 		void operator -= (array5D<C>temp){
		 	if( temp.Ne != Ne){
				cout << "Add:Conflicting Dims " << Ne << " != " << temp.Ne << endl;
			}
			for(int t=0; t<Ne; t++){
				vals[t]-=temp.vals[t];					
			}		 
		 }
		 
		 double energy( void){
		 	double temp=0;
			for(int t=0; t<Ne; t++){
				temp+= vals[t].energy();  					
			}
			return(temp);		 
		 }
		 
		 C sum( void){
		 	C temp=0;
			for(int t=0; t<Ne; t++){
				temp+= vals[t].sum();  					
			}
			return(temp);		 
		 }
		 
		 C max( void){
		 	C m=0;
			for(int t=0; t<Ne; t++){
				C temp = vals[t].max();
				m = ( abs(m) > abs(temp) ) ? ( m ) : ( temp );					
			}
			return(m);		 
		 }
		 
		 C min( void){
		 	C m=0;
			for(int t=0; t<Ne; t++){
				C temp = vals[t].max();
				m = ( abs(m) < abs(temp) ) ? ( m ) : ( temp );					
			}
			return(m);		 
		 }
		 
		void operator += (array5D<C>temp){
		 	if( temp.Ne != Ne){
				cout << "Add:Conflicting Dims " << Ne << " != " << temp.Ne << endl;
			}
			for(int t=0; t<Ne; t++){
				vals[t]+=temp.vals[t];					
			}		 
		 }
		 
		 void operator *= (array5D<C>temp){
		 	if( temp.Ne != Ne){
				cout << "Add:Conflicting Dims " << Ne << " != " << temp.Ne << endl;
			}
			for(int t=0; t<Ne; t++){
				vals[t]*=temp.vals[t];					
			}		 
		 }
		 
		 void operator *= (C temp){
		 	for(int t=0; t<Ne; t++){
				vals[t]*=temp;					
			}		 
		 }
		
	private:
	
};


#endif


