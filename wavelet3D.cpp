#include "wavelet3D.h"

using namespace std;
using namespace NDarray;

WAVELET3D::WAVELET3D(){
 
}


WAVELET3D::~WAVELET3D(){
 	// Not here yet
}

/**
 * Constructor for 5D Wavelet. Will only do 3D wavelets.
 * @see setup()
 * @param temp Array (or array of same size) to be transformed
 * @param l Array listing what type of transform along each dim
 * @param type The actual wavelet to be used.
 */
WAVELET3D::WAVELET3D( Array< complex<float>, 5 >& temp, int *l, int type){
	Array< complex<float>,3>temp3 = temp(Range::all(),Range::all(),Range::all(),0,0);
	setup(  temp3, l, type);		
}

/**
 * Constructor for 4D Wavelet. Will only do 3D wavelets.
 * @see setup()
 * @param temp Array (or array of same size) to be transformed
 * @param l Array listing what type of transform along each dim
 * @param type The actual wavelet to be used.
 */
WAVELET3D::WAVELET3D( Array< complex<float>, 4 >& temp, int *l, int type){
	Array< complex<float>,3>temp3 = temp(Range::all(),Range::all(),Range::all(),0);
	setup(  temp3, l, type);		
}

/**
 * Constructor for 3D Wavelet.
 * @see setup()
 * @param temp Array (or array of same size) to be transformed
 * @param l Array listing what type of transform along each dim
 * @param type The actual wavelet to be used.
 */
WAVELET3D::WAVELET3D( Array< complex<float>, 3 >& temp, int *l, int type){
	setup(  temp, l, type);
}


/**
 * Constructor for 3D Wavelet.
 * @see setup()
 * @param nn size of array to be transformd
 * @param ll levels at each dimension
 * @param type The actual wavelet to be used.
 */
WAVELET3D::WAVELET3D( TinyVector<int,3> nn, TinyVector<int,3>  ll, int type){
	N[2] = nn(2);
	N[1] = nn(1);
	N[0] = nn(0);
	
	L[2] = ll(2);
	L[1] = ll(1);
	L[0] = ll(0);

	wType = type;
			
	get_filter_banks();
	setup_wavelet_directions();
}



/**
 * Code to setup a wavelet transform based on 3D Array
 * @see get_filter_banks()
 * @see setup_wavelet_directions()
 * @param temp 3D Array to set allowable levels
 * @param l Array listing what level of transform along each dim
 * @param type The actual wavelet to be used.
 */
void WAVELET3D::setup( Array< complex<float>, 3 >& temp, int *l, int type){
	
	/// cout << "Init Wavelet " << endl;
	N[2] = temp.length(thirdDim);
	N[1] = temp.length(secondDim);
	N[0] = temp.length(firstDim);
	
	L[2] = l[2];
	L[1] = l[1];
	L[0] = l[0];

	wType = type;
			
	get_filter_banks();
	setup_wavelet_directions();
}

/**
 * Code to setup a wavelet transform based on 3D Array
 * @see setup()
 * This code check to ensure the number of levels is not greater than that allowed by the wavelet basis.
 */
void WAVELET3D::setup_wavelet_directions(void){
	cout << "Setting up dimensions to do wavelet transform" << endl;
	max_level = 0;
	for( int dir=0; dir<3; dir++){
		if(floor(N[dir]/pow(2.0,L[dir])) != (N[dir]/pow(2.0,L[dir]))){
			int max_levels = (int)log2( (float)N[dir]/(float)(wN-1));
			cout << "\tError: Can't perform " << L[dir] << "level wavelet for dir" << dir << " N=" << N[dir] << endl;
			L[dir] = ( max_levels <= 0 ) ? ( 0 ) : ( max_levels);
			cout << "\tUsing " << L[dir] << "levels instead" << endl;
		}
		max_level = ( max_level > L[dir] ) ? ( max_level ) : ( L[dir] );
	}
	cout << "Max Level = " << max_level << endl;
}


/**
 * Code to setup a wavelet basis
 * @see setup()
 * This copies discrete wavelet basis made in matlab to structures in class
 */
void WAVELET3D::get_filter_banks( void ){
	
	//  Get the filter banks for wavelet (Coef are generated from matlab script)
	S[0]=0;
	S[1]=0;
	S[2]=0;
	
	switch(wType){
		case(WAVE_DB2):{
			wN = 2; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
			
			lpf[0]=0.7071067811865475700000;
			lpf[1]=0.7071067811865475700000;
			hpf[0]=-0.7071067811865475700000;
			hpf[1]=0.7071067811865475700000;
			Slpf[0]=0.7071067811865475700000;
			Slpf[1]=0.7071067811865475700000;
			Shpf[0]=0.7071067811865475700000;
			Shpf[1]=-0.7071067811865475700000;
			
		}break;
		
		default:
		case(WAVE_DB4):{
			wN = 4; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
			
			lpf[3]=-0.1294095225509214500000;
			lpf[2]=0.2241438680418573500000;
			lpf[1]=0.8365163037374689900000;
			lpf[0]=0.4829629131446902500000;
			
			hpf[3]= 0.4829629131446902500000;
			hpf[2]=-0.8365163037374689900000;
			hpf[1]= 0.2241438680418573500000;
			hpf[0]= 0.1294095225509214500000;
			
			Slpf[3]= 0.4829629131446902500000;
			Slpf[2]= 0.8365163037374689900000;
			Slpf[1]= 0.2241438680418573500000;
			Slpf[0]=-0.1294095225509214500000;
			
			Shpf[3]=0.1294095225509214500000;
			Shpf[2]=0.2241438680418573500000;
			Shpf[1]=-0.8365163037374689900000;
			Shpf[0]=0.4829629131446902500000;
			
		}break;
			
		case(WAVE_DB6):{
			wN = 6; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
	
			lpf[0]=0.0352262918821006560000;
			lpf[1]=-0.0854412738822414860000;
			lpf[2]=-0.1350110200103908400000;
			lpf[3]=0.4598775021193313200000;
			lpf[4]=0.8068915093133387500000;
			lpf[5]=0.3326705529509568800000;
			hpf[0]=-0.3326705529509568800000;
			hpf[1]=0.8068915093133387500000;
			hpf[2]=-0.4598775021193313200000;
			hpf[3]=-0.1350110200103908400000;
			hpf[4]=0.0854412738822414860000;
			hpf[5]=0.0352262918821006560000;
			Slpf[0]=0.3326705529509568800000;
			Slpf[1]=0.8068915093133387500000;
			Slpf[2]=0.4598775021193313200000;
			Slpf[3]=-0.1350110200103908400000;
			Slpf[4]=-0.0854412738822414860000;
			Slpf[5]=0.0352262918821006560000;
			Shpf[0]=0.0352262918821006560000;
			Shpf[1]=0.0854412738822414860000;
			Shpf[2]=-0.1350110200103908400000;
			Shpf[3]=-0.4598775021193313200000;
			Shpf[4]=0.8068915093133387500000;
			Shpf[5]=-0.3326705529509568800000;
	
		}break;
		
		case(WAVE_DB8):{
			wN = 8; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
			
			lpf[0]=-0.0105974017849972780000;
			lpf[1]=0.0328830116669829450000;
			lpf[2]=0.0308413818359869650000;
			lpf[3]=-0.1870348117188811400000;
			lpf[4]=-0.0279837694169838490000;
			lpf[5]=0.6308807679295903600000;
			lpf[6]=0.7148465705525415300000;
			lpf[7]=0.2303778133088552300000;
			hpf[0]=-0.2303778133088552300000;
			hpf[1]=0.7148465705525415300000;
			hpf[2]=-0.6308807679295903600000;
			hpf[3]=-0.0279837694169838490000;
			hpf[4]=0.1870348117188811400000;
			hpf[5]=0.0308413818359869650000;
			hpf[6]=-0.0328830116669829450000;
			hpf[7]=-0.0105974017849972780000;
			Slpf[0]=0.2303778133088552300000;
			Slpf[1]=0.7148465705525415300000;
			Slpf[2]=0.6308807679295903600000;
			Slpf[3]=-0.0279837694169838490000;
			Slpf[4]=-0.1870348117188811400000;
			Slpf[5]=0.0308413818359869650000;
			Slpf[6]=0.0328830116669829450000;
			Slpf[7]=-0.0105974017849972780000;
			Shpf[0]=-0.0105974017849972780000;
			Shpf[1]=-0.0328830116669829450000;
			Shpf[2]=0.0308413818359869650000;
			Shpf[3]=0.1870348117188811400000;
			Shpf[4]=-0.0279837694169838490000;
			Shpf[5]=-0.6308807679295903600000;
			Shpf[6]=0.7148465705525415300000;
			Shpf[7]=-0.2303778133088552300000;
			
		}break;
		
		case(WAVE_SYM2):{
			wN = 2; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
			
			lpf[0]=0.7071067811865475700000;
			lpf[1]=0.7071067811865475700000;
			hpf[0]=-0.7071067811865475700000;
			hpf[1]=0.7071067811865475700000;
			Slpf[0]=0.7071067811865475700000;
			Slpf[1]=0.7071067811865475700000;
			Shpf[0]=0.7071067811865475700000;
			Shpf[1]=-0.7071067811865475700000;
			
		}break;
		
		case(WAVE_SYM4):{
			wN = 4; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
			
			lpf[0]=-0.1294095225509214500000;
			lpf[1]=0.2241438680418573500000;
			lpf[2]=0.8365163037374689900000;
			lpf[3]=0.4829629131446902500000;
			hpf[0]=-0.4829629131446902500000;
			hpf[1]=0.8365163037374689900000;
			hpf[2]=-0.2241438680418573500000;
			hpf[3]=-0.1294095225509214500000;
			Slpf[0]=0.4829629131446902500000;
			Slpf[1]=0.8365163037374689900000;
			Slpf[2]=0.2241438680418573500000;
			Slpf[3]=-0.1294095225509214500000;
			Shpf[0]=-0.1294095225509214500000;
			Shpf[1]=-0.2241438680418573500000;
			Shpf[2]=0.8365163037374689900000;
			Shpf[3]=-0.4829629131446902500000;
			
		}break;
			
		case(WAVE_SYM6):{
			wN = 6; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
	
			lpf[0]=0.0352262918821006560000;
			lpf[1]=-0.0854412738822414860000;
			lpf[2]=-0.1350110200103908400000;
			lpf[3]=0.4598775021193313200000;
			lpf[4]=0.8068915093133387500000;
			lpf[5]=0.3326705529509568800000;
			hpf[0]=-0.3326705529509568800000;
			hpf[1]=0.8068915093133387500000;
			hpf[2]=-0.4598775021193313200000;
			hpf[3]=-0.1350110200103908400000;
			hpf[4]=0.0854412738822414860000;
			hpf[5]=0.0352262918821006560000;
			Slpf[0]=0.3326705529509568800000;
			Slpf[1]=0.8068915093133387500000;
			Slpf[2]=0.4598775021193313200000;
			Slpf[3]=-0.1350110200103908400000;
			Slpf[4]=-0.0854412738822414860000;
			Slpf[5]=0.0352262918821006560000;
			Shpf[0]=0.0352262918821006560000;
			Shpf[1]=0.0854412738822414860000;
			Shpf[2]=-0.1350110200103908400000;
			Shpf[3]=-0.4598775021193313200000;
			Shpf[4]=0.8068915093133387500000;
			Shpf[5]=-0.3326705529509568800000;
	
		}break;
		
		case(WAVE_SYM8):{
			wN = 8; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
			
			lpf[0]=-0.0757657147892733250000;
			lpf[1]=-0.0296355276459985100000;
			lpf[2]=0.4976186676320154500000;
			lpf[3]=0.8037387518059161400000;
			lpf[4]=0.2978577956052773600000;
			lpf[5]=-0.0992195435768472160000;
			lpf[6]=-0.0126039672620378330000;
			lpf[7]=0.0322231006040427020000;
			
			hpf[0]=-0.0322231006040427020000;
			hpf[1]=-0.0126039672620378330000;
			hpf[2]=0.0992195435768472160000;
			hpf[3]=0.2978577956052773600000;
			hpf[4]=-0.8037387518059161400000;
			hpf[5]=0.4976186676320154500000;
			hpf[6]=0.0296355276459985100000;
			hpf[7]=-0.0757657147892733250000;
			
			Slpf[0]=0.0322231006040427020000;
			Slpf[1]=-0.0126039672620378330000;
			Slpf[2]=-0.0992195435768472160000;
			Slpf[3]=0.2978577956052773600000;
			Slpf[4]=0.8037387518059161400000;
			Slpf[5]=0.4976186676320154500000;
			Slpf[6]=-0.0296355276459985100000;
			Slpf[7]=-0.0757657147892733250000;
			
			Shpf[0]=-0.0757657147892733250000;
			Shpf[1]=0.0296355276459985100000;
			Shpf[2]=0.4976186676320154500000;
			Shpf[3]=-0.8037387518059161400000;
			Shpf[4]=0.2978577956052773600000;
			Shpf[5]=0.0992195435768472160000;
			Shpf[6]=-0.0126039672620378330000;
			Shpf[7]=-0.0322231006040427020000;
			
		}break;
		
		case(WAVE_BO33):{
			wN = 8; 
			lpf=new float[wN];
			hpf=new float[wN];
			Slpf=new float[wN];
			Shpf=new float[wN];
			
			lpf[0]=0.0662912607362388380000;
			lpf[1]=-0.1988737822087165200000;
			lpf[2]=-0.1546796083845572700000;
			lpf[3]=0.9943689110435824900000;
			lpf[4]=0.9943689110435824900000;
			lpf[5]=-0.1546796083845572700000;
			lpf[6]=-0.1988737822087165200000;
			lpf[7]=0.0662912607362388380000;
			
			hpf[0]=-0.0000000000000000000000;
			hpf[1]=0.0000000000000000000000;
			hpf[2]=-0.1767766952966368900000;
			hpf[3]=0.5303300858899107100000;
			hpf[4]=-0.5303300858899107100000;
			hpf[5]=0.1767766952966368900000;
			hpf[6]=-0.0000000000000000000000;
			hpf[7]=0.0000000000000000000000;
			
			// Reverse biorthogonal33 for transpose not inverse		
			Slpf[0]=0.0662912607362388380000;
			Slpf[1]=-0.1988737822087165200000;
			Slpf[2]=-0.1546796083845572700000;
			Slpf[3]=0.9943689110435824900000;
			Slpf[4]=0.9943689110435824900000;
			Slpf[5]=-0.1546796083845572700000;
			Slpf[6]=-0.1988737822087165200000;
			Slpf[7]=0.0662912607362388380000;
			
			Shpf[0]=0.0000000000000000000000;
			Shpf[1]=-0.0000000000000000000000;
			Shpf[2]=0.1767766952966368900000;
			Shpf[3]=-0.5303300858899107100000;
			Shpf[4]=0.5303300858899107100000;
			Shpf[5]=-0.1767766952966368900000;
			Shpf[6]=0.0000000000000000000000;
			Shpf[7]=-0.0000000000000000000000;
		}break;
	}
}

/**
 * Code to randomize shifts to prevent coherencies in L1
 */
void WAVELET3D::random_shift( void){
	for(int dir=0; dir< 3; dir++){
		S[dir] = rand() % wN;
		cout << "Dir = " << dir << "Shifting by " << S[dir] << endl;	
	}
}

/**
 * Call for forward transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::forward( Array<complex<float>,5>&Coef){
	for(int e=0;e<Coef.extent(fifthDim);e++){
	for(int t=0;t<Coef.extent(fourthDim);t++){
		Array<complex<float>,3>CoefRef=Coef(Range::all(),Range::all(),Range::all(),t,e);
		forward3D(CoefRef);
	}}
}

/**
 * Call for backward transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::backward( Array<complex<float>,5>&Coef){
	for(int e=0;e<Coef.extent(fifthDim);e++){
	for(int t=0;t<Coef.extent(fourthDim);t++){
		Array<complex<float>,3>CoefRef=Coef(Range::all(),Range::all(),Range::all(),t,e);
		backward3D(CoefRef);
	}}
}

/**
 * Call for forward transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::forward( Array<complex<float>,4>&Coef){
	for(int t=0;t<Coef.extent(fourthDim);t++){
		Array<complex<float>,3>CoefRef=Coef(Range::all(),Range::all(),Range::all(),t);
		forward3D(CoefRef);
	}
}

/**
 * Call for backward transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::backward( Array<complex<float>,4>&Coef){
	for(int t=0;t<Coef.extent(fourthDim);t++){
		Array<complex<float>,3>CoefRef=Coef(Range::all(),Range::all(),Range::all(),t);
		backward3D(CoefRef);
	}
}

/**
 * Call for forward transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::forward( Array<complex<float>,3>&Coef){
	forward3D(Coef);
}

/**
 * Call for backward transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::backward( Array<complex<float>,3>&Coef){
	backward3D(Coef);
}

/**
 * Private function for actually performing forward 3D Wavelet Transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::forward3D( Array<complex<float>,3>&Coef){
	int sz = N[2];
	int sy = N[1];
	int sx = N[0];
	for(int l=0; l<max_level; l++){
		if( l < L[0] ){
			wave_x(Coef,sz,sy,sx,0);
		}
		
		if( l < L[1] ){
			wave_y(Coef,sz,sy,sx,0);
		}
		
		if( l < L[2] ){
			wave_z(Coef,sz,sy,sx,0);
		}
		
		if( l < L[2] ){
			sz = N[2]/(int)pow(2.0,l+1);
		}
		if( l < L[1] ){
			sy = N[1]/(int)pow(2.0,l+1);
		}
		if( l < L[0] ){
			sx = N[0]/(int)pow(2.0,l+1);
		}
	}
}

/**
 * Private function for actually performing backward 3D Wavelet Transform
 * @param Coef Array to be transformed (in-place)
 */
void WAVELET3D::backward3D( Array<complex<float>,3>&Coef){
	
	int sz = N[2]/(int)pow(2.0,L[2]);
	int sy = N[1]/(int)pow(2.0,L[1]);
	int sx = N[0]/(int)pow(2.0,L[0]);
	
	for(int l=max_level-1; l>=0; l--){
		
		if( l < L[2] ){
			sz = N[2]/(int)pow(2.0,l);
		}
		if( l < L[1] ){
			sy = N[1]/(int)pow(2.0,l);
		}
		if( l < L[0] ){
			sx = N[0]/(int)pow(2.0,l);
		}
			
		if( l < L[2] ){
			wave_z(Coef,sz,sy,sx,1);
		}
		
		if( l < L[1] ){
			wave_y(Coef,sz,sy,sx,1);
		}
		
		if( l < L[0] ){
			wave_x(Coef,sz,sy,sx,1);
		}
	}
}



/**
 * Transform in X
 * @param Coef Array to be transformed (in-place)
 * @param Lx  Size of subarray to be transformed
 * @param Ly  Size of subarray to be transformed
 * @param Lz  Size of subarray to be transformed
 */
void WAVELET3D::wave_x(Array<complex<float>,3>&Coef,int Lz,int Ly,int Lx, int dir){
	// cout << "Wave X = " << Le << "," << Lt << "," << Lz << "," << Ly << "," << Lx << "  Dir = " << dir << endl;
		
	#pragma omp parallel for 
	for(int k=0; k<Lz; k++){
		
		complex<float> *buffer = new complex<float>[Lx];
		complex<float> *buffer2 = new complex<float>[Lx];
	
		for(int j=0; j<Ly; j++){
		
			//Copy
			for(int i=0; i<Lx; i++){	
				buffer[i]= Coef(i,j,k);
			}
			
			//Transform
			if(dir==0){
				wave1D(buffer,buffer2,Lx,0);
			}else{
				iwave1D(buffer,buffer2,Lx,0);
			}
			
			//Copy Back
			for(int i=0; i<Lx; i++){
				Coef(i,j,k) = buffer2[i];
			}
			
		}
		
		delete [] buffer;
		delete [] buffer2;
	}
}


/**
 * Transform in Y
 * @param Coef Array to be transformed (in-place)
 * @param Lx  Size of subarray to be transformed
 * @param Ly  Size of subarray to be transformed
 * @param Lz  Size of subarray to be transformed
 */
void WAVELET3D::wave_y(Array<complex<float>,3>&Coef,int Lz,int Ly,int Lx, int dir){
	// cout << "Wave Y = " << Le << "," << Lt << "," << Lz << "," << Ly << "," << Lx << "  Dir = " << dir << endl;
	
	#pragma omp parallel for 
	for(int k=0; k<Lz; k++){
		complex<float> *buffer = new complex<float>[Ly];
		complex<float> *buffer2 = new complex<float>[Ly];
	
		for(int i=0; i<Lx; i++){	
			
			//Copy
			for(int j=0; j<Ly; j++){
				buffer[j]= Coef(i,j,k);
			}
			
			//Transform
			if(dir==0){
				wave1D(buffer,buffer2,Ly,1);
			}else{
				iwave1D(buffer,buffer2,Ly,1);
			}
			
			//Copy Back
			for(int j=0; j<Ly; j++){
				Coef(i,j,k) = buffer2[j];
			}
			
		}
		delete [] buffer;
		delete [] buffer2;
	}
}

/**
 * Transform in Z
 * @param Coef Array to be transformed (in-place)
 * @param Lx  Size of subarray to be transformed
 * @param Ly  Size of subarray to be transformed
 * @param Lz  Size of subarray to be transformed
 */
void WAVELET3D::wave_z(Array<complex<float>,3>&Coef,int Lz,int Ly,int Lx, int dir){
	// cout << "Wave Z = " << Le << "," << Lt << "," << Lz << "," << Ly << "," << Lx << "  Dir = " << dir << endl;
	
	#pragma omp parallel for 
	for(int j=0; j<Ly; j++){
		complex<float> *buffer = new complex<float>[Lz];
		complex<float> *buffer2 = new complex<float>[Lz];
		
		for(int i=0; i<Lx; i++){	
			
			//Copy
			for(int k=0; k<Lz; k++){
				buffer[k]= Coef(i,j,k);
			}
			
			//Transform
			if(dir==0){
				wave1D(buffer,buffer2,Lz,2);
			}else{
				iwave1D(buffer,buffer2,Lz,2);
			}
			
			//Copy Back
			for(int k=0; k<Lz; k++){
				Coef(i,j,k) = buffer2[k];
			}
			
		}
		delete [] buffer;
		delete [] buffer2;
		
	}
}


/**
 * Actual 1D Transform for 1D Arrays
 * @param in data to transfrom
 * @param out output location
 * @param pts  size of in / out
 * @param dir  direction(x,y,z) for shifts
 */
void WAVELET3D::wave1D( complex<float> in[], complex<float> out[],int pts,int dir){
	
	/*Scaling Coef*/
	for(int pos=0;pos<pts/2; pos++){
		out[pos]=complex<float>(0,0);
		out[pos+pts/2]=complex<float>(0,0);
		for(int fp =0; fp< wN; fp++){
			int cpos = ( (2*pos - fp + S[dir]) + pts)% pts;
			out[pos]+= lpf[fp]*in[cpos];
			out[pos+pts/2]+= hpf[fp]*in[cpos];
		}
	}
}

/**
 * Actual 1D Synthesis for 1D Arrays
 * @param in data to transfrom
 * @param out output location
 * @param pts  size of in / out
 * @param dir  direction(x,y,z) for shifts
 */
void WAVELET3D::iwave1D( complex<float> in[],complex<float> out[],int pts,int dir){
	for(int pos=0;pos<pts; pos++){
		out[pos]= complex<float>(0,0);
	}
	
	/*Scaling Coef*/
	for(int pos=0;pos<pts/2; pos++){
		//Add to nearby points
		for(int fp =0; fp< wN; fp++){
			int cpos = ( (2*pos + fp-wN+1+S[dir]) + pts)% pts; // Three is for shift correction
			out[cpos]+= in[pos]*Slpf[fp];
			out[cpos]+= in[pos+pts/2]*Shpf[fp];
		}
	}	
}


