/************************************************
3D Wavlet Libraries for pcvipr.e

Initial Author: Kevin M. Johnson
Description: This code contains functions for 3D wavelet transform. I tried blitzwave 
  but ended up not liking the lack of updates for newer gcc compilers and lack of threads.

*************************************************/

#include "wavelet3D.h"

WAVELET3D::~WAVELET3D(){
 	// Not here yet
}

#ifdef XXX
WAVELET3D::WAVELET3D( int sz,int sy, int sx, int l, int wavelet_type){
	
	/// cout << "Init Wavelet " << endl;
	levels= l;
	
	// Check dimensions to make sure we don't have to zero pad
	N[0] = sx;
	N[1] = sy;
	N[2] = sz;
	N[3] = 1;
	N[4] = 1;
	Coef.alloc(1,1,N[2],N[1],N[0]);
	get_filter_banks( wavelet_type );
	setup_wavelet_directions();
}
#endif

WAVELET3D::WAVELET3D( array3D< complex<float> >*temp, int *l, int type){
	
	/// cout << "Init Wavelet " << endl;
	N[4] = 1;
	N[3] = 1;
	N[2] = temp->Nz;
	N[1] = temp->Ny;
	N[0] = temp->Nx;
	
	L[2] = l[2];
	L[1] = l[1];
	L[0] = l[0];

	wType = type;
			
	get_filter_banks();
	Coef.point_to_3D(temp);
	setup_wavelet_directions();
}



WAVELET3D::WAVELET3D( array4D< complex<float> >*temp, int l, int type){
	
	/// cout << "Init Wavelet " << endl;
	N[4] = 1;
	N[3] = temp->Nt;
	N[2] = temp->Nz;
	N[1] = temp->Ny;
	N[0] = temp->Nx;
	
	L[2] = l;
	L[1] = l;
	L[0] = l;
	
	wType = type;
	
	get_filter_banks();
	Coef.point_to_4D(temp);
	setup_wavelet_directions();
	
}


WAVELET3D::WAVELET3D( array4D< complex<float> >*temp, int *l, int type){
	
	/// cout << "Init Wavelet " << endl;
	N[4] = 1;
	N[3] = temp->Nt;
	N[2] = temp->Nz;
	N[1] = temp->Ny;
	N[0] = temp->Nx;
	
	L[2] = l[2];
	L[1] = l[1];
	L[0] = l[0];
	
	wType = type;
			
	get_filter_banks();
	Coef.point_to_4D(temp);
	setup_wavelet_directions();
}

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


WAVELET3D::WAVELET3D( array5D< complex<float> >*temp, int *l, int type){
	
	/// cout << "Init Wavelet " << endl;
	N[4] = temp->Ne;
	N[3] = temp->Nt;
	N[2] = temp->Nz;
	N[1] = temp->Ny;
	N[0] = temp->Nx;
	
	L[4] = l[4];
	L[3] = l[3];
	L[2] = l[2];
	L[1] = l[1];
	L[0] = l[0];
	
	wType = type;		

	Coef.Nx=temp->Nx;
	Coef.Ny=temp->Ny;
	Coef.Nz=temp->Nz;
	Coef.Nt=temp->Nt;
	Coef.Ne=temp->Ne;
	Coef.vals=temp->vals;

	get_filter_banks();
	setup_wavelet_directions();
	
	cout << "Wavelet Parameters" << endl;
	for(int dir=0; dir<3; dir++){
		cout << "\tDir = " << dir << " N = " << N[dir] << " Levels= " << L[dir]  << endl;
	}
		
}

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

void WAVELET3D::random_shift( void){
	for(int dir=0; dir< 3; dir++){
		S[dir] = rand() % wN;
		cout << "Dir = " << dir << "Shifting by " << S[dir] << endl;	
	}
}


void WAVELET3D::forward( void){
	int sz = N[2];
	int sy = N[1];
	int sx = N[0];
	for(int l=0; l<max_level; l++){
		if( l < L[0] ){
			wave_x(N[4],N[3],sz,sy,sx,0);
		}
		
		if( l < L[1] ){
			wave_y(N[4],N[3],sz,sy,sx,0);
		}
		
		if( l < L[2] ){
			wave_z(N[4],N[3],sz,sy,sx,0);
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

void WAVELET3D::backward( void){
	
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
			wave_z(N[4],N[3],sz,sy,sx,1);
		}
		
		if( l < L[1] ){
			wave_y(N[4],N[3],sz,sy,sx,1);
		}
		
		if( l < L[0] ){
			wave_x(N[4],N[3],sz,sy,sx,1);
		}
	}
}


/*----------------------------------------------
     Transform in X - Threaded version
 *----------------------------------------------*/ 
void WAVELET3D::wave_x(int Le,int Lt,int Lz,int Ly,int Lx, int dir){
	// cout << "Wave X = " << Le << "," << Lt << "," << Lz << "," << Ly << "," << Lx << "  Dir = " << dir << endl;
	for(int e=0; e<Le; e++){
	for(int t=0; t<Lt; t++){
		
	#pragma omp parallel for 
	for(int k=0; k<Lz; k++){
		
		complex<float> *buffer = new complex<float>[Lx];
		complex<float> *buffer2 = new complex<float>[Lx];
	
		for(int j=0; j<Ly; j++){
		
			//Copy
			for(int i=0; i<Lx; i++){	
				buffer[i]= Coef[e][t][k][j][i];
			}
			
			//Transform
			if(dir==0){
				wave1D(buffer,buffer2,Lx,0);
			}else{
				iwave1D(buffer,buffer2,Lx,0);
			}
			
			//Copy Back
			for(int i=0; i<Lx; i++){
				Coef[e][t][k][j][i] = buffer2[i];
			}
			
		}
		
		delete [] buffer;
		delete [] buffer2;
	}
	}
	}
}


/*----------------------------------------------
     Transform in Y - Pthread version
 *----------------------------------------------*/ 
void WAVELET3D::wave_y(int Le,int Lt,int Lz,int Ly,int Lx, int dir){
	// cout << "Wave Y = " << Le << "," << Lt << "," << Lz << "," << Ly << "," << Lx << "  Dir = " << dir << endl;
	
	for(int e=0; e<Le; e++){
	for(int t=0; t<Lt; t++){
	
	#pragma omp parallel for 
	for(int k=0; k<Lz; k++){
		complex<float> *buffer = new complex<float>[Ly];
		complex<float> *buffer2 = new complex<float>[Ly];
	
		for(int i=0; i<Lx; i++){	
			
			//Copy
			for(int j=0; j<Ly; j++){
				buffer[j]= Coef[e][t][k][j][i];
			}
			
			//Transform
			if(dir==0){
				wave1D(buffer,buffer2,Ly,1);
			}else{
				iwave1D(buffer,buffer2,Ly,1);
			}
			
			//Copy Back
			for(int j=0; j<Ly; j++){
				Coef[e][t][k][j][i] = buffer2[j];
			}
			
		}
		delete [] buffer;
		delete [] buffer2;
	}
	}
	}
}

/*----------------------------------------------
     Transform in Z - Pthread version
 *----------------------------------------------*/ 
void WAVELET3D::wave_z(int Le,int Lt,int Lz,int Ly,int Lx, int dir){
	// cout << "Wave Z = " << Le << "," << Lt << "," << Lz << "," << Ly << "," << Lx << "  Dir = " << dir << endl;
	
	for(int e=0; e<Le; e++){
	for(int t=0; t<Lt; t++){
		
	#pragma omp parallel for 
	for(int j=0; j<Ly; j++){
		complex<float> *buffer = new complex<float>[Lz];
		complex<float> *buffer2 = new complex<float>[Lz];
		
		for(int i=0; i<Lx; i++){	
			
			//Copy
			for(int k=0; k<Lz; k++){
				buffer[k]= Coef[e][t][k][j][i];
			}
			
			//Transform
			if(dir==0){
				wave1D(buffer,buffer2,Lz,2);
			}else{
				iwave1D(buffer,buffer2,Lz,2);
			}
			
			//Copy Back
			for(int k=0; k<Lz; k++){
				Coef[e][t][k][j][i] = buffer2[k];
			}
			
		}
		delete [] buffer;
		delete [] buffer2;
		
	}
	}
	}
}


/*----------------------------------------------
     1D Wavelet transform
*----------------------------------------------*/ 
void WAVELET3D::wave1D( complex<float> in[],complex<float> out[],int pts,int dir){
	
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

/*----------------------------------------------
     1D Wavelet inverse transform
*----------------------------------------------*/ 
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


