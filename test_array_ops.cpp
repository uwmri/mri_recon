#include "recon_lib.h"

using namespace NDarray;
using namespace std;

int main(int argc, char **argv){
	
	int N = 512;
	Array< complex<float>,3> X(N,N,N);
	
	
	for(int i=0; i<N; i++){
	for(int j=0; j<N; j++){
	for(int k=0; k<N; k++){
		X(k,j,i) = complex<float>(i+j,k);
	}}}
	
	
	tictoc T; 
	
	
	T.tic();
	ArrayWriteMag(X,"X_OF.dat");
	cout << "Tool " << T << endl;
	
	T.tic();
	ofstream ofs("Xbuff.dat", ios_base::binary);
	for( Array<complex<float>,3>::iterator miter=X.begin(); miter!=X.end(); miter++){
		float val=abs( *miter);
		ofs.write( (char *)&val,sizeof(T));
   	}
	cout << "Tool " << T << endl;
	
	
	
	
	
	T.tic();
	ArrayWriteMag(X,"X_OF.dat");
	cout << "Tool " << T << endl;
	
	for(int buffsize = 8; buffsize < 8192; buffsize*=2){
		char *buffer = new char[buffsize];
	
	
		T.tic();
		ofstream ofs("Xbuff.dat", ios_base::binary);
		ofs.rdbuf()->pubsetbuf(buffer, buffsize);
		for( Array<complex<float>,3>::iterator miter=X.begin(); miter!=X.end(); miter++){
			float val=abs( *miter);
			ofs.write( (char *)&val,sizeof(T));
    	}
		cout << "Buffer size = " << buffsize << " took " << T << endl;
		delete [] buffer;
	}	
	T.tic();
	ArrayWriteMag(X,"X_OF.dat");
	cout << "Tool " << T << endl;
	
	
	
	T.tic();
	Array< float,3>Xmag(N,N,N);
	Xmag = abs(X);
	FILE *fid = fopen("Xcont.dat","w");
	fwrite(Xmag.data(),sizeof(float),N*N*N,fid);
	fclose(fid);
	 cout << "Tool " << T << endl;
	
	
	return 0;
}
	



