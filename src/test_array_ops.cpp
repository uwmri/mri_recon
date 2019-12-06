#include "recon_lib.h"

using namespace NDarray;
using namespace std;

int main(void){
	
	
	tictoc T; 
	
	{
	cout << "Check SVD Performance" << endl;
	int N = 256*256*8;
	int Np = 385;

	// preallocate some arrays
	arma::cx_mat A = arma::randn< arma::cx_mat >(N,Np);
    
    arma::cx_mat U;
	arma::cx_mat V;
	arma::vec s;

	cout << "working on SVD ( " << N << " , " << Np << ")" << endl << std::flush;
	T.tic();
	arma::svd_econ(U,s,V,A,"both");
	cout << "it took " << T << " s per iteraiton" << endl << flush;

	cout << "working on rotation" << endl;
	T.tic();
	A = U*arma::diagmat(s)*V.t();
	cout << "it took " << T << " s per iteraiton" << endl << flush;
	
	cout << "Took " << T << endl;
	}
		
	
	
	T.tic();
	{
	cout << "Check BART Export" << endl;
	Array< Array<complex<float>,3>, 2> TEMP = Alloc5DContainer< complex<float> >( 3, 4, 5, 6, 7);
	//WriteCFL( TEMP, "BART_CFL");
	cout << "Took " << T << endl;
	}
	
	T.tic();
	{
	cout << "Check BART Export" << endl;
	Array< float, 3> TEMP( 3, 4, 5, ColumnMajorArray<3>());
	WriteCFL( TEMP, "BART_CFL_FLOAT");
	cout << "Took " << T << endl;
	}
	
	{
	cout << "Check BART Export" << endl;
	Array< Array<complex<float>,3>, 2> TEMP1 = Alloc5DContainer< complex<float> >( 3, 4, 1, 1, 1);
	Array< Array<complex<float>,3>, 2> TEMP2 = Alloc5DContainer< complex<float> >( 3, 4, 1, 1, 1);
	Array< Array<complex<float>,3>, 2> TEMP3 = Alloc5DContainer< complex<float> >( 3, 4, 1, 1, 1);
	//WriteCFL_triplet( TEMP1,TEMP2,TEMP3, "BART_TRIPLE_CFL_FLOAT");
	cout << "Took " << T << endl;
	}
	
	MRI_DATA data;
	data.Num_Encodings = 1;
	data.Num_Coils = 1;
	data.init_memory(256,1000,1);
	data.write_bart_data("MRI_Raw_Bart");
	
	
	T.tic();
	int N= 256;
	Array< complex<float>,3> X(N,N,N);
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
	



