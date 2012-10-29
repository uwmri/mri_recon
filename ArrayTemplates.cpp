#include "ArrayTemplates.hpp"

void fftshift( Array< complex<float>,3>& temp){
	for(int k=0; k<temp.extent(thirdDim);k++){
	for(int j=0; j<temp.extent(secondDim);j++){
	for(int i=0; i<temp.extent(firstDim);i++){
		float mod = ((float)( 2*(( i+j+k)%2) - 1));
        temp(i,j,k) *= mod;
	}}}	
}

void ifft( Array< complex<float>,3>& temp){
	// Shift to Center of K-space
	fftshift( temp);
	
	// Do FFT 
	fftwf_complex *ptr=NULL;
	complex<float>*data=NULL;
	if(temp.isStorageContiguous()){ 
	 	ptr = reinterpret_cast<fftwf_complex*>(temp.data());
	}else{
		cout << "Warning:: Doing FFT of Non-Contiguous Data. Will be slow!" << endl;
		data = new complex<float>[temp.numElements()];
		int pos=0;
		for(int k=0;k<temp.extent(thirdDim);k++){
		for(int j=0;j<temp.extent(secondDim);j++){
		for(int i=0;i<temp.extent(firstDim);i++){
			data[pos]=temp(i,j,k);
			pos++;
		}}}
		ptr = reinterpret_cast<fftwf_complex*>(data);
	}	
	
	fftwf_plan fft_plan = fftwf_plan_dft_3d(temp.length(thirdDim),temp.length(secondDim),temp.length(firstDim),ptr,ptr,FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwf_execute(fft_plan);	
	
	if(!temp.isStorageContiguous()){ 
		int pos=0;
		for(int k=0;k<temp.extent(thirdDim);k++){
		for(int j=0;j<temp.extent(secondDim);j++){
		for(int i=0;i<temp.extent(firstDim);i++){
			temp(i,j,k)=data[pos];
			pos++;
		}}}
		free(ptr);
	}	
	
	// Fix Phase from Shift
	fftshift( temp);
	
	// Cleanup	
	fftwf_destroy_plan(fft_plan);
}


void fft( Array< complex<float>,3>& temp){
	// Shift to Center of K-space
	fftshift( temp);
	
	// Do FFT 
	fftwf_complex *ptr=NULL;
	complex<float>*data=NULL;
	if(temp.isStorageContiguous()){ 
	 	ptr = reinterpret_cast<fftwf_complex*>(temp.data());
	}else{
		cout << "Warning:: Doing FFT of Non-Contiguous Data. Will be slow!" << endl;
		data = new complex<float>[temp.numElements()];
		int pos=0;
		for(int k=0;k<temp.extent(thirdDim);k++){
		for(int j=0;j<temp.extent(secondDim);j++){
		for(int i=0;i<temp.extent(firstDim);i++){
			data[pos]=temp(i,j,k);
			pos++;
		}}}
		ptr = reinterpret_cast<fftwf_complex*>(data);
	}	
	
	fftwf_plan fft_plan = fftwf_plan_dft_3d(temp.length(thirdDim),temp.length(secondDim),temp.length(firstDim),ptr,ptr,FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(fft_plan);	
	
	if(!temp.isStorageContiguous()){ 
		int pos=0;
		for(int k=0;k<temp.extent(thirdDim);k++){
		for(int j=0;j<temp.extent(secondDim);j++){
		for(int i=0;i<temp.extent(firstDim);i++){
			temp(i,j,k)=data[pos];
			pos++;
		}}}
		free(ptr);
	}	
	
	// Fix Phase from Shift
	fftshift( temp);
	
	// Cleanup	
	fftwf_destroy_plan(fft_plan);
}


