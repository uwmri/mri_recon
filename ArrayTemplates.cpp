#include "ArrayTemplates.hpp"

void NDarray::fftshift( Array< complex<float>,3>& temp){
	for(int k=0; k<temp.extent(thirdDim);k++){
		for(int j=0; j<temp.extent(secondDim);j++){
			for(int i=0; i<temp.extent(firstDim);i++){
				float mod = ((float)( 2*(( i+j+k)%2) - 1));
				temp(i,j,k) *= mod;
			}}}	
}

void NDarray::ifft( Array< complex<float>,3>& temp){
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

	// Scale
	float scale = 1.0/sqrt( (float)temp.extent(thirdDim) * (float)temp.extent(secondDim) * (float)temp.extent(firstDim));
	temp *= scale;

	// Cleanup	
	fftwf_destroy_plan(fft_plan);
}


void NDarray::fft( Array< complex<float>,3>& temp){
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

	// Scale
	float scale = 1.0/sqrt( (float)temp.extent(thirdDim) * (float)temp.extent(secondDim) * (float)temp.extent(firstDim));
	temp *= scale;

	// Cleanup	
	fftwf_destroy_plan(fft_plan);
}



// FFT in only one dimension
void NDarray::fft( Array< complex<float>,3>& temp, int dim){
	fft3( temp, dim, FFTW_FORWARD,1);
}

void NDarray::ifft( Array< complex<float>,3>& temp, int dim){
	fft3( temp, dim, FFTW_BACKWARD,1);
}

void NDarray::fft3( Array< complex<float>,3>& temp, int dim, int direction, bool chop){

	fftwf_init_threads();
	fftwf_plan_with_nthreads(1);

	// Get Size		
	int N;
	switch(dim){
		case(0):{ N=temp.length(firstDim); }break;
		case(1):{ N=temp.length(secondDim); }break;
		case(2):{ N=temp.length(thirdDim); }break;
		default:{
			cout << "Error trying to FFT a dimension that doesn't exist" << endl;
			exit(1);
		}
	}
	
	
	// Get a plan but never use
	complex<float> *data_temp = new complex<float>[N];
	fftwf_complex *ptr = reinterpret_cast<fftwf_complex*>(data_temp);
	fftwf_plan plan = fftwf_plan_dft_1d(N,ptr,ptr,direction, FFTW_MEASURE);
	delete [] data_temp;
	
	float scale = 1./sqrt(N);
	switch(dim){
		case(0):{
			#pragma omp parallel for
			for(int k=0;k<temp.extent(thirdDim);k++){
		
			complex<float> *data = new complex<float>[N];
			fftwf_complex *data_ptr = reinterpret_cast<fftwf_complex*>(data);
			
			for(int j=0;j<temp.extent(secondDim);j++){
				// Copy
				if(chop==1){
					for(int i=0;i<temp.extent(firstDim);i++){
						data[i] = temp(i,j,k)*((float)( 2*(i%2)-1));
					}
				}else{
					for(int i=0;i<temp.extent(firstDim);i++){
						data[i] = temp(i,j,k);
					}
				}	
							
				// FFT
				fftwf_execute_dft(plan,data_ptr,data_ptr);
				
				// Copy Back
				if(chop){
					for(int i=0;i<temp.extent(firstDim);i++){
						temp(i,j,k)=data[i]*((float)( 2*(i%2)-1))*scale;
					}
				}else{
					for(int i=0;i<temp.extent(firstDim);i++){
						temp(i,j,k)=data[i];
					}
				}
			}
			delete [] data;
			}
		}break;
	
		case(1):{
			#pragma omp parallel for
			for(int k=0;k<temp.extent(thirdDim);k++){
		
			complex<float> *data = new complex<float>[N];
			fftwf_complex *data_ptr = reinterpret_cast<fftwf_complex*>(data);
			
			for(int i=0;i<temp.extent(firstDim);i++){
				
				// Copy
				if(chop==1){
					for(int j=0;j<temp.extent(secondDim);j++){
						data[j] = temp(i,j,k)*((float)( 2*(j%2)-1));
					}
				}else{
					for(int j=0;j<temp.extent(secondDim);j++){
						data[j] = temp(i,j,k);
					}
				}
				// FFT
				fftwf_execute_dft(plan,data_ptr,data_ptr);
				
				// Copy Back
				if(chop){
					for(int j=0;j<temp.extent(secondDim);j++){
						temp(i,j,k) = data[j]*((float)( 2*(j%2)-1))*scale;
					}
				}else{
					for(int j=0;j<temp.extent(secondDim);j++){
						temp(i,j,k) = data[j];
					}
				}
			}
			delete [] data;
			}
		
		}break;
	
		case(2):{
			#pragma omp parallel for
			for(int j=0;j<temp.extent(secondDim);j++){
				
			complex<float> *data = new complex<float>[N];
			fftwf_complex *data_ptr = reinterpret_cast<fftwf_complex*>(data);
			
			for(int i=0;i<temp.extent(firstDim);i++){
				
				// Copy
				if(chop){
					for(int k=0;k<temp.extent(thirdDim);k++){
						data[k] = temp(i,j,k)*((float)( 2*(k%2)-1));
					}
				}else{
					for(int k=0;k<temp.extent(thirdDim);k++){
						data[k] = temp(i,j,k);
					}
				}
				
				// FFT
				fftwf_execute_dft(plan,data_ptr,data_ptr);
				
				// Copy Back
				if(chop){
					for(int k=0;k<temp.extent(thirdDim);k++){
						temp(i,j,k)= data[k]*((float)( 2*(k%2)-1))*scale;
					}
				}else{
					for(int k=0;k<temp.extent(thirdDim);k++){
						temp(i,j,k)= data[k];
					}
				}
				
			}
			
			delete [] data;
			}
			
		}break;
	}// Switch
	
	// Cleanup	
	fftwf_destroy_plan(plan);
}

