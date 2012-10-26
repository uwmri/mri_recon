#ifndef hTIMEARRAY
#define hTIMEARRAY

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <omp.h>

class tictoc{
  public:
    double start;
    double stop;
    
    void tic(){
      start = omp_get_wtime();
    }
    
    void toc(){
      stop = omp_get_wtime();
>>>>>>>>>>>>>>>>>>>> File 1
      printf("Elapsed time: %.6fs\n",stop-start);
>>>>>>>>>>>>>>>>>>>> File 2
<<<<<<<<<<<<<<<<<<<<
    }
    
	tictoc& operator=(const tictoc&);
	
	friend ostream& operator<<(ostream& out, tictoc& v){
        v.toc();
		return out << (v.stop-v.start);
    }
	
	
};


#endif
