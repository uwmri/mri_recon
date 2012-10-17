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
    }
    
	tictoc& operator=(const tictoc&);
	
	friend ostream& operator<<(ostream& out, tictoc& v){
        v.toc();
		return out << (v.stop-v.start);
    }
	
	
};


#endif