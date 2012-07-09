#ifndef hTICTOC
#define hTICTOC

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
      printf("Elapsed time: %.6fs\n",stop-start);
    }
    
    void toc(char* message){
      stop = omp_get_wtime();
      printf("%s %.6fs\n", message, stop-start);
    }
    
};

#endif
