/************************************************
Soft Thresholding Library

Description: This code contains functions to perform thresholding operations


*************************************************/
#include "pcvipr_lib.h"

/*----------------------------------------------
     Constructor  - read command line
 *----------------------------------------------*/ 

PCVIPR::PCVIPR( int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if (strcmp("-h", pstring[pos] ) == 0) {
	  	printf("\n*********************************************\n");
	  	printf("PCVIPR Control:\n");
	  	printf("*********************************************\n");
	  	
	
	}
   }    

}

void PCVIPR::interpret_pfile_header(PFILE pfile){


}

void PCVIPR::get_kspace_trajectory(PFILE pfile){



}

