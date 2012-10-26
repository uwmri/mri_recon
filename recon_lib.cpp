#include "recon_lib.h"
#include "io_templates.cpp"

// ----------------------
//  Basic constructor (no args)
// ----------------------
RECON::RECON(void){
	set_defaults();
}

// ----------------------
//  Sets Default Recon Parameters
// ----------------------
void RECON::set_defaults( void){
	// Help Message for recon
	recon_type = RECON_SOS;
	data_type = RECON_EXTERNAL;
	
	numrecv = 1;
	zero_fill = 1.0;
	zoom = 1.0;
	zoom_x = 1.0;
	zoom_y = 1.0;
	zoom_z = 1.0;
	  
	rcxres=-1;
	rcyres=-1;
	rczres=-1;
	rcframes=1;
	rcencodes=1;
	num_slices =1;    
	lp_frac=1.0;
	smap_res=16;
	
	acc = 1;
	compress_coils = 0.0;
	max_iter = 50;		
}

// ----------------------
//  Constructor with Command Line Read
// ----------------------

RECON::RECON(int numarg, char **pstring){
	set_defaults();	
	parse_commandline(numarg,pstring);
}

// ----------------------
// Help Message
// ----------------------
void RECON::help_message(void){
	cout << "----------------------------------------------" << endl;
	cout << "   Basic Recon Control " << endl;
	cout << "----------------------------------------------" << endl;
	cout << "Usage:" << endl;
	cout << "   recon_binary -f header.txt [flags]" << endl;
	
	cout << "Recon Size:" << endl;
	help_flag("-rcxres []","matrix size in x");
	help_flag("-rcyres []","matrix size in y");
	help_flag("-rczres []","matrix size in z");
	help_flag("-rcframes []","reconstructed temporal frames");
	help_flag("-zoom_x []","zoom factor in x");
	help_flag("-zoom_y []","zoom factor in x");
	help_flag("-zoom_z []","zoom factor in x");
	
	cout << "Recon Types:" << endl;
	help_flag("-sos","sum of squares");
	help_flag("-pils","pils (coil combine with low resolution images)");
	help_flag("-ist","iterative soft thresholding");
	help_flag("-fista","fast iterative soft thresholding");
	
	cout << "Recon Control:" << endl;
	help_flag("-max_iter []","max iterations for iterative recons");

}

// --------------------
//  Read command line and set variables
// --------------------
void RECON::parse_commandline(int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]);
  	  
  for(int pos=0; pos < numarg; pos++){
 	
  	if (strcmp("-h", pstring[pos] ) == 0) {
		char_flag("-f",filename);
		
		// Reconstruction Geometry
		int_flag("-rcxres",rcxres);
		int_flag("-rcyres",rcyres);
		int_flag("-rczres",rczres);
		int_flag("-rcframes",rcframes);
		
		float_flag("-zoom",zoom);
		float_flag("-zoom_x",zoom_x);
		float_flag("-zoom_y",zoom_y);
		float_flag("-zoom_z",zoom_z);
		
		// Type of Recons		
		trig_flag(RECON_SOS,"-sos",recon_type);
		trig_flag(RECON_CG,"-isense",recon_type);
		trig_flag(RECON_PILS,"-pils",recon_type);
		trig_flag(RECON_IST,"-ist",recon_type);
		trig_flag(RECON_FISTA,"-fista",recon_type);
		
		// Source of data
		trig_flag(RECON_EXTERNAL,"-external_data",data_type);
		trig_flag(RECON_PFILE,"-pfile",data_type);
		
		// Data modification
		int_flag("-acc",acc);
		float_flag("-compress_coils",compress_coils);
		
		// Coil Combination + Resolution		
		float_flag("-lp_frac",lp_frac);
		float_flag("-smap_res",smap_res);
		
		// Iterations for IST
		int_flag("-max_iter",max_iter);
		
	}
  }
} 

//--------------------------------------------------
//  Read external header
//--------------------------------------------------

void RECON::parse_external_header(void){
	
	char parameter[80];
	float value;
	float value2;
	float value3;
	float value4;
	FILE *fid;
	char line[200];
	
	cout << "Reading External Header: " << endl;
	
	fid = fopen(filename,"r");
	while( fgets(line, sizeof(line),fid) != NULL ) {
    	if(sscanf(line,"%s\t%f\t%f\t%f\t%f",parameter,&value,&value2,&value3,&value4) < 2){
		}else if(strcmp("acq_bw",parameter) == 0){ acq_bw = value;
		}else if(strcmp("xres",parameter) == 0){ xres = (int)value;
		}else if(strcmp("numrecv",parameter) == 0){ num_coils = (int)value;
		}else if(strcmp("slices",parameter) == 0){ num_slices = (int)value;
		}else if(strcmp("2d_flag",parameter) == 0){ ss_2d = (int)value;
		}else if(strcmp("nproj",parameter) == 0){ num_readouts = (int)value;
		}else if(strcmp("rcxres",parameter) == 0){ 
			rcxres     = (rcxres == -1 ) ?  ( (int)value ) : ( rcxres);
		}else if(strcmp("rcyres",parameter) == 0){ 
			rcyres     = (rcyres == -1 ) ?  ( (int)value ) : ( rcyres);
		}else if(strcmp("rczres",parameter) == 0){ 
			rczres     = (rczres == -1 ) ?  ( (int)value ) : ( rczres);
		}else if(strcmp("multi_echo",parameter) == 0){ multi_echo = (int)value;
		}else if(strcmp("num_encodes",parameter) == 0){ rcencodes = (int)value;
		}
	}
	fclose(fid);
}




