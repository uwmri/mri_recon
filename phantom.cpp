/*

Title: Phantom Library   

Summary: This code exists to create simulated phantom data 

Usage: Currently 

	PHANTOM phantom;						
	phantom.read_commandline(argc,argv);  
	phantom.init(recon.rcxres,recon.rcyres,recon.rczres);
	phantom.update_smap_biotsavart(coil,recon.num_coils);	
	
	// Which will populate:
	phantom.IMAGE // Image (density)
	phantom.SMAP // Sensitivity Map 
	phantom.TOA  // Time of arrival

*/

#include "phantom.h"
using namespace arma;
using namespace std;

// Default Constructor
//	
//
PHANTOM::PHANTOM( void){
	over_res=2;
	phantom_type = FRACTAL;
	phantom_noise= 0.0;	
	debug = 0;
}

// Initialization:
//	-This function sets up a digital phantom of size (Nx,Ny,Nz)
//
void PHANTOM::init(int Sx, int Sy, int Sz){
	
	// Don't change resolution
	if(phantom_type==EXTERNAL){
		over_res =1.0;
	}	
	
	// Copy Size
	Nx = (int)(over_res*Sx);
	Ny = (int)(over_res*Sy);
	Nz = (int)(over_res*Sz);
	
	// Phantom Density
	cout << "Init IMAGE: " << Nx << " by " << Ny << " by " << Nz << endl << flush;
	IMAGE.setStorage( ColumnMajorArray<3>());
  	IMAGE.resize(Nx,Ny,Nz);
  	IMAGE = 0;
	
	// Sensitivity Map	
	cout << "Init Smap" << endl << flush;
	SMAP.setStorage( ColumnMajorArray<3>());
  	SMAP.resize(Nx,Ny,Nz);
  	SMAP = 1;
		
	// Switch For Phantom
	switch(EXTERNAL){
		case(FRACTAL):{
		
			
			cout << "Init TOA" << endl << flush;
			TOA.setStorage( ColumnMajorArray<3>());
  			TOA.resize(Nx,Ny,Nz);
  			TOA = 0;
			
			cout << "Build Tree" << endl << flush;
			fractal3D_new(Nx,Ny,Nz);
		}break;
		
		case(SHEPP):{
			// Doesn't exist yet (but can be k-space based, Koay et al.)
		}break;
		
		case(PSF):{
			// Just Get PSF 
		}break;
		
		case(EXTERNAL):{
			
			
			cout << "Reading External Phantom" << external_phantom_name << endl;
			ArrayRead(IMAGE,external_phantom_name);
		}break;
	}
}


// Get a sensitivity based on coil geometry
// with Biot-Savart simulation of field.
//
void  PHANTOM::update_smap_biotsavart(int coil, int Ncoils){
	
	// Doesn't Make sense for one coil
	if(Ncoils==1){
		SMAP = 1;
		cout << "No Coil Mapping! (one coil)" << endl;	
		return;
	}
			
	// Made up Coil Geometry
	float coil_length = 1.25*Nx/2;	
	float coil_radius = 1.25*Nx/2; 
	float overlap_angle = 2.0*PI/(float)Ncoils/4;
	float overlap_length= coil_length/3;
	
	// Wire Positions
	int RES = 100;
	arma::fmat loop_pos(3,2*5*RES); // RES pts a side
	arma::fmat loop_current(3,2*5*RES); // RES pts a side
	int count=0;
	
	//----------------------------------------------
	//  Create Wire Loop
	//----------------------------------------------	
		
	/*Coil element looks like this
	
	/------\
	|      |
	|      |
	|      |
	\------/
	
	Elements are placed in the Surface of a cylinder with beveled edges for overlap.
	*/
	
	// First Rail Section ( left )
	//cout << "First Rail" << endl << flush;
	for(int pos=0;pos< RES; pos++){
		loop_pos(0,count) = coil_radius*cos( -overlap_angle);
		loop_pos(1,count) = coil_radius*sin( -overlap_angle);
		loop_pos(2,count) = (coil_length -overlap_length)*( (float)pos/(float)RES);
		count++;
	}
		
	// First Overlap Section (top left )
	//cout << "Overlap 1" << endl << flush;
	for(int pos=0;pos< RES; pos++){
		loop_pos(0,count) = coil_radius*cos( -overlap_angle + (2*overlap_angle)*( (float)pos/(float)RES));
		loop_pos(1,count) = coil_radius*sin( -overlap_angle + (2*overlap_angle)*( (float)pos/(float)RES));
		loop_pos(2,count) = (coil_length -overlap_length) + overlap_length*( (float)pos/(float)RES);
		count++;
	}
	
	// Loop (top line)
	//cout << "Top Loop" << endl << flush;
	for(int pos=0;pos< RES; pos++){
		loop_pos(0,count) = coil_radius*cos(overlap_angle+ (2*PI/(float)Ncoils-2.0*overlap_angle)*( (float)pos/(float)RES));
		loop_pos(1,count) = coil_radius*sin(overlap_angle+ (2*PI/(float)Ncoils-2.0*overlap_angle)*( (float)pos/(float)RES));
		loop_pos(2,count) = coil_length;
		count++;
	}
	
	// Second Overlap Section (top right)
	//cout << "Overlap 2" << endl << flush;
	for(int pos=0;pos< RES; pos++){
		loop_pos(0,count) = coil_radius*cos( 2*PI/(float)Ncoils -overlap_angle+ (2*overlap_angle)*( (float)pos/(float)RES));
		loop_pos(1,count) = coil_radius*sin( 2*PI/(float)Ncoils -overlap_angle+ (2*overlap_angle)*( (float)pos/(float)RES));
		loop_pos(2,count) = coil_length - overlap_length*( (float)(pos+1)/(float)RES);
		count++;
	}
	
	// Rail right side
	//cout << "Second Rail" << endl << flush;
	for(int pos=0;pos< RES; pos++){
		loop_pos(0,count) = coil_radius*cos( overlap_angle+2*PI/(float)Ncoils);
		loop_pos(1,count) = coil_radius*sin( overlap_angle+2*PI/(float)Ncoils);
		loop_pos(2,count) = ( coil_length - overlap_length ) - (coil_length-overlap_length)*( (float)pos/(float)RES);
		count++;
	}
	int backtrace_count=count;
	
	//cout << "Backtrace" << endl << flush;
	// Backtrace to get symetric other side
	for(int pos=0; pos< backtrace_count; pos++){
		loop_pos(0,count)= loop_pos(0,backtrace_count-pos-1);
		loop_pos(1,count)= loop_pos(1,backtrace_count-pos-1);
		loop_pos(2,count)=-loop_pos(2,backtrace_count-pos-1);	
		count++;
	}	
	
	//----------------------------------------------
	//  Rotate the Coil 
	//----------------------------------------------	
	arma::fmat R;
	float a = 2*PI/(float)Ncoils * (float)coil;
	
	R << cos(a) << sin(a) << 0.0 << endr
  	  <<-sin(a) << cos(a) << 0.0 << endr
	  << 0.0    << 0.0    << 1.0 << endr;
	
	loop_pos = R*loop_pos;
	
	//----------------------------------------------
	//  Get Current Vector (unscaled)
	//----------------------------------------------	
	
	for(unsigned int pos=0; pos< loop_pos.n_cols; pos++){
		loop_current(0,pos)= loop_pos(0,(pos+1)%loop_pos.n_cols)-loop_pos(0,pos);
		loop_current(1,pos)= loop_pos(1,(pos+1)%loop_pos.n_cols)-loop_pos(1,pos);
		loop_current(2,pos)= loop_pos(2,(pos+1)%loop_pos.n_cols)-loop_pos(2,pos);
	}
	
	if(debug){
		loop_pos.save("LoopPos.dat",arma::raw_ascii);
		loop_current.save("LoopCur.dat",arma::raw_ascii);
	}
	
	//----------------------------------------------
	//   Now do Biot-Savart to Solve for Field
	//----------------------------------------------
	
	float cx = Nx/2.0;
	float cy = Ny/2.0;
	float cz = Nz/2.0;
	
	int mapshrink = 4;
	SMAP= complex<float>(0.0,0.0);
	
	#pragma omp parallel for
	for(int k=0; k< (int)Nz; k+=mapshrink){
	for(int j=0; j< (int)Ny; j+=mapshrink){
	for(int i=0; i< (int)Nx; i+=mapshrink){
		complex<float>temp_smap(0.0,0.0); 
		
		// Array Position
		float x= ( (float)i - cx);
		float y= ( (float)j - cy);
		float z= ( (float)k - cz);
		 
		 // Objects have to be in the coil!	
		if( (x*x + y*y) > (0.64*coil_radius*coil_radius) ){
			continue;
		}
		
		// Store Position
		fvec S(3);
		S(0)=x;
		S(1)=y;
		S(2)=z;
		
		for(unsigned int pos=0; pos< loop_pos.n_cols; pos++){	
			fvec r = loop_pos.col(pos);
			fvec I = loop_current.col(pos);
			fvec diff =  S - r;
			float dist = (float)pow( (float)(diff(0)*diff(0) + diff(1)*diff(1) + diff(2)*diff(2)),(float)(3.0/2.0));			
			
			fvec temp = cross(I,diff);  
			temp_smap +=complex<float>( temp(0)/dist,temp(1)/dist);
		}
		
		// Now add to higher resolution smap
		int sz = k - mapshrink + 1;
		int ez = k + mapshrink - 1;
		sz = ( sz < 0 ) ? ( 0 ) : ( sz);
		ez = ( ez >=Nz ) ? ( Nz-1 ) : ( ez); 
		for(int kk=sz; kk<=ez; kk++){
			float wz = 1 - abs((float)kk-(float)k)/(float)mapshrink;
			
						
			int sy = j - mapshrink + 1;
			int ey = j + mapshrink - 1;
			sy = ( sy < 0 ) ? ( 0 ) : ( sy);
			ey = ( ey >=Ny ) ? ( Ny-1 ) : ( ey); 
			for(int jj=sy; jj<=ey; jj++){
				float wy = wz*( 1 - abs((float)jj-(float)j)/(float)mapshrink);
				
				int sx = i - mapshrink + 1;
				int ex = i + mapshrink - 1;
				sx = ( sx < 0 ) ? ( 0 ) : ( sx);
				ex = ( ex >=Nx ) ? ( Nx-1 ) : ( ex); 
				for(int ii=sx; ii<=ex; ii++){
					float wx = wy*( 1 - abs((float)ii-(float)i)/(float)mapshrink);
			
				
					complex<float>temp2 = wx*temp_smap;
					float RD = real(temp2);
					float ID = imag(temp2);
					float *I = reinterpret_cast<float *>(&SMAP(ii,jj,kk));
					float *R = I++;
					
					// Prevent Race conditions in multi-threaded
					#pragma omp atomic
					*R+=RD;
					
					#pragma omp atomic
					*I+=ID;
		}}}
			
				
	}}}

  
  	if(debug==1){ 
		cout << "Write";
		ArrayWriteMag(SMAP,"Smap_Resolution.dat");
	}
}


// 
// Add some phase to the phantom this could be better
//
void PHANTOM::add_phase( void ){	

	for(int k=0; k<SMAP.length(thirdDim); k++){
	for(int j=0; j<SMAP.length(secondDim); j++){
	for(int i=0; i<SMAP.length(firstDim); i++){
		float temp_phase = i*PI/Nx;
		IMAGE(i,j,k) = IMAGE(i,j,k)*polar((float)1.0, temp_phase);

	}}}

}

// 
//	 Add Noise to the phantom
//
void PHANTOM::add_noise( Array<complex<float>,5>&kdata){
	
	// Estimate 
	float mean_kdata = sum(abs(kdata)) / kdata.size();
	cout << "Mean kdata = "<< mean_kdata << endl;
	
	// Add Complex Noise  
	float noise_level = phantom_noise/sqrt(2.0)*mean_kdata;
	cout << "Noise Level = "<< noise_level << endl;
	
	for(int c=0; c<kdata.extent(fifthDim); c++){
		for(int ee=0; ee< kdata.extent(fourthDim); ee++){
			for(int k=0; k< kdata.extent(thirdDim); k++){
				for(int j=0; j< kdata.extent(secondDim); j++){
					for(int i=0; i< kdata.extent(firstDim); i++){
						arma::vec N = arma::randn<vec>(2);
						kdata(i,j,k,ee,c)+= complex<float>(N(0)*noise_level,N(1)*noise_level);
	}}}}}
}


// 
//	 Class to keep track of fractal tree
//
class TFRACT_RAND{
  public:
	fvec start;
	fvec stop;
	int root;
	int children[2];
	float Q;
	float l;
	float cost;
	float time;
    
	void update_cost( float qPower){
		l = norm( stop - start,2);
		cost = l*pow(Q,qPower);
	}
	
};

Array< int,3 > PHANTOM::synthetic_perfusion(int xs, int ys, int zs, PerfType ptype){
				
	
	// Input Perfusion Image
	Array< int,3>perf;
	perf.setStorage( ColumnMajorArray<3>());
	
	// Now Create Synthetic Image to Start
	switch(ptype){
		case(ELLIPSE):{
			
  			perf.resize(xs,ys,zs);
  			perf = 0;
			
			for(int k=0; k<perf.extent(thirdDim); k++){
			for(int j=0; j<perf.extent(secondDim); j++){
			for(int i=0; i<perf.extent(firstDim); i++){
			
			float x = (i - perf.extent(firstDim)*0.5)/perf.extent(firstDim);
			float y = (j - perf.extent(secondDim)*0.5)/perf.extent(secondDim);
			float z = (k - perf.extent(thirdDim)*0.5)/perf.extent(thirdDim);
	
			if( ( (x*x)/(0.5*0.5) + (y*y)/(0.35*0.35) + (z*z)/(0.25*0.25) ) < 0.75){
				perf(i,j,k)=1;
			}else if( ( (x*x)/(0.5*0.5) + (y*y)/(0.35*0.35) + (z*z)/(0.25*0.25) ) < 1){
				perf(i,j,k)=4;
			}  
			
			}}}
			
		}break;
			
		case(ASL):{
			// Input Perfusion Image
			perf.resize(128,128,128);
  			
			ArrayRead(perf,"ASL_perf.dat");
			int max_perf = max(perf);
			cout << "Max of Perf = " << max_perf << endl;
		}break;
	
		
	}
	
	if(debug==1){
		ArrayWrite(perf,"Perf.dat");	
	}
	
	return(perf);

}

void update_children(field<TFRACT_RAND>&tree, int pos){
	if( tree(pos).children[0] != -1){
		int cpos = tree(pos).children[0];
		float rsquare = pow( tree(pos).Q,2.0f/3.0f);
		float vel = tree(pos).Q/rsquare;
		tree(cpos).time =tree(pos).time + tree(pos).l/vel;
		update_children(tree,cpos);	
	}
	
	if( tree(pos).children[0] != -1){
		int cpos = tree(pos).children[1];
		float rsquare = pow( tree(pos).Q,2.0f/3.0f);
		float vel = tree(pos).Q/rsquare;
		tree(cpos).time =tree(pos).time + tree(pos).l/vel;
		update_children(tree,cpos);		
	}
} 


field<TFRACT_RAND> create_tree( fmat seeds_start, fmat seeds_stop,fmat X){ 
	
	// Variables						
	int Nseeds = seeds_start.n_cols;
	int terminal_pts = X.n_cols;
	float qPower = (2.0/3.0);
	 
	field<TFRACT_RAND>tree(terminal_pts*2+Nseeds);
	fvec all_costs(terminal_pts*2+Nseeds);
	int count = Nseeds;

	// Generate Initial Tree
	for( int s =0; s<Nseeds; s++){	
		tree(s).start= seeds_start.col(s);
		tree(s).stop= seeds_stop.col(s);
		tree(s).children[0] =-1; /*No Children*/
		tree(s).children[1] =-1; 
		tree(s).root = -1;
		tree(s).Q = 1;
		tree(s).update_cost(qPower);
	}
	
	// Loop to generate points
	for(int n=0; n < terminal_pts; n++){
		
		if(n%100==0){
			cout << "Starting tree : " << n << " of " << terminal_pts << endl;
		}		
		// Location of Perfusion Point
		fvec XT = X.col(n);
			
		// Search Loop for best connection
		float best_cost=0;
		int best_idx=0;
		
		#pragma omp parallel for
		for(int pos=0; pos < count; pos++){
			
			float delta_cost = -tree(pos).cost;
						
			float wTrunk = pow( tree(pos).Q + 1.0f, qPower);
			float wBranch = pow( tree(pos).Q, qPower);
			fvec centroid = ( XT + wTrunk*tree(pos).start + wBranch*tree(pos).stop);
			centroid *= 1.0f/(1.0f + wTrunk + wBranch);
			
			// Old Segment Branch
			float old_branch_l =  norm(centroid - tree(pos).stop,2);
			delta_cost += old_branch_l*pow( tree(pos).Q,qPower);
			
			// New Segment Branch
			float new_branch_l  = norm( centroid - XT,2);
			delta_cost += new_branch_l;
		
			//Old Segment (parent)
			float old_root_l  = norm( centroid - tree(pos).start,2);
			delta_cost += old_root_l*pow(tree(pos).Q+1,qPower); 
						
			// Add new flow from perfusion point
			int to_initial =0;
			int loc = pos;
			while( to_initial == 0){
				int root = tree(loc).root;
				if( root == -1){
					to_initial =1;
				}else{
					//tree2(root).Q+=1.0;
					//tree2(root).update_cost(qPower);
					delta_cost += tree(root).l*( pow(tree(root).Q+1,qPower) - pow(tree(root).Q,qPower));		
					loc = root;
				}
			}
			
			/// Sum Cost
			all_costs(pos)=delta_cost;
			
		}// Count
		
		// Find Best Position
		for(int pos=0; pos < count; pos++){
			if( (all_costs(pos)<best_cost) || (pos==0)){
				best_idx = pos;
				best_cost =all_costs(pos);
			}
		}
		
		int pos = best_idx;
		
		// Now Update to tree
		float wTrunk = pow( tree(pos).Q + 1.0f, qPower);
		float wBranch = pow( tree(pos).Q, qPower);
		fvec centroid = ( XT + wTrunk*tree(pos).start + wBranch*tree(pos).stop);
		centroid *= 1.0f/(1.0f + wTrunk + wBranch);
			
		// Old Segment Branch (Same children as old parent)
		tree(count).start = centroid;
		tree(count).stop = tree(pos).stop;
		tree(count).Q = tree(pos).Q;
		tree(count).root = pos;
		tree(count).update_cost(qPower);
		tree(count).children[0] = tree(pos).children[0];
		tree(count).children[1] = tree(pos).children[1];
									
		// Update Parent for Downstream Branches
		for(int jj=0; jj<count; jj++){
			if( tree(jj).root ==pos){
				tree(jj).root = count;	
			}
		}
			
		// New Segment Branch			
		tree(count+1).start = centroid;
		tree(count+1).stop = XT;
		tree(count+1).Q = 1;
		tree(count+1).root = pos;
		tree(count+1).update_cost(qPower);
		tree(count+1).children[0] = -1;
		tree(count+1).children[1] = -1;
			
		//Old Segment (parent)
		tree(pos).stop = centroid;
		tree(pos).children[0] = count;
		tree(pos).children[1] = count+1;
				
		// Add new flow from perfusion point
		int to_initial =0;
		int loc = count +1;
		while( to_initial == 0){
			int root = tree(loc).root;
			if( root== -1){
				to_initial =1;
			}else{
				tree(root).Q+=1.0;
				tree(root).update_cost(qPower);
				loc = root;
			}
		}
		count +=2;
		
	}// Terminal Pts
	
	/*Get the time*/	
	for(int pos=0; pos < count; pos++){
		if(tree(pos).root==-1){
			tree(pos).time = 0;
			
			//Iterate Through Children (recursive)
			update_children(tree,pos);			
		}
	}
	
	return(tree);
}



// 
//	 Generate Fractal Based on Random Generation within Volume
//     Similar to Karch et al. Computers in Biology and Medicine 29(1999)19-38
//
void  PHANTOM::fractal3D_new(int Nx, int Ny, int Nz){
	
	// Set number of threads 
	if( omp_get_max_threads() > 8){
		omp_set_num_threads(omp_get_max_threads()-2);
	}
		
	// Inputs
	int terminal_pts = fractal_pts; // Number of Endpoints
		
	// Input Perfusion Image
	Array< int,3>perf = synthetic_perfusion(128,128,128,ASL);
		
	// Convert The Perf to Random Ordered Point
	int perf_max = max(perf);
	cout << "Max perf : " << perf_max << endl;
	
	// Convert to points	
	fmat X= zeros<fmat>(3,terminal_pts);
	for(int pos=0; pos < terminal_pts; pos++){
		int found =0;
		while(found==0){
			int k =  (int)(( (float)perf.extent(thirdDim)*(float)rand() ) / ((float)RAND_MAX));
			int j =  (int)(( (float)perf.extent(secondDim)*(float)rand() ) / ((float)RAND_MAX));
			int i =  (int)(( (float)perf.extent(firstDim)*(float)rand() ) / ((float)RAND_MAX));
			int val = (int)(( (float)perf_max*(float)rand() ) / ((float)RAND_MAX));
			if(val < perf(i,j,k)){
				X(0,pos) = ((float)i + (float)rand()/(float)RAND_MAX - 0.5f - (float)perf.extent(firstDim)*0.5)/(float)perf.extent(firstDim);
				X(1,pos) = ((float)j + (float)rand()/(float)RAND_MAX - 0.5f- (float)perf.extent(secondDim)*0.5)/(float)perf.extent(secondDim);
				X(2,pos) = ((float)k + (float)rand()/(float)RAND_MAX - 0.5f- (float)perf.extent(thirdDim)*0.5)/(float)perf.extent(thirdDim);
				found =1;
			}
		}
	} 	
	
	// ---------------------------------------------
	//    Arterial Tree
	// ---------------------------------------------
	
	cout << "Arterial Tree" << endl << flush;
	// Get Seed Points	
	int Nseeds =3;
	fmat seeds_start= zeros<fmat>(3,Nseeds);
	fmat seeds_stop= zeros<fmat>(3,Nseeds);
	seeds_start(0,0) =  -0.05;
	seeds_start(1,0) =  0.1;
	seeds_start(2,0) = -0.3;
	
	seeds_stop(0,0) =  -0.05;
	seeds_stop(1,0) =  0.1;
	seeds_stop(2,0) =  0.0;
	
	
	seeds_start(0,1) = -0.05;
	seeds_start(1,1) = -0.1;
	seeds_start(2,1) = -0.3;
	seeds_stop(0,1) =  -0.05;
	seeds_stop(1,1) = -0.1;
	seeds_stop(2,1) =  0.0;
	
	seeds_start(0,2) =  0.05;
	seeds_start(1,2) =  0.0;
	seeds_start(2,2) =  -0.3;
	seeds_stop(0,2) =  0.05;
	seeds_stop(1,2) =  0.0;
	seeds_stop(2,2) =  -0.28;
	
	// Arterial Tree 
	field<TFRACT_RAND>arterial_tree = create_tree(seeds_start,seeds_stop,X);
	
	// ---------------------------------------------
	//    Venous Tree
	// ---------------------------------------------
	
	cout << "Venous Tree" << endl << flush;
	// Get Seed Points	
	fmat Vseeds_start= zeros<fmat>(3,1);
	fmat Vseeds_stop= zeros<fmat>(3,1);
	
	// Venous Seeds
	Vseeds_start(0,0) =  0.25;
	Vseeds_start(1,0) =  0.0;
	Vseeds_start(2,0) =  -0.3;
	Vseeds_stop(0,0) =  seeds_stop(0,0);
	Vseeds_stop(1,0) =  seeds_stop(1,0);
	Vseeds_stop(2,0) =  seeds_stop(2,0);
	
	// Need to add arterial seed points to X
	X.resize(3,terminal_pts+Nseeds-1); 
	for(int t=1; t<Nseeds; t++){
		X(0,terminal_pts+t-1)=seeds_stop(0,t);
		X(1,terminal_pts+t-1)=seeds_stop(1,t);
		X(2,terminal_pts+t-1)=seeds_stop(2,t);
	}
	
	// Venous Tree
	field<TFRACT_RAND>venous_tree = create_tree(Vseeds_start,Vseeds_stop,X);
	
	// ---------------------------------------------
	//    Tissue Tree
	// ---------------------------------------------
	cout << "Tissue Tree" << endl << flush;
	field<TFRACT_RAND>tissue_tree(terminal_pts+Nseeds);
	int cpos =0;
	for( int t=0; t< (int)arterial_tree.n_elem; t++){
		if(arterial_tree(t).Q==1){
			tissue_tree(cpos).stop= arterial_tree(t).stop;
			tissue_tree(cpos).start= arterial_tree(t).stop;
			tissue_tree(cpos).time= arterial_tree(t).time;
			cpos++;
		}
	}	
			
	// ---------------------------------------------
	//    Now Grid The Data
	// ---------------------------------------------

	float tissue_radius = 0.05;
	// Phantom Density
	
	// Arterial, Perfusion, and Venous Compartments
	FUZZY.setStorage( ColumnMajorArray<4>());
  	FUZZY.resize(Nx,Ny,Nz,3);
  	FUZZY = 0;
	
	// Arterial, Perfusion, and Venous Compartments
	FUZZYT.setStorage( ColumnMajorArray<4>());
  	FUZZYT.resize(Nx,Ny,Nz,3);
  	FUZZYT = 0;
		
	
	// Find Min/Max
	float SX = 10000;
	float SY = 10000;
	float SZ = 10000;
	float EX = 0;
	float EY = 0;
	float EZ = 0;
	float Qmax = 0.0;
	for(unsigned int t=0; t< arterial_tree.n_elem; t++){
		SX = min( arterial_tree(t).start(0),SX);
		SY = min( arterial_tree(t).start(1),SY);
		SZ = min( arterial_tree(t).start(2),SZ);
		EX = max( arterial_tree(t).start(0),EX);
		EY = max( arterial_tree(t).start(1),EY);
		EZ = max( arterial_tree(t).start(2),EZ);
		SX = min( arterial_tree(t).stop(0),SX);
		SY = min( arterial_tree(t).stop(1),SY);
		SZ = min( arterial_tree(t).stop(2),SZ);
		EX = max( arterial_tree(t).stop(0),EX);
		EY = max( arterial_tree(t).stop(1),EY);
		EZ = max( arterial_tree(t).stop(2),EZ);
		Qmax = max(arterial_tree[t].Q,Qmax);
	}
	float Qscale = 0.010 / pow( Qmax,1.0f/3.0f);
	float Rmin =0;
	float Rmax = max(0.010f,tissue_radius);
		
	cout << "Range is" << endl;
	cout << "\t X:" << SX << "to " << EX << endl;
	cout << "\t Y:" << SY << "to " << EY << endl;
	cout << "\t Z:" << SZ << "to " << EZ << endl;
	cout << "\t R:" << Rmin << "to " << Rmax << endl;
	
	float scaleX = 0.95*((float)Nx)/( 2*tissue_radius + EX - SX );
	float scaleY = 0.95*((float)Ny)/( 2*tissue_radius + EY - SY );
	float scaleZ = 0.95*((float)Nz)/( 2*tissue_radius + EZ - SZ );
	float scale = min(min(scaleX,scaleY),scaleZ);
		
	float CX = scale*(EX + SX)/2 - Nx/2;
	float CY = scale*(EY + SY)/2 - Ny/2;
	float CZ = scale*(EZ + SZ)/2 - Nz/2;
			
	cout << "Gridding Tree" << endl;
	cout << " scale = " <<  scale << endl;
	
	// ---------------------------------------------
	//    Arterial
	// ---------------------------------------------

	#pragma omp parallel for	
	for(int t=0; t< (int)arterial_tree.n_elem; t++){
  								
		// Get  Length	
		fvec vessel_dir = scale*(arterial_tree(t).stop-arterial_tree(t).start);
		float vessel_length = norm(vessel_dir,2);
		int n_vessel = (int)round(vessel_length /0.1);
		float R= scale*Qscale*pow( arterial_tree(t).Q,1.0f/3.0f);
		
		float rsquare = pow( arterial_tree(t).Q,2.0f/3.0f);
		float vel = arterial_tree(t).Q/rsquare;
		float total_time =arterial_tree(t).l/vel;					
		
		// Vessels
		if( n_vessel > 1){
		for(int ii=0; ii<n_vessel; ii++){
			
			fvec vessel_pos = scale*arterial_tree(t).start + vessel_dir*( (float)ii/(n_vessel-1));
			
			float time = arterial_tree(t).time + total_time*( (float)ii/(n_vessel-1));
			
			int xx = (int)round(vessel_pos(0) - CX);
        	int yy = (int)round(vessel_pos(1) - CY);
        	int zz = (int)round(vessel_pos(2) - CZ);
			float dx = vessel_pos(0) - CX - (float)xx;
			float dy = vessel_pos(1) - CY - (float)yy;
			float dz = vessel_pos(2) - CZ - (float)zz;
			int rr = (int)ceil(R);
        	
			// cout << "Vessel Position = " << xx << "," << yy << "," << zz << " R = " << rr << endl << flush;
		
			// cout << "Position = " << xx << "," << yy << "," << zz << " R = " << rr << endl << flush;
			for(int k=-rr; k<=rr; k++){
			for(int j=-rr; j<=rr; j++){
			for(int i=-rr; i<=rr; i++){
				
				float iact = (float)i - dx;
				float jact = (float)j - dy;
				float kact = (float)k - dz;
				
				float radius = sqrt( iact*iact + jact*jact + kact*kact);
				float s;
				if(R < 1.5){
				   // For small vessels below resolution of matrix need to normalize size
				   s=(float)(1.0/( 1.0 + exp(2.0*(radius-(float)R))) )*( R*R /2.25);
				}else{
				   s=(float)(1.0/( 1.0 + exp(2.0*(radius-(float)R))) );
				}
				
				if( s > FUZZY( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,0) ){
					FUZZY( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,0) = s;
					FUZZYT( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,0)= s*time;
				}
				
			}}}
		}
		}
	}// Arterial tree
	

	// ---------------------------------------------
	//    Venous 
	// ---------------------------------------------

	#pragma omp parallel for	
	for( int t=0; t< (int)venous_tree.n_elem; t++){
  								
		// Get  Length	
		fvec vessel_dir = scale*(venous_tree(t).stop-venous_tree(t).start);
		float vessel_length = norm(vessel_dir,2);
		int n_vessel = (int)round(vessel_length /0.1);
		float R= scale*Qscale*pow( venous_tree(t).Q,1.0f/3.0f);
		
		float rsquare = pow( venous_tree(t).Q,2.0f/3.0f);
		float vel = venous_tree(t).Q/rsquare;
		float total_time =venous_tree(t).l/vel;					
		
		// Vessels
		if( n_vessel > 1){
		for(int ii=0; ii<n_vessel; ii++){
			
			fvec vessel_pos = scale*venous_tree(t).start + vessel_dir*( (float)ii/(n_vessel-1));
			
			float time = venous_tree(t).time + total_time*( (float)ii/(n_vessel-1));
			
			int xx = (int)round(vessel_pos(0) - CX);
        	int yy = (int)round(vessel_pos(1) - CY);
        	int zz = (int)round(vessel_pos(2) - CZ);
			float dx = vessel_pos(0) - CX - (float)xx;
			float dy = vessel_pos(1) - CY - (float)yy;
			float dz = vessel_pos(2) - CZ - (float)zz;
			int rr = (int)ceil(R);
        	
			// cout << "Vessel Position = " << xx << "," << yy << "," << zz << " R = " << rr << endl << flush;
		
			// cout << "Position = " << xx << "," << yy << "," << zz << " R = " << rr << endl << flush;
			for(int k=-rr; k<=rr; k++){
			for(int j=-rr; j<=rr; j++){
			for(int i=-rr; i<=rr; i++){
				
				float iact = (float)i - dx;
				float jact = (float)j - dy;
				float kact = (float)k - dz;
				
				float radius = sqrt( iact*iact + jact*jact + kact*kact);
				float s;
				if(R < 1.5){
				   // For small vessels below resolution of matrix need to normalize size
				   s=(float)(1.0/( 1.0 + exp(2.0*(radius-(float)R))) )*( R*R /2.25);
				}else{
				   s=(float)(1.0/( 1.0 + exp(2.0*(radius-(float)R))) );
				}
				
				if( s > FUZZY( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,2) ){
					FUZZY( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,2) = s;
					FUZZYT( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,2)= s*time;
				}
				
			}}}
		}
		}
	}// Venous_tree
	
	
	// ---------------------------------------------
	//    Tissue
	// ---------------------------------------------

	#pragma omp parallel for	
	for( int t=0; t< (int)tissue_tree.n_elem; t++){
		
			fvec vessel_pos = scale*tissue_tree(t).start;
			
			float R= scale*tissue_radius;
						
			int xx = (int)round(vessel_pos(0) - CX);
        	int yy = (int)round(vessel_pos(1) - CY);
        	int zz = (int)round(vessel_pos(2) - CZ);
			float dx = vessel_pos(0) - CX - (float)xx;
			float dy = vessel_pos(1) - CY - (float)yy;
			float dz = vessel_pos(2) - CZ - (float)zz;
			float time = tissue_tree(t).time;
			int rr = (int)ceil(R);
        	// cout << "Tissue Position = " << xx << "," << yy << "," << zz << " R = " << rr << endl << flush;
								
			for(int k=-rr; k<=rr; k++){
			for(int j=-rr; j<=rr; j++){
			for(int i=-rr; i<=rr; i++){
				
				float iact = (float)i - dx;
				float jact = (float)j - dy;
				float kact = (float)k - dz;
						
				float radius = sqrt( iact*iact + jact*jact + kact*kact);
				float s=(float)(1.0/( 1.0 + exp(8.0*(radius-0.5*(float)R))) );
				FUZZY( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,1) += s;
				FUZZYT( (i+xx)%Nx,(j+yy)%Ny,(k+zz)%Nz,1) += s*time;
				
			}}}
	}
	
	// Times Are Weighted in Averaging Calcs
	FUZZYT/=(FUZZY+1e-6*max(FUZZY));
		
	// Adjust Perfusion by Max
	// CBV is about 10%
	Array< float, 3>Tissue = FUZZY(Range::all(),Range::all(),Range::all(),1);
	Array< float, 3>Arterial = FUZZY(Range::all(),Range::all(),Range::all(),0);
	Tissue*= 0.1*max(Arterial)/max(Tissue);
	
	// Normalize times to 0.5
	FUZZYT/= 2.0*max(FUZZYT);
	
	// Make venous go backwards and rescale to 1
	FUZZYT(Range::all(),Range::all(),Range::all(),2 )= 1.0-FUZZYT(Range::all(),Range::all(),Range::all(),2 );
	
		
	cout << "Writing to files " << endl;
	if(debug==1){
		ArrayWrite(FUZZY,"Phantom.dat"); 
		ArrayWrite(FUZZYT,"PhantomT.dat"); 
	}
	
	// Export Magnitude
	if(debug==1){
		ArrayWriteMag(IMAGE,"Phantom.dat"); 
		ArrayWrite(TOA,"TOA.dat"); 
	}
	
	}



inline float Factorial(float x) {
  return (x == 1.0 ? x : x * Factorial(x - 1.0));
}

inline float gamma_variate(float x,float beta,float alpha){
 	if(x < 0.0){
		return(0.0);
	}else{
		return(  (1.0/Factorial(alpha)/pow(beta,alpha))*pow(x,alpha)*exp(-x/beta) );
	}
}

void PHANTOM::write_matlab_truth_script( char *folder){
	
	char fname[1024];
	sprintf(fname,"%s%s",folder,"get_pixel_timecourse.m");
	
	ofstream script;
	script.open( fname);
	script << "function st = get_pixel_timecourse( i,j,k ) " << endl;
	script << "Nx = " << Nx << ";" << endl;
	script << "Ny = " << Ny << ";" << endl;
	script << "Nz = " << Nz << ";" << endl;
	script << "Nt = " << Nt << ";" << endl;
	
	// Read in a script (density)
	script << "fid = fopen('Phantom.dat','r'); " << endl;
	script << "fseek(fid,( Nx*Ny*(k-1)+Nx*(j-1)+(i-1))*4,'bof'); " << endl;
	script << "dens_a = fread(fid,1,'float'); " << endl;
	script << "fseek(fid,(Nx*Ny*Nz + Nx*Ny*(k-1)+Nx*(j-1)+(i-1))*4,'bof'); " << endl;
	script << "dens_v = fread(fid,1,'float'); " << endl;
	script << "fseek(fid,(2*Nx*Ny*Nz + Nx*Ny*(k-1)+Nx*(j-1)+(i-1))*4,'bof'); " << endl;
	script << "dens_t = fread(fid,1,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	// Read in a script (TOA)
	script << "fid = fopen('PhantomT.dat','r'); " << endl;
	script << "fseek(fid,( Nx*Ny*(k-1)+Nx*(j-1)+(i-1))*4,'bof'); " << endl;
	script << "toa_a = fread(fid,1,'float'); " << endl;
	script << "fseek(fid,(Nx*Ny*Nz + Nx*Ny*(k-1)+Nx*(j-1)+(i-1))*4,'bof'); " << endl;
	script << "toa_v = fread(fid,1,'float'); " << endl;
	script << "fseek(fid,(2*Nx*Ny*Nz + Nx*Ny*(k-1)+Nx*(j-1)+(i-1))*4,'bof'); " << endl;
	script << "toa_t = fread(fid,1,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	// Calc signal
	script << "t = 3.0* (0:Nt-1)/Nt;" << endl;
	script << "st =t*0;" << endl;
	script << "beta = " << beta << endl;
	script << "alpha = " << alpha << endl;
	script << "st = st + dens_a*gamma_variate((t-toa_a),beta*(1+2*toa_a),alpha);" << endl;		
	script << "st = st + dens_v*gamma_variate((t-toa_v),beta*(1+2*toa_v),alpha);" << endl;		
	script << "st = st + dens_t*gamma_variate((t-toa_t),beta*(1+2*toa_t),alpha);" << endl;		
	
	script << "function s = gamma_variate( x,beta,alpha) " << endl;
	script << "s = double(x>0)*(  (1.0./factorial(alpha)/(beta^alpha)).*(x.^alpha).*exp(-x/beta));" << endl;
	script.close();
	
		
	sprintf(fname,"%s%s",folder,"get_truth_image.m");
	script.open(fname);
	script << "function s = get_truth_image( tt ) " << endl;
	script << "%t=position in waveform (0 < tt < 1)" << endl;
	script << "Nx = " << Nx << ";" << endl;
	script << "Ny = " << Ny << ";" << endl;
	script << "Nz = " << Nz << ";" << endl;
	script << "Nt = " << Nt << ";" << endl;
	script << "beta = " << beta << endl;
	script << "alpha = " << alpha << endl;
	script << "t = tt*3.0;" << endl;
	
	// Arterial
	script << "fid = fopen('Phantom.dat','r'); " << endl;
	script << "dens = fread(fid,Nx*Ny*Nz,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	script << "fid = fopen('PhantomT.dat','r'); " << endl;
	script << "toa = fread(fid,Nx*Ny*Nz,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	script << "s = dens.*gamma_variate((t-toa),beta*(1+2*toa),alpha);" << endl;
		
	// Venous
	script << "fid = fopen('Phantom.dat','r'); " << endl;
	script << "fseek(fid,Nx*Ny*Nz*4,'bof'); " << endl;
	script << "dens = fread(fid,Nx*Ny*Nz,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	script << "fid = fopen('PhantomT.dat','r'); " << endl;
	script << "fseek(fid,Nx*Ny*Nz*4,'bof'); " << endl;
	script << "toa = fread(fid,Nx*Ny*Nz,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	script << "s = s + dens.*gamma_variate((t-toa),beta*(1+2*toa),alpha);" << endl;
		
	// Tissue
	script << "fid = fopen('Phantom.dat','r'); " << endl;
	script << "fseek(fid,2*Nx*Ny*Nz*4,'bof'); " << endl;
	script << "dens = fread(fid,Nx*Ny*Nz,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	script << "fid = fopen('PhantomT.dat','r'); " << endl;
	script << "fseek(fid,2*Nx*Ny*Nz*4,'bof'); " << endl;
	script << "toa = fread(fid,Nx*Ny*Nz,'float'); " << endl;
	script << "fclose(fid);" << endl;
	
	script << "s = s + dens.*gamma_variate((t-toa),beta*(1+2*toa),alpha);" << endl;
	script << "s = reshape(s,[" << Nx << " " << Ny << " " << Nz << "]);" << endl;
	
	script << "function s = gamma_variate( x,beta,alpha) " << endl;
	script << "s = double(x>0).*(  (1.0./factorial(alpha)./(beta.^alpha)).*(x.^alpha).*exp(-x./beta));" << endl;
	
	script.close();
	
	
	
}

	
// ----------------------
// Create Image (simple for test)
// ----------------------
void  PHANTOM::calc_image(int t, int frames){
	switch(phantom_type){
	
	case(FRACTAL):{	
	float current_time = (float)t* 3.0 / (float)frames;
	IMAGE = complex<float>(0,0);
	
	beta = 0.25;
	alpha = 2;
	
	for(int compartment =0; compartment < FUZZYT.length(fourthDim); compartment++){
	#pragma omp parallel for
	for(int k= 0; k< FUZZYT.length(thirdDim); k++){
	for(int j= 0; j< FUZZYT.length(secondDim); j++){
	for(int i= 0; i< FUZZYT.length(firstDim); i++){
		float time = FUZZYT(i,j,k,compartment);
		float dens = FUZZY(i,j,k,compartment);
		float beta_t = ( (1+2*time)*beta); // Curve Smooths in time
		
		if( current_time > time){
			float tdiff = (current_time - time);
			IMAGE(i,j,k) += dens*gamma_variate(tdiff,beta_t,alpha); // Bolus
			IMAGE(i,j,k) += 0.25*dens*gamma_variate(tdiff,4*beta_t,alpha); // Persistent signal 
		}			
	}}}}
	}break;
	
	default:
	case(EXTERNAL):{
		// No changes
	}break;
	
	}
	
	if(debug){
		char fname[1024];
		sprintf(fname,"Frame_%d.dat",t);
		ArrayWriteMag(IMAGE,fname);
	}	
	// Multiply by Sensitivity
	IMAGE *= SMAP;

}
// ----------------------
// Help Message
// ----------------------
#include "io_templates.cpp"
void PHANTOM::help_message(void){
	cout << "----------------------------------------------" << endl;
	cout << "   Phantom Control " << endl;
	cout << "----------------------------------------------" << endl;
	
	cout<<"Control" << endl;
	help_flag("-phantom_noise []","Noise as fraction of mean of abs(kdata)");
	help_flag("-fractal_pts []","Terminal Points for Fractal");
	help_flag("-external_phantom []","Use phantom of name []");
}


//----------------------------------------
// Phantom Read Command Line
//----------------------------------------
void PHANTOM::read_commandline(int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]); 

  for(int pos=0; pos < numarg; pos++){
  
  	if (strcmp("-h", pstring[pos] ) == 0) {
		float_flag("-phantom_noise",phantom_noise);
		int_flag("-fractal_pts",fractal_pts);
	}else if(strcmp("-external_phantom",pstring[pos]) == 0) {
		pos++;
		phantom_type = EXTERNAL;
		if( pos==numarg){
			cout << "Please provide external phantom name" << endl;
		}else{
			strcpy(external_phantom_name,pstring[pos]);
		}
  	}
  }  
}  
