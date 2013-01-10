/*

Title: Phantom Library   

Summary: This code exists to create simulated phantom data 

Usage: Currently 

	PHANTOM phantom;						
	phantom.read_commandline(argc,argv);  
	phantom.init(recon.rcxres,recon.rcyres,recon.rczres);
	phantom.update_smap_biotsavart(coil,recon.num_coils);	
	
	// Which will populate:
	phantom.IMAGE // Image (density)
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
void PHANTOM::init(int Nx, int Ny, int Nz){
	
	// Phantom Density
	IMAGE.setStorage( ColumnMajorArray<3>());
  	IMAGE.resize((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
  	IMAGE = 0;
	
	// Sensitivity Map	
	SMAP.setStorage( ColumnMajorArray<3>());
  	SMAP.resize((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
  	SMAP = 1;
		
	// Switch For Phantom
	switch(phantom_type){
		case(FRACTAL):{
			TOA.setStorage( ColumnMajorArray<3>());
  			TOA.resize((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
  			TOA = 0;
			
			fractal3D((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
		}break;
		
		case(SHEPP):{
			// Doesn't exist yet (but can be k-space based, Koay et al.)
		}break;
		
		case(PSF):{
			// Just Get PSF 
		}
	}
}


// Get a sensitivity based on coil geometry
// with Biot-Savart simulation of field.
//
void  PHANTOM::update_smap_biotsavart(int coil, int Ncoils){
	
	// Doesn't Make sense for one coil
	if(Ncoils==1){
		SMAP = 1;	
		return;
	}
	
	// Grab size for Short Hand
	float Nx = (float)SMAP.length(firstDim);
	float Ny = (float)SMAP.length(secondDim);
	float Nz = (float)SMAP.length(thirdDim);
		
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
		
	#pragma omp parallel for
	for(int k=0; k< (int)Nz; k++){
	for(int j=0; j< (int)Ny; j++){
	for(int i=0; i< (int)Nx; i++){
		SMAP(i,j,k)=complex<float>(0.0,0.0);
		
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
			SMAP(i,j,k)+=complex<float>( temp(0)/dist,temp(1)/dist);
		}
	}}}
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
class TFRACT{
  public:
	double r[3];
	double v[3];
	double vt[3];
	double l;
	int count;
	float rad;
	float blength;
	int tactive;
	int tn;
	
	vector<float> btime;
	vector<float> xlist;
	vector<float> ylist;
	vector<float> zlist;
	vector<float> rlist;
  	
	void resize(int n){
		xlist.resize(n);
		ylist.resize(n);
		zlist.resize(n);
		rlist.resize(n);
		btime.resize(n);
	};
  
};

// 
//	 Generate Base factal Tree 
//
void  PHANTOM::fractal3D(int Nx, int Ny, int Nz){

	int n = 13; // Maximum Branches
	float branch_angle= PI/5;
	float branch_rotation = PI/4; //(3-sqrt(5));
	float branch_fraction = 0.8;
	float branch_length = 60;
	float tissue_radius = 60;
	
	//Stores all the data
	std::vector<TFRACT>tree(1);
	tree[0].resize(1);
			
	//Initialize with starting tree
	cout << "Initializing Tree" << endl;
	int trees=1;
	int active_trees = 1;
	int time = 0;
	
	tree[0].r[0]=0.0;
	tree[0].r[1]=0.0;
	tree[0].r[2]=0.0;
		
	tree[0].v[0]=0.0;
	tree[0].v[1]=0.0;
	tree[0].v[2]=0.05;
	
	tree[0].vt[0]=1.0;
	tree[0].vt[1]=0.0;
	tree[0].vt[2]=0.0;
	
	tree[0].l=0.0;
	tree[0].count = 1;
	tree[0].rad = 16.0;	
	tree[0].blength = branch_length;	
	tree[0].tactive=1;
	tree[0].tn=1;
	
	tree[0].btime[0]=time;
	tree[0].xlist[0]=tree[0].r[0];
	tree[0].ylist[0]=tree[0].r[1];
	tree[0].zlist[0]=tree[0].r[2];
	tree[0].rlist[0]=tree[0].rad;
	
	while( active_trees > 0){
	
    	for(int t=0; t< trees; t++){
        	if(tree[t].tactive==0){
            	continue;
        	}
        	
			tree[t].count++;
			
			//Step for Tree
			tree[t].r[0]+=tree[t].v[0];
			tree[t].r[1]+=tree[t].v[1];
			tree[t].r[2]+=tree[t].v[2];
			
			// Length
			tree[t].l += sqrt( pow(tree[t].v[0],2.0) + pow(tree[t].v[1],2.0) + pow(tree[t].v[2],2.0));
							
			//Update List
        	tree[t].resize(tree[t].count);
			
			int ii = tree[t].count -1;
			tree[t].btime[ii]=time;
			tree[t].xlist[ii]=tree[t].r[0];
			tree[t].ylist[ii]=tree[t].r[1];
			tree[t].zlist[ii]=tree[t].r[2];
			tree[t].rlist[ii]=tree[t].rad;
			
			if(tree[t].l > tree[t].blength){
				if( tree[t].tn==n){
					tree[t].tactive=0;
					continue;
				}
				
				// cout << "Branching:" << tree[t].tn << endl;
				trees=trees+1;
				tree.resize(trees);
				
				//----------------------------------
				// Rotation About Velocity Axis
				//----------------------------------
				
				int nt = trees-1;
            	tree[t].tn++;
				tree[nt].tn = tree[t].tn;
				tree[nt].tactive=1;
				
				
				tree[nt].resize(1);
				tree[t].blength *=branch_fraction;
				tree[nt].blength=tree[t].blength;
				
				tree[t].rad/= pow(2.0,1.0/3.0);
				tree[nt].rad=tree[t].rad;
				
				tree[t].l=0;
				tree[nt].l=0;
				                       
            	tree[nt].count=1;
            	tree[nt].r[0]=tree[t].r[0];
				tree[nt].r[1]=tree[t].r[1];
				tree[nt].r[2]=tree[t].r[2];
				           					
				tree[nt].btime[0]=time;
				tree[nt].xlist[0]=tree[t].r[0];
				tree[nt].ylist[0]=tree[t].r[1];
				tree[nt].zlist[0]=tree[t].r[2];
				tree[nt].rlist[0]=tree[t].rad;
				
				//----------------------------------
				// Rotation About Velocity Axis
				//----------------------------------
            	float vmag = sqrt( pow(tree[t].v[0],2.0) + pow(tree[t].v[1],2.0) + pow(tree[t].v[2],2.0));
				float ux = tree[t].v[0]/vmag;
				float uy = tree[t].v[1]/vmag;
				float uz = tree[t].v[2]/vmag;
				
				mat utu,usk;
				utu << ux*ux << ux*uy << ux*uz << endr
				    << ux*uy << uy*uy << uy*uz << endr
					<< ux*uz << uy*uz << uz*uz << endr;
				
				usk << 0.0 << -uz <<  uy << endr 
					<< uz  << 0.0 << -ux << endr
					<< -uy << ux  << 0.0 << endr;
                
            
            	float rot = branch_rotation;
            	mat I;
				I.eye(3,3);
				mat Raxis = utu + cos(rot)*(I-utu) + sin(rot)*usk;
            	
				vec vt(3);
				vt(0)=tree[t].vt[0];
				vt(1)=tree[t].vt[1];
				vt(2)=tree[t].vt[2];
				vt = Raxis*vt;
							
				tree[t].vt[0]=vt(0);
				tree[t].vt[1]=vt(1);
				tree[t].vt[2]=vt(2);
				tree[nt].vt[0]=vt(0);
				tree[nt].vt[1]=vt(1);
				tree[nt].vt[2]=vt(2);
				
				//----------------------------------
				// Rotation About Tangential Axis
				//----------------------------------
				 ux = tree[t].vt[0];
				 uy = tree[t].vt[1];
				 uz = tree[t].vt[2];
				
				utu << ux*ux << ux*uy << ux*uz << endr
				    << ux*uy << uy*uy << uy*uz << endr
					<< ux*uz << uy*uz << uz*uz << endr;
				
				usk << 0.0 << -uz <<  uy << endr 
					<< uz  << 0.0 << -ux << endr
					<< -uy << ux  << 0.0 << endr;
                
            
            	rot = branch_angle;
            				
				// Branch 1
				Raxis = utu + cos(rot)*(I-utu) + sin(rot)*usk;
            	vec v(3);
				v(0)=tree[t].v[0];
				v(1)=tree[t].v[1];
				v(2)=tree[t].v[2];
				v = Raxis*v;
				tree[nt].v[0]=v(0);
				tree[nt].v[1]=v(1);
				tree[nt].v[2]=v(2);
				
				// Branch 2
				v(0)=tree[t].v[0];
				v(1)=tree[t].v[1];
				v(2)=tree[t].v[2];
				Raxis = utu + cos(-rot)*(I-utu) + sin(-rot)*usk;
				v = Raxis*v;
				tree[t].v[0]=v(0);
				tree[t].v[1]=v(1);
				tree[t].v[2]=v(2);
			}//if branch
		}//tree loop
		
		
		// Count Active trees
		active_trees =0;
		for(int t=0; t< trees; t++){
			active_trees += tree[t].tactive;
		}
		//cout << "Active trees = " << active_trees << endl;				
	}

	
	// Count the total points
	int total_points=0;
	for(int t=0; t< trees; t++){
		total_points+= tree[t].count;
	}
	cout << "Total Points= " << total_points << endl;
	
	// Find Min/Max
	float SX = 10000;
	float SY = 10000;
	float SZ = 10000;
	float EX = 0;
	float EY = 0;
	float EZ = 0;
	float Rmax = 0;
	float Rmin =1000;
	for(int t=0; t< trees; t++){
		for(int ii=0; ii<tree[t].count; ii++){
			SX = min( tree[t].xlist[ii],SX);
			SY = min( tree[t].ylist[ii],SY);
			SZ = min( tree[t].zlist[ii],SZ);
			EX = max( tree[t].xlist[ii],EX);
			EY = max( tree[t].ylist[ii],EY);
			EZ = max( tree[t].zlist[ii],EZ);
			Rmax = max( tree[t].rlist[ii],Rmax);
			Rmin = min( tree[t].rlist[ii],Rmin);
		}
	}
	cout << "Range is" << endl;
	cout << "\t X:" << SX << "to " << EX << endl;
	cout << "\t Y:" << SY << "to " << EY << endl;
	cout << "\t Z:" << SZ << "to " << EZ << endl;
	cout << "\t R:" << Rmin << "to " << Rmax << endl;
	float res = (float)Nx;
	
	float scaleX = 0.95*res/(tissue_radius + EX - SX);
	float scaleY = 0.95*res/(tissue_radius + EY - SY);
	float scaleZ = 0.95*res/(tissue_radius + EZ - SZ);
	float scale = min(min(scaleX,scaleY),scaleZ);

	float CX = scale*(EX + SX)/2 - res/2;
	float CY = scale*(EY + SY)/2 - res/2;
	float CZ = scale*(EZ + SZ)/2 - res/2;
			
	cout << "Gridding Tree" << endl;
	int count =0;
	
	for(int t=0; t< trees; t++){
    	ostringstream temp;
		temp << count << " of " << trees << endl;
		cout << temp.str() << flush;
		
		#pragma omp atomic
		count++;
		
		// Vessels
		for(int ii=0; ii<tree[t].count; ii++){
        	float X=scale*tree[t].xlist[ii];
    		float Y=scale*tree[t].ylist[ii];
    		float Z=scale*tree[t].zlist[ii];
    		float R=scale*tree[t].rlist[ii];
			float T=tree[t].btime[ii];
			
			int xx = (int)round(X - CX);
        	int yy = (int)round(Y - CY);
        	int zz = (int)round(Z - CZ);
			int rr = (int)ceil(R);
        							
			for(int k=-rr; k<=rr; k++){
			for(int j=-rr; j<=rr; j++){
			for(int i=-rr; i<=rr; i++){
				float radius = sqrt( (float)(i*i) + (float)(j*j) + (float)(k*k)); 
				float s=(float)(1.0/( 1.0 + exp(2.0*(radius-0.5*(float)R))) );
				
				if( s > real(IMAGE(i+xx,j+yy,k+zz)) ){
					IMAGE(i+xx,j+yy,k+zz) = complex<float>( s,0.0); 
					TOA(i+xx,j+yy,k+zz) = T;
				}
			}}}
		}
		
		// Tissue 
		int ii=tree[t].count-1;
        float X=scale*tree[t].xlist[ii];
    	float Y=scale*tree[t].ylist[ii];
    	float Z=scale*tree[t].zlist[ii];
    	float R=scale*tissue_radius;
		float T=tree[t].btime[ii];
			
		int xx = (int)floor(X - CX);
        int yy = (int)floor(Y - CY);
        int zz = (int)floor(Z - CZ);
		int rr = (int)ceil(R);
       	
		for(int k=-rr; k<=rr; k++){
		for(int j=-rr; j<=rr; j++){
		for(int i=-rr; i<=rr; i++){
			float radius = sqrt( (float)(i*i) + (float)(j*j) + (float)(k*k)); 
			float s=(float)(0.1/( 1.0 + exp(2.0*(radius-0.5*(float)R))) );
			
			if( s > real(IMAGE(i+xx,j+yy,k+zz)) ){
				IMAGE(i+xx,j+yy,k+zz) = complex<float>( s,0.0); 
				TOA(i+xx,j+yy,k+zz) = T;
			}

		}}}
	}
	
	// Export Magnitude
	if(debug==1){
		ArrayWriteMag(IMAGE,"Phantom.dat"); 
		ArrayWrite(TOA,"TOA.dat"); 
	}
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
}


//----------------------------------------
// Phantom Read Command Line
//----------------------------------------
void PHANTOM::read_commandline(int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if (strcmp("-h", pstring[pos] ) == 0) {
		float_flag("-phantom_noise",phantom_noise);
	}
  }
}    
