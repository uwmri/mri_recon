/************************************************

Library to create data out of nothing!

*************************************************/

#include "phantom.h"

using namespace arma;
using namespace std;

PHANTOM::PHANTOM( void){
	over_res=2;
	phantom_type = FRACTAL_PHANTOM;
	phantom_noise= 1.0;	
}

void PHANTOM::init(int Nx, int Ny, int Nz){
	

	IMAGE.setStorage( ColumnMajorArray<3>());
  	IMAGE.resize(Nx*over_res,Ny*over_res,Nz*over_res);
  	IMAGE = 0;
	
	SMAP.setStorage( ColumnMajorArray<3>());
  	SMAP.resize(Nx*over_res,Ny*over_res,Nz*over_res);
  	SMAP = 1;
		
	// Array for Storage of Image
	int Nt=1;
	switch(phantom_type){
		case(FRACTAL_PHANTOM):{
			fractal3D(over_res*Nx,over_res*Ny,over_res*Nz);
		}break;
		
		case(SHEPP_PHANTOM):{
			// Doesn't exist yet (but can be k-space based, Koay et al.)
		}break;
	}
}

void  PHANTOM::update_smap(int coil, int Ncoils){
	float Nx = (float)SMAP.length(firstDim);
	float Ny = (float)SMAP.length(secondDim);
	float Nz = (float)SMAP.length(thirdDim);
	float scale = 1./(Nx*Ny*Nz);
	
	switch(Ncoils){
		case(1):{
			SMAP = 1;
		}break;
		
		case(2):{
			for(int k=0; k<SMAP.length(thirdDim); k++){
			for(int j=0; j<SMAP.length(secondDim); j++){
			for(int i=0; i<SMAP.length(firstDim); i++){
				switch(coil){
					case(0):{ SMAP(i,j,k)=scale*(Nx - (float)i); }break;
					case(1):{ SMAP(i,j,k)=scale*(     (float)i); }break;
				}
			}}}
		}break;
		
		case(4):{
			for(int k=0; k<SMAP.length(thirdDim); k++){
			for(int j=0; j<SMAP.length(secondDim); j++){
			for(int i=0; i<SMAP.length(firstDim); i++){
				switch(coil){
					case(0):{ SMAP(i,j,k)=scale*(Nx - (float)i)*(Ny - (float)j); }break;
					case(1):{ SMAP(i,j,k)=scale*(     (float)i)*(Ny - (float)j); }break;
					case(2):{ SMAP(i,j,k)=scale*(Nx - (float)i)*(     (float)j); }break;
					case(3):{ SMAP(i,j,k)=scale*(     (float)i)*(     (float)j); }break;
				}
			}}}
		}break;
		
		case(8):{
			for(int k=0; k<SMAP.length(thirdDim); k++){
			for(int j=0; j<SMAP.length(secondDim); j++){
			for(int i=0; i<SMAP.length(firstDim); i++){
				switch(coil){
					case(0):{ SMAP(i,j,k)=scale*(Nx - (float)i)*(Ny - (float)j)*(Nz - (float)k); }break;
					case(1):{ SMAP(i,j,k)=scale*((float)i)*(Ny - (float)j)*(Nz - (float)k); }break;
					case(2):{ SMAP(i,j,k)=scale*(Nx - (float)i)*((float)j)*(Nz - (float)k); }break;
					case(3):{ SMAP(i,j,k)=scale*((float)i)*((float)j)*(Nz - (float)k); }break;
					case(4):{ SMAP(i,j,k)=scale*(Nx - (float)i)*(Ny - (float)j)*((float)k); }break;
					case(5):{ SMAP(i,j,k)=scale*((float)i)*(Ny - (float)j)*((float)k); }break;
					case(6):{ SMAP(i,j,k)=scale*(Nx - (float)i)*((float)j)*((float)k); }break;
					case(7):{ SMAP(i,j,k)=scale*((float)i)*((float)j)*((float)k); }break;
				}
			}}}
		}break;
		
	
	
	}
}


// Add Noise
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


void  PHANTOM::fractal3D(int Nx, int Ny, int Nz){

	int n = 10; // Maximum Branches
	float branch_angle= PI/5;
	float branch_rotation = PI/4; //(3-sqrt(5));
	float branch_fraction = 0.8;
	float branch_length = 60;
	
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
	for(int t=0; t< trees; t++){
		for(int ii=0; ii<tree[t].count; ii++){
			SX = min( tree[t].xlist[ii],SX);
			SY = min( tree[t].ylist[ii],SY);
			SZ = min( tree[t].zlist[ii],SZ);
			EX = max( tree[t].xlist[ii],EX);
			EY = max( tree[t].ylist[ii],EY);
			EZ = max( tree[t].zlist[ii],EZ);
		}
	}
	cout << "Range is" << endl;
	cout << "\t X:" << SX << "to " << EX << endl;
	cout << "\t Y:" << SY << "to " << EY << endl;
	cout << "\t Z:" << SZ << "to " << EZ << endl;
	
	float res = (float)Nx;
	
	float scaleX = res/(60 + EX - SX);
	float scaleY = res/(60 + EY - SY);
	float scaleZ = res/(60 + EZ - SZ);
	float scale = min(min(scaleX,scaleY),scaleZ);

	float CX = scale*(EX + SX)/2 - res/2;
	float CY = scale*(EY + SY)/2 - res/2;
	float CZ = scale*(EZ + SZ)/2 - res/2;
	
	for(int t=0; t< trees; t++){
    	// Vessels
		for(int ii=0; ii<tree[t].count; ii++){
        	float X=scale*tree[t].xlist[ii];
    		float Y=scale*tree[t].ylist[ii];
    		float Z=scale*tree[t].zlist[ii];
    		float R=scale*tree[t].rlist[ii];
			
			int xx = (int)floor(X - CX);
        	int yy = (int)floor(Y - CY);
        	int zz = (int)floor(Z - CZ);
        	int rr = (int)ceil(R);
        	
			for(int k=-rr; k<=rr; k++){
			for(int j=-rr; j<=rr; j++){
			for(int i=-rr; i<=rr; i++){
				float radius = sqrt( (float)(i*i) + (float)(j*j) + (float)(k*k)); 
				IMAGE(i+xx,j+yy,k+zz)= complex<float>( max(real(IMAGE(i+xx,j+yy,k+zz)), (float)(1.0/( 1.0 + exp(2.0*(radius-0.5*(float)R))) )),0.0); 
			}}}
		}
		
		// Tissue
		int ii=tree[t].count-1;
        float X=scale*tree[t].xlist[ii];
    	float Y=scale*tree[t].ylist[ii];
    	float Z=scale*tree[t].zlist[ii];
    	float R=scale*30.0;
			
		int xx = (int)floor(X - CX);
        int yy = (int)floor(Y - CY);
        int zz = (int)floor(Z - CZ);
        int rr = (int)ceil(R);
        	
		for(int k=-rr; k<=rr; k++){
		for(int j=-rr; j<=rr; j++){
		for(int i=-rr; i<=rr; i++){
			float radius = sqrt( (float)(i*i) + (float)(j*j) + (float)(k*k)); 
			IMAGE(i+xx,j+yy,k+zz)= complex<float>( max(real(IMAGE(i+xx,j+yy,k+zz)), (float)(0.1/( 1.0 + exp(2.0*(radius-0.5*(float)R))) )),0.0); 
		}}}
	}
	ArrayWriteMag(IMAGE,"Phantom.dat"); 
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
// Phantome Read Command Line
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
