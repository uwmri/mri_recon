#include "phantom.h"

using namespace arma;
using namespace std;

// Constructor:
//	
//
PHANTOM::PHANTOM( void){
	over_res=2;
	phantom_type = FRACTAL_PHANTOM;
	phantom_noise= 1.0;	
}

// Initialization:
//	
//
void PHANTOM::init(int Nx, int Ny, int Nz){
	
	// Phantom Density
	IMAGE.setStorage( ColumnMajorArray<3>());
  	IMAGE.resize((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
  	IMAGE = 0;
	
	//Sensitivity Map	
	SMAP.setStorage( ColumnMajorArray<3>());
  	SMAP.resize((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
  	SMAP = 1;
		
	// Switch For Phantom
	switch(phantom_type){
		case(FRACTAL_PHANTOM):{
			TOA.setStorage( ColumnMajorArray<3>());
  			TOA.resize((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
  			TOA = 0;
			
			fractal3D((int)((float)Nx*over_res),(int)((float)Ny*over_res),(int)((float)Nz*over_res));
		}break;
		
		case(SHEPP_PHANTOM):{
			// Doesn't exist yet (but can be k-space based, Koay et al.)
		}break;
	}
}


// Get a made up sensitivity map
//	 Need to update to Biot-Savart of loop coil
//
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
//	 Add some phase to the phantom
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
	
	// Fermi-Lookup
	float edge_width = 1.0;
	
	
	cout << "Gridding Tree" << endl;
	int count =0;
//	#pragma omp parallel for schedule(static,1) 
	for(int t=0; t< trees; t++){
    	cout << count << " of " << trees << endl;
		
//		#pragma omp atomic
		count++;
		
		// Vessels
		for(int ii=0; ii<tree[t].count; ii++){
        	float X=scale*tree[t].xlist[ii];
    		float Y=scale*tree[t].ylist[ii];
    		float Z=scale*tree[t].zlist[ii];
    		float R=scale*tree[t].rlist[ii];
			float T=tree[t].btime[ii];
			if(R<1.0){
				continue;
			}
			int xx = (int)round(X - CX);
        	int yy = (int)round(Y - CY);
        	int zz = (int)round(Z - CZ);
			int rr = (int)ceil(R);
        							
			for(int k=-rr; k<=rr; k++){
			for(int j=-rr; j<=rr; j++){
			for(int i=-rr; i<=rr; i++){
				float radius = sqrt( (float)(i*i) + (float)(j*j) + (float)(k*k)); 
				//if(radius < 0.5*R){
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
			float s=(float)(1.0/( 1.0 + exp(2.0*(radius-0.5*(float)R))) );
			if( s > real(IMAGE(i+xx,j+yy,k+zz))){
				IMAGE(i+xx,j+yy,k+zz) = complex<float>( s,0.0); 
				TOA(i+xx,j+yy,k+zz) = T;
			}
		}}}
	}
	
	// Export Magnitude
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
