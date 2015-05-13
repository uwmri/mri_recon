#include "voronoi_dcf.h"

using namespace std;
using namespace NDarray;
using namespace voro;

void VORONOI_DCF::dcf_3D( Array<float,3> &kw,Array<float,3> &kx,Array<float,3> &ky,Array<float,3> &kz){

	double x_min = min(kx);
	double y_min = min(ky);
	double z_min = min(kz);
	
	double x_max = max(kx);
	double y_max = max(ky);
	double z_max = max(kz);
	
	int n_x = ceil((x_max - x_min)/4);
	int n_y = ceil((y_max - y_min)/4);
	int n_z = ceil((z_max - z_min)/4);
	
	cout << "Starting Vor DCF " << endl;
	cout << "  X: " << x_min << " to " << x_max << endl;
	cout << "  Y: " << y_min << " to " << y_max << endl;
	cout << "  Z: " << z_min << " to " << z_max << endl;
	cout << "  Grid Size: " << n_x << " , " << n_y << " , " << n_z  << endl;
	
	// Create Container
	container con(x_min,x_max, \
				  y_min,y_max, \
				  z_min,z_max, \
				  n_x,n_y,n_z, \
				  false, \
				  false, \
				  false, \
				  50);
	
	// Add a spherical wall to the container
    wall_sphere kmax(0,0,0,x_max);
    con.add_wall(kmax);
	
	float fermi_r = x_max - 6;
	float fermi_w = 1;
	
	
	// Add Points
	particle_order order(kx.numElements() );
	Array< int,1> index_order( kx.numElements() );
	index_order = 0;
	
	int count = 0;
	int tot_count = 0;
	for(int k =0; k< kx.length(thirdDim); k++){
		for(int j =0; j< kx.length(secondDim); j++){
			for(int i =0; i< kx.length(firstDim); i++){
				
				double x = (double)kx(i,j,k) + 0.00001*(double)( (double)rand()/(double)RAND_MAX - 0.5);
				double y = (double)ky(i,j,k) + 0.00001*(double)( (double)rand()/(double)RAND_MAX - 0.5);
				double z = (double)kz(i,j,k) + 0.00001*(double)( (double)rand()/(double)RAND_MAX - 0.5);
						
				if(con.point_inside(x,y,z)){ 
        			con.put(order,count,x,y,z);
					index_order(count) = tot_count;
					count++;
				}
				tot_count++;
			}
		}	
	}

	
	// Now calculate and copy back
	voronoicell c;
  	c_loop_order vl(con,order);
  	
	count=0;
	if(vl.start()) do{
		
		double vv = 0.0;
		
		if(count%10000==0){
			cout << "Counted " << count/10000 << " of " << kx.numElements()/10000 << endl;
		}
				
		
		if(con.compute_cell(c,vl)){
			vv = c.volume();
        }else{
			vv = 0.0;
		}
		
		int t = index_order( count);
		
		// Convert count to index
		int i = t % kx.length(firstDim);
		int tempjj = t / kx.length(firstDim);
		int j = tempjj % kx.length(secondDim);
		int k = tempjj / kx.length(secondDim);
		
		float r = sqrt( kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) + kz(i,j,k)*kz(i,j,k));
		
		kw(i,j,k) = vv*(1 + exp( (r-fermi_r)/fermi_w));
		
		count++;
		
	}while( vl.inc());

}

void VORONOI_DCF::dcf_2D( Array<float,3> &kw,Array<float,3> &kx,Array<float,3> &ky){
	
	double x_min = min(kx);
	double y_min = min(ky);
	double z_min = -0.5;
	
	double x_max = max(kx);
	double y_max = max(ky);
	double z_max = 0.5;
	
	int n_x = ceil((x_max - x_min)/4);
	int n_y = ceil((y_max - y_min)/4);
	int n_z = 1;
	
	cout << "Starting Vor DCF " << endl;
	cout << "  X: " << x_min << " to " << x_max << endl;
	cout << "  Y: " << y_min << " to " << y_max << endl;
	cout << "  Z: " << z_min << " to " << z_max << endl;
	cout << "  Grid Size: " << n_x << " , " << n_y << " , " << n_z  << endl;
	
	// Create Container
	container con(x_min,x_max, \
				  y_min,y_max, \
				  z_min,z_max, \
				  n_x,n_y,n_z, \
				  false, \
				  false, \
				  false, \
				  64);
	
	// Add a spherical wall to the container
    wall_cylinder kmax(0,0,-0.5,0,0,0.5,x_max);
    con.add_wall(kmax);
	
	float fermi_r = x_max - 6;
	float fermi_w = 1;
	
	
	// Add Points
	particle_order order(kx.numElements() );
	Array< int,1> index_order( kx.numElements() );
	index_order = 0;
	
	int count = 0;
	int tot_count = 0;
	for(int k =0; k< kx.length(thirdDim); k++){
		for(int j =0; j< kx.length(secondDim); j++){
			for(int i =0; i< kx.length(firstDim); i++){
				
				double x = (double)kx(i,j,k);
				double y = (double)ky(i,j,k);
				double z = 0.0;
							
				if(con.point_inside(x,y,z)){ 
        			con.put(order,count,x,y,z);
					index_order(count) = tot_count;
					count++;
				}
				tot_count++;
			}
		}	
	}

	
	// Now calculate and copy back
	voronoicell c;
  	c_loop_order vl(con,order);
  	vl.start();
	count=0;


	if(vl.start()) do{
		
		double vv = 0.0;
		
		if(count%10000==0){
			cout << "Counted " << count/10000 << " of " << kx.numElements()/10000 << endl;
		}
				
		
		if(con.compute_cell(c,vl)){
			vv = c.volume();
        }else{
			vv = 0.0;
		}
		
		int t = index_order( count);
		
		// Convert count to index
		int i = t % kx.length(firstDim);
		int tempjj = t / kx.length(firstDim);
		int j = tempjj % kx.length(secondDim);
		int k = tempjj / kx.length(secondDim);
		
		float r = sqrt( kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) );
		
		kw(i,j,k) = vv/(1 + exp( (r-fermi_r)/fermi_w) );
		
		count++;
		
	}while( vl.inc());
	

}

