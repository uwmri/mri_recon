#include "voronoi_dcf.h"

using namespace std;
using namespace NDarray;
using namespace voro;
using namespace arma;

void VORONOI_DCF::vor_dcf( Array<float,3> &kw,Array<float,3> &kx,Array<float,3> &ky,Array<float,3> &kz,VORONOI_DCF::KShape shape ){

	double x_min = min(kx);
	double y_min = min(ky);
	double z_min = min(kz);
	
	double x_max = max(kx);
	double y_max = max(ky);
	double z_max = max(kz);
	
	int n_x = ceil((x_max - x_min)/4);
	int n_y = ceil((y_max - y_min)/4);
	int n_z = max( 1, (int)ceil((z_max - z_min)/4) );
	
	cout << "Starting Vor DCF " << endl;
	cout << "  X: " << x_min << " to " << x_max << endl;
	cout << "  Y: " << y_min << " to " << y_max << endl;
	cout << "  Z: " << z_min << " to " << z_max << endl;
	cout << "  Grid Size: " << n_x << " , " << n_y << " , " << n_z  << endl;
	
	// Create Container
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,512);

	// Add wall
	wall_sphere ksphere(0,0,0,x_max);
	wall_cylinder kcylinder(0,0,-0.5+z_min,0,0,0.5+z_max,x_max);
	switch(shape){
		case(CYLINDER):{
			// Add a cylinder to the container
    		cout << "Add Cylinder ( " << x_max << " ) " << endl;
			con.add_wall(kcylinder);
		}break;
	
		case(SPHERE):{
			// Add a spherical wall to the container
   			cout << "Add Sphere" << endl;
			con.add_wall(ksphere);
		}break;
	}
	


	float fermi_r = x_max - 6;
	float fermi_w = 1;
	
	int Npts = kx.numElements();
	
	// Add Points
	particle_order order(Npts);
		
	// Need to remove points that are duplicates to some precision
	cout << "Assign to Arma pts = " << Npts << endl;
	arma::fvec kxx (Npts);
	arma::fvec kyy (Npts);
	arma::fvec kzz (Npts);
	int count= 0;
	for(int k =0; k< kx.length(thirdDim); k++){
		for(int j =0; j< kx.length(secondDim); j++){
			for(int i =0; i< kx.length(firstDim); i++){
				kxx(count) = kx(i,j,k);
				kyy(count) = ky(i,j,k);
				kzz(count) = kz(i,j,k);
				count++;				
	}}}
				
	cout << "Sorting in x" << endl;
	uvec index = sort_index( kxx ); // Sort in xx
	
	cout << "Asigning " << endl;	
	arma::fvec kxx_sorted (Npts);
	arma::fvec kyy_sorted (Npts);
	arma::fvec kzz_sorted (Npts);
	for(int pos=0; pos< (int)kxx.n_elem; pos++){
		kxx_sorted( pos )=kxx(index(pos));
		kyy_sorted( pos )=kyy(index(pos));
		kzz_sorted( pos )=kzz(index(pos));
	}
	
	
	arma::fvec kn( Npts); // Number of points
	arma::vec kpos( Npts); // Position in container
	arma::vec counted( Npts); // Sets if pint belongs
	arma::uvec set_idx( Npts); // Sets if pint belongs
	counted.zeros();
	
	cout << "Finding Unique and adding" << endl;
	int unique_points=0;
	for(unsigned int pos=0; pos< kxx.n_elem; pos++ ){
			
			if( (pos % 10000) ==0){
				cout << "Pos = " << pos << " found " << unique_points << endl;
			}
			
			// Check to see if point has been collected
			if(counted(pos) ){
				continue;
			}
						
			// Now search over points
			float kxt = kxx_sorted(pos); 
			float kyt = kyy_sorted(pos); 
			float kzt = kzz_sorted(pos); 
			
			// Gather all non-unique points and average
			float kxavg = kxt;
			float kyavg = kyt;
			float kzavg = kzt;
			int n = 1;
			set_idx(n) = pos;
			
			// Gather unique points
			int forward_pos=pos+1;
			while( forward_pos < (int)kxx.n_elem){
				// Difference in X (determines search range)
				float xdiff = abs(kxx_sorted(forward_pos) - kxt);
				
				// Difference in radius				
				float diff = pow(xdiff,2.0);
				diff      += pow(kyy_sorted(forward_pos) - kyt,2.0);
				diff      += pow(kzz_sorted(forward_pos) - kzt,2.0);
				diff = sqrt(diff);
				
				
				if(xdiff > 0.05){
					// out of search range
					break;
				}else if(diff > 0.5){
					// Just don't add point
				}else{
					// Count this shot
					kxavg += kxx_sorted(forward_pos);
					kyavg += kyy_sorted(forward_pos);
					kzavg += kzz_sorted(forward_pos);
					set_idx(n) = forward_pos;
					n++;
				}
				forward_pos++;
			}
									
			// Compute actual average
			kxavg /= (float)n;
			kyavg /= (float)n;
			kzavg /= (float)n;
			
			// cout << "Pt = (" << kxavg << "," << kyavg << "," << kzavg << ")" << endl;
			
			// Check to see if it's in the container. If it is add it
			if(con.point_inside(kxavg,kyavg,kzavg)){
				con.put(order,unique_points,kxavg,kyavg,kzavg);
				
				// Update decoding array
				for(int dpos= 0; dpos < n; dpos++){
					kpos(set_idx(dpos)) = unique_points;
					kn(set_idx(dpos))   = n;
					counted(set_idx(dpos)) = 1;
				}
				unique_points++;
			}else{
				// Update decoding array
				// Update decoding array
				for(int dpos= 0; dpos < n; dpos++){
					kpos(set_idx(dpos)) = -1;
					kn(set_idx(dpos))   = n;
					counted(set_idx(dpos)) = 1;
				}
			}
	}
	cout << "Done combining, found " << unique_points << " points " << endl << flush;
	

	// Now calculate and copy back
	arma::fvec kw_calculated( unique_points );
	voronoicell c;
  	c_loop_order vl(con,order);
  	
	count=0;
	if(vl.start()) do{
		if(count%10000==0){
			cout << "Counted " << count/10000 << " of " << unique_points/10000 << endl;
		}
				
		if(con.compute_cell(c,vl)){
			kw_calculated(count) = c.volume();
        }else{
			kw_calculated(count) = 0.0;
		}
		count++;
	}while( vl.inc());

	// Convert to sorted index
	arma::fvec kw_sorted( Npts );
	for(int pos=0; pos < Npts; pos++){
		if( kpos(pos) == -1){
			kw_sorted(pos) = 0.0;
		}else{
			kw_sorted(pos) = kw_calculated( kpos(pos))/kn(pos);
		}
	}
	
	// Unsort
	cout << " Unsort " << endl;
	arma::fvec kw_unsorted( Npts );
	for(int pos=0; pos< (int)kxx.n_elem; pos++){
		kw_unsorted(index(pos) )=kw_sorted(pos);
	}
		
	// Unsort and put in array 
	cout << "Copy back" << endl;
	count = 0;
	for(int k =0; k< kx.length(thirdDim); k++){
		for(int j =0; j< kx.length(secondDim); j++){
			for(int i =0; i< kx.length(firstDim); i++){
				float kww = kw_unsorted(count);
				
				float r = sqrt( kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) + kz(i,j,k)*kz(i,j,k));
				kw(i,j,k) = kww*(1 + exp( (r-fermi_r)/fermi_w));
				count++;				
	}}}
}


