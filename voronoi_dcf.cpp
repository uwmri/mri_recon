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
	
	int n_x = 1+ceil(1.0*(x_max - x_min));
	int n_y = 1+ceil(1.0*(y_max - y_min));
	int n_z = 1+ceil(1.0*(z_max - z_min));
	
	cout << "Starting Vor DCF " << endl;
	cout << "  X: " << x_min << " to " << x_max << endl;
	cout << "  Y: " << y_min << " to " << y_max << endl;
	cout << "  Z: " << z_min << " to " << z_max << endl;
	cout << "  Grid Size: " << n_x << " , " << n_y << " , " << n_z  << endl;
	
	// Create Container
	container con(x_min-0.5,x_max+0.5,y_min-0.5,y_max+0.5,z_min-0.5,z_max+0.5,n_x,n_y,n_z,false,false,false,64);
	
	double kmax = x_max*0.98;
	
	// Add wall
	wall_sphere ksphere(0,0,0,kmax+0.5);
	wall_cylinder kcylinder(0,0,0,0,0,1.0,kmax+0.5);
	
	switch(shape){
		case(CYLINDER):{
			// Add a cylinder to the container
    		cout << "Add Cylinder ( " << kmax+0.5 << " ) " << endl;
			con.add_wall(kcylinder);
		}break;
	
		case(SPHERE):{
			// Add a spherical wall to the container
   			cout << "Add Sphere" << endl;
			con.add_wall(ksphere);
		}break;
	}
	


	float fermi_r = kmax - 6;
	float fermi_w = 1;
	
	int Npts = kx.numElements();
	
	
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
	int non_unique_points=0;
	particle_order order(Npts);
	
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
			set_idx(0) = pos;
			
			// Gather unique points
			int forward_pos=pos+1;
			while( forward_pos < (int)kxx.n_elem){
				// Difference in X (determines search range)
				float xdiff = abs(kxx_sorted(forward_pos) - kxt);
				
				// Difference in radius				
				float ydiff = abs(kyy_sorted(forward_pos) - kyt);
				float zdiff = abs(kzz_sorted(forward_pos) - kzt);
								
				if(xdiff > 0.01){
					// out of search range
					break;
				}else if( (ydiff>0.01) || ( zdiff > 0.01)){
					// Just don't add point
				}else{
					// Add this
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
			if(con.point_inside(kxavg,kyavg,kzavg) ){
				
				// Actually add points
				con.put(order,unique_points,kxavg,kyavg,kzavg);
				
				// Update decoding array
				for(int dpos= 0; dpos < n; dpos++){
					kpos(set_idx(dpos)) = unique_points;
					kn(set_idx(dpos))   = (float)n;
					counted(set_idx(dpos)) = 1;
				}
				unique_points++;
			}else{
				// Update decoding array
				// Update decoding array
				for(int dpos= 0; dpos < n; dpos++){
					kpos(set_idx(dpos)) = -1;
					kn(set_idx(dpos))   = (float)n;
					counted(set_idx(dpos)) = 1;
				}
				non_unique_points += n;
			}
	}
	
	cout << "Min counted" <<  min(counted) << endl;
	cout << "Done combining, found " << unique_points << " unique points " << endl << flush;
	
	// Add a ring of points at edge
	switch(shape){
	
	case(CYLINDER):{
		int radial_pts = (int)( 20*2*3.14156*kmax);
		cout << "Adding " << radial_pts << " radial points to kmax" << endl;
		int extra_points=0;
		for(int pos =0; pos<radial_pts; pos++){
			double z = 0.0;
			double y = ( kmax)*sin( 2*3.14156*(double)pos/ (double)radial_pts );
			double x = ( kmax)*cos( 2*3.14156*(double)pos/ (double)radial_pts );
			if(con.point_inside(x,y,z) ){
				con.put(order,(unique_points+extra_points),x,y,z);
				extra_points++;
			}
		
		}
		cout << "Actually added " << extra_points << endl;
	}break;
	
	case(SPHERE):{
		int N = (int)( 4*3.14156*kmax*kmax);
		cout << "Adding " << N << " radial points to kmax" << endl;
		
		int extra_points=0;
		for(int pos =0; pos< N; pos++){
			double z =  (2*pos - (double)N + 1)/(double)N;
			double y = sin( sqrt((double)N*3.15156)*asin(z) )*sqrt( 1 - z*z);
			double x = cos( sqrt((double)N*3.15156)*asin(z) )*sqrt( 1 - z*z);
			z *= (0.25+kmax);
			y *= (0.25+kmax);
			x *= (0.25+kmax);
			
			
			if(con.point_inside(x,y,z) ){
				con.put(order,(unique_points+extra_points),x,y,z);
				extra_points++;
			}
		
		}
		cout << "Actually added " << extra_points << endl;
	}break;
	
	}
	
	
	//con.draw_cells_pov("cells.pov");
	//con.draw_particles_pov("particles.pov");
	

	// Now calculate and copy back
	arma::fvec kw_calculated( unique_points );
	
	voronoicell cell;
  	c_loop_order vl(con,order);
  	
	cout << "vl start = " << vl.start() << endl;
	count=0;
	if( vl.start() ) do{
		if(count%1000==0){
			
			cout << "\rProgress " << ( 100*count/unique_points ) <<  "%" << flush;
		}
		
		if(con.compute_cell(cell,vl)){
			kw_calculated(count) = cell.volume();
        }else{
			kw_calculated(count) = 0.0;
		}
		// cout << kw_calculated(count) << endl;
		count++;
	}while( vl.inc() && (count< unique_points) );
	cout << "Counted " << count << " points in dcf calc" << endl;	


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
				kw(i,j,k) = kww; //*(1 + exp( (r-fermi_r)/fermi_w));
				count++;				
	}}}
}


