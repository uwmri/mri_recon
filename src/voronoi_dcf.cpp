#include "voronoi_dcf.h"

using namespace std;
using namespace NDarray;
using namespace voro;
using namespace arma;

/** \brief A class representing a spherical wall object.
 *
 * This class represents a spherical wall object. */
struct wall_sphere2 : public wall {
 public:
  /** Constructs a spherical wall object.
   * \param[in] w_id_ an ID number to associate with the wall for
   *		    neighbor tracking.
   * \param[in] (xc_,yc_,zc_) a position vector for the sphere's
   * 			    center.
   * \param[in] rc_ the radius of the sphere. */
  wall_sphere2(double xc_, double yc_, double zc_, double rc_, int w_id_ = -99)
      : w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), rc(rc_) {}
  bool point_inside(double x, double y, double z);
  template <class v_cell>
  bool cut_cell_base(v_cell &c, double x, double y, double z);
  bool cut_cell(voronoicell &c, double x, double y, double z) {
    return cut_cell_base(c, x, y, z);
  }
  bool cut_cell(voronoicell_neighbor &c, double x, double y, double z) {
    return cut_cell_base(c, x, y, z);
  }

 private:
  const int w_id;
  const double xc, yc, zc, rc;
};

/** Tests to see whether a point is inside the sphere wall object.
 * \param[in,out] (x,y,z) the vector to test.
 * \return True if the point is inside, false if the point is outside. */
bool wall_sphere2::point_inside(double x, double y, double z) {
  return (x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc) >
         rc * rc;
}

/** Cuts a cell by the sphere wall object. The spherical wall is approximated by
 * a single plane applied at the point on the sphere which is closest to the
 * center of the cell. This works well for particle arrangements that are packed
 * against the wall, but loses accuracy for sparse particle distributions.
 * \param[in,out] c the Voronoi cell to be cut.
 * \param[in] (x,y,z) the location of the Voronoi cell.
 * \return True if the cell still exists, false if the cell is deleted. */
template <class v_cell>
bool wall_sphere2::cut_cell_base(v_cell &c, double x, double y, double z) {
  double xd = x - xc;
  double yd = y - yc;
  double zd = z - zc;
  double dq = xd * xd + yd * yd + zd * zd;  // Distance ^2
  // double rq=sqrt(dq); // Distance

  if (dq > 1e-5) {
    dq = 2 * (sqrt(dq) * rc - dq);
    // dq=2*(sqrt(dq)*rc-rc*rc);
    return c.nplane(xd, yd, zd, dq, w_id);
  }
  return true;
}

void VORONOI_DCF::vor_dcf(Array<float, 3> &kw, Array<float, 3> &kx,
                          Array<float, 3> &ky, Array<float, 3> &kz,
                          VORONOI_DCF::KShape shape) {
  double x_min = min(kx);
  double y_min = min(ky);
  double z_min = min(kz);

  double x_max = max(kx);
  double y_max = max(ky);
  double z_max = max(kz);

  int n_x = 1 + ceil(1.0 * (x_max - x_min));
  int n_y = 1 + ceil(1.0 * (y_max - y_min));
  int n_z = 1 + ceil(1.0 * (z_max - z_min));

  cout << "Starting Vor DCF " << endl;
  cout << "  X: " << x_min << " to " << x_max << endl;
  cout << "  Y: " << y_min << " to " << y_max << endl;
  cout << "  Z: " << z_min << " to " << z_max << endl;
  cout << "  Grid Size: " << n_x << " , " << n_y << " , " << n_z << endl;

  // Create Container
  container con(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5, z_min - 0.5,
                z_max + 0.5, n_x, n_y, n_z, false, false, false, 64);

  double kmax = x_max * 0.98;

  // Add wall
  wall_sphere ksphere(0, 0, 0, kmax + 0.5);
  wall_cylinder kcylinder(0, 0, 0, 0, 0, 1.0, kmax + 0.5);

  switch (shape) {
    case (CYLINDER): {
      // Add a cylinder to the container
      cout << "Add Cylinder ( " << kmax + 0.5 << " ) " << endl;
      con.add_wall(kcylinder);
    } break;

    case (SPHERE): {
      // Add a spherical wall to the container
      cout << "Add Sphere" << endl;
      con.add_wall(ksphere);
    } break;

    case (CUBE): {
      // Cube
    } break;
  }

  cout << "Sizes of Input Array" << endl;
  cout << "Kx = " << kx.length(firstDim) << " x " << kx.length(secondDim) << "x"
       << kx.length(thirdDim) << endl;
  cout << "Ky = " << ky.length(firstDim) << " x " << ky.length(secondDim) << "x"
       << ky.length(thirdDim) << endl;
  cout << "Kz = " << kz.length(firstDim) << " x " << kz.length(secondDim) << "x"
       << kz.length(thirdDim) << endl;
  cout << "Kw = " << kw.length(firstDim) << " x " << kw.length(secondDim) << "x"
       << kw.length(thirdDim) << endl;

  // float fermi_r = kmax - 6;
  // float fermi_w = 1;

  int Npts = kx.numElements();

  // Need to remove points that are duplicates to some precision
  cout << "Assign to Arma pts = " << Npts << endl;
  arma::fvec kxx(Npts);
  arma::fvec kyy(Npts);
  arma::fvec kzz(Npts);
  int count = 0;
  for (int k = 0; k < kx.length(thirdDim); k++) {
    for (int j = 0; j < kx.length(secondDim); j++) {
      for (int i = 0; i < kx.length(firstDim); i++) {
        kxx(count) = kx(i, j, k);
        kyy(count) = ky(i, j, k);
        kzz(count) = kz(i, j, k);
        count++;
      }
    }
  }

  cout << "Sorting in x" << endl;
  uvec index = sort_index(kxx);  // Sort in xx

  cout << "Asigning " << endl;
  arma::fvec kxx_sorted(Npts);
  arma::fvec kyy_sorted(Npts);
  arma::fvec kzz_sorted(Npts);
  for (int pos = 0; pos < (int)kxx.n_elem; pos++) {
    kxx_sorted(pos) = kxx(index(pos));
    kyy_sorted(pos) = kyy(index(pos));
    kzz_sorted(pos) = kzz(index(pos));
  }

  arma::fvec kn(Npts);       // Number of points
  arma::vec kpos(Npts);      // Position in container
  arma::vec counted(Npts);   // Sets if pint belongs
  arma::uvec set_idx(Npts);  // Sets if pint belongs
  counted.zeros();

  cout << "Finding Unique and adding" << endl;
  int unique_points = 0;
  int non_unique_points = 0;
  particle_order order(Npts);

  for (unsigned int pos = 0; pos < kxx.n_elem; pos++) {
    if ((pos % 10000) == 0) {
      cout << "Pos = " << pos << " found " << unique_points << endl;
    }

    // Check to see if point has been collected
    if (counted(pos)) {
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
    int forward_pos = pos + 1;
    while (forward_pos < (int)kxx.n_elem) {
      // Difference in X (determines search range)
      float xdiff = abs(kxx_sorted(forward_pos) - kxt);

      // Difference in radius
      float ydiff = abs(kyy_sorted(forward_pos) - kyt);
      float zdiff = abs(kzz_sorted(forward_pos) - kzt);

      if (xdiff > 0.01) {
        // out of search range
        break;
      } else if ((ydiff > 0.01) || (zdiff > 0.01)) {
        // Just don't add point
      } else {
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
    if (con.point_inside(kxavg, kyavg, kzavg)) {
      // Actually add points
      con.put(order, unique_points, kxavg, kyavg, kzavg);

      // Update decoding array
      for (int dpos = 0; dpos < n; dpos++) {
        kpos(set_idx(dpos)) = unique_points;
        kn(set_idx(dpos)) = (float)n;
        counted(set_idx(dpos)) = 1;
      }
      unique_points++;
    } else {
      // Update decoding array
      // Update decoding array
      for (int dpos = 0; dpos < n; dpos++) {
        kpos(set_idx(dpos)) = -1;
        kn(set_idx(dpos)) = (float)n;
        counted(set_idx(dpos)) = 1;
      }
      non_unique_points += n;
    }
  }

  cout << "Min counted" << min(counted) << endl;
  cout << "Done combining, found " << unique_points << " unique points " << endl
       << flush;

  // Add a ring of points at edge
  switch (shape) {
    case (CYLINDER): {
      int radial_pts = (int)(20 * 2 * 3.14156 * kmax);
      cout << "Adding " << radial_pts << " radial points to kmax" << endl;
      int extra_points = 0;
      for (int pos = 0; pos < radial_pts; pos++) {
        double z = 0.0;
        double y = (kmax)*sin(2 * 3.14156 * (double)pos / (double)radial_pts);
        double x = (kmax)*cos(2 * 3.14156 * (double)pos / (double)radial_pts);
        if (con.point_inside(x, y, z)) {
          con.put(order, (unique_points + extra_points), x, y, z);
          extra_points++;
        }
      }
      cout << "Actually added " << extra_points << endl;
    } break;

    case (SPHERE): {
      int N = (int)(4 * 3.14156 * kmax * kmax);
      cout << "Adding " << N << " radial points to kmax" << endl;

      int extra_points = 0;
      for (int pos = 0; pos < N; pos++) {
        double z = (2 * pos - (double)N + 1) / (double)N;
        double y = sin(sqrt((double)N * 3.15156) * asin(z)) * sqrt(1 - z * z);
        double x = cos(sqrt((double)N * 3.15156) * asin(z)) * sqrt(1 - z * z);
        z *= (0.25 + kmax);
        y *= (0.25 + kmax);
        x *= (0.25 + kmax);

        if (con.point_inside(x, y, z)) {
          con.put(order, (unique_points + extra_points), x, y, z);
          extra_points++;
        }
      }
      cout << "Actually added " << extra_points << endl;
    } break;

    case (CUBE): {
      // Don't add extra points
    } break;
  }

  // con.draw_cells_pov("cells.pov");
  // con.draw_particles_pov("particles.pov");

  // Now calculate and copy back
  arma::fvec kw_calculated(unique_points);

  voronoicell cell;
  c_loop_order vl(con, order);

  cout << "vl start = " << vl.start() << endl;
  count = 0;
  if (vl.start()) do {
      if (count % 1000 == 0) {
        cout << "\rProgress " << (100 * count / unique_points) << "%" << flush;
      }

      if (con.compute_cell(cell, vl)) {
        kw_calculated(count) = cell.volume();
      } else {
        kw_calculated(count) = 0.0;
      }
      // cout << kw_calculated(count) << endl;
      count++;
    } while (vl.inc() && (count < unique_points));
  cout << endl << "Counted " << count << " points in dcf calc" << endl;

  // Convert to sorted index
  cout << "Convert to Sorted" << endl;
  arma::fvec kw_sorted(Npts);
  for (int pos = 0; pos < Npts; pos++) {
    if (kpos(pos) == -1) {
      kw_sorted(pos) = 0.0;
    } else {
      kw_sorted(pos) = kw_calculated(kpos(pos)) / kn(pos);
    }
  }

  // Unsort
  cout << " Unsort " << endl;
  arma::fvec kw_unsorted(Npts);
  for (int pos = 0; pos < (int)Npts; pos++) {
    kw_unsorted(index(pos)) = kw_sorted(pos);
  }

  // Unsort and put in array
  cout << "Copy back" << endl;
  count = 0;
  for (int k = 0; k < kx.length(thirdDim); k++) {
    for (int j = 0; j < kx.length(secondDim); j++) {
      for (int i = 0; i < kx.length(firstDim); i++) {
        float kww = kw_unsorted(count);

        // float r = sqrt( kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) +
        // kz(i,j,k)*kz(i,j,k));
        kw(i, j, k) = kww;  //*(1 + exp( (r-fermi_r)/fermi_w));
        count++;
      }
    }
  }

  cout << "Done with Copy" << endl;
}

void VORONOI_DCF::vor_sphere(Array<float, 1> &kw, Array<float, 1> &kx,
                             Array<float, 1> &ky, Array<float, 1> &kz) {
  cout << "Starting Spherical Vor DCF " << endl;

  // Use a fixed radius
  double kmax = 16;

  double x_min = -kmax;
  double y_min = -kmax;
  double z_min = -kmax;
  double x_max = kmax;
  double y_max = kmax;
  double z_max = kmax;

  int n_x = 1 + ceil(1.0 * (x_max - x_min));
  int n_y = 1 + ceil(1.0 * (y_max - y_min));
  int n_z = 1 + ceil(1.0 * (z_max - z_min));

  // Create Container
  container con(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5, z_min - 0.5,
                z_max + 0.5, n_x, n_y, n_z, false, false, false, 64);

  // Add a spherical wall to the container
  wall_sphere ksphereOuter(0, 0, 0, kmax + 0.5);
  wall_sphere2 ksphereInner(0, 0, 0, kmax - 0.5);

  cout << "Add Sphere" << endl << flush;
  con.add_wall(ksphereOuter);
  con.add_wall(ksphereInner);

  int Npts = kx.numElements();

  // Need to remove points that are duplicates to some precision
  cout << "Assign to Arma pts = " << Npts << endl;
  arma::fvec kxx(Npts);
  arma::fvec kyy(Npts);
  arma::fvec kzz(Npts);
  arma::fvec krr(Npts);
  int count = 0;
  for (int i = 0; i < (int)kx.numElements(); i++) {
    double radius = sqrt(kx(i) * kx(i) + ky(i) * ky(i) + kz(i) * kz(i));

    krr(count) = radius;
    kxx(count) = kmax * kx(i) / radius;
    kyy(count) = kmax * ky(i) / radius;
    kzz(count) = kmax * kz(i) / radius;
    count++;
  }

  cout << "Sorting in x" << endl;
  uvec index = sort_index(kxx);  // Sort in xx

  cout << "Asigning " << endl;
  arma::fvec kxx_sorted(Npts);
  arma::fvec kyy_sorted(Npts);
  arma::fvec kzz_sorted(Npts);
  for (int pos = 0; pos < (int)kxx.n_elem; pos++) {
    kxx_sorted(pos) = kxx(index(pos));
    kyy_sorted(pos) = kyy(index(pos));
    kzz_sorted(pos) = kzz(index(pos));
  }

  arma::fvec kn(Npts);       // Number of points
  arma::vec kpos(Npts);      // Position in container
  arma::vec counted(Npts);   // Sets if pint belongs
  arma::uvec set_idx(Npts);  // Sets if pint belongs
  counted.zeros();

  cout << "Finding Unique and adding" << endl;
  int unique_points = 0;
  int non_unique_points = 0;
  particle_order order(Npts);

  for (unsigned int pos = 0; pos < kxx.n_elem; pos++) {
    if ((pos % 10000) == 0) {
      cout << "Pos = " << pos << " found " << unique_points << endl;
    }

    // Check to see if point has been collected
    if (counted(pos)) {
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
    int forward_pos = pos + 1;
    while (forward_pos < (int)kxx.n_elem) {
      // Difference in X (determines search range)
      float xdiff = abs(kxx_sorted(forward_pos) - kxt);

      // Difference in radius
      float ydiff = abs(kyy_sorted(forward_pos) - kyt);
      float zdiff = abs(kzz_sorted(forward_pos) - kzt);

      if (xdiff > 0.01) {
        // out of search range
        break;
      } else if ((ydiff > 0.01) || (zdiff > 0.01)) {
        // Just don't add point
      } else {
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
    if (con.point_inside(kxavg, kyavg, kzavg)) {
      // Actually add points
      con.put(order, unique_points, kxavg, kyavg, kzavg);

      // Update decoding array
      for (int dpos = 0; dpos < n; dpos++) {
        kpos(set_idx(dpos)) = unique_points;
        kn(set_idx(dpos)) = (float)n;
        counted(set_idx(dpos)) = 1;
      }
      unique_points++;
    } else {
      // Update decoding array
      // Update decoding array
      for (int dpos = 0; dpos < n; dpos++) {
        kpos(set_idx(dpos)) = -1;
        kn(set_idx(dpos)) = (float)n;
        counted(set_idx(dpos)) = 1;
      }
      non_unique_points += n;
    }
  }

  cout << "Min counted" << min(counted) << endl;
  cout << "Done combining, found " << unique_points << " unique points " << endl
       << flush;

  // Now calculate and copy back
  arma::fvec kw_calculated(unique_points);

  voronoicell cell;
  c_loop_order vl(con, order);

  cout << "vl start = " << vl.start() << endl;
  count = 0;
  if (vl.start()) do {
      if (count % 10 == 0) {
        cout << "\rProgress " << (100 * count / unique_points) << "%" << flush;
      }

      if (con.compute_cell(cell, vl)) {
        kw_calculated(count) = cell.volume();
      } else {
        kw_calculated(count) = 0.0;
      }
      // cout << kw_calculated(count) << endl;
      count++;
    } while (vl.inc() && (count < unique_points));
  cout << endl << "Counted " << count << " points in dcf calc" << endl;

  // Convert to sorted index
  cout << "Convert to Sorted" << endl;
  arma::fvec kw_sorted(Npts);
  for (int pos = 0; pos < Npts; pos++) {
    if (kpos(pos) == -1) {
      kw_sorted(pos) = 0.0;
    } else {
      kw_sorted(pos) = kw_calculated(kpos(pos)) / kn(pos);
    }
  }

  // Unsort
  cout << " Unsort " << endl;
  arma::fvec kw_unsorted(Npts);
  for (int pos = 0; pos < (int)Npts; pos++) {
    kw_unsorted(index(pos)) = kw_sorted(pos);
  }

  // Unsort and put in array
  cout << "Copy back" << endl;
  count = 0;
  for (int i = 0; i < kx.length(firstDim); i++) {
    float kww = kw_unsorted(count);
    kw(i) = kww * krr(i);  //*krr(i);
    count++;
  }

  cout << "Done with Copy" << endl;
}
