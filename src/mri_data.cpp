#include "mri_data.h"

using arma::cx_fmat;
using arma::fvec;
using namespace NDarray;

void MRI_DATA::dump_stats(const string name, const Array<Array<float, 3>, 1> &in) {
  cout << name << endl;
  cout << "\tContainer Size = " << in.length(firstDim) << endl;
  cout << "\tElement size = " << in(0).length(firstDim) << " x " << in(0).length(secondDim) << " x " << in(0).length(thirdDim) << endl;
  float high = 0;
  float low = 9999;
  for (Array<Array<float, 3>, 1>::const_iterator miter = in.begin(); miter != in.end(); miter++) {
    high = max(high, max(*miter));
    low = min(low, min(*miter));
  }
  cout << "\tRange = " << low << " to " << high << endl;
}

void MRI_DATA::dump_stats(const string name, const Array<Array<complex<float>, 3>, 2> &in) {
  cout << name << endl;
  cout << "\tContainer Size = " << in.length(firstDim) << " x " << in.length(secondDim) << endl;
  cout << "\tElement size = " << in(0).length(firstDim) << " x " << in(0).length(secondDim) << " x " << in(0).length(thirdDim) << endl;
  float high = 0;
  float low = 9999;
  for (Array<Array<complex<float>, 3>, 2>::const_iterator miter = in.begin(); miter != in.end(); miter++) {
    float temp = max(abs(*miter));
    high = max(high, temp);

    temp = min(abs(*miter));
    low = min(low, temp);
  }
  cout << "\tRange = " << low << " to " << high << endl;
}

void MRI_DATA::scale_fov(float scale_x, float scale_y, float scale_z) {
  if ((scale_x == 1.0) && (scale_y == 1.0) && (scale_z == 1.0)) {
    return;
  }

  cout << "Scaling the fov by" << scale_x << "," << scale_y << "," << scale_z << endl;
  float scale_kx = (dft_needed(0)) ? (scale_x) : (scale_x);
  float scale_ky = (dft_needed(1)) ? (scale_y) : (scale_y);
  float scale_kz = (dft_needed(2)) ? (scale_z) : (scale_z);
  cout << "Scaling the kspace by" << scale_kx << "," << scale_ky << "," << scale_kz << endl;

  // Multiply kspace by inverse
  for (Array<Array<float, 3>, 1>::iterator miter = kx.begin(); miter != kx.end(); miter++) {
    (*miter) *= scale_kx;
  }

  for (Array<Array<float, 3>, 1>::iterator miter = ky.begin(); miter != ky.end(); miter++) {
    (*miter) *= scale_ky;
  }

  for (Array<Array<float, 3>, 1>::iterator miter = kz.begin(); miter != kz.end(); miter++) {
    (*miter) *= scale_kz;
  }

  // Multiple image by scale
  zfov *= std::abs(scale_z);
  yfov *= std::abs(scale_y);
  xfov *= std::abs(scale_x);

  if (sms_type == SMSon) {
    for (Array<Array<float, 3>, 2>::iterator miter = z.begin(); miter != z.end(); miter++) {
      (*miter) /= scale_z;
    }
  }
}

void MRI_DATA::convert_encodes_to_coils(int passes) {
  // This function treats encodes as sperate coils
  // Only works if encodes are identical
  int new_Num_Coils = (Num_Encodings / passes) * Num_Coils;
  // int new_Num_Encodings =  passes;
  /*
  this->kx.resizeAndPreserve(1);
  this->ky.resizeAndPreserve(1);
  this->kz.resizeAndPreserve(1);
  this->kw.resizeAndPreserve(1);
  this->kt.resizeAndPreserve(1);
  */

  cout << "Create temp array" << endl;
  Array<Array<complex<float>, 3>, 2> kdata_temp;
  kdata_temp.setStorage(ColumnMajorArray<2>());
  kdata_temp.resize(1, new_Num_Coils);

  cout << "Copy array" << endl;
  int count = 0;
  for (int ee = 0; ee < this->Num_Encodings; ee++) {
    for (int coil = 0; coil < this->Num_Coils; coil++) {
      kdata_temp(0, count).setStorage(ColumnMajorArray<3>());
      kdata_temp(0, count).resize(this->kdata(ee, coil).shape());
      kdata_temp(0, count) = this->kdata(ee, coil);
      count++;
    }
  }

  cout << "Copy Back" << endl;
  this->kdata.resize(1, new_Num_Coils);
  for (int pos = 0; pos < new_Num_Coils; pos++) {
    this->kdata(0, pos).setStorage(ColumnMajorArray<3>());
    this->kdata(0, pos).resize(kdata_temp(0, pos).shape());
    this->kdata(0, pos) = kdata_temp(0, pos);
  }

  // Reset size
  this->Num_Encodings = 1;
  this->Num_Coils = new_Num_Coils;
}

MRI_DATA MRI_DATA::subframe(int eStart, int eStop, int eStride) {
  MRI_DATA data2;
  data2.clone_attributes(*this);

  // Copy info
  data2.Num_Encodings = floor((1 + eStop - eStart) / eStride);
  data2.init_memory();

  int count = 0;
  for (int ee = eStart; ee <= eStop; ee += eStride) {
    data2.kx(count) = this->kx(ee);
    data2.ky(count) = this->ky(ee);
    data2.kz(count) = this->kz(ee);
    data2.kw(count) = this->kw(ee);
    data2.kt(count) = this->kt(ee);

    // SMS need z positions
    if (data2.sms_type == SMSon) {
      for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
        data2.z(sms_pos, count) = this->z(sms_pos, ee);
      }
    }

    // Copy the data
    for (int coil = 0; coil < data2.Num_Coils; coil++) {
      data2.kdata(count, coil) = this->kdata(ee, coil);
    }
    count = count + 1;
  }

  return (data2);
}

void MRI_DATA::stats(void) {
  if (kx.numElements() == 0) {
    cout << "Kspace does not exist" << endl;
  } else {
    dump_stats("Kx", kx);
    dump_stats("Ky", ky);
    dump_stats("Kz", kz);
    dump_stats("Kw", kw);
  }

  if (kdata.numElements() == 0) {
    cout << "Kdata does not exist yet" << endl;
  } else {
    dump_stats("Kdata", kdata);
  }

  // gating
  if (ecg.numElements() == 0) {
    cout << "Physiologic data does not exist yet" << endl;
  } else {
    cout << "Range ECG = " << Dmin(ecg) << " to " << Dmax(ecg) << endl;
    cout << "Range RESP = " << Dmin(resp) << " to " << Dmax(resp) << endl;
    cout << "Range TIME = " << Dmin(time) << " to " << Dmax(time) << endl;
    cout << "Range PREP = " << Dmin(prep) << " to " << Dmax(prep) << endl;
  }
}

void MRI_DATA::gating_stats(void) {
  if (ecg.numElements() == 0) {
    cout << "Physiologic data does not exist yet" << endl;
  } else {
    std::cout << "---------Gating Stats--------------" << std::endl;
    std::cout << "Range ECG = " << Dmin(ecg) << " to " << Dmax(ecg) << std::endl;
    std::cout << "Range Resp = " << Dmin(resp) << " to " << Dmax(resp) << std::endl;
    std::cout << "Range Time = " << Dmin(time) << " to " << Dmax(time) << std::endl;
    std::cout << "Range Prep = " << Dmin(prep) << " to " << Dmax(prep) << std::endl;
  }
}

void MRI_DATA::mod_time(double mod_time) {
  if (time.numElements() == 0) {
    cout << "Mod time - timing data does not exist yet" << endl;
  } else {
    NDarray::Array<NDarray::Array<double, 2>, 1> time;

    Array<Array<double, 2>, 1>::iterator miter = this->time.begin();
    Array<Array<double, 2>, 1>::iterator miter_end = this->time.end();

    for (; (miter != miter_end); miter++) {
      Array<double, 2>::iterator titer = (*miter).begin();
      Array<double, 2>::iterator titer_end = (*miter).end();
      for (; (titer != titer_end); titer++) {
        *titer = std::fmod(*titer, mod_time);
      }  // Inner time array
    }  // Encodes
  }  // If data exists
}

// Data for Whitening
void MRI_DATA::init_noise_samples(int total_samples) {
  noise_samples.setStorage(ColumnMajorArray<2>());
  noise_samples.resize(total_samples, Num_Coils);
  noise_samples = complex<float>(0, 0);
}

void MRI_DATA::clone_attributes(MRI_DATA &data) {
  Num_Encodings = data.Num_Encodings;

  Num_Coils = data.Num_Coils;
  Num_Frames = data.Num_Frames;

  dft_needed = data.dft_needed;
  trajectory_type = data.trajectory_type;

  xres = data.xres;
  yres = data.yres;
  zres = data.zres;
  tres = data.tres;

  sms_factor = data.sms_factor;
  sms_type = data.sms_type;

  zfov = data.xfov;
  yfov = data.yfov;
  zfov = data.zfov;
}

// Constructer for MRI data type
void MRI_DATA::init_memory(void) {
  cout << "Container Size = " << Num_Encodings << " x " << Num_Coils << endl;
  kx.resize(Num_Encodings);
  ky.resize(Num_Encodings);
  kz.resize(Num_Encodings);
  kw.resize(Num_Encodings);
  kt.resize(Num_Encodings);

  if (sms_type == SMSon) {
    z.resize(Num_Encodings, sms_factor);
  }

  kdata.setStorage(ColumnMajorArray<2>());
  kdata.resize(Num_Encodings, Num_Coils);

  kdata_gating.setStorage(ColumnMajorArray<2>());
  kdata_gating.resize(Num_Encodings, Num_Coils);

  // Times
  time.resize(Num_Encodings);
  ecg.resize(Num_Encodings);
  resp.resize(Num_Encodings);
  prep.resize(Num_Encodings);
}

void MRI_DATA::init_memory(int readouts, int shots, int slices) {
  cout << "Initializing Container for " << readouts << " x " << shots << " x "
       << slices << endl;
  init_memory();
  for (int e = 0; e < Num_Encodings; e++) {
    init_encode(e, readouts, shots, slices);
  }
}

void MRI_DATA::init_encode(int e, int readouts, int shots, int slices) {
  kx(e).setStorage(ColumnMajorArray<3>());
  kx(e).resize(readouts, shots, slices);

  ky(e).setStorage(ColumnMajorArray<3>());
  ky(e).resize(readouts, shots, slices);

  kz(e).setStorage(ColumnMajorArray<3>());
  kz(e).resize(readouts, shots, slices);

  if (sms_type == SMSon) {
    for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
      z(e, sms_pos).setStorage(ColumnMajorArray<3>());
      z(e, sms_pos).resize(readouts, shots, slices);
    }
  }

  kw(e).setStorage(ColumnMajorArray<3>());
  kw(e).resize(readouts, shots, slices);

  kt(e).setStorage(ColumnMajorArray<3>());
  kt(e).resize(readouts, shots, slices);

  for (int coil = 0; coil < Num_Coils; coil++) {
    kdata(e, coil).setStorage(ColumnMajorArray<3>());
    kdata(e, coil).resize(readouts, shots, slices);

    kdata_gating(e, coil).setStorage(ColumnMajorArray<2>());
    kdata_gating(e, coil).resize(shots, slices);
  }

  time(e).setStorage(ColumnMajorArray<2>());
  time(e).resize(shots, slices);

  ecg(e).setStorage(ColumnMajorArray<2>());
  ecg(e).resize(shots, slices);

  resp(e).setStorage(ColumnMajorArray<2>());
  resp(e).resize(shots, slices);

  prep(e).setStorage(ColumnMajorArray<2>());
  prep(e).resize(shots, slices);
}

void MRI_DATA::init_gating_kdata(int gating_samples) {
  cout << "Gating samples = " << gating_samples << endl;
  /*kdata_gating.setStorage( ColumnMajorArray<5>());
  kdata_gating.resize(gating_samples,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils);
  */
}

MRI_DATA::MRI_DATA(void) {
  // Initial Values
  Num_Encodings = -1;
  Num_Coils = -1;
  Num_Frames = 1;
  for (int dir = 0; dir < 3; dir++) {
    trajectory_type(dir) = NONCARTESIAN;
    dft_needed(dir) = true;
  }
  sms_type = SMSoff;
  sms_factor = 1;
  recon_fov = 0;
  recon_res = 0;
}

//---------------------------------------------------
//    This function allocates and reads all data into memory
//---------------------------------------------------

void MRI_DATA::demod_kdata(float demod) {
  for (int e = 0; e < Num_Encodings; e++) {
    for (int c = 0; c < Num_Coils; c++) {
      // Each Dataset
      for (int slice = 0; slice < kdata(e, c).length(thirdDim); slice++) {
#pragma omp parallel for
        for (int readout = 0; readout < kdata(e, c).length(secondDim);
             readout++) {
          for (int i = 0; i < kdata(e, c).length(firstDim); i++) {
            kdata(e, c)(i, readout, slice) *= polar<float>(
                1.0, demod * 2.0 * arma::datum::pi * kt(e)(i, readout, slice));
          }
        }
      }
    }
  }
}

//
//  Export for Berkely Advanced Reconstruction Tools
//
void MRI_DATA::write_bart_data(string fname_in) {
  const char *fname = fname_in.c_str();

  // Now determine the size of the ND array
  int total_shots = 0;
  int total_elements = 0;
  int xres = 1;
  // this->kx(0).length(firstDim);
  for (int pos = 0; pos < this->kx.length(firstDim); pos++) {
    total_shots = max((int)this->kx(pos).numElements(), total_shots);
    cout << "Encode " << pos << " elements = " << total_shots
         << " Shape = " << this->kx(pos).shape() << endl;
  }
  total_elements = total_shots;
  total_shots /= xres;
  cout << "Bart export setting total shots to " << total_shots << endl;

  // Need a new export order for pregated data
  // Time should be the last dimension
  Array<int, 1> export_order(this->Num_Encodings);
  {
    int count = 0;
    for (int t = 0; t < (this->Num_Frames); t++) {
      for (int e = 0; e < (this->Num_Encodings / this->Num_Frames); e++) {
        export_order(count) = e * this->Num_Frames + t;
        count++;
      }
    }
  }

  {
    // Export kdata
    Array<int, 1> Dims(11);
    Dims = 1;
    Dims(0) = 1;
    Dims(1) = xres;
    Dims(2) = total_shots;
    Dims(3) = this->Num_Coils;
    Dims(4) = this->Num_Encodings / this->Num_Frames;
    Dims(10) = this->Num_Frames;  // Encodes? Encodes*Frames

    // Header for kdata
    char name_hdr[1024];
    sprintf(name_hdr, "%s_data.hdr", fname);

    // Write the header
    FILE *fid;
    fid = fopen(name_hdr, "w");
    fprintf(fid, "# Dimensions\n");
    for (int i = 0; i < Dims.length(firstDim); i++) {
      fprintf(fid, "%d ", Dims(i));
    }
    fclose(fid);

    // Write the binary data
    char name_bin[1024];
    sprintf(name_bin, "%s_data.cfl", fname);

    remove(name_bin);
    ofstream ofs(name_bin, ios_base::binary | ios_base::app);
    for (int count = 0; count < this->Num_Encodings; count++) {
      // For reordered export
      int encode = export_order(count);

      for (int coil = 0; coil < this->Num_Coils; coil++) {
        // int subframe_shots = this->kdata(encode,coil).numElements();
        Array<complex<float>, 3>::iterator miter = this->kdata(encode, coil).begin();
        Array<complex<float>, 3>::iterator miter_stop = this->kdata(encode, coil).end();
        for (; (miter != miter_stop); miter++) {
          complex<float> val = (*miter);
          ofs.write((char *)&val, sizeof(complex<float>));
        }

        // If the subframe is smaller, pad with zeros to make bart happy
        if ((int)this->kdata(encode, coil).numElements() < total_elements) {
          complex<float> val(0.0, 0.0);
          for (int i = 0; i < (total_elements - (int)this->kdata(encode, coil).numElements()); i++) {
            ofs.write((char *)&val, sizeof(complex<float>));
          }
        }
      }
    }
    ofs.close();
  }

  {
    // Export kx,ky,kz
    Array<int, 1> Dims(11);
    Dims = 1;
    Dims(0) = 3;
    Dims(1) = xres;
    Dims(2) = total_shots;
    Dims(3) = 1;
    Dims(4) = this->Num_Encodings / this->Num_Frames;
    Dims(10) = this->Num_Frames;  // Encodes? Encodes*Frames

    // Header for kdata
    char name_hdr[1024];
    sprintf(name_hdr, "%s_traj.hdr", fname);

    // Write the header
    FILE *fid;
    fid = fopen(name_hdr, "w");
    fprintf(fid, "# Dimensions\n");
    for (int i = 0; i < Dims.length(firstDim); i++) {
      fprintf(fid, "%d ", Dims(i));
    }
    fclose(fid);

    // Write the binary data
    char name_bin[1024];
    sprintf(name_bin, "%s_traj.cfl", fname);

    remove(name_bin);
    ofstream ofs(name_bin, ios_base::binary | ios_base::app);
    for (int count = 0; count < this->Num_Encodings; count++) {
      // For reordered export
      int encode = export_order(count);

      // int subframe_shots = this->kdata(encode,coil).numElements();

      Array<float, 3>::iterator kx_iter = this->kx(encode).begin();
      Array<float, 3>::iterator ky_iter = this->ky(encode).begin();
      Array<float, 3>::iterator kz_iter = this->kz(encode).begin();

      Array<float, 3>::iterator kx_iter_stop = this->kx(encode).end();
      Array<float, 3>::iterator ky_iter_stop = this->ky(encode).end();
      Array<float, 3>::iterator kz_iter_stop = this->kz(encode).end();

      for (; (kx_iter != kx_iter_stop); kx_iter++, ky_iter++, kz_iter++) {
        complex<float> val = complex<float>(*kx_iter, 0);
        ofs.write((char *)&val, sizeof(complex<float>));
        val = complex<float>(*ky_iter, 0);
        ofs.write((char *)&val, sizeof(complex<float>));
        val = complex<float>(*kz_iter, 0);
        ofs.write((char *)&val, sizeof(complex<float>));
      }

      // If the subframe is smaller, pad with zeros to make bart happy
      if ((int)this->kx(encode).numElements() < total_elements) {
        complex<float> val(0.0, 0.0);
        for (int i = 0;
             i < 3 * (total_elements - (int)this->kx(encode).numElements());
             i++) {
          ofs.write((char *)&val, sizeof(complex<float>));
        }
      }
    }
    ofs.close();
  }

  {
    // Export kw
    Array<int, 1> Dims(11);
    Dims = 1;
    Dims(0) = 1;
    Dims(1) = xres;
    Dims(2) = total_shots;
    Dims(3) = 1;
    Dims(4) = this->Num_Encodings / this->Num_Frames;
    Dims(10) = this->Num_Frames;  // Encodes? Encodes*Frames

    // Header for kdata
    char name_hdr[1024];
    sprintf(name_hdr, "%s_dcf.hdr", fname);

    // Write the header
    FILE *fid;
    fid = fopen(name_hdr, "w");
    fprintf(fid, "# Dimensions\n");
    for (int i = 0; i < Dims.length(firstDim); i++) {
      fprintf(fid, "%d ", Dims(i));
    }
    fclose(fid);

    // Write the binary data
    char name_bin[1024];
    sprintf(name_bin, "%s_dcf.cfl", fname);

    remove(name_bin);
    ofstream ofs(name_bin, ios_base::binary | ios_base::app);
    for (int count = 0; count < this->Num_Encodings; count++) {
      // For reordered export
      int encode = export_order(count);

      // int subframe_shots = this->kdata(encode,coil).numElements();

      Array<float, 3>::iterator kw_iter = this->kw(encode).begin();
      Array<float, 3>::iterator kw_iter_stop = this->kw(encode).end();

      for (; (kw_iter != kw_iter_stop); kw_iter++) {
        complex<float> val = complex<float>(*kw_iter, 0);
        ofs.write((char *)&val, sizeof(complex<float>));
      }

      // If the subframe is smaller, pad with zeros to make bart happy
      if ((int)this->kw(encode).numElements() < total_elements) {
        complex<float> val(0.0, 0.0);
        for (int i = 0;
             i < (total_elements - (int)this->kw(encode).numElements()); i++) {
          ofs.write((char *)&val, sizeof(complex<float>));
        }
      }
    }
    ofs.close();
  }
}

//---------------------------------------------------
//  Temporary Function to Write Data ( will be replaced by ismrmd )
//---------------------------------------------------

void MRI_DATA::write_external_data(string fname) {
  HDF5 file = HDF5(fname, "w");

  cout << "Exporting Attributes" << endl;

  // Add dimensions
  file.AddH5Scaler("Kdata", "Num_Encodings", Num_Encodings);
  file.AddH5Scaler("Kdata", "Num_Coils", Num_Coils);
  file.AddH5Scaler("Kdata", "Num_Frames", Num_Frames);

  // 2D/3D Cartesian/Non-Cartesian
  file.AddH5Scaler("Kdata", "trajectory_typeX", (int)trajectory_type(0));
  file.AddH5Scaler("Kdata", "trajectory_typeY", (int)trajectory_type(1));
  file.AddH5Scaler("Kdata", "trajectory_typeZ", (int)trajectory_type(2));

  file.AddH5Scaler("Kdata", "dft_neededX", (int)dft_needed(0));
  file.AddH5Scaler("Kdata", "dft_neededY", (int)dft_needed(1));
  file.AddH5Scaler("Kdata", "dft_neededZ", (int)dft_needed(2));

  for (int encode = 0; encode < kdata.length(firstDim); encode++) {
    cout << "Exporting " << encode << endl;
    {
      try {
        stringstream ss;
        ss << "KT_E" << encode;
        string s = ss.str();
        file.AddH5Array("Kdata", s.c_str(), kt(encode));
      } catch (...) {
        cout << "Can't export KT for encode " << encode << endl;
      }
    }

    {
      try {
        stringstream ss;
        ss << "KX_E" << encode;
        string s = ss.str();
        file.AddH5Array("Kdata", s.c_str(), kx(encode));
      } catch (...) {
        cout << "Can't export KX for encode " << encode << endl;
      }
    }

    {
      try {
        stringstream ss;
        ss << "KY_E" << encode;
        string s = ss.str();
        file.AddH5Array("Kdata", s.c_str(), ky(encode));
      } catch (...) {
        cout << "Can't export KY for encode " << encode << endl;
      }
    }

    {
      try {
        stringstream ss;
        ss << "KZ_E" << encode;
        string s = ss.str();
        file.AddH5Array("Kdata", s.c_str(), kz(encode));
      } catch (...) {
        cout << "Can't export KZ for encode " << encode << endl;
      }
    }

    if (sms_type == SMSon) {
      cout << "Export SMS Z2" << endl
           << flush;
      for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
        stringstream ss;
        ss << "Z_E" << encode << "_S" << sms_pos;
        string s = ss.str();
        file.AddH5Array("Kdata", s.c_str(), z(encode, sms_pos));
      }
    }

    {
      try {
        stringstream ss;
        ss << "KW_E" << encode;
        string s = ss.str();
        file.AddH5Array("Kdata", s.c_str(), kw(encode));
      } catch (...) {
        cout << "Can't export KW for encode " << encode << endl;
      }
    }

    cout << "Exporting data" << endl
         << flush;
    for (int coil = 0; coil < kdata.length(secondDim); coil++) {
      {
        try {
          stringstream ss;
          ss << "KData_E" << encode << "_C" << coil;
          string s = ss.str();
          file.AddH5Array("Kdata", s.c_str(), kdata(encode, coil));
        } catch (...) {
          cout << "Can't export Kdata for encode " << encode << ", coil "
               << coil << endl;
        }
      }
    }

    // Gating
    if (ecg.numElements() != 0) {
      cout << "Exporting Gating " << endl
           << flush;

      {
        try {
          stringstream ss;
          ss << "ECG_E" << encode;
          string s = ss.str();
          file.AddH5Array("Gating", s.c_str(), ecg(encode));
        } catch (...) {
          cout << "Can't export ECG data" << endl;
        }
      }

      {
        try {
          stringstream ss;
          ss << "RESP_E" << encode;
          string s = ss.str();
          file.AddH5Array("Gating", s.c_str(), resp(encode));
        } catch (...) {
          cout << "Can't export Resp data" << endl;
        }
      }

      {
        try {
          stringstream ss;
          ss << "PREP_E" << encode;
          string s = ss.str();
          file.AddH5Array("Gating", s.c_str(), prep(encode));
        } catch (...) {
          cout << "Can't export PREP data" << endl;
        }
      }

      {
        try {
          stringstream ss;
          ss << "TIME_E" << encode;
          string s = ss.str();
          file.AddH5Array("Gating", s.c_str(), time(encode));
        } catch (...) {
          cout << "Can't export TIME data" << endl;
        }
      }

      for (int coil = 0; coil < kdata.length(secondDim); coil++) {
        try {
          stringstream ss;
          ss << "K0_E" << encode << "_C" << coil;
          string s = ss.str();
          file.AddH5Array("Gating", s.c_str(), kdata_gating(encode, coil));
        } catch (...) {
          cout << "Can't export k0 data" << endl;
        }
      }
    }
  }

  // Noise Samples
  if (noise_samples.numElements() != 0) {
    file.AddH5Array("Kdata", "Noise", noise_samples);
  }
}

void MRI_DATA::read_external_data(string fname) {
  HDF5 file = HDF5(fname, "r");

  cout << "Reading External File " << fname << endl;

  // Read atrributes common to all
  file.ReadH5Scaler("Kdata", "Num_Encodings", &Num_Encodings);
  file.ReadH5Scaler("Kdata", "Num_Coils", &Num_Coils);

  init_memory();

  // 2D/3D Cartesian/Non-Cartesian
  int temp1;
  file.ReadH5Scaler("Kdata", "trajectory_typeX", &temp1);
  trajectory_type(0) = static_cast<MRI_DATA::TrajType>(temp1);

  file.ReadH5Scaler("Kdata", "trajectory_typeY", &temp1);
  trajectory_type(1) = static_cast<MRI_DATA::TrajType>(temp1);

  file.ReadH5Scaler("Kdata", "trajectory_typeZ", &temp1);
  trajectory_type(2) = static_cast<MRI_DATA::TrajType>(temp1);

  file.ReadH5Scaler("Kdata", "dft_neededX", &temp1);
  dft_needed(0) = static_cast<bool>(temp1);

  file.ReadH5Scaler("Kdata", "dft_neededY", &temp1);
  dft_needed(1) = static_cast<bool>(temp1);

  file.ReadH5Scaler("Kdata", "dft_neededZ", &temp1);
  dft_needed(2) = static_cast<bool>(temp1);

  for (int encode = 0; encode < kdata.length(firstDim); encode++) {
    cout << "Importing " << encode << endl;

    {
      stringstream ss;
      ss << "KT_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Kdata", s, kt(encode));
    }

    {
      stringstream ss;
      ss << "KX_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Kdata", s, kx(encode));
    }

    {
      stringstream ss;
      ss << "KY_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Kdata", s, ky(encode));
    }

    {
      stringstream ss;
      ss << "KZ_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Kdata", s, kz(encode));
    }

    if (sms_type == SMSon) {
      for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
        stringstream ss;
        ss << "Z_E" << encode << "_S" << sms_pos;
        string s = ss.str();
        file.ReadH5Array("Kdata", s, z(encode, sms_pos));
      }
    }

    {
      stringstream ss;
      ss << "KW_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Kdata", s, kw(encode));
    }

    for (int coil = 0; coil < kdata.length(secondDim); coil++) {
      {
        stringstream ss;
        ss << "KData_E" << encode << "_C" << coil;
        string s = ss.str();
        file.ReadH5Array("Kdata", s, kdata(encode, coil));
      }
    }

    cout << "Read Gating " << endl
         << flush;
    {
      stringstream ss;
      ss << "ECG_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Gating", s.c_str(), ecg(encode));
    }

    {
      stringstream ss;
      ss << "RESP_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Gating", s.c_str(), resp(encode));
    }

    {
      stringstream ss;
      ss << "PREP_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Gating", s.c_str(), prep(encode));
    }

    {
      stringstream ss;
      ss << "TIME_E" << encode;
      string s = ss.str();
      file.ReadH5Array("Gating", s.c_str(), time(encode));
    }
  }

  // Noise Samples
  try {
    file.ReadH5Array("Kdata", "Noise", noise_samples);
  } catch (...) {
    cout << "Can't import noise samples " << endl;
  }

  // file.ReadH5Array( "Gating","kdata_gating",kdata_gating);
}

/** Coil compress data with set number of coils
 *
 */
void MRI_DATA::coilcompress(float Num_VCoils, float kr_thresh) {
  cout << "about to compress coils" << endl
       << flush;

  if (Num_VCoils > this->Num_Coils) {
    std::cout
        << "Target virtual coils is greater than actual coils, doing nothing"
        << std::endl;
    return;
  }

  tictoc ctimer;

  // calculate kr
  cout << "Using kr = 0  to " << kr_thresh << " for SVD in coil compression "
       << endl
       << flush;

  // Get the values of kspace points less than kr_thresh
  int Num_Pixels = 0;
  float kr_thresh_squared = kr_thresh * kr_thresh;

#pragma omp parallel for reduction(+ : Num_Pixels)
  for (int encode = 0; encode < kx.length(firstDim); encode++) {
    // Assign iterators to go over the data
    Array<float, 3>::const_iterator kx_iter = this->kx(encode).begin();
    Array<float, 3>::const_iterator ky_iter = this->ky(encode).begin();
    Array<float, 3>::const_iterator kz_iter = this->kz(encode).begin();
    Array<float, 3>::const_iterator kx_iter_end = this->kx(encode).end();
    Array<float, 3>::const_iterator ky_iter_end = this->ky(encode).end();
    Array<float, 3>::const_iterator kz_iter_end = this->kz(encode).end();

    // Iterate through and check to see if less than the kr_thresh
    for (; (kx_iter != kx_iter_end) && (ky_iter != ky_iter_end) &&
           (kz_iter != kz_iter_end);
         kx_iter++, ky_iter++, kz_iter++) {
      float kr =
          pow((*kx_iter), 2.0) + pow((*ky_iter), 2.0) + pow((*kz_iter), 2.0);
      if (kr < kr_thresh_squared) {
        Num_Pixels++;
      }
    }
  }

  // Allocate memory for the data
  arma::cx_fmat all_data;
  all_data.zeros(Num_Pixels, Num_Coils);

  cout << "Collect Data" << endl
       << flush;
  ctimer.tic();
  int idx = 0;
  for (int encode = 0; encode < Num_Encodings; encode++) {
    for (int slice = 0; slice < kx(encode).length(thirdDim); slice++) {
      for (int view = 0; view < kx(encode).length(secondDim); view++) {
        for (int pos = 0; pos < kx(encode).length(firstDim); pos++) {
          float kr = pow(kx(encode)(pos, view, slice), 2.0) +
                     pow(ky(encode)(pos, view, slice), 2.0) +
                     pow(kz(encode)(pos, view, slice), 2.0);
          if (kr < kr_thresh_squared) {
            for (int coil = 0; coil < this->Num_Coils; coil++) {
              all_data(idx, coil) = kdata(encode, coil)(pos, view, slice);
            }
            idx++;
          }
        }  // pos
      }  // view
    }  // slice
  }  // encode
  cout << "Copied pixels = " << idx << endl
       << flush;
  cout << "took " << ctimer << " s to copy data" << endl;

  cout << "SVD " << endl
       << flush;
  arma::fvec s;
  arma::cx_fmat U;
  arma::cx_fmat V;
  ctimer.tic();
  arma::svd_econ(U, s, V, all_data);
  cout << "took " << ctimer << " s to perform SVD" << endl;
  s = s / s(0);

  s.print("S");

  arma::cx_fmat VV = V.cols(0, (int)Num_VCoils - 1);
  VV.print("V");

  cout << "Rotate to " << Num_VCoils << " coils " << endl
       << flush;

  ctimer.tic();
  for (int encode = 0; encode < Num_Encodings; encode++) {
    for (int slice = 0; slice < kx(encode).length(thirdDim); slice++) {
#pragma omp parallel for
      for (int view = 0; view < kx(encode).length(secondDim); view++) {
        Array<complex<float>, 1> tmp_data(Num_Coils);
        for (int pos = 0; pos < kx(encode).length(firstDim); pos++) {
          for (int coil = 0; coil < Num_Coils; coil++) {
            tmp_data(coil) = kdata(encode, coil)(pos, view, slice);
          }  // coil

          for (int vcoil = 0; vcoil < Num_VCoils; vcoil++) {
            complex<float> tmp(0.0, 0.0);
            for (int coil = 0; coil < Num_Coils; coil++) {
              tmp += tmp_data(coil) * V(coil, vcoil);
            }  // coil
            kdata(encode, vcoil)(pos, view, slice) = tmp;
          }  // vcoil
        }  // pos
      }  // view
    }  // slice
  }  // encode
  cout << "took " << ctimer << " s to rotate data" << endl;

  /* Temp data structure
  Array< Array<complex<float>,3>, 2>
  kdata2(Num_Encodings,Num_VCoils,ColumnMajorArray<2>()); for( int encode = 0;
  encode < Num_Encodings; encode++){ for( int coil =0; coil < Num_VCoils;
  coil++){ kdata2(encode,coil) = this->kdata(encode,coil);
          }
  }
  //this->kdata.resize(Num_Encodings, Num_VCoils);
  this->kdata.reference(kdata2);
  */
  this->Num_Coils = Num_VCoils;

  cout << "done with coil compression" << endl;
}

template <typename T>
arma::Mat<complex<float> > whitening_matrix_calc(const Array<complex<T>, 2> &noise_samples, int reference_coil) {
  std::cout << "WHITEN::Noise Samples : " << noise_samples.length(firstDim) << std::endl;
  std::cout << "WHITEN::Reference coil : " << reference_coil << std::endl;
  std::cout << "WHITEN::Noise Pre-Whitening" << std::endl;

  // Copy into matrix (coil x samples)
  int Num_Coils = noise_samples.length(secondDim);
  arma::cx_mat NoiseData = arma::randu<arma::cx_mat>(noise_samples.length(secondDim), noise_samples.length(firstDim));
  for (int coil = 0; coil < Num_Coils; coil++) {
    for (int i = 0; i < noise_samples.length(firstDim); i++) {
      NoiseData(coil, i) = noise_samples(i, coil);
    }
  }

  // Swap reference coil to first coil to avoid rotation
  arma::cx_mat CSWAP = arma::eye<arma::cx_mat>(Num_Coils, Num_Coils);
  if (reference_coil != 0) {
    std::cout << "WHITEN::Swap " << reference_coil << " to 0 for referencing" << std::endl;
    CSWAP.swap_rows(0, reference_coil);
  }
  NoiseData = CSWAP * NoiseData;

  std::cout << "WHITEN::Calc Cov" << std::endl;
  arma::cx_mat CV = NoiseData * NoiseData.t();
  arma::mat CV_ABS = arma::abs(CV);
  cout << "WHITEN::Noise Covariance" << endl;
  cout << "  size = " << CV_ABS.n_rows << " x " << CV_ABS.n_cols << endl;
  // cout << CV_ABS << endl;

  std::cout << "WHITEN::Calc whitening matrix" << std::endl;
  arma::cx_mat V = chol(CV);
  arma::cx_mat VT = V.t();
  arma::cx_mat Decorr = VT.i();

  // Swap the rotation matrix back
  if (reference_coil != 0) {
    Decorr = CSWAP * Decorr * CSWAP;

    cout << "WHITEN::Noise Decorrelation matrix" << endl;
    cout << "  size = " << Decorr.n_rows << " x " << Decorr.n_cols << endl;
    // cout << Decorr << endl;
  }

  // Test Whitening
  arma::cx_mat W = NoiseData;
  arma::cx_mat temp = arma::randu<arma::cx_mat>(Num_Coils);
  for (int i = 0; i < noise_samples.length(firstDim); i++) {
    // Copy coil data
    for (int coil = 0; coil < Num_Coils; coil++) {
      temp(coil, 0) = noise_samples(i, coil);
    }

    // Whiten
    arma::cx_mat temp2 = Decorr * temp;

    // Copy back
    for (int coil = 0; coil < Num_Coils; coil++) {
      W(coil, i) = temp2(coil, 0);
    }
  }

  arma::cx_mat CV_POST = W * W.t();
  arma::mat CV_POST_ABS = arma::abs(CV_POST);
  cout << "WHITEN::Noise Covariance post whiten" << endl;
  cout << "  size = " << W.n_rows << " x " << W.n_cols << endl;
  // cout << CV_POST_ABS << endl;

  arma::cx_fmat whitening_matrix;
  whitening_matrix.copy_size(Decorr);
  for (unsigned int i = 0; i < Decorr.n_rows; i++) {
    for (unsigned int j = 0; j < Decorr.n_cols; j++) {
      whitening_matrix(i, j) = Decorr(i, j);
    }
  }

  return whitening_matrix;
}

arma::cx_fmat MRI_DATA::get_whitening_matrix(const Array<complex<float>, 2> &noise_samples, int reference_coil) {
  arma::cx_fmat decorr = whitening_matrix_calc<float>(noise_samples, reference_coil);
  return (decorr);
}

arma::cx_fmat MRI_DATA::get_whitening_matrix(const Array<complex<double>, 2> &noise_samples, int reference_coil) {
  arma::cx_fmat decorr = whitening_matrix_calc<double>(noise_samples, reference_coil);
  return (decorr);
}

//--------------------------------------------------
//  Whiten Data
//--------------------------------------------------

void MRI_DATA::whiten(void) {
  if (noise_samples.numElements() == 0) {
    cout << "Noise Samples do not exist:: Can't whiten data" << endl;
    return;
  }

  // Get the whitening matrix
  arma::cx_fmat Decorr = MRI_DATA::get_whitening_matrix(noise_samples, 0);

  // Now Whiten Actual Data
  cout << "Whiten all data" << endl;
  for (int e = 0; e < Num_Encodings; e++) {
    for (int k = 0; k < kdata(e, 0).length(thirdDim); k++) {
#pragma omp parallel for
      for (int j = 0; j < kdata(e, 0).length(secondDim); j++) {
        arma::cx_fmat AA;
        AA.zeros(Num_Coils, kdata(e, 0).length(firstDim));  // Working memory

        // Copy into memory
        for (int coil = 0; coil < Num_Coils; coil++) {
          for (int i = 0; i < kdata(e, 0).length(firstDim); i++) {
            AA(coil, i) = kdata(e, coil)(i, j, k);
          }
        }

        // Transform to matrix
        arma::cx_fmat temp2 = Decorr * AA;

        // Copy Back
        for (int coil = 0; coil < Num_Coils; coil++) {
          for (int i = 0; i < kdata(e, 0).length(firstDim); i++) {
            kdata(e, coil)(i, j, k) = temp2(coil, i);
          }
        }
      }
    }
  }
}

void MRI_DATA::add_noise(float noise_factor) {
  cout << "Adding " << (100.0 * (noise_factor - 1.0)) << "% more noise." << endl;
  float noise_std = sqrt(noise_factor * noise_factor - 1.0) / sqrt(2.0);  // variance is additive but standard deviation is not
  for (int e = 0; e < Num_Encodings; e++) {
    for (int coil = 0; coil < Num_Coils; coil++) {
      arma::cube noise_real = arma::randn<arma::cube>(
          kdata(e, coil).length(firstDim), kdata(e, coil).length(secondDim),
          kdata(e, coil).length(thirdDim));
      arma::cube noise_imag = arma::randn<arma::cube>(
          kdata(e, coil).length(firstDim), kdata(e, coil).length(secondDim),
          kdata(e, coil).length(thirdDim));
      for (int k = 0; k < kdata(e, coil).length(thirdDim); k++) {
        for (int j = 0; j < kdata(e, coil).length(secondDim); j++) {
          for (int i = 0; i < kdata(e, coil).length(firstDim); i++) {
            complex<float> noise(noise_real(i, j, k) * noise_std,
                                 noise_imag(i, j, k) * noise_std);
            kdata(e, coil)(i, j, k) = kdata(e, coil)(i, j, k) + noise;
          }
        }
      }
    }
  }
}

arma::cx_fmat covariance(arma::cx_fmat X, int Nx, int Ny) {
  // Initialize covariance matrix
  arma::cx_fmat cov;
  cov.zeros(Nx, Nx);

  // Loop Through Coils
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Nx; j++) {
      std::complex<float> meani = 0;
      std::complex<float> meanj = 0;
      for (int k = 0; k < Ny; k++) {
        meani += X(i, k);
        meanj += X(j, k);
      }
      meani /= (float)Ny;
      meanj /= (float)Ny;
      for (int k = 0; k < Ny; k++) {
        cov(i, j) += std::conj(X(i, k) - meani) * (X(j, k) - meanj);
      }
      cov(i, j) /= ((float)Ny - 1);
    }
  }

  return cov;
}
