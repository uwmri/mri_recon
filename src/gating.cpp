/************************************************
View sharing techniques

The class uses information from Times.dat file to create
a mask used during reconstruction of each frame in recon.cxx

Initial Author:
        Grzegorz Bauman (gbauman@wisc.edu)

Changelog:
        Stan Kruger (sjkruger@wisc.edu) 130107
        tornado filter should now be more robust.  Bins should be much closer to
the same number of projections, and can now overlap as desired depending on <
"frames," "vs_a," and "vs_b" > init 2013-03-08  KMJ: Major changes. Renamed
many variables to practical names. Fixed tornado filer for non equidistant
spacing. Added respiratory gating. Better commenting,etc.

Init:
        recon_binary -f data_header.txt -rcframes 32 -vs tornado -vs_a 1 -vs_b 5
-vs_shape 2 GATING vs(argc,**argv); vs.createmask(TimeWeight,timesE,t);


*************************************************/

#include "gating.h"
#include "io_templates.hpp"
using namespace NDarray;

arma::vec array_to_vec(const Array<Array<double, 2>, 1> &A) {
  int total_count = 0;
  for (Array<Array<double, 2>, 1>::const_iterator miter = A.begin();
       miter != A.end(); miter++) {
    total_count += (*miter).numElements();
  }
  arma::vec out(total_count);

  int count = 0;
  for (Array<Array<double, 2>, 1>::const_iterator miter1 = A.begin();
       miter1 != A.end(); miter1++) {
    for (Array<double, 2>::const_iterator miter2 = (*miter1).begin();
         miter2 != (*miter1).end(); miter2++) {
      out(count) = (*miter2);
      count++;
    }
  }
  return (out);
}

void vec_to_array(Array<Array<double, 2>, 1> &A, arma::vec &out) {
  int count = 0;
  for (Array<Array<double, 2>, 1>::iterator miter1 = A.begin();
       miter1 != A.end(); miter1++) {
    for (Array<double, 2>::iterator miter2 = (*miter1).begin();
         miter2 != (*miter1).end(); miter2++) {
      (*miter2) = out(count);
      count++;
    }
  }
}

GATING::GATING(void) {
}

// Setup of
GATING::GATING(int numarg, const char **pstring) {
  // Setting default values, configurable
  wdth_low = 1;
  wdth_high = 4;

  vs_type = NONE;
  tornado_shape = VIPR;  // Kr^2 shape
  kmax = 128;            // TEMP
  gate_type = GATE_NONE;

  // Retrospective scan time control (s)
  start_proj = 0;
  end_proj = 100;

  // Respiratory Efficiency
  correct_resp_drift = 0;
  adaptive_resp_window = 5.0;
  resp_gate_efficiency = 0.5;
  resp_phase_lower = 0.0;
  resp_phase_upper = 0.5;
  resp_phase_type = 0;
  resp_gate_weight = 3;  // This controls exponential decay of resp weights for soft-thresholding.
  resp_gate_type = RESP_NONE;
  resp_gate_signal = BELLOWS;
  gate_min_quantile = 0.0;
  gate_max_quantile = 1.0;
  resp_sign = 1.0;

// Catch command line switches
#define trig_flag(num, name, val)             \
  }                                           \
  else if (strcmp(name, pstring[pos]) == 0) { \
    val = num;
#define float_flag(name, val)                 \
  }                                           \
  else if (strcmp(name, pstring[pos]) == 0) { \
    pos++;                                    \
    val = atof(pstring[pos]);
#define int_flag(name, val)                   \
  }                                           \
  else if (strcmp(name, pstring[pos]) == 0) { \
    pos++;                                    \
    val = atoi(pstring[pos]);
#define char_flag(name, val)                  \
  }                                           \
  else if (strcmp(name, pstring[pos]) == 0) { \
    pos++;                                    \
    strcpy(val, pstring[pos]);

  for (int pos = 0; pos < numarg; pos++) {
    if (strcmp("-viewshare_type", pstring[pos]) == 0) {
      pos++;
      if (pos == numarg) {
        cout << "Please provide viewshare type (-h for usage)" << endl;
        exit(1);
        trig_flag(TORNADO, "tornado", vs_type);
        trig_flag(HIST_MODE, "hist", vs_type);
        trig_flag(NONE, "none", vs_type);

      } else {
        cout << "Please provide viewshare type (-h for usage)" << endl;
        exit(1);
      }
    } else if (strcmp("-gating_type", pstring[pos]) == 0) {
      pos++;
      if (pos == numarg) {
        cout << "Please provide gating type (-h for usage)" << endl;
        exit(1);
        trig_flag(RETRO_ECG, "retro_ecg", gate_type);
        trig_flag(RESP, "resp", gate_type);
        trig_flag(ECG, "ecg", gate_type);
        trig_flag(TIME, "time", gate_type);
        trig_flag(PREP, "prep", gate_type);

      } else {
        cout << "Please provide gating type..none/dft/diff/pca" << endl;
        exit(1);
      }
    } else if (strcmp("-resp_gate", pstring[pos]) == 0) {
      pos++;
      if (pos == numarg) {
        cout << "Please provide respiratory gating type..thresh/weight (-h for "
                "usage)"
             << endl;
        exit(1);
        trig_flag(RESP_THRESH, "thresh", resp_gate_type);
        trig_flag(RESP_PHASE, "phase", resp_gate_type);
        trig_flag(RESP_WEIGHT, "weight", resp_gate_type);
        trig_flag(RESP_HARD, "hard", resp_gate_type);

      } else {
        cout << "Please provide respiratory gating type..thresh/weight/hard "
                "(-h for usage)"
             << endl;
        exit(1);
      }
    } else if (strcmp("-resp_gate_signal", pstring[pos]) == 0) {
      pos++;
      if (pos == numarg) {
        cout << "Please provide a data source for estimation of respiratory "
                "phase..bellows/internal (-h for usage)"
             << endl;
        exit(1);
        trig_flag(BELLOWS, "bellows", resp_gate_signal);
        trig_flag(DC_DATA, "internal", resp_gate_signal);

      } else {
        cout << "Please provide a data source for estimation of respiratory "
                "phase..bellows/internal (-h for usage)"
             << endl;
        exit(1);
      }

      int_flag("-vs_wdth_low", wdth_low);
      int_flag("-vs_wdth_high", wdth_high);
      trig_flag(1, "-correct_resp_drift", correct_resp_drift);
      float_flag("-resp_gate_efficiency", resp_gate_efficiency);
      float_flag("-adaptive_resp_window", adaptive_resp_window);
      float_flag("-resp_phase_lower", resp_phase_lower);
      float_flag("-resp_phase_upper", resp_phase_upper);
      float_flag("-resp_phase_type", resp_phase_type);
      float_flag("-resp_gate_weight", resp_gate_weight);
      float_flag("-resp_sign", resp_sign);
      float_flag("-gate_min_quantile", gate_min_quantile);
      float_flag("-gate_max_quantile", gate_max_quantile);

      trig_flag(1, "-retro_scan_time", retro_scan_time);
      float_flag("-start_proj", start_proj);
      float_flag("-end_proj", end_proj);
    }
  }

  if (resp_gate_type == RESP_THRESH) {
    cout << "Using moving average threshold based respiratory gating" << endl;
  } else if (resp_gate_type == RESP_PHASE) {
    cout << "Using thresholding with arbitrary limits for respiratory gating" << endl;
  } else if (resp_gate_type == RESP_WEIGHT) {
    cout << "Using (fuzzy) weight based respiratory gating" << endl;
  } else if (resp_gate_type == RESP_HARD) {
    cout << "Using hard threshold respiratory gating" << endl;
  }
}

void GATING::help_message() {
  cout << "----------------------------------------------" << endl;
  cout << "   View Sharing Control " << endl;
  cout << "----------------------------------------------" << endl;
  cout << "Usage:" << endl;
  help_flag("-viewshare_type []", "view sharing method tornado/none/hist (defult=none)");
  help_flag("", "  tornado = variable width filter in kr");
  help_flag("", "  none = images at equal time intervals, no sharing between frames");
  help_flag("", "  hist = images with equal data points");

  help_flag("-gating_type []", "how to gate images");
  help_flag("", "  resp = respiratory phases");
  help_flag("", "  ecg = gate by cardiac");
  help_flag("", "  retro_ecg = retrospective gate by cardiac");
  help_flag("", "  time = bin by acquisition time");
  help_flag("", "  prep = bin by time from prep pulses");

  help_flag("-resp_gate []", "In addition to other gating, perform respiratory gating");
  help_flag("", "  thresh = threshold respiratory values");
  help_flag("", "  phase = threshold with upper and lower bounds");
  help_flag("", "  weight = downweight bad values (see Johnson et al. MRM 67(6):1600");

  help_flag("-resp_gate_signal", "Specify source for the data used to estimate respiratory phase");
  help_flag("",
            "  bellows = signal from respiratory bellows belt in gating file "
            "(default)");
  help_flag("",
            "  internal = use dc/low spatial frequency information extracted "
            "from acquired data");

  cout << "Filter parameters for tornado:" << endl;
  help_flag("-vs_wdth_low []", "width in the center of k-space in frames");
  help_flag("-vs_wdth_high []", "width in the periphery of k-space in frames");
  help_flag("-vs_vipr_tornado", "3D radial tornado (default)");
  help_flag("-vs_radial_tornado", "2D radial tornado");

  cout << "Control for Resp Data" << endl;
  help_flag("-correct_resp_drift", "Median filter with 10s interval");
  help_flag("-resp_gate_efficiency", "Fraction of data to accept");
  help_flag("-resp_sign", "Flips respiratory waveform");
  help_flag("-resp_phase_lower", "Lower bound of threshold");
  help_flag("-resp_phase_upper", "Upper bound of threshold");
  help_flag("-resp_phase_type", "Phase of respiratory waveform to reconstruct (1=expiration, 2=inspiration, 0=both)");
  help_flag("-resp_gate_weight", "Soft-Gating Decay Constant (Exponential)");
  help_flag("-adaptive_resp_window", "Length of window to use for thresholding (default=5.0s)");

  cout << "Control for ECG Data" << endl;
  help_flag("-bad_ecg_filter", "Filter Bad ECG Vals (>10,000ms)");

  cout << "Control for retrospective scan time adjustment" << endl;
  help_flag("-retro_scan_time", "Use a subset of acquired projections (time sorted) to reconstruct");
  help_flag("-start_proj", "Projection range start");
  help_flag("-end_proj", "Projection range end");
}

/*----------------------------------------------
     Smooths the Resp and subtracts off to correct drift
 *----------------------------------------------*/
void GATING::filter_resp(const MRI_DATA &data) {
  // Assume all the data is contigous
  double min_time = Dmin(data.time);
  double max_time = Dmax(data.time);
  cout << "Time range [ " << min_time << " to " << max_time << " ] span = " << (max_time - min_time) << endl;

  // Time in seconds to grab median from
  double fsize = 5.0;

  // Gate times
  cout << "Sorting Gate Data by Acquisition Time" << endl;

  // Use Aradillo Sort function
  arma::vec time = array_to_vec(data.time);
  arma::vec resp = array_to_vec(data.resp);
  arma::vec time_linear_resp(resp);

  // Sort
  arma::uvec idx = arma::sort_index(time);

  // Copy Resp
  idx.save("Sorted.dat", arma::raw_ascii);
  for (int i = 0; i < (int)this->gate_times.numElements(); i++) {
    time_linear_resp(i) = resp(idx(i));
  }
  time_linear_resp.save("RSorted.dat", arma::raw_ascii);

  // Now Filter
  cout << "Filtering Resp Data by" << fsize << endl;
  for (int i = 0; i < (int)resp.n_elem; i++) {
    int start = i - fsize;
    int stop = i + fsize;
    if (start < 0) {
      stop = 2 * fsize;
      start = 0;
    }

    if (stop >= (int)resp.n_elem) {
      stop = resp.n_elem - 1;
      start = resp.n_elem - 1 - 2 * fsize;
    }

    double thresh = median(time_linear_resp.rows(start, stop));
    resp(idx(i)) -= thresh;
  }
  resp.save("RFiltered.dat", arma::raw_ascii);

  // Put back into an array
  vec_to_array(this->gate_times, resp);
}

void GATING::init(const MRI_DATA &data, int target_frames) {
  init_resp_gating(data);
  init_time_resolved(data, target_frames);
}

int GATING::number_of_frames(void) {
  return (this->recon_frames);
}

float GATING::temporal_resolution(void) {
  return (actual_temporal_resolution);
}

/* Assuming multiple readouts and channels - combine the data*/
NDarray::Array<NDarray::Array<complex<float>, 2>, 1> GATING::combine_kspace_channels(
    const NDarray::Array<NDarray::Array<complex<float>, 2>, 2> &kdata_gating) {
  cout << "Combining k-space gating data" << endl;

  int Num_Channels = kdata_gating.length(secondDim);

  int Num_Times = 0;
  for (int e = 0; e < kdata_gating.length(firstDim); e++) {
    Num_Times += kdata_gating(e, 0).numElements();
  }

  arma::cx_fmat full_data;
  full_data.zeros(Num_Times, Num_Channels);

  cout << "Collect Data" << endl;
  for (int coil = 0; coil < kdata_gating.length(secondDim); coil++) {
    int pos = 0;
    for (int e = 0; e < kdata_gating.length(firstDim); e++) {
      // Now grab the values
      for (Array<complex<float>, 2>::const_iterator miter = kdata_gating(e, coil).begin(); miter != kdata_gating(e, coil).end(); miter++) {
        full_data(pos, coil) = (*miter);
        pos++;
      }
    }
  }

  cout << "SVD " << endl
       << flush;
  arma::fvec s;
  arma::cx_fmat U;
  arma::cx_fmat V;
  arma::svd_econ(U, s, V, full_data);

  arma::cx_fmat VV = V.cols(0, 0);
  cout << "Rotate " << endl
       << flush;
  full_data = full_data * VV;

  cout << "Copy Back" << endl;
  int nencodes = kdata_gating.length(firstDim);
  Array<Array<complex<float>, 2>, 1> combined_kdata(nencodes);
  for (int e = 0; e < nencodes; e++) {
    combined_kdata(e).setStorage(ColumnMajorArray<2>());
    combined_kdata(e).resize(kdata_gating(e, 0).shape());
  }

  cout << "Collect Data" << endl;
  int pos = 0;
  for (int e = 0; e < kdata_gating.length(firstDim); e++) {
    // Now grab the values
    for (Array<complex<float>, 2>::iterator miter = combined_kdata(e).begin();
         miter != combined_kdata(e).end(); miter++) {
      (*miter) = full_data(pos, 0);
    }
  }

  return (combined_kdata);
}

void GATING::init_resp_gating(const MRI_DATA &data) {
  std::cout << "Initializing Respiratory Gating" << endl;
  this->resp_weight.resize(data.Num_Encodings);
  for (int e = 0; e < this->resp_weight.length(firstDim); e++) {
    this->resp_weight(e).setStorage(ColumnMajorArray<2>());
    this->resp_weight(e).resize(1, 1);
  }

  if (resp_gate_type != RESP_NONE) {
    this->resp_weight.resize(data.resp.shape());
    for (int e = 0; e < this->resp_weight.length(firstDim); e++) {
      this->resp_weight(e).setStorage(ColumnMajorArray<2>());
      this->resp_weight(e).resize(data.resp(e).shape());
    }

    switch (resp_gate_signal) {
      case (BELLOWS): {
        cout << "Respiratory gating using bellows belt waveform" << endl;
        this->resp_weight = resp_sign * data.resp;
      } break;

      case (DC_DATA): {
        cout << "Respiratory gating using internal DC waveform(s)" << endl;

        /* Do a SVD - coil/channel combine - might need more complex for slice
         * based*/
        Array<Array<complex<float>, 2>, 1> kdc = combine_kspace_channels(data.kdata_gating);

        cout << "Sizes of SVD combined K0" << endl;
        cout << kdc(0).shape() << endl;
        cout << kdc.shape() << endl;

        /*  The question is now how to process/convert this data for the gating
         * algorithm below? For now, just use the vector magnitude over all
         * coils of k = 0, then apply a moving average to upsample to full
         * temporal resolution					*/

        int views_per_grid = 64;
        cout << "views_per_grid = " << views_per_grid << endl;

        // Use Armadillo Sort function
        arma::vec time(this->resp_weight(0).numElements());
        arma::cx_fvec signal(this->resp_weight(0).numElements());

        // Put into Matrix for Armadillo
        {
          int count = 0;
          for (int e = 0; e < data.resp.length(firstDim); e++) {
            for (Array<double, 2>::const_iterator miter = data.time(e).begin(); miter != data.time(e).end(); miter++) {
              time(count) = (*miter);
              count++;
            }
          }
        }

        // Sort
        arma::uvec idx = arma::sort_index(time);
        arma::cx_fvec time_linear_signal(this->resp_weight.numElements());
        for (int i = 0; i < (int)this->resp_weight.numElements(); i++) {
          time_linear_signal(i) = signal(idx(i));
        }

        cout << "Moving average filter" << endl;

        /*  Now compute moving average, with a window 2x the views per grid,
         * store the result in the original order */
        int total_samples = this->resp_weight.numElements();
        arma::cx_fvec filtered_time_linear_signal(this->resp_weight.numElements());

        for (int pos = 0; (int)pos < total_samples; pos++) {
          int start = max(pos - views_per_grid, 0);
          int stop = min(pos + views_per_grid, total_samples - 1);

          filtered_time_linear_signal(pos) = 0;
          for (int i = start; i <= stop; i++) {
            filtered_time_linear_signal(pos) += time_linear_signal(i);
          }
          filtered_time_linear_signal(pos) /= (double)(stop - start);
        }

        // Put back into
        arma::vec filtered_signal(resp_weight.numElements());
        for (int i = 0; i < (int)resp_weight.numElements(); i++) {
          filtered_signal(idx(i)) = abs(filtered_time_linear_signal(i));
        }

        // Sort back
        {
          int count = 0;
          for (int e = 0; e < this->resp_weight.length(thirdDim); e++) {
            for (int slice = 0; slice < this->resp_weight.length(secondDim); slice++) {
              for (int view = 0; view < this->resp_weight.length(firstDim); view++) {
                this->resp_weight(view, slice, e) = resp_sign * filtered_signal(count);
                count++;
              }
            }
          }
        }

        cout << "DC signal processing complete" << endl;
      } break;

        /*No default - should lead to enum error*/
    }

    switch (resp_gate_type) {
      case (RESP_THRESH): {
        cout << "Time Sorting Data" << endl;

        // Use Armadillo Sort function
        arma::vec time = array_to_vec(data.time);
        arma::vec resp = array_to_vec(this->resp_weight);
        time.save("Time.txt", arma::raw_ascii);
        resp.save("Resp.txt", arma::raw_ascii);

        // Vectors for working on data
        int N = resp.n_elem;
        arma::vec arma_resp_weight = arma::zeros<arma::vec>(N);
        arma::vec time_linear_resp = arma::zeros<arma::vec>(N);
        arma::vec time_sort_resp_weight = arma::zeros<arma::vec>(N);

        // Sort
        arma::uvec time_idx = arma::sort_index(time);
        time_idx.save("Sorted.txt", arma::raw_ascii);

        // Copy Resp
        for (int i = 0; i < N; i++) {
          time_linear_resp(i) = resp(time_idx(i));
        }
        time_linear_resp.save("TimeResp.txt", arma::raw_ascii);

        // Size of histogram
        cout << "Time range = " << (Dmax(data.time) - Dmin(data.time)) << endl;
        int fsize = (int)(adaptive_resp_window / ((Dmax(data.time) - Dmin(data.time)) / resp.n_elem));  // 5s filter / delta time

        float min_proj;
        float max_proj;
        if (retro_scan_time) {
          min_proj = start_proj;
          max_proj = end_proj;
          cout << "Retrospective projection range = " << min_proj << " to " << max_proj << endl;
        } else {
          min_proj = 0;
          max_proj = N;
        }

        // Now Filter
        cout << "Thresholding Data Frame Size = " << fsize << endl;
        for (int i = 0; i < N; i++) {
          int start = i - fsize;
          int stop = i + fsize;
          if (start < 0) {
            stop = 2 * fsize;
            start = 0;
          }

          if (stop >= N) {
            stop = N - 1;
            start = time_linear_resp.n_elem - 1 - 2 * fsize;
          }

          arma::vec temp = time_linear_resp.rows(start, stop);
          arma::vec temp2 = sort(temp);
          double thresh = temp2((int)((double)temp2.n_elem * (1.0 - resp_gate_efficiency)));

          arma_resp_weight(time_idx(i)) = (time_linear_resp(i) >= thresh) && (i >= min_proj) && (i < max_proj) ? (1.0) : (0.0);
          time_sort_resp_weight(i) = arma_resp_weight(time_idx(i));
        }
        time_sort_resp_weight.save("TimeWeight.txt", arma::raw_ascii);
        arma_resp_weight.save("Weight.txt", arma::raw_ascii);

        // Copy Back
        vec_to_array(this->resp_weight, arma_resp_weight);

      } break;

      case (RESP_PHASE): {
        cout << "Time Sorting Data with Upper and Lower Bounds" << endl;

        // Use Aradillo Sort function
        arma::vec time = array_to_vec(data.time);
        arma::vec resp = array_to_vec(this->resp_weight);
        time.save("Time.txt", arma::raw_ascii);
        resp.save("Resp.txt", arma::raw_ascii);

        // Vectors for working on data
        int N = resp.n_elem;
        arma::vec arma_resp_weight = arma::zeros<arma::vec>(N);
        arma::vec time_linear_resp = arma::zeros<arma::vec>(N);
        arma::vec time_sort_resp_weight = arma::zeros<arma::vec>(N);
        arma::vec derivatives = arma::zeros<arma::vec>(N);

        // Sort
        arma::uvec time_idx = arma::sort_index(time);
        time_idx.save("Sorted.txt", arma::raw_ascii);

        // Copy Resp
        for (int i = 0; i < N; i++) {
          time_linear_resp(i) = resp(time_idx(i));
        }
        time_linear_resp.save("TimeResp.txt", arma::raw_ascii);

        // Size of histogram
        cout << "Time range = " << (Dmax(data.time) - Dmin(data.time)) << endl;
        int fsize = (int)(adaptive_resp_window / ((Dmax(data.time) - Dmin(data.time)) / resp.n_elem));  // 5s filter / delta time
        int slope_window = (int) fsize/(adaptive_resp_window*10);                                       // 100 ms window for determining slope
        
        float min_proj;
        float max_proj;
        if (retro_scan_time) {
          min_proj = start_proj;
          max_proj = end_proj;
          cout << "Retrospective projection range = " << min_proj << " to " << max_proj << endl;
        } else {
          min_proj = 0;
          max_proj = N;
        }

        // Now Filter
        cout << "Thresholding Window Size = " << fsize << " (" << adaptive_resp_window << "s filter)" << endl;
        cout << "Phase Window Size (0.1s filter) = " << slope_window << endl;
        for (int i = 0; i < N; i++) {
          int start = i - fsize;
          int stop = i + fsize;
          int slope_start = i - slope_window;
          int slope_stop = i + slope_window;
          if (start < 0) {
            stop = 2 * fsize;
            start = 0;
          }
          if (stop >= N) {
            stop = N - 1;
            start = time_linear_resp.n_elem - 1 - 2 * fsize;
          }

          if (slope_start < 0) {
            slope_stop = 2 * slope_window;
            slope_start = 0;
          }
          if (slope_stop >= N) {
            slope_stop = N - 1;
            slope_start = time_linear_resp.n_elem - 1 - 2 * slope_window;
          }

          arma::vec temp = time_linear_resp.rows(start, stop);
          arma::vec temp2 = sort(temp);

          arma::vec temp3 = time_linear_resp.rows(slope_start, slope_stop);
          double deriv = (arma::mean(arma::diff(temp3)) >= 0) ? (1.0) : (-1.0);
          derivatives(i) = deriv;

          double thresh1 = temp2((int)((double)temp2.n_elem * (1.0 - resp_phase_upper)));
          double thresh2;
          if (resp_phase_lower == 0.0) {
            thresh2 = temp2((int)(temp2.n_elem-1));
          } else {
            thresh2 = temp2((int)((double)temp2.n_elem * (1.0 - resp_phase_lower)));
          }
          
          if (resp_phase_type == 0) {
            arma_resp_weight(time_idx(i)) = (time_linear_resp(i) >= thresh1) && (time_linear_resp(i) < thresh2) && (i >= min_proj) && (i < max_proj) ? (1.0) : (0.0);
          } else if (resp_phase_type == 1) {
            arma_resp_weight(time_idx(i)) = (time_linear_resp(i) >= thresh1) && (time_linear_resp(i) < thresh2) && (deriv == 1.0) && (i >= min_proj) && (i < max_proj) ? (1.0) : (0.0);
          } else if (resp_phase_type == 2) {
            arma_resp_weight(time_idx(i)) = (time_linear_resp(i) >= thresh1) && (time_linear_resp(i) < thresh2) && (deriv == -1.0) && (i >= min_proj) && (i < max_proj) ? (1.0) : (0.0);
          } else {
            cout << "invalid argument for resp_phase_type" << endl;
          }
          time_sort_resp_weight(i) = arma_resp_weight(time_idx(i));
        }
        time_sort_resp_weight.save("TimeWeight.txt", arma::raw_ascii);
        arma_resp_weight.save("Weight.txt", arma::raw_ascii);
        derivatives.save("RespDerivatives.txt", arma::raw_ascii);

        // Copy Back
        vec_to_array(this->resp_weight, arma_resp_weight);

      } break;

      case (RESP_HARD): {
        std::cout << "Setting hard threshold without moving average"
                  << std::endl;

        // Use Aradillo Sort function
        arma::vec resp = array_to_vec(this->resp_weight);
        resp.save("Resp.txt", arma::raw_ascii);
        int N = resp.n_elem;

        // Sort
        arma::vec resp_sorted = arma::sort(resp);

        // Get threshold
        int thresh_idx = (int)((float)resp_sorted.n_elem * (1.0 - this->resp_gate_efficiency));
        float resp_thresh = resp_sorted(thresh_idx);
        cout << "Hard resp thresh = " << resp_thresh << endl;

        float min_proj;
        float max_proj;
        if (retro_scan_time) {
          min_proj = start_proj;
          max_proj = end_proj;
          cout << "Retrospective projection range = " << min_proj << " to " << max_proj << endl;
        } else {
          min_proj = 0;
          max_proj = N;
        }

        // Copy to weight
        arma::vec arma_resp_weight = arma::zeros<arma::vec>(resp.n_elem);
        for (int i = 0; i < (int)resp.n_elem; i++) {
          arma_resp_weight(i) = (resp(i) > resp_thresh) && (i >= min_proj) && (i < max_proj) ? (1.0) : (0.0);
        }

        // Copy Back
        vec_to_array(this->resp_weight, arma_resp_weight);

      } break;

      case (RESP_WEIGHT): {
        if (correct_resp_drift == 1) {
          cout << "Correcting Drift" << endl;
          filter_resp(data);
        }
        cout << "Time Sorting Data" << endl;

        // Use Aradillo Sort function
        arma::vec time = array_to_vec(data.time);
        arma::vec resp = array_to_vec(data.resp);
        time.save("Time.txt", arma::raw_ascii);
        resp.save("Resp.txt", arma::raw_ascii);

        // Vectors for working on data
        int N = resp.n_elem;
        arma::vec arma_resp_weight = arma::zeros<arma::vec>(N);
        arma::vec time_linear_resp = arma::zeros<arma::vec>(N);
        arma::vec time_sort_resp_weight = arma::zeros<arma::vec>(N);

        // Sort
        arma::uvec time_idx = arma::sort_index(time);
        time_idx.save("Sorted.txt", arma::raw_ascii);

        // Copy Resp
        for (int i = 0; i < N; i++) {
          time_linear_resp(i) = resp(time_idx(i));
        }
        time_linear_resp.save("TimeResp.txt", arma::raw_ascii);

        // This version uses a 10 second window to define the threshold and
        // weights. Size of histogram
        cout << "Time range = " << (Dmax(data.time) - Dmin(data.time)) << endl;
        int fsize = (int)(5.0 / ((Dmax(data.time) - Dmin(data.time)) / resp.n_elem));  // 10s filter / delta time

        arma::vec temp = time_linear_resp;
        arma::vec temp2 = sort(temp);
        double thresh = temp2((int)((double)temp2.n_elem * (1.0 - resp_gate_efficiency)));
        double decay_const = resp_gate_weight / arma::max(temp);
        // cout << "Decay Const: " << decay_const << endl;

        float min_proj;
        float max_proj;
        if (retro_scan_time) {
          min_proj = start_proj;
          max_proj = end_proj;
        } else {
          min_proj = 0;
          max_proj = N;
        }

        // Now Filter
        cout << "Stencil Radius = " << fsize << endl;
        for (int i = 0; i < N; i++) {
          // For exponentially decaying weights. Implemented Soft-Gating as in:
          // https://doi.org/10.1002/mrm.26958
          arma_resp_weight(time_idx(i)) = (time_linear_resp(i) >= thresh) && (i >= min_proj) && (i < max_proj) ? (1.0) : (exp(-decay_const * ((thresh - time_linear_resp(i)))));

          time_sort_resp_weight(i) = arma_resp_weight(time_idx(i));
        }

        // Save
        time_sort_resp_weight.save("TimeWeight.txt", arma::raw_ascii);
        arma_resp_weight.save("Weight.txt", arma::raw_ascii);

        // Copy Back
        vec_to_array(this->resp_weight, arma_resp_weight);

      } break;

      case (RESP_NONE):
      default: {
        return;
      }
    }
  }
}

void GATING::init_time_resolved(const MRI_DATA &data, int target_frames) {
  cout << "Initializing Time resolved for " << (target_frames) << " target frames" << endl;

  // Create Array and Fill with Base
  this->gate_times.resize(data.Num_Encodings);
  for (int e = 0; e < data.Num_Encodings; e++) {
    this->gate_times(e).setStorage(ColumnMajorArray<2>());
    this->gate_times(e).resize(1, 1);
  }

  switch (gate_type) {
    case (RESP): {
      cout << "Using Resp gate" << endl;
      for (int e = 0; e < data.Num_Encodings; e++) {
        this->gate_times(e).resize(data.resp(e).shape());
        this->gate_times(e) = data.resp(e);
      }

      if (correct_resp_drift == 1) {
        cout << "Correcting Drift" << endl;
        filter_resp(data);
      }
    } break;

    case (RETRO_ECG):
    case (ECG): {
      cout << "Using ECG " << endl;
      for (int e = 0; e < data.Num_Encodings; e++) {
        this->gate_times(e).resize(data.ecg(e).shape());
        this->gate_times(e) = data.ecg(e);
      }

    } break;

    case (TIME): {
      cout << "Using Time Resolved" << endl;
      for (int e = 0; e < data.Num_Encodings; e++) {
        this->gate_times(e).resize(data.time(e).shape());
        if (!data.dft_needed(2)) {
          for (int shot = 0; shot < data.time(e).length(0); shot++) {
            for (int slice = 0; slice < data.time(e).length(1); slice++) {
              int projPerSlice = data.time(e).length(0) / data.zres;
              int firstShotOfSlice = floor(shot / projPerSlice) * projPerSlice;
              gate_times(e)(shot, slice) = data.time(e)(shot, slice) - data.time(e)(firstShotOfSlice, slice);
            }
          }
        } else {
          gate_times(e) = data.time(e);
        }
      }
    } break;

    case (PREP): {
      cout << "Using Prep Timer" << endl;
      for (int e = 0; e < data.Num_Encodings; e++) {
        this->gate_times(e).resize(data.prep(e).shape());
        this->gate_times(e) = data.prep(e);
      }
    } break;

    case (GATE_NONE): {
      this->recon_frames = 1;
      std::cout << "Data is time averaged" << std::endl;
      return;
    } break;
  }

  /* Bounds check*/
  if (target_frames == -1) {
    this->recon_frames = data.tres;
    cout << "Using native number of time frames:" << this->recon_frames << endl;
  } else {
    this->recon_frames = target_frames;
  }

  // Get Range
  double max_time = 0;
  double min_time = 0;

  switch (gate_type) {
    case (RETRO_ECG): {
      int total_elements = 0;
      for (int e = 0; e < gate_times.length(firstDim); e++) {
        this->gate_times(e) -= min_time;
        total_elements += this->gate_times(e).numElements();
      }

      min_time = 0;

      // Use Median to set value
      arma::vec temp(total_elements);
      int count = 0;
      for (int e = 0; e < this->gate_times.length(firstDim); e++) {
        for (Array<double, 2>::iterator miter = this->gate_times(e).begin(); miter != this->gate_times(e).end(); miter++) {
          temp(count) = *miter;
          count++;
        }
      }
      max_time = 2.0 * median(temp);
      cout << "Retro ECG::RR is estimated to be " << max_time << endl;
    } break;

    case (RESP): {
      // For respiratory weighting it is useful to reject outliers in the data
      int total_elements = 0;
      for (int e = 0; e < gate_times.length(firstDim); e++) {
        this->gate_times(e) -= min_time;
        total_elements += this->gate_times(e).numElements();
      }

      min_time = 0;
      max_time = 0;

      // Create histogram
      arma::vec temp(total_elements);
      int count = 0;
      for (int e = 0; e < this->gate_times.length(firstDim); e++) {
        for (Array<double, 2>::iterator miter = this->gate_times(e).begin(); miter != this->gate_times(e).end(); miter++) {
          temp(count) = *miter;
          count++;
        }
      }
      arma::vec P(2);
      P(0) = this->gate_min_quantile;
      P(1) = this->gate_max_quantile;
      arma::vec V = arma::quantile(temp, P);
      cout << "Respiratory Quantiles = " << endl
           << V << endl;
      double frame_width = (V(1) - V(0)) / (this->recon_frames - 1.0);

      // Add half frame on ends
      max_time = V(1) + 0.5 * frame_width;
      min_time = V(0) - 0.5 * frame_width;

    } break;

    default: {
      max_time = Dmax(this->gate_times);
      min_time = Dmin(this->gate_times);
    } break;
  }

  // Rescale to Frames
  scale_time = (float)this->recon_frames / (max_time - min_time) * (1 - 1e-9);  // Extra factor is to map last point to < frames
  offset_time = min_time;

  // Temporal resolution
  actual_temporal_resolution = (max_time - min_time) / ((float)this->recon_frames);

  cout << "Time Range :: " << min_time << " to " << max_time << endl;
  cout << "Actual temporal resolution = " << actual_temporal_resolution << endl;
  cout << " Gate offset = " << offset_time << endl;
  cout << " Gate scale = " << scale_time << endl;

  for (int e = 0; e < this->gate_times.length(firstDim); e++) {
    this->gate_times(e) -= offset_time;
    this->gate_times(e) *= scale_time;
  }

  cout << "Min Time - Post scale = " << Dmin(this->gate_times) << endl;
  cout << "Max Time - Post scale = " << Dmax(this->gate_times) << endl;

  /* Histogram*/
  {
    arma::vec temp(this->recon_frames);
    temp.fill(0);
    for (int e = 0; e < this->gate_times.length(firstDim); e++) {
      for (Array<double, 2>::iterator miter = this->gate_times(e).begin(); miter != this->gate_times(e).end(); miter++) {
        int pos = floor(*miter);
        if ((pos < (this->recon_frames)) && (pos >= 0)) {
          temp(pos)++;
        }
      }
    }

    // Export
    cout << "Values per frames" << endl;
    for (int i = 0; i < (this->recon_frames); i++) {
      cout << " Frame " << i << " ,count = " << temp(i) << endl;
    }
  }

  cout << "Setting up view share" << endl;
  switch (vs_type) {
    case (TORNADO): {
      // Set time points
      gate_frames.resize(this->recon_frames);
      for (int i = 0; i < (this->recon_frames); i++) {
        gate_frames(i) = 0.5 + (double)i;
      }
      switch (tornado_shape) {
        case (VIPR): {
          kmax = max(data.kx(0) * data.kx(0) + data.ky(0) * data.ky(0) + data.kz(0) * data.kz(0));
          kmax = sqrt(kmax);
        } break;

        case (RADIAL): {
          kmax = max(data.kx(0) * data.kx(0) + data.ky(0) * data.ky(0));
          kmax = sqrt(kmax);
        } break;

        case (FLAT): {
          kmax = 999;
        }
      }
      cout << "Kmax = " << kmax << endl;
    } break;

    case (NONE): {
      this->gate_times = floor(this->gate_times);
      cout << "Max gate time = " << Dmin(gate_times) << endl;
      cout << "Min gate time = " << Dmax(gate_times) << endl;
    } break;

    case (HIST_MODE): {
      cout << "gating::histmode::Sorting Data into Histogram" << endl;

      // Use Aradillo Sort function
      int total_views = 0;
      for (int e = 0; e < this->gate_times.length(firstDim); e++) {
        total_views += this->gate_times(e).numElements();
      }
      cout << "gating::histmode::Total views = " << total_views << endl;

      arma::vec time_sort(total_views);

      // Copy into array
      int count = 0;
      for (int e = 0; e < this->gate_times.length(firstDim); e++) {
        for (Array<double, 2>::iterator miter = this->gate_times(e).begin(); miter != this->gate_times(e).end(); miter++) {
          time_sort(count) = (*miter);
          count++;
        }
      }
      int Ncount = count;

      // Sort
      arma::uvec sort_temp = sort_index(time_sort);
      arma::uvec sort_idx = sort_index(sort_temp);

      // Now Split into frames
      count = 0;
      for (int e = 0; e < this->gate_times.length(firstDim); e++) {
        for (Array<double, 2>::iterator miter = this->gate_times(e).begin(); miter != this->gate_times(e).end(); miter++) {
          int t_frame = (int)((double)(sort_idx(count) * (this->recon_frames)) / Ncount);
          *miter = t_frame;
          count++;
        }
      }

      cout << "gating::histmode::Max gate time = " << Dmin(gate_times) << endl;
      cout << "gating::histmode::Min gate time = " << Dmax(gate_times) << endl;

    } break;
  }  // Switch vs_type
}

void GATING::weight_data(Array<float, 3> &Tw, int e, const Array<float, 3> &kx,
                         const Array<float, 3> &ky, const Array<float, 3> &kz,
                         int t, WeightType w_type, FrameType comp_type) {
  switch (resp_gate_type) {
    case (RESP_HARD):
    case (RESP_THRESH):
    case (RESP_PHASE):
    case (RESP_WEIGHT): {
      for (int k = 0; k < Tw.length(thirdDim); k++) {
        for (int j = 0; j < Tw.length(secondDim); j++) {
          for (int i = 0; i < Tw.length(firstDim); i++) {
            Tw(i, j, k) *= this->resp_weight(e)(j, k);
          }
        }
      }
      cout << "Resp weighting done" << endl
           << flush;
    } break;

    default: {
      cout << "No Additional Resp Gating" << endl
           << flush;
    }
  }

  if ((gate_type != GATE_NONE) && (comp_type != COMPOSITE)) {
    switch (vs_type) {
      case (TORNADO): {
        tornado_weight(Tw, e, kx, ky, kz, t, w_type);
      } break;

      case (HIST_MODE):
      case (NONE): {
        hist_weight(Tw, e, t);
      } break;
    }
  }

  // Normalize Weighting
  // float sum_Tw = sum(Tw);
  // cout << "Sum Time Weight = " << sum_Tw << endl;
  // Tw /= sum_Tw;
}

/*
Simple 1 to 1 frames. No sharing
*/

void GATING::hist_weight(Array<float, 3> &Tw, int e, int t) {
  for (int k = 0; k < Tw.length(thirdDim); k++) {
    for (int j = 0; j < Tw.length(secondDim); j++) {
      for (int i = 0; i < Tw.length(firstDim); i++) {
        // Get K-space Radius
        Tw(i, j, k) *= (floor(this->gate_times(e)(j, k)) == t) ? (1.0) : (0.0);
      }
    }
  }
}

/* Tornado-like filter in k-space in rcframe units
  ________________________________________________
              /wdth_low\
            /           \
          /              \
        /                 \
     /                     \
    <-----wdth_high-------->

        KMJ: Rewrote entirely. Otherwise would be wrong for all but center out
  without ramp sampling/variable density!
*/
void GATING::tornado_weight(Array<float, 3> &Tw, int e,
                            const Array<float, 3> &kx,
                            const Array<float, 3> &ky,
                            const Array<float, 3> &kz, int t,
                            WeightType w_type) {
  double current_time = gate_frames(t);

  for (int k = 0; k < Tw.length(thirdDim); k++) {
    for (int j = 0; j < Tw.length(secondDim); j++) {
      for (int i = 0; i < Tw.length(firstDim); i++) {
        double t_diff = abs(this->gate_times(e)(j, k) - current_time);

        // Get K-space Radius
        float kr = 0;
        float k_power = 0;
        switch (tornado_shape) {
          case (FLAT): {
            kr = 0.0;
            k_power = 0.0;
          } break;

          case (RADIAL): {
            kr = sqrt(kx(i, j, k) * kx(i, j, k) + ky(i, j, k) * ky(i, j, k));
            k_power = 1.0;
          } break;

          case (VIPR): {
            kr = sqrt(kx(i, j, k) * kx(i, j, k) + ky(i, j, k) * ky(i, j, k) + kz(i, j, k) * kz(i, j, k));
            k_power = 2.0;
          } break;
        }

        double wdth = 0.5 * ((wdth_high - wdth_low) * pow(kr / kmax, k_power) + wdth_low);
        if (w_type == ITERATIVE) {
          // Don't Divide (i.e. perform density compensation)
          Tw(i, j, k) *= (t_diff < wdth) ? (1.0) : (0.0);
        } else {
          Tw(i, j, k) *= (t_diff < wdth) ? (1. / wdth) : (0.0);
        }
      }
    }
  }
  // ArrayWrite(Tw,"TimeWeight.dat");
}
