/*
   This contains Commandline Interface to the Recon
   */
#include "recon_lib.h"

using namespace NDarray;
using namespace std;

int main(int argc, char **argv){

    // ------------------------------------
    // Initialize Recon - reading the command line
    // ------------------------------------
    RECON recon(argc,argv);
    if(	recon.threads != -1){
        omp_set_num_threads( recon.threads );
    }

    // ------------------------------------
    // Read Data -
    // ------------------------------------

    MRI_DATA data;
    cout << "----Read Data-----" << endl;
    switch(recon.data_type){
        case(RECON::BENCHMARK):{
                                   // This runs a system benchmark


                               }break;

        case(RECON::SIMULATE):{
                                  // Use completely made up data

                              }break;

        case(RECON::PSF):{
                             // Use external Kx,Ky,Kz
                             data.read_external_data(recon.filename);
                             for(Array< Array<complex<float>,3>,2>::iterator miter=data.kdata.begin(); miter!=data.kdata.end(); miter++){
                                 *miter = complex<float>(1.0,0.0);
                             }
                         }break;

        case(RECON::PHANTOM):{
                                 // Use external Kx,Ky,Kz
                                 data.read_external_data(recon.filename);
                                 for(Array< Array<complex<float>,3>,2>::iterator miter=data.kdata.begin(); miter!=data.kdata.end(); miter++){
                                     *miter = complex<float>(0.0,0.0);
                                 }

                                 // Initialize Phantom
                                 cout << "Phantom " << endl;
                                 PHANTOM phantom;
                                 phantom.read_commandline(argc,argv);
                                 phantom.init(recon.rcxres,recon.rcyres,recon.rczres,recon.rcframes);


                                 // More accurate gridding for Phantom
                                 cout << "Grid " << endl;
                                 gridFFT phantom_gridding;
                                 phantom_gridding.kernel_type = KAISER_KERNEL;
                                 phantom_gridding.overgrid = 1.5;
                                 phantom_gridding.dwinX = 6;
                                 phantom_gridding.dwinY = 6;
                                 phantom_gridding.dwinZ = 6;
                                 phantom_gridding.precalc_gridding(phantom.IMAGE.length(firstDim),phantom.IMAGE.length(secondDim),phantom.IMAGE.length(thirdDim),data);


                                 GATING gate(argc,argv);
                                 // Weighting Array for Time coding
                                 Array< float, 3 >TimeWeight(data.kx(0).shape(),ColumnMajorArray<3>());
                                 gate.init(data,&(recon.rcframes));


                                 /*-----------------------------
                                   Collect Data
                                   ------------------------------*/
                                 int e= 0;
                                 Array< float,3 >kxE = data.kx(e);
                                 Array< float,3 >kyE = data.ky(e);
                                 Array< float,3 >kzE = data.kz(e);
                                 Array< float,3 >kwE = data.kw(e);

                                 cout << "Num coils = " << data.Num_Coils << endl;
                                 for(int coil =0; coil < data.Num_Coils; coil++){
                                     cout << "Getting phantom" << coil << ":" << flush;
                                     cout << "Smap " << endl << flush;
                                     phantom.update_smap_biotsavart(coil,data.Num_Coils);

                                     // Update Image (Entire Filled Array)
                                     Array<complex<float>,3>kdataE = data.kdata(e,coil); // Get one encoding, one coil
                                     kdataE=0;	// Zero data

                                     //Temp array to zero and fill.
                                     Array<complex<float>,3>temp(kdataE.shape(),ColumnMajorArray<3>());
                                     temp = complex<float>(0.0,0.0);

                                     for(int t =0; t < recon.rcframes; t++){
                                         // Get Image
                                         phantom.calc_image(t,recon.rcframes);

                                         // Weight Image
                                         TimeWeight = 1;
                                         //TimeWeight = kwE;
                                         gate.weight_data( TimeWeight,e, kxE, kyE,kzE,t,GATING::NON_ITERATIVE,GATING::TIME_FRAME);

                                         // Now Inverse Grid
                                         cout << " Inverse Grid :: " << t << endl;
                                         phantom_gridding.backward(phantom.IMAGE,temp,kxE,kyE,kzE,TimeWeight);
                                         kdataE += temp;
                                     }
                                 }

                                 // Add Noise
                                 phantom.add_noise( data.kdata );
                                 data.write_external_data("PhantomData/Phantom.h5");
                                 phantom.write_matlab_truth_script("PhantomData/");

                                 cout << "Only generating data" << endl;
                                 exit(1);


                             }break;

        default:
        case(RECON::EXTERNAL):{
                                  // Read in External Data Format
                                  data.read_external_data(recon.filename);
                                  data.scale_fov(1./recon.zoom_x,1./recon.zoom_y,1./recon.zoom_z);
                              }break;
    }

    // --------------------------------------------------
    // Code for recon (no PSD specific data/structures)
    // --------------------------------------------------
    recon.init_recon( argc,argv,data);
    Array< Array<complex<float>,3>,2 >X = recon.reconstruct_all_frames(data);

    // ------------------------------------
    // Post Processing + Export
    // ------------------------------------

    HDF5 clearRAW("Images.h5","w");
    for(int ee=0; ee<recon.rcencodes; ee++){
        for(int tt=0; tt<recon.rcframes; tt++){
            char fname[80];
            sprintf(fname,"X_%03d_%03d.dat.complex",ee,tt);
            Array< float, 3> IMAGE;
            IMAGE.setStorage( ColumnMajorArray<3>());
            IMAGE.resize(X(0,0).shape());
            IMAGE = abs( X(tt,ee) );
            clearRAW.AddH5Array("IMAGES",fname,IMAGE);
        }}

    for(int ee=0; ee<recon.rcencodes; ee++){
        for(int tt=0; tt<recon.rcframes; tt++){
            char fname[80];
            sprintf(fname,"X_%03d_%03d.dat.complex",ee,tt);
            ArrayWrite( X(tt,ee),fname);
            sprintf(fname,"X_%03d_%03d.dat",ee,tt);
            ArrayWriteMag( X(tt,ee),fname);
        }}


    return(0);
}
