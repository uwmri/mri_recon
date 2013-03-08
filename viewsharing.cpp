/************************************************
View sharing techniques

The class uses information from Times.dat file to create
a mask used during reconstruction of each frame in recon.cxx

Initial Author:
        Grzegorz Bauman (gbauman@wisc.edu)

Changelog: 
        Stan Kruger (sjkruger@wisc.edu) 130107
        tornado filter should now be more robust.  Bins should be much closer to the same number of projections, and can now overlap as desired depending on < "frames," "vs_a," and "vs_b" >

Init:
        recon_binary -f data_header.txt -rcframes 32 -vs tornado -vs_a 1 -vs_b 5 -vs_shape 2
        VIEWSHARING vs(argc,**argv);
        vs.createmask(TimeWeight,timesE,t);


*************************************************/

#include "viewsharing.h"

// Setup of VIEWSHARING
VIEWSHARING::VIEWSHARING(int numarg, char **pstring) {

        // Setting default values, configurable
        a = 1;
        b = 4;
        vs_type = 0;
        vs_shape = 3;
        // Catch command line switches
#define trig_flag(num,name,val) if(strcmp(name,pstring[pos]) == 0){ val=num; }
#define float_flag(name,val) if(strcmp(name,pstring[pos]) == 0){ pos++; val=atof(pstring[pos]); }
#define int_flag(name,val) if(strcmp(name,pstring[pos]) == 0){ pos++; val=atoi(pstring[pos]); }
#define char_flag(name,val)   if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]); }

        for(int pos=0; pos<numarg; pos++) {
	        char_flag("-vs",cvs_type);
                int_flag("-vs_a",a);
                int_flag("-vs_b",b);
    		int_flag("-vs_shape",vs_shape);
        }
    
    if(strcmp("tornado",cvs_type) == 0) {
        vs_type = VS_TORNADO;
        cout<<"View sharing: tornado"<<endl;
    } else if(strcmp("sliding",cvs_type) == 0) {
        vs_type = VS_SLIDING;
        cout<<"View sharing: sliding"<<endl;
    } else {
        vs_type = VS_NONE;
        cout<<"View sharing: none"<<endl;
    }

        
}

void VIEWSHARING::help_message() {
        cout << "----------------------------------------------" << endl;
        cout << "   View Sharing Control " << endl;
        cout << "----------------------------------------------" << endl;
        cout << "Usage:" << endl;
        help_flag("-vs []","view sharing method none/sliding/tornado");
        cout << "Filter parameters for tornado:" << endl;
        help_flag("-vs_a []","width in the center of k-space in frames");
        help_flag("-vs_b []","width in the periphery of k-space in frames");
        help_flag("-vs_shape []","shape parameter (default = 3)");
        cout << "Filter parameters for sliding window:" << endl;
        help_flag("-vs_a []","width in frames");
        help_flag("-vs_b []","0 - centered, -1 - backward, 1 - forward sliding");
}

// Assign chosen view sharing filter and create mask
void VIEWSHARING::createmask(Array<float,3>&Tw,Array<float,3>&times,int t,int frames) { // sjk 121213 recon now needs to know frames
                
        sk = Tw.length(2);
        sj = Tw.length(1);
        si = Tw.length(0);
                
        switch(vs_type) {
                
                case VS_TORNADO:{
                  tornado(Tw,times,t,frames); // sjk 121213 recon now needs to know frames
                }break;
        case VS_SLIDING: {
            slidingwindow(Tw,times,t);            
        }break;
                case VS_NONE:
                default:{
                        // Put zeros to other frames than current frame t
                        for(int k=0;k<sk;k++) {
                                for(int j=0;j<sj;j++) {
                                        for(int i=0;i<si;i++) {
                                                idif = (int)floor(times(i,j,k)) - t;
                                                if(idif) {
                                                        Tw(i,j,k)=0.0;
                                                }
                                        }
                                }
                        }
                }
        }
        
}

/* Tornado-like filter in k-space in rcframe units
  ________________________________________________
              /a\
            /     \
          /   | |   \
        /     | |     \
     /    c   | |   c    \
    <----------b--------->
*/
void VIEWSHARING::tornado(Array<float,3>&Tw,Array<float,3>&times,int t,int frames) {
        
    //int til,c; sjk 121214 c needs to be float
    int til; // sjk
    float c; // sjk

    int zcntl = 0; // sjk 121213 to count how many times Tw is not zeroed for each frame (left)
    int zcntr = 0; // on the right
    int zcnto = 0; // outside fc

    if(a==b) {b++;}
    else if(a>b) {c=(float)b;b=a;a=(int)c;} // swap incorrect ranges
    if(b>frames) {b = frames;} // sjk 121213 ensure upper range of krad tracks correctly
    c = (float)(b-a)/2; // sjk make sure float takes
    //cout<<"c = "<<c<<endl; //sjk
    float fa,fb,fc;
    fa=(float)a;
    fb=(float)b;
    fc=(float)c;
    
    /**** sjk 121213 *****/
    float cf1 = fa/2; // center of 1st frame
    float cff = frames-fa/2; //center of last frame
    float dc = (cff-cf1) / (float)(frames-1); // spacing between frames
    float cc = cf1 + dc*(float)t; // center of current frame
  
    float fcl = fc; // if upper krad reaches of tornado is cut off, make up for it on the other side
    float fcr = fc;
    if (cc-fa/2-fc <0) {
      fcl -= fabs(cc-fa/2-fc);
      fcr += fabs(cc-fa/2-fc);
    } 
    if (cc+fa/2+fc > frames) {
      fcr -= fabs(cc+fa/2+fc - frames);
      fcl += fabs(cc+fa/2+fc - frames);
    }
    //    cout<<"cc= "<<cc<<endl; //output for debug
    //cout<<"dc= "<<dc<<endl;
    //cout<<"fcl= "<<fcl<<endl;
    //cout<<"fcr= "<<fcr<<endl;
    /**** end sjk ****/

    for(int k=0;k<sk;k++) {
        for(int j=0;j<sj;j++) {
            
          fdif = times(0,j,k) - cc; // sjk 121213 difference is now measured from the center of the frame, not the left
            
          if((fabs(fdif)<(fcr+fa/2)) && (fdif)>fa/2) { // sjk 121213 divided fa by 2 because we're checking from the center to the right
                til=(int)floor(si-si*pow((fcr-fdif+fa/2)/fcr,vs_shape)); // sjk 121213 divided fa by 2 because we're checking from the center out
                if(til<0 || til>si) til=si;
                for(int i=0;i<til;i++) {        // for [0, k] spokes, could be generalized to [-k,k] trajectory
                    Tw(i,j,k)=0.0;
                    zcntr++; // sjk
                }
                
           } else if((fabs(fdif)<fcl+fa/2) && (fdif<-fa/2)) { // sjk 121213 same as above, but checking left of cc now.
                til=(int)floor(si-si*pow((fcl+fdif+fa/2)/fcl,vs_shape)); // sjk 121213 divided fa by 2 because we're checking from the center out
                if(til<0 || til>si) til=si;
                for(int i=0;i<til;i++) {
                    Tw(i,j,k)=0.0;
                    zcntl++; // sjk
                    } 
          } else if(fabs(fdif)<fa/2) { // divided fa by 2 because we're checking from center out.  
                    continue;
          } else {
                for(int i=0;i<si;i++) {
                    Tw(i,j,k)=0.0;
                    zcnto++; // sjk
                }
          } // close tornado if
        } // close j
    } // close k

    //cout<<"times zeroed = "<<endl; // output for debug sjk 
    //cout<<"left:    "<<zcntl<<endl;
    //cout<<"right:   "<<zcntr<<endl;
    //cout<<"outside: "<<zcnto<<endl;
    //cout<<"total:   "<<zcntl+zcntr+zcnto<<endl;
}

// Sliding window thru k-space containing 'a' frames
// For sliding window set zeroes when:
// abs(idif)>a          b=0 centered sliding window
// idif>0 || idif<a     b=-1  backward sliding window
// idif<0 || idif>a     b=1 forward sliding window
void VIEWSHARING::slidingwindow(Array<float,3>&Tw,Array<float,3>&times,int t) {

    if(b>1) b=1;
    if(b<-1) b=-1;
    
    for(int k=0;k<sk;k++) {
        for(int j=0;j<sj;j++) {
            for(int i=0;i<si;i++) {
                idif = (int)floor(times(i,j,k)) - t;
                switch(b) {
                    case -1:{
                        if((idif>0) || (idif<a)) {
                            Tw(i,j,k)=0.0;
                        }
                    }break;
                    case 1:{
                        if((idif<0) || (idif>a)) {
                            Tw(i,j,k)=0.0;
                        }
                    }break;
                    case 0:
                    default:{
                        if(abs(idif)>a) {
                            Tw(i,j,k)=0.0;
                        }
                    }
                }
            }
        }
    }
    
}
