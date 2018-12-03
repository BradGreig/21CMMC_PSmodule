#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <complex.h>
#include <fftw3.h>

// A module for constructing the 21cm power spectrum using C code directly accessible from Python
// Written by: Brad Greig (3rd December 2018)

#define PI (double) (3.14159265358979323846264338327)
#define TWOPI (double) (2.0*PI)

// Taken directly from 21cmFAST. For indexing boxes etc.

/* INDEXING MACROS */
#define DIM_C_INDEX(x,y,z)((unsigned long long)((z)+(DIM_MID+1llu)*((y)+DIM_BOX*(x))))// for 3D complex array
#define DIM_R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(DIM_MID+1llu)*((y)+DIM_BOX*(x)))) // for 3D real array with the FFT padding
#define DIM_R_INDEX(x,y,z)((unsigned long long)((z)+DIM_BOX*((y)+DIM_BOX*(x)))) // for 3D real array with no padding

// For initialising/constructing the power spectrum data
int NUM_BINS;

double *p_box, *k_ave;
unsigned long long *in_bin_ct;

float k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor, Delta_k, Volume;

void init_21cmPS_arrays(int DIM_BOX, float L_BOX);

// Data structure to be returned to Python using Ctypes
struct ReturnData{
    int PSbins;
    float *PSData_k;
    float *PSData;
    float *PSData_error;
};

// Data passed from Python to C
struct ReturnData Compute_21cmPS(int LightCone, int DIMENSIONAL_T_POWER_SPEC, int DIM_BOX, float L_BOX, float *box) {
    
    // (int) Lightcone: Whether the cubic box is co-eval (0) or light-cone (1)
    // (int) DIMENSIONAL_T_POWER_SPEC: Whether the 21cm power spectrum is dimensionless (0) or dimensional (i.e. mK^2) (1)
    // (int) DIM_BOX: Number of voxels per side length
    // (float) L_BOX: Length of box in Mpc (comoving)
    // (float) *box: Cubic box of data from which the 21cm PS is to be constructed
    
    struct ReturnData DataToBeReturned;
    
    unsigned long long DIM_TOT_NUM_PIXELS, DIM_TOT_FFT_NUM_PIXELS, DIM_KSPACE_NUM_PIXELS, DIM_MID;
    
    DIM_MID = DIM_BOX/2;
    
    DIM_TOT_NUM_PIXELS = (unsigned long long)(DIM_BOX*DIM_BOX*DIM_BOX);
    DIM_TOT_FFT_NUM_PIXELS = ((unsigned long long)(DIM_BOX*DIM_BOX*2llu*(DIM_MID+1llu)));
    DIM_KSPACE_NUM_PIXELS = ((unsigned long long)(DIM_BOX*DIM_BOX*(DIM_MID+1llu)));
    
    int i,j,k;
    
    double ave;
    int n_x, n_y, n_z, counter;
    float k_x, k_y, k_z, k_mag;
    unsigned long long ct;
    
    float *delta_box = (float *) malloc(sizeof(float)*DIM_TOT_NUM_PIXELS);
    
    // Average up the box to be used to construct the power spectrum
    ave = 0.0;
    for (i=0; i<DIM_BOX; i++){
        for (j=0; j<DIM_BOX; j++){
            for (k=0; k<DIM_BOX; k++){
                ave += box[DIM_R_INDEX(i,j,k)];
            }
        }
    }
    ave /= DIM_TOT_NUM_PIXELS;
    
    // Initialise the power spectrum arrays
    init_21cmPS_arrays(DIM_BOX,L_BOX);
    
    DataToBeReturned.PSData_k = calloc(NUM_BINS,sizeof(float));
    DataToBeReturned.PSData = calloc(NUM_BINS,sizeof(float));
    DataToBeReturned.PSData_error = calloc(NUM_BINS,sizeof(float));
    
    // Initialise FFT arrays
    fftwf_complex *deldel_T;
    fftwf_plan plan;
    
    deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*DIM_KSPACE_NUM_PIXELS);

    // Fill FFT-array with the fluctuation data
    for (i=0; i<DIM_BOX; i++){
        for (j=0; j<DIM_BOX; j++){
            for (k=0; k<DIM_BOX; k++){
                *((float *)deldel_T + DIM_R_FFT_INDEX(i,j,k)) = (box[DIM_R_INDEX(i,j,k)]/ave - 1)*Volume/(DIM_TOT_NUM_PIXELS+0.0);
                if (DIMENSIONAL_T_POWER_SPEC){
                    *((float *)deldel_T + DIM_R_FFT_INDEX(i,j,k)) *= ave;
                }
                // Note: we include the V/N factor for the scaling after the fft
            }
        }
    }

    // Perform FFT
    plan = fftwf_plan_dft_r2c_3d(DIM_BOX, DIM_BOX, DIM_BOX, (float *)deldel_T, (fftwf_complex *)deldel_T, FFTW_ESTIMATE);
    fftwf_execute(plan);
    
    // Compute 21cm power spectrum
    // Logic for light-cone or co-eval box
    if(LightCone) {
        
        for (n_x=0; n_x<DIM_BOX; n_x++){
            if (n_x>DIM_MID)
                k_x =(n_x-(int)DIM_BOX) * Delta_k;  // wrap around for FFT convention
            else
                k_x = n_x * Delta_k;
            
            for (n_y=0; n_y<DIM_BOX; n_y++){
                // avoid the k(k_x = 0, k_y = 0, k_z) modes (pure line-of-sight modes)
                if(n_x != 0 && n_y != 0) {

                    if (n_y>DIM_MID) {
                        k_y =(n_y-(int)DIM_BOX) * Delta_k;
                    }
                    else
                        k_y = n_y * Delta_k;
                
                    for (n_z=0; n_z<=DIM_MID; n_z++){
                        k_z = n_z * Delta_k;
                    
                        k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
                    
                        // now go through the k bins and update
                        ct = 0;
                        k_floor = 0;
                        k_ceil = k_first_bin_ceil;
                        while (k_ceil < k_max){
                            // check if we fal in this bin
                            if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                                in_bin_ct[ct]++;
                                p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T[DIM_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*Volume);
                                // note the 1/VOLUME factor, which turns this into a power density in k-space
                            
                                k_ave[ct] += k_mag;
                                break;
                            }
                        
                            ct++;
                            k_floor=k_ceil;
                            k_ceil*=k_factor;
                        }
                    }
                }
            }
        }
    }
    else {
        for (n_x=0; n_x<DIM_BOX; n_x++){
            if (n_x>DIM_MID)
                k_x =(n_x-(int)DIM_BOX) * Delta_k;  // wrap around for FFT convention
            else
                k_x = n_x * Delta_k;
        
            for (n_y=0; n_y<DIM_BOX; n_y++){
                if (n_y>DIM_MID) {
                    k_y =(n_y-(int)DIM_BOX) * Delta_k;
                }
                else
                    k_y = n_y * (TWOPI/L_BOX);
            
                for (n_z=0; n_z<=DIM_MID; n_z++){
                    k_z = n_z * Delta_k;
                
                    k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
                
                    // now go through the k bins and update
                    ct = 0;
                    k_floor = 0;
                    k_ceil = k_first_bin_ceil;
                    while (k_ceil < k_max){
                        // check if we fal in this bin
                        if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                            in_bin_ct[ct]++;
                            p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T[DIM_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*Volume);
                            // note the 1/VOLUME factor, which turns this into a power density in k-space
                        
                            k_ave[ct] += k_mag;
                            break;
                        }
                    
                        ct++;
                        k_floor=k_ceil;
                        k_ceil*=k_factor;
                    }
                }
            }
        }
    }
    
    counter = 0;
    
    // Bin the data
    for (ct=1; ct<NUM_BINS; ct++){
        if (in_bin_ct[ct]>0) {
            DataToBeReturned.PSData_k[counter] = k_ave[ct]/(in_bin_ct[ct]+0.0);
            DataToBeReturned.PSData[counter] = p_box[ct]/(in_bin_ct[ct]+0.0);
            DataToBeReturned.PSData_error[counter] = p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0);
            
            counter += 1;
        }
    }

    DataToBeReturned.PSbins = counter;
    
    // return the data
    return DataToBeReturned;
}

void init_21cmPS_arrays(int DIM_BOX, float L_BOX) {
    
    // initialise binning of the power spectrum
    
    Delta_k = (TWOPI/L_BOX);
    Volume = L_BOX*L_BOX*L_BOX;
    
    k_factor = 1.35;
    k_first_bin_ceil = Delta_k;
    k_max = Delta_k*(int)DIM_BOX;

    // initialize arrays
    // ghetto counting (lookup how to do logs of arbitrary bases in c...)
    NUM_BINS = 0;
    k_floor = 0;
    k_ceil = k_first_bin_ceil;
    while (k_ceil < k_max){
        NUM_BINS++;
        k_floor=k_ceil;
        k_ceil*=k_factor;
    }
    
    p_box = calloc(NUM_BINS,sizeof(double));
    k_ave = calloc(NUM_BINS,sizeof(double));
    in_bin_ct = (unsigned long long *)calloc(NUM_BINS,sizeof(unsigned long long));
}


