#include <stdio.h>
#include <math.h>
#include <string.h>

#define CROSS_VBINS 1000000
#define VMAX 1.0E5
#define VMIN 1.0E-2
#define LIGHTSPEED 299972.0 // in km/s
#define CM_CONVERT 100000 // converts km/s to cm/s
#define V_NORM 10000000 // normalizes velocity
#define MAX_FNAME_LENGTH 30
#define NO_ELASTIC_FNAMES 4
#define NO_INELASTIC_FNAMES 8


/*
 * This code generates 2cDM cross section tables for
 * user-provided sigma0, scattering power law,
 * and conversion power law.
 *
 */


int main(){

  printf("Welcome!\n");

  /* Default values */
  double sigma0 = 0.1;
  int ps = 0; // power for scattering
  int pc = 0; // power for conversion; momentum ratios subtract 1 from this

  double Dvlog = log( VMAX/VMIN ) / CROSS_VBINS;
  
  char elastic_filenames[NO_ELASTIC_FNAMES][MAX_FNAME_LENGTH] = {
    "sidm_cross_reaction_0.txt",
    "sidm_cross_reaction_4.txt",
    "sidm_cross_reaction_7.txt",
    "sidm_cross_reaction_11.txt"
  };

  char inelastic_filenames[NO_INELASTIC_FNAMES][MAX_FNAME_LENGTH] = {
    "sidm_cross_reaction_1.txt",
    "sidm_cross_reaction_2.txt",
    "sidm_cross_reaction_3.txt",
    "sidm_cross_reaction_5.txt",
    "sidm_cross_reaction_6.txt",
    "sidm_cross_reaction_8.txt",
    "sidm_cross_reaction_9.txt",
    "sidm_cross_reaction_10.txt"
  };



  printf("This will create a cross section table for sigma: %.2f \n", sigma0);
  printf("and power laws of (%d, %d - 1).\n", ps, pc);
  
  /* ELASTIC SCATTERING BLOCK */

  static double velocity_array[CROSS_VBINS];
  static double cross_section_array[CROSS_VBINS];

  double velocity = 0;
  
  for( int i = 0; i < CROSS_VBINS; i++ ) {
    
    velocity = exp( Dvlog * (i + 0.5) + log(VMIN) ); // is in km/s

    double vel_cm = velocity * CM_CONVERT; // converts to cm/s

    double norm_vel = vel_cm / V_NORM; // normalizes velocity

    double cross_sec = pow( norm_vel, ps ) * sigma0; // cross section at this velocity

    /* Storing the result */

    velocity_array[i] = velocity;

    cross_section_array[i] = cross_sec;

  }

  
  for( int i = 0; i < NO_ELASTIC_FNAMES; i++) {

    FILE *fp;

    fp = fopen( elastic_filenames[i], "w+");

    for( int j = 0; j < CROSS_VBINS; j++ ) {

      fprintf(fp, "%.18e\t", velocity_array[j]);
      fprintf(fp, "%.18e\n", cross_section_array[j]);
    }

    fclose(fp);
  }


  /* INELASTIC SCATTERING BLOCK */

  velocity = 0;

  for( int i = 0; i < CROSS_VBINS; i++ ) {
    
    velocity = exp( Dvlog * (i + 0.5) + log(VMIN) ); // is in km/s

    double vel_cm = velocity * CM_CONVERT; // converts to cm/s

    double norm_vel = vel_cm / V_NORM; // normalizes velocity

    double cross_sec = 0.5 * pow( norm_vel, pc - 1 ) * sigma0; // cross section at this velocity

    /* Storing the result */

    velocity_array[i] = velocity;

    cross_section_array[i] = cross_sec;

  }

  
  for( int i = 0; i < NO_INELASTIC_FNAMES; i++) {

    FILE *fp;

    fp = fopen( inelastic_filenames[i], "w+");

    for( int j = 0; j < CROSS_VBINS; j++ ) {

      fprintf(fp, "%.18e\t", velocity_array[j]);
      fprintf(fp, "%.18e\n", cross_section_array[j]);
    }

    fclose(fp);
  }

  printf("Finished!\n");

  return 0;
}

