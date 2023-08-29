#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h> //PATH_MAX

#define CROSS_VBINS 1000000
#define VMAX 1.0E5
#define VMIN 1.0E-2
#define LIGHTSPEED 299972.0 // in km/s
#define CM_CONVERT 100000 // converts km/s to cm/s
#define V_NORM 10000000 // normalizes velocity
#define MAX_FNAME_LENGTH 30
#define NO_ELASTIC_FNAMES 4
#define NO_INELASTIC_FNAMES 8
#define M2 -2 // For (-2, -2) case


/*
 * This code generates 2cDM cross section tables for
 * user-provided sigma0, scattering power law,
 * and conversion power law.
 *
 */

/* Got this from Stackexchange */
int build_path(char *dest, size_t size, char *subdir, char *filename)
{
    // Don't assume the input is valid
    if ( (dest == NULL) || (size < 1) || (filename == NULL) ) {
        return 0;
    }

    // Make no subdir work  (sane default behavior)
    if ( subdir == NULL ) {
        subdir = "";
    }

    // Prevent buffer overruns by using a safe formatting function
    size_t stored = snprintf(dest, size, "/%s/%s", subdir, filename);

    // Returns success (true) if the path fit the buffer
    return ((stored+1) <= size);
}


int main(int argc, char* argv[]){

  printf("Welcome!\n");
  
  int printhelp = 0;
  char symlinkpath[PATH_MAX] = ".";
  char actualpath [PATH_MAX];
  char *ptr;
  char *sptr;
  
  /* Default values */
  double sigma0 = 0.1;
  int ps = 0; // power for scattering
  int pc = 0; // power for conversion; momentum ratios subtract 1 from this

  int pst = 0; // input variables, will copy them later
  int pct = 0;

  for( int i = 0; i < argc; i++ ){

    if (strncmp(argv[i], "--help", 6) == 0){
      printhelp = 1;
    }
    if (strncmp(argv[i], "-h", 2) == 0){
      printhelp = 1;
    }
    if (strncmp(argv[i], "-o", 2) == 0){
      i++;
      strcpy(symlinkpath, argv[i]);
    }
    if (strncmp(argv[i], "--sigma", 7) == 0){
      i++;
      sigma0 = strtod(argv[i], &sptr);
    }
    if (strncmp(argv[i], "-ps", 3) == 0) {
      i++;
      pst = strtol(argv[i], &sptr, 10);
    }
    if (strncmp(argv[i], "-pc", 3) == 0) {
      i++;
      pct = strtol(argv[i], &sptr, 10);
    }

  }

  if(printhelp){

    printf("Usage: %s [options]\n", argv[0]);
    printf("  options:\n");
    printf("\t --help (-h)\tprint options\n");
    printf("\t --sigma\tset sigma0\n");
    printf("\t -ps\tset scattering power law (default 0)\n");
    printf("\t -pc\tset conversion power law (default 0)\n");
    printf("\t -o\tset output path\n");
    return 0;
  }

  printf("This will create a cross section table for sigma: %.2f. \n", sigma0);

  if ( (pst == pct) && (pct == M2) ){
    printf("We are in the power law (%d, %d) case. Setting pc to (%d) exactly.\n", pst, pct, M2);
    
    // Doesn't hurt to be too careful.
    ps = M2;
    pc = M2;

    printf("Power laws are (%d, %d).\n", pst, pct);
  }
  else {
    //Do the minus 1 here
    pc = pct-1;
    printf("Power laws are (%d, %d - 1).\n", pst, pct);
  }


  ptr = realpath(symlinkpath, actualpath);
  
  printf("Cross section files will be output to:\n");
  printf("%s\n", actualpath);
  
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

  char fpath [PATH_MAX];

  for( int i = 0; i < NO_ELASTIC_FNAMES; i++) {

    if( ! build_path( fpath, sizeof(fpath), actualpath, elastic_filenames[i] ) ){
      fprintf( stderr, "ERROR: Invalid Path\n");
      return 1;
    }


    FILE *fp;

    fp = fopen( fpath, "w+");

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

    double cross_sec = 0.5 * pow( norm_vel, pc ) * sigma0; // cross section at this velocity

    /* Storing the result */

    velocity_array[i] = velocity;

    cross_section_array[i] = cross_sec;

  }

  
  for( int i = 0; i < NO_INELASTIC_FNAMES; i++) {

    if( ! build_path( fpath, sizeof(fpath), actualpath, inelastic_filenames[i] ) ){
      fprintf( stderr, "ERROR: Invalid Path\n");
      return 1;
    }

    FILE *fp;

    fp = fopen( fpath, "w+");

    for( int j = 0; j < CROSS_VBINS; j++ ) {

      fprintf(fp, "%.18e\t", velocity_array[j]);
      fprintf(fp, "%.18e\n", cross_section_array[j]);
    }

    fclose(fp);
  }

  printf("Finished!\n");

  return 0;
}

