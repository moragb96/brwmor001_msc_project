/******************************************************************************
 calc_r2c_fftw.c
 Author: Morag Brown
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main(int argc, char *argv[])
{
// argv[1] is input file name 
// argv[2] is output file name
// argv[3] is complex output file name
size_t NFFT = atoi(argv[4]);     // FFT length
size_t N = NFFT;                 // input sig == FFT length
size_t i=0;
char line[128];

if (argc!=5){
        printf("usage: %s <input file> <output file> <complex output file>  <FFT len> \n", argv[0]);
        return 0;
}

 FILE * file;
 FILE * file1;

 double *in;
 double complex z;
 double c_out;

 fftw_complex *out;
 fftw_plan p;
 in = (double*) malloc(sizeof(double)*N);
 out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(N/2 +1));

 p = fftw_plan_dft_r2c_1d(NFFT, in, out, FFTW_ESTIMATE);

 file = fopen(argv[1],"r");
 if(file == NULL) {
     perror("Error opening file");
     return(-1);
   }

int count = 0;

while( fgets (line, sizeof(line), file)) {
    if (count < N){
        in[i] = atof(line);
        i++;
    }
    count++;
}

 fclose(file);

 fftw_execute(p);
 
 
 file = fopen(argv[2],"w");
 file1 = fopen(argv[3],"w");
 for (i = 0; i < N/2; i++){
         z = creal(out[i]) + I*cimag(out[i]);
         c_out = cabs(z);
         fprintf(file, "%lf\n", c_out);
         fprintf(file1, "%lf + %lfj\n", creal(out[i]), cimag(out[i]));
 }

 fclose(file);
 fclose(file1);

 fftw_destroy_plan(p);
 fftw_free(out);

 return 0;

 }

