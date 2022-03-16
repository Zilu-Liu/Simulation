//QBL_V5: Normal priors for PCs. Number of PCs is arbitrary.

//To create dynamic library .so file:
// on Mac: PKG_LIBS="-lgsl -lblas" R CMD SHLIB QBL_V5.c
// on Linux: PKG_CPPFLAGS="-I/usr/include" R CMD SHLIB QBL_V5.c -lgsl -lgslcblas -lm

// The QBL_V5.so file provided in this folder was compiled on Mac, it may not work on other platforms.

//prior for beta0 is N(0, 1)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "libnew.h"
#include "R.h"
#include <Rmath.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

//function gsl_ran_multivariate_gaussian is in gsl library
int gsl_ran_multivariate_gaussian (const gsl_rng * r, const gsl_vector * mu, const gsl_matrix * L, gsl_vector * result);

double GXE_calc_a(double *freq, int *per_freq, double D);

void GXE_update_z(int **per_map, int *per_hz, double *beta, double *freq, double D, double sigmas,double per_y_new, double *per_xz, double **per_xhaplo, int num_haplo_per, int x_length, double *R_i, double *gamma, int q);

double GXE_update_lambda(double *beta, double a, double b, int x_length);

double GXE_update_sigma(double *beta, double asig, double bsig,double **xz, int N, int q, double *GXE_y_new,int x_length, double **R, double *gamma);

//double randn (double mu, double sigma);

//double GXE_calc_den_post(double *beta, double *freq, int h_length, double D, int tot_hap, int **xz, int n_cov, int N);

void GXE_update_beta(double sigma,double **xz, double *beta, int which_beta, double lambda, double *freq, double D, int x_length, int h_length, double *GXE_y_new, int N, int q, int n_cov, double **R, double *gamma);
double GXE_update_gamma_quadraticPart(double *gamma, int q);
double GXE_update_randomAsigma(double asigRandomA, double bsigRandomA, int q, double quadratic_sum);
double GXE_gen_double_exp(double mean, double SD);

void GXE_update_freq(double *freq, double *beta, double D, int h_length, double **xz, int N, int **hz, int n_cov, int Con);

void GXE_update_gamma(double *gamma, double *mean_vector, double **cov_matrix, double **R, double sigmas, double randomA_sigmas, double *y, int N, int q, double *beta, double **xz, int x_length);//new function

double GXE_update_D(double *freq, double *beta, double D, int h_length, double **xz, int N, int **hz, int n_cov);

void GXE_dirichlet(double *param, int dim, double *gen_sample);

double GXE_sum(double *x, int n);

double GXE_find_min(double *arr, int n);

int accept_count=0;

int **GXE_uniq_map, GXE_total_freq_update=0; //define unique matrix x, response y
double *GXE_y_new;

void GXE_mcmc(double *x, int *tot_hap, double *y,
int *N, int *q, int *num_haplo_id, int *x_length, int *h_length,
int *freq_num, double *freq, double *D,
double *sigmas, double *beta, double *gamma, double *a, double *b,
double *asig, double *bsig, double *lambda,
int *NUM_IT, int *BURN_IN, double *beta_out, double *gamma_out, double *lambda_out,
double *sigmas_out, double *randomA_sigmas_out, double *freq_out, double *D_out, int *n_cov, double *up_xz, int *CC, double *BDIC_out, double *randomA_sigmas, double *randomA_sigmas_a, double *randomA_sigmas_b, double *R_vec)

/*
 *x:        design matrix as a vector(col1',col2'...)

 *tot_hap:  # of rows in design matrix

 *y:        case/control(affected)

 *N:        sum of weight(total # of subject)
 
 *q:        number of PCs
 
 *R:        n by q matrix, columns as pcs

 *num_haplo_id: the # of rows for each original ID

 *x_length: # of hyplotypes (m)

 *freq:     freq vector, the freq of the most common hyplotype is the last

 *D:        D=0

 *sigmas:   sigmas=1

 *beta:     all beta=start #(0.01)
 
 *gamma:    parameter for PCs with Normal priors

 *a:        gamma prior of lamda a=20

 *b:        gamma prior of lamda b=20

 *asig:     inv-gamma prior of sigma squared asig=

 *bsig:     inv-gamma prior of sigma squared bsig=
 
 *randomA_sigmas_a: inv-gamma prior for randomA_sigmas
 
 *randomA_sigmas_b: inv-gamma prior for randomA_sigmas

 *uniq_x:   get all unique rows of design matrix, then do the transpose, then                   read it a vector(row1, row2,...)

 *tot_uniq_x: # of unique rows

 *lambda:   initial of lamda(prior of beta) =1

 *NUM_IT:   # of iteration = 50000

 *BURN_IN:  start burning iteration

 *beta_out: # of sampled beta's

 *lambda_out: # num of sampled lambda's(prior of beta)

 *freq_out: #num of sampled f's

 *D_out:    # number of 'd'
 
 *Rgamma:   the product of R and gamma for ith individual
 */

{


printf("\n Check values of arguments before put into MCMC but after put into C program\n");
printf("CC is %d\n",*CC);

/*
printf("\n Check values of arguments before put into MCMC but after put into C program\n");

  printf("First x is %lf\n",x[1]);
  printf("Number of people in design matrix x is %d\n",*tot_hap);
  printf("first several responses %lf\n",y[1]);
  printf("Total number of people in original design matrix is %d\n",*N);
  printf("Total number of people in original design matrix is %d\n",*num_haplo_id);
  printf("x.length is equal to %d\n",*x_length);
  printf("total number of possible haplotypes is %d\n",*h_length);
  printf("Head of freq.num is %d\n",*freq_num);
  printf("freq.new is %lf\n",freq[0]);

  printf("The initial value of D is %lf\n",*D);
  printf("The initial value of sigmas is %lf\n",*sigmas);
  printf("The initial value of beta is %lf\n",*beta);
  printf("The initial value of a is %lf\n",*a);
  printf("The initial value of b is %lf\n",*b);
  printf("The initial value of asig is %lf\n",*asig);
  printf("The initial value of bsig is %lf\n",*bsig);
  printf("Value of num.it is %d\n",*NUM_IT);
  printf("Value of burn.in is %d\n", *BURN_IN);
  printf("number of number of covariate is %d\n", *n_cov);


  printf("Initial value of lambda is %lf\n",*lambda);

  printf("Initial values beat is as below\n");
  for (int k=0;k<*x_length+1;k++){
  printf("%d th initial beta is %lf\n",k, beta[k]);
  }


  printf("First lambda.out is %lf\n", lambda_out[0]);
  printf("First sigmas.out is %lf\n", sigmas_out[0]);


  printf("Initial freq.out is as below\n");

  for (int k=0;k<*h_length;k++){
  printf("%d th Initial freq.out is %lf\n",k, freq_out[k]);
  }

  printf("Initial D.out is %lf\n", D_out[0]);

printf("End of Check values of arguments before put into MCMC but after put into C program\n\n");

*/

//printf("%lf\n",x[1]);
//printf("%lf\n",x[0]);
//printf("Total number of haplotype%d\n",*tot_hap);




    int i, j, k, l, m, n,  ***freq_map, which_beta, **hz, h_mat[*tot_hap][2], it=0, it1=0, it2=0, it3=0, heat_it;
    double  **xz, x_mat[*tot_hap][*x_length],***xhaplo,**R,*mean_vector, **cov_matrix;

    //allocate memory to pointers!!!
    //New: mean_vector and cov_matrix for gamma
    mean_vector=calloc(*q, sizeof(double));
    cov_matrix=calloc(*q, sizeof(double));
    
    for (i =0; i<*q; ++i)
    {
        cov_matrix[i]=calloc(*q, sizeof(double));
    }
    
    

    GXE_y_new = calloc(*N, sizeof(double));
    xhaplo = calloc(*N, sizeof(double));
    R=calloc(*N, sizeof(double));
    freq_map = calloc(*N, sizeof(int*));


    for (i =0; i<*N; ++i)

    {

        freq_map[i] = calloc(num_haplo_id[i], sizeof(int*));
        R[i]=calloc(*q, sizeof(double)); //*N not num_haplo_id[i]
        xhaplo[i] = calloc(num_haplo_id[i], sizeof(double));

    }

    for (i =0; i<*N; ++i)

    {

        for (j = 0; j < num_haplo_id[i]; ++j)

        {

            xhaplo[i][j] = calloc(*x_length, sizeof(double));

            freq_map[i][j] = calloc(2, sizeof(int));

        }

    }

    /* separating x vector from R as x matrix */

    l=0;

    for (j=0; j<*x_length; ++j)

    {

        for (i=0; i<*tot_hap; ++i)

        {

            x_mat[i][j] = x[l]; //ZL: x_mat row is design matrix i, column is design matrix j

            ++l;

        }

    }

    //printf("Kinship matrix first entry: %g\n",kin_mat_inv[0][0]);


    /* separating haplo map vector from R as h matrix */

    l=0;

    for (j=0; j<2; ++j)

    {

        for (i=0; i<*tot_hap; ++i)

        {

            h_mat[i][j] = freq_num[l];  //ZL: h_mat row is 6 haps, 2 columns

            ++l;

        }

    }



    xz = calloc(*N, sizeof(double));

    hz = calloc(*N, sizeof(int*));

    for (i=0; i<*N; ++i)

    {

        xz[i] = calloc(*x_length, sizeof(double));

        hz[i] = calloc(2, sizeof(int));

    }
    
    l=0;
    for (i=0;i<*N;++i)
    {
    	for (j=0;j<*q;++j)
    	{
    		R[i][j]=R_vec[l];
    		++l;
    	}
    }
    //printf("%lf\n",x[0]);
    printf("R matrix first entry %lf\n",R[0][0]);
    printf("R matrix second entry %lf\n",R[0][1]);
    // x_mat is the design matrix, xz is the matrix with orignal # of subject.
    /* separating haplotypes per person from the x matrix and assigning z (missing haplotype) for persons with only one compatible haplotype */

    l = 0;
    for (i =0; i<*N; ++i)
    {
        GXE_y_new[i] = y[l];
       // n_cases += GXE_y_new[i];
        for (j = 0; j < num_haplo_id[i]; ++j)
        {
            m = 0;
            for (k = 0; k < *x_length; ++k)
            {
                xhaplo[i][j][k] = x_mat[l][m];
                if (num_haplo_id[i]==1)
                    xz[i][k] = x_mat[l][m];
                else
                    xz[i][k] = 0;
                ++m;
            }
            ++l;
        }
    }
    l = 0;
    for (i =0; i<*N; ++i)
    {
        for (j = 0; j < num_haplo_id[i]; ++j)
        {
            m = 0;
            for (k = 0; k < 2; ++k)
            {
                freq_map[i][j][k] = h_mat[l][m];
                if (num_haplo_id[i]==1)
                    hz[i][k] = h_mat[l][m];
                else
                    hz[i][k] = 0;
                ++m;
            }
            ++l;
        }
    }
//ZL: xz,hz for individuals whose haps are fixed
    /* separating unique x vector (haplotypes) from R as unique x matrix (to be used in denominator calculation) */






    /********************** start MCMC here ***************************/

    heat_it = 0; /* tracks # of heating iterations */

      // printf("Begin Design matrix for beta1=%lf\n",xz[1][0]);
      //	printf("Begin Design matrix for beta2=%lf\n",xz[1][1]);
	//printf("Begin Design matrix for beta3=%lf\n",xz[1][2]);
	//printf("Begin Design matrix for beta4=%lf\n",xz[1][3]);
	//printf("Begin Design matrix for beta5=%lf\n",xz[1][4]);
    //   printf("Begin Updated beta0=%lf\n",beta[0]);
     //  printf("Begin Updated beta1=%lf\n",beta[1]);
//	printf("Begin Updated beta2=%lf\n",beta[2]);
//	printf("Begin Updated beta3=%lf\n",beta[3]);
//	printf("Begin Updated beta4=%lf\n",beta[4]);
//	printf("Begin Updated beta5=%lf\n",beta[5]);

    for (n=0; n<*NUM_IT; ++n)

    {
    	double quadratic_sum=0;
        heat_it = heat_it+1;
        //printf("Current Iteration is %d\n",n+1);
        //printf("Random_A first entry is: %g\n",random_A[0]);
        //printf("Variance of random A is: %g\n",*randomA_sigmas);
        //printf("Variance of Error is: %g\n",*sigmas);
        if (heat_it == 101) heat_it = 1;

        int h=0;

        /* assigning or updating z (missing haplotype) for persons with more than one compatible haplotype */

        for (i =0; i<*N; ++i)

        {
            if (num_haplo_id[i]>1)
            {
			//   printf("which z update %d\n",i);
			//   for (k = 0; k < *x_length; ++k)
         //   {
         //      printf("xz before update is %lf\n",xz[i][k]);
         //   }
			  //     printf("sigmas is %lf\n",*sigmas);

                GXE_update_z(freq_map[i], hz[i], beta, freq, *D,  *sigmas, GXE_y_new[i], xz[i], xhaplo[i], num_haplo_id[i], *x_length, R[i], gamma, *q); //individual i

			// for (k = 0; k < *x_length; ++k)
           // {
          //     printf("xz after update is %lf\n",xz[i][k]);
          //  }
	        }
      // find out updated xz to test the algorithm
           for (k = 0; k < *x_length; ++k)
           {
            up_xz[h]=xz[i][k];
           //printf("%d th row, %d th column of xz is %lf\n",i+1,k+1,xz[i][k]);
          // printf("%d th up_xz is %lf\n",h+1,up_xz[h]);
            ++h;
           }
        }
        
        
        
        /* update beta parameters */
       // printf("sigmas is  =%lf\n",*sigmas);
        for (i=0; i<*x_length+1; ++i)

        {
            which_beta=i;

            GXE_update_beta(*sigmas, xz, beta, which_beta, *lambda, freq, *D, *x_length, *h_length, GXE_y_new, *N, *q, *n_cov, R, gamma);

        }
        /*
        for (i=0; i<*x_length+1; ++i)

           {
               printf("updated beta i is  =%lf\n",beta[i]);

           }
        */
        //update quadratic form
        quadratic_sum=GXE_update_gamma_quadraticPart(gamma, *q);
        //printf("Quadratic sum is  =%lf\n",quadratic_sum);
        *randomA_sigmas=GXE_update_randomAsigma(*randomA_sigmas_a, *randomA_sigmas_b, *q, quadratic_sum);
        //printf("updated randomA_sigmas is  =%lf\n",*randomA_sigmas);
        
        /* update sigma parameters */
        *sigmas = GXE_update_sigma(beta,*asig,*bsig,xz, *N, *q, GXE_y_new, *x_length, R, gamma);
        //printf("updated sigmas is  =%lf\n",*sigmas);

        /* update lambda */

        *lambda = GXE_update_lambda(beta, *a, *b, *x_length);

        /* update frequencies and D */
        GXE_update_freq(freq, beta, *D, *h_length, xz,  *N, hz,*n_cov, *CC);
        /* update D parameter */
        *D = GXE_update_D(freq, beta, *D, *h_length, xz,  *N, hz,*n_cov);
        

        //update gamma
        GXE_update_gamma(gamma,mean_vector,cov_matrix, R, *sigmas, *randomA_sigmas, GXE_y_new, *N, *q, beta, xz, *x_length);
        
        /*
           for (i=0;i<*q;i++)
           {
                printf("elements in mean vector for gamma is %lf\n",mean_vector[i]);
               for (j=0;j<*q;j++)
               {
                   printf("elements in cov matrix for gamma row i column j is %lf\n",cov_matrix[i][j]);
               }
               printf("updated gamma i is %lf\n",gamma[i]);
               
           }
        */
        // printf("Printing D: %lf", D);
        // n greater than burn in, store information in .._out
        if (n >= *BURN_IN)
        {
            for (i=0; i<*x_length+1; ++i)
            {
                beta_out[it] = beta[i];
                ++it;
            }
            for (i=0; i<*q; ++i)
            {
                gamma_out[it3] = gamma[i];
                ++it3;
            }
            lambda_out[it2] = *lambda;
            sigmas_out[it2] = *sigmas;
            randomA_sigmas_out[it2] = *randomA_sigmas;
            //printf("Updated sigma squared is %lf",*sigmas);

            for (i=0; i<*h_length; ++i)
            {
                freq_out[it1] = freq[i];
                ++it1;
            }
            D_out[it2] = *D;
			/* Calcuate BDIC*/
		 double fit=0,xxbeta[*N], resid=0;
		  for (i=0; i<*N; ++i)
		  {
			  xxbeta[i]=beta[0];
		    for (j=0; j<*x_length; ++j)
            {
                xxbeta[i]=xxbeta[i]+xz[i][j]*beta[j+1];
            }
		 resid=	GXE_y_new[i]-xxbeta[i];
         fit=fit+pow(resid,2)/(*sigmas*2);
		//printf("resid=%lf\n",resid);
        //printf("resid to the power of 2=%lf\n",pow(resid,2));
		  }
        BDIC_out[it2] = -*N*log(sqrt(2*M_PI*pow(*sigmas,2)))-fit;
	    /*printf("After xxbeta[0]=%lf\n",xxbeta[0]);
        printf("fit=%lf\n",fit);
		*/
            ++it2;
        }
        
    }
    
// end pf MCMC
    //free memory
    
       for (i=0; i<*N; ++i)

       {
           free(hz[i]);
           free(xz[i]);
       }
       
       free(xz);
       free(hz);
       
       
       
       
       for (i =0; i<*N; ++i)

       {

           for (j = 0; j < num_haplo_id[i]; ++j)

           {

               free(freq_map[i][j]);
               free(xhaplo[i][j]);

           }

       }
       
       for (i =0; i<*N; ++i)

       {

           free(xhaplo[i]);
           free(freq_map[i] );

       }
       
       
       free(freq_map);
       free(xhaplo);
       free(GXE_y_new);
    
    
}




// functions
void GXE_update_z(int **per_map, int *per_hz, double *beta, double *freq, double D,double sigmas, double per_y_new, double *per_xz, double **per_xhaplo, int num_haplo_per, int x_length,
		double *R_i, double *gamma, int q) // for one individual
{

    int i, k;

    double prob[num_haplo_per], cum_prob[num_haplo_per], sum_prob, x, a[num_haplo_per], xbeta[num_haplo_per], Rgamma;
    Rgamma=0;
    for (k=0; k < q; ++k)
    {
        Rgamma=Rgamma+R_i[k]*gamma[k];
    }
    for (i=0; i<num_haplo_per; ++i)

    {

        a[i] = GXE_calc_a(freq, per_map[i], D);
        prob[i] = a[i];
        xbeta[i]=beta[0];
        /*if (per_y_new == 1)
         * {
         */
            for (k=0; k < x_length; ++k)
            {
                xbeta[i]=xbeta[i]+per_xhaplo[i][k]*beta[k+1];
            }

        /*}
    */
        
       prob[i] = prob[i]*exp(-(per_y_new-xbeta[i]-Rgamma)*(per_y_new-xbeta[i]-Rgamma)/(2*sigmas)); //random

    }
    sum_prob = GXE_sum(prob, num_haplo_per);
    for (i=0; i<num_haplo_per; ++i)

    {

        prob[i] = prob[i]/sum_prob;

    }
    cum_prob[0]= prob[0];

    for (i=1; i<num_haplo_per; ++i)

    {

        cum_prob[i] = cum_prob[i-1]+prob[i];

    }
    GetRNGstate();

    x=runif(0,1);

    PutRNGstate();
    if (x < cum_prob[0])

    {

        for (k=0; k<x_length; ++k)

        {

            per_xz[k] = per_xhaplo[0][k];

        }

        for (k=0; k<2; ++k)

        {
            //printf("old hap 1st is %d\n",per_hz[k]);
            per_hz[k] = per_map[0][k];
            //printf("updated hap 1st is %d\n",per_hz[k]);
        }

    }
    for (i=1; i<num_haplo_per; ++i)

    {

        if (x > cum_prob[i-1] && x < cum_prob[i])

        {



            for (k=0; k<x_length; ++k)

            {

                per_xz[k] = per_xhaplo[i][k];

            }

            for (k=0; k<2; ++k)

            {
                // printf("old hap 1st is %d\n",per_hz[k]);
            	per_hz[k] = per_map[i][k];
               // printf("updated hap 1st is %d\n",per_hz[k]);
            }




        }



    }


}

void GXE_update_beta(double sigma, double **xz, double *beta, int which_beta, double lambda,
 double *freq, double D, int x_length, int h_length, double *GXE_y_new, int N, int q, int n_cov, double **R, double *gamma) // add mu_0, std_0 (for the intercept)

{

    double beta_new,g_1,g_2,f_old,f_new,beta_new_vec[x_length+1], SD, accept_prob, xbeta[N], xbetan[N], Rgamma;
    double mu_0=0, std_0;
    if (which_beta==0){ //update intercept
        
    int i, k;

    //printf("which beta to updata is %d\n",which_beta);
    //printf("old beta is %lf\n",beta[which_beta]);


     for (i=0; i<N; ++i)

     {
         Rgamma=0;
         for (k=0; k < q; ++k)
         {
             Rgamma=Rgamma+R[i][k]*gamma[k];
         }
         
		 //y_mean=y_mean+GXE_y_new[i]/N;

         xbetan[i]=GXE_y_new[i];

         for (k=0; k < x_length; ++k)

         {

            xbetan[i]=xbetan[i]-xz[i][k]*beta[k+1];

         }

         mu_0=mu_0+xbetan[i]-Rgamma;


     }
        
        
     //for (i=0; i<N; ++i)
     //{

         //y_var=y_var+(GXE_y_new[i]-y_mean)*(GXE_y_new[i]-y_mean);

     //}
        
        //y_var=y_var/(N-1);
     
        
     mu_0=(mu_0/sigma)/(1+N/sigma);  //prior for beta0 is N(0, 1)
     std_0=sqrt(1/(1+N/sigma));
     beta_new = rnorm(mu_0, std_0);

     beta[which_beta] = beta_new;
        
     /*
     printf("mean for y is %lf\n",y_mean);
     printf("var for y is %lf\n",y_var);
     printf("mean for beta_0 is %lf\n",mu_0);
     printf("std for beta_0 is %lf\n",std_0);
     printf("updated beta_0 is %lf\n",beta[which_beta]);
     */
    
    }else{
    int i, x, k;
    //printf("which beta to updata is %d\n",which_beta);
    //printf("old beta is %lf\n",beta[which_beta]);
    //printf("mean based on old beta is %lf\n",beta[which_beta]);
    //printf("std based on old beta is %lf\n",sqrt(fabs(beta[which_beta])));
    beta_new = GXE_gen_double_exp(beta[which_beta], sqrt(fabs(beta[which_beta])));
    //printf("new beta before test is %lf\n",beta_new);
    //printf("sigma squared is %lf\n",sigma);
    g_1 = -lambda*fabs(beta[which_beta]) ;
    g_2 =-lambda*fabs(beta_new);

    for (i=0; i<x_length+1; ++i){
		beta_new_vec[i] = beta[i];
		}
    beta_new_vec[which_beta] = beta_new;
    for (i=0; i<N; ++i){   //individual i
	        xbeta[i]=beta[0];
            xbetan[i]=beta_new_vec[0];
            for (k=0; k < x_length; ++k){
                xbeta[i]=xbeta[i]+xz[i][k]*beta[k+1];
                xbetan[i]=xbetan[i]+xz[i][k]*beta_new_vec[k+1];
            }
     Rgamma=0;
     for (k=0; k < q; ++k)
     {
         Rgamma=Rgamma+R[i][k]*gamma[k];
     }
     g_1 = g_1-((GXE_y_new[i]-xbeta[i]-Rgamma)*(GXE_y_new[i]-xbeta[i]-Rgamma))/(2*sigma); //random effect
     g_2 = g_2-((GXE_y_new[i]-xbetan[i]-Rgamma)*(GXE_y_new[i]-xbetan[i]-Rgamma))/(2*sigma);
    }

    SD = sqrt(fabs(beta_new));
    f_old = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
    SD = sqrt(fabs(beta[which_beta]));
    f_new = exp(-sqrt(2)*fabs(beta[which_beta]-beta_new)/SD)/(sqrt(2)*SD);
	accept_prob = exp(g_2-g_1)*f_old/f_new;

    //printf("Accepted probability is %lf\n",accept_prob);

    //int abc=accept_prob > 1;
    //printf("if accept_prob bigger than 1 %d\n",abc);

    if (accept_prob > 1)

        beta[which_beta] = beta_new;

    else
    {

        GetRNGstate();

        x = rbinom(1,accept_prob);

       // printf("The bernulli r.v. is %d\n",x);

        PutRNGstate();

        if (x == 1) {beta[which_beta] = beta_new;accept_count++;}

    }

        //printf("updated beta_1 is %lf\n",beta[which_beta]);
 }


}


/*
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}
*/

double GXE_update_lambda(double *beta, double a, double b, int x_length)

{

    double lambda, beta_abs[x_length];

    int i;



    for (i=1; i<x_length+1; ++i)

        beta_abs[i-1] = fabs(beta[i]);

    GetRNGstate();

    lambda=rgamma((double) a+x_length, 1/(GXE_sum(beta_abs, x_length)+b));

    PutRNGstate();



    return lambda;

}




double GXE_update_sigma(double *beta, double asig, double bsig,double **xz, int N, int q, double *yall,int x_length,double **R, double *gamma)
{

    double sigma, ele[N], xbeta[N], Rgamma;
    int i, k;
    for (i=0; i<N; ++i)
    {
        
      Rgamma=0;
      for (k=0; k < q; ++k)
      {
          Rgamma=Rgamma+R[i][k]*gamma[k];
      }

      xbeta[i]=beta[0];


      for (k=0; k < x_length; ++k)

      {
          xbeta[i]=xbeta[i]+xz[i][k]*beta[k+1];
          
      }


        ele[i] = (yall[i]-xbeta[i]-Rgamma)*(yall[i]-xbeta[i]-Rgamma);
        
    }

    GetRNGstate();

    sigma=rgamma((double) asig+N/2,1/(bsig+GXE_sum(ele,N)/2));
    
    
    PutRNGstate();

    sigma = 1/sigma;




    return sigma;

}












void GXE_update_freq(double *freq, double *beta, double D, int h_length, double **xz, int N, int **hz, int n_cov, int Con)

{

    int i, C=Con, update=0;
    //printf("C is %d",C);
    double prop_freq_with_last[h_length], g_old=0, g_new=0, accept_prob=0, f_old, f_new, b_old[h_length], b_new[h_length], min_f_old, min_f_new;
    //double g_new_GXE_calc;
    GetRNGstate();

    for (i=0; i<h_length; ++i)

    {

        b_old[i] = freq[i]*C;
    }
    //printf("number of freq is %d",h_length);
    // see the trend of updated frequency

 /*  for (i=0; i<h_length; ++i)

    {

      printf("%e ", freq[i]);
      printf("%e ", b_old[i]);
    }
    printf("\n");
*/

    GXE_dirichlet(b_old, h_length, prop_freq_with_last);

	//for (i=0; i<h_length; ++i)
//	{
//	printf("this is the %d th freq \n",i);
//	printf("prop freq is %lf\n",prop_freq_with_last[i]);
	//printf("old freq is %lf\n",freq[i]);
//	}



    min_f_old = GXE_find_min(freq, h_length);

    /* check if the constraint -min fk/(1-min fk) < d is satisfied */

    min_f_new = GXE_find_min(prop_freq_with_last, h_length);

    if (-min_f_new/(1-min_f_new)  < D)

    {

        /* needed in acceptance prob. computation */
        for (i=0; i<h_length; ++i)

        {

            b_new[i] = prop_freq_with_last[i]*C;

            if (b_new[i] <= 0)

            {

                printf("Warning: prop_freq_with_last=%e, b_new = %lf\n", prop_freq_with_last[i], b_new[i]);

                for (i=0; i<h_length; ++i)
                {

                    printf("%e ", prop_freq_with_last[i]);

                }

            }

        }

        /* calculate g(f^(t)) and g(f*) */

        /*g_old=log(1-min_f_old)-GXE_calc_den_post(beta, freq, h_length, D, tot_uniq_x, xz, n_cov, N);*/

        /*g_new=log(1-min_f_new)-GXE_calc_den_post(beta, prop_freq_with_last, h_length, D, tot_uniq_x, xz, n_cov, N);*/

        g_old=log(1-min_f_old);

        g_new=log(1-min_f_new);

        //g_new_GXE_calc=log(1-min_f_new)-g_new;

        for (i=0; i<N; ++i)
        {
            g_old = g_old+log(GXE_calc_a(freq, hz[i], D));
            g_new = g_new+log(GXE_calc_a(prop_freq_with_last, hz[i], D));
          //  if (g_new_GXE_calc < -pow(10,10))
         //       printf("Warning: calc_den_new = %e, g_old = %lf, g_new = %lf\n", g_new_GXE_calc,g_old, g_new);


        }

	//	printf("g_new is %lf\n",g_new);
    //    printf("g_old is %lf\n",g_old);
        /* calculate f(f*|f^(t)) = f_new and f(f^(t)|f*) = f_old */

        f_old = lgammafn(C);

        f_new = f_old;

        for (i=0; i<h_length; ++i)

        {

            f_old += (b_new[i]-1)*log(freq[i]) - lgammafn(b_new[i]);
            f_new += (b_old[i]-1)*log(prop_freq_with_last[i]) - lgammafn(b_old[i]);

        }
	//	printf("f_old is %lf\n",f_old);
     //   printf("f_new is %lf\n",f_new);

        accept_prob = exp(g_new-g_old+f_old-f_new);



      //  printf("accept probability is %lf\n",accept_prob);

        if (-min_f_new/(1-min_f_new)  >= D) accept_prob = 0;

        if (accept_prob > 1) update = 1;

        else update = rbinom(1, accept_prob);

        if (update ==1)

        {

            for (i=0; i<h_length; ++i)
                  {
                freq[i] = prop_freq_with_last[i];
                  }
        }
    }

    if (-min_f_old/(1-min_f_old)  > D)

    {

        printf("error in updating f: Min f_old and D constraint violated min_f_old=%f\n", min_f_old);

        for (i=0; i<h_length; ++i)

           printf("%f\n", freq[i]);

        printf("%f\n", D);

        exit(1);

    }

    if (-min_f_new/(1-min_f_new)  > D &&  accept_prob>0)

    {
        printf("error in updating f: Min f_new and D constraint violated ");
        exit(1);
    }

    PutRNGstate();

}






double GXE_update_D(double *freq, double *beta, double D,int h_length, double **xz, int N, int **hz, int n_cov)

{

    int i, update = 0;

    double prop_D, accept_prob, g_old=0, g_new=0, min_f, delta=0.05, lower, upper, f_old, f_new;

    GetRNGstate();

    min_f = GXE_find_min(freq, h_length);

    lower = D-delta;

    upper = D+delta;


    if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);

    if (upper > 1) upper = 1;

    prop_D = runif(lower, upper);
    g_old=0;
    g_new=0;

    for (i=0; i<N; ++i)

	{
        g_old = g_old+log(GXE_calc_a(freq, hz[i], D));
        g_new = g_new+log(GXE_calc_a(freq, hz[i], prop_D));
	}

    f_new = 1/(upper-lower);

    lower = prop_D-delta;

    upper = prop_D+delta;

    if (lower < -min_f/(1-min_f)) lower = -min_f/(1-min_f);

    if (upper > 1) upper = 1;

    f_old = 1/(upper-lower);

    accept_prob = exp(g_new-g_old)*f_old/f_new;

    if (-min_f/(1-min_f) > D || -min_f/(1-min_f) > prop_D)

    {

        printf("error in updating D:  Min f and D constraint violated ");

        exit(1);

    }

    if (accept_prob > 1) update = 1;

    else update = rbinom(1, accept_prob);



    if (update == 1)

        return prop_D;

    else

        return D;



    PutRNGstate();

}

/* Calculate log of denominator of the posterior distribution which involves both beta and a(F) parameters; used in updating beta, f, and d

double GXE_calc_den_post(double *beta, double *freq, int h_length,
double D, int tot_hap, int **xz, int n_cov,  int N)
{
  int i, j, k;
  double  term[tot_hap], sum_term[N],sum_temp,sum;
  sum_temp=0;
  for (j=0; j<tot_hap; ++j)
    {

      term[j] =GXE_calc_a(freq, GXE_uniq_map[j], D)*exp(beta[0]);

      for (k=0; k < h_length-1; ++k)
	{
	  term[j] = term[j]*exp(GXE_uniq_x_mat[j][k]*beta[k+1]);
	}

      sum_temp =sum_temp+ term[j];
    }

   for (i=0;i<N;++i)
    {
     sum_term[i]=sum_temp;

      for (k=0; k < n_cov; ++k)
       {
      sum_term[i] =sum_term[i]* exp(xz[i][h_length+k-1]*beta[h_length+k]);
       }
     sum_term[i] = sum_term[i]+1;
     sum_term[i]=log(sum_term[i]);
    }
  sum=GXE_sum(sum_term,N);
  return sum;
}

*/






/* Calculate a(F) that is in the denominator of the likelihood*/

double GXE_calc_a(double *freq, int *per_freq, double D)

{

    int i, j, k;
    double a;
    i=per_freq[0];
     j=per_freq[1];

    if (i==j)
	{
        k=1;
    }
    else
    {
        k=0;
    }
    a = k*D*freq[i-1]+(2-k)*(1-D)*freq[i-1]*freq[j-1];
    return a;
}


/* function to find sum of real numbers */
double GXE_sum(double *x, int n)

{

    double sum=0.0;

    int i;
    for (i=0; i<n ; ++i)

        sum = sum + x[i];
    return sum;
}

/* function to calculate min. of an array of numbers of length n */

double GXE_find_min(double *arr, int n)

{

    int i;

    double min=arr[0];

    for(i=1;i<n; ++i)

    {

        if(min > arr[i])

            min = arr[i];

    }

    return min;

}

/* function to generate from double exponential distribution */

double GXE_gen_double_exp(double mean, double SD)

{

    double x, gen_exp;

    GetRNGstate();

    x = runif(0,1);

//    printf("uniform is %lf\n",x);

    gen_exp = rexp(SD/sqrt(2));

//    printf("Exp is %lf\n",gen_exp);

    PutRNGstate();

    if(x > 0.5)

        return gen_exp+mean;

    else

        return -gen_exp+mean;

}

/* function to generate from Dirichet distribution */

void GXE_dirichlet(double *param, int dim, double *gen_sample)

{

    int i,j=0;

    double gen_gamma[dim], sum_gamma;

    GetRNGstate();

    for (i=0; i<dim; ++i)

    {

        //if (param[i]<=0) //printf("Warning: param=%lf\n", param[i]);

        gen_gamma[i] = rgamma(param[i], 1);
        if (gen_gamma[i]<0.000001)
	{
		//printf("WarningOld: gen_gamma=%lf\n", gen_gamma[i]);
		gen_gamma[i]=0.000001;
		//printf("WarningNew: gen_gamma=%lf\n", gen_gamma[i]);


	}

     j=j+1;


      //  if (gen_gamma[i]<=0) //printf("Warning: gen_gamma=%lf, param=%lf\n", gen_gamma[i], param[i]);

    }



    sum_gamma = GXE_sum(gen_gamma, dim);


    for (i=0; i<dim; ++i)

    {

        gen_sample[i] = gen_gamma[i]/sum_gamma;

    }

    PutRNGstate();

}

double GXE_update_gamma_quadraticPart(double *gamma, int q)
{
	double result=0;
	int i;
	for (i=0; i<q; ++i)
	{
        result=result+gamma[i]*gamma[i];
	}
	return result;
}


double GXE_update_randomAsigma(double asigRandomA, double bsigRandomA, int q, double quadratic_sum)
{
	double randomA_sigma;

	GetRNGstate();
	randomA_sigma=rgamma((double) asigRandomA+q/2,1/(bsigRandomA+quadratic_sum/2));
	PutRNGstate();
	randomA_sigma = 1/randomA_sigma;
	return randomA_sigma;
}






void GXE_update_gamma(double *gamma, double *mean_vector, double **cov_matrix, double **R, double sigmas, double randomA_sigmas, double *y, int N, int q, double *beta, double **xz, int x_length) //new!
{
	int i2,i4,s,k;
	double xbeta[N];

	gsl_matrix *m1=gsl_matrix_alloc(q,q);
	gsl_matrix *m2=gsl_matrix_alloc(q,q);
	gsl_matrix *mnorm_cov_mat=gsl_matrix_alloc(q,q);
	gsl_matrix *y_gsl=gsl_matrix_alloc(N,1); //matrix
	gsl_matrix *y_hap_only=gsl_matrix_alloc(N,1);
    gsl_matrix *temp=gsl_matrix_alloc(q,1);
	gsl_matrix *mean_mat=gsl_matrix_alloc(q,1);
	gsl_matrix *inv_m=gsl_matrix_alloc(q,q); //cov matrix of gamma
	gsl_permutation * p = gsl_permutation_alloc (q);
	gsl_vector *mean_vec=gsl_vector_alloc(q);
	gsl_vector *gamma_result=gsl_vector_alloc(q);
    gsl_matrix *R_gsl=gsl_matrix_alloc(N,q);
    

	gsl_matrix_set_identity (m2);

    for (i2=0;i2<N;i2++)
    {
        xbeta[i2]=beta[0];
        for (k=0; k<x_length; k++)
        {
            xbeta[i2]=xbeta[i2]+xz[i2][k]*beta[k+1]; //x*beta for individual i2 in a pedigree
        }
    }

    for (i2=0;i2<N;i2++)
    {
        gsl_matrix_set(y_gsl,i2,0,y[i2]);
        gsl_matrix_set(y_hap_only,i2,0,xbeta[i2]);
    }
    
    
    for (i2=0;i2<N;i2++)
    {
        for (i4=0; i4<q;i4++)
        {
            gsl_matrix_set(R_gsl,i2,i4,R[i2][i4]); //m2[i2][i4]: kinship_inv[i2][i4]
        }
        
    }

    //printf("R ",R_gsl[0][0],R_gsl[0][1]);
    //t(R)%*%R
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                    1.0, R_gsl, R_gsl, //matrix calculation
                    0.0, m1); //m1=t(R)%*%R
    
    gsl_matrix_scale (m1, 1/sigmas);
    gsl_matrix_scale (m2, 1/randomA_sigmas);
    gsl_matrix_add (m1, m2); //sum becomes m1, m1 is the inverse of cov matrix
    //printf("inverse of cov ",m1[0][0],R_gsl[0][1]);
    
    //find the inverse
    gsl_linalg_LU_decomp (m1, p, &s); //has to be &s
    gsl_linalg_LU_invert (m1, p, inv_m);
    gsl_matrix_memcpy (mnorm_cov_mat, inv_m);


    gsl_matrix_scale (inv_m, 1/sigmas); // inv_m is the cov matrix of posterior of A
    gsl_matrix_sub (y_gsl, y_hap_only);// subtraction
    gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                    1.0, R_gsl, y_gsl, //matrix calculation
                    0.0, temp); //new: temp=t(R)Y !!!!!!!!! R_gsl and y_gsl have to be constant

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
					1.0, inv_m, temp, //matrix calculation
					0.0, mean_mat); //understand the parameters here

    for (i2=0;i2<q;i2++)
    {
        gsl_vector_set(mean_vec,i2,gsl_matrix_get(mean_mat,i2,0));
    }

    
    //new !!!!!!!!!! store and then print out the mean_vector and cov_matrix of gamma
      for (i2=0;i2<q;i2++)
      {
          mean_vector[i2]=gsl_vector_get(mean_vec,i2);
          for(i4=0;i4<q;i4++)
          {
              cov_matrix[i2][i4]=gsl_matrix_get(mnorm_cov_mat,i2,i4);
          }
          
      }
    
    
    
    

    const gsl_rng_type * T;
          gsl_rng *r;

          gsl_rng_env_setup();
          T=gsl_rng_default;
          r=gsl_rng_alloc(T);

    gsl_ran_multivariate_gaussian(r,mean_vec,mnorm_cov_mat,gamma_result);

    for (i2=0;i2<q;i2++)
    {
        gamma[i2]=gsl_vector_get(gamma_result,i2);
    }



    gsl_rng_free (r);
    gsl_matrix_free (m1);
    gsl_matrix_free (m2);
    gsl_matrix_free (mnorm_cov_mat);
    gsl_matrix_free (y_gsl);
    gsl_matrix_free (y_hap_only);
    gsl_matrix_free (temp);
    gsl_matrix_free (mean_mat);
    gsl_matrix_free (inv_m);
    gsl_permutation_free (p);
    gsl_vector_free (mean_vec);
    gsl_vector_free (gamma_result);
    gsl_matrix_free (R_gsl);


}

/* Generate a random vector from a multivariate Gaussian distribution using
 * the Cholesky decomposition of the variance-covariance matrix, following
 * "Computational Statistics" from Gentle (2009), section 7.4.
 *
 * mu      mean vector (dimension d)
 * L       matrix resulting from the Cholesky decomposition of
 *         variance-covariance matrix Sigma = L L^T (dimension d x d)
 * result  output vector (dimension d)
 */
int gsl_ran_multivariate_gaussian (const gsl_rng * r, const gsl_vector * mu, const gsl_matrix * L, gsl_vector * result)
{
  const size_t M = L->size1;
  const size_t N = L->size2;

  if (M != N)
    {
      GSL_ERROR("requires square matrix", GSL_ENOTSQR);
    }
  else if (mu->size != M)
    {
      GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
    }
  else if (result->size != M)
    {
      GSL_ERROR("incompatible dimension of result vector", GSL_EBADLEN);
    }
  else
    {
      size_t i;

      for (i = 0; i < M; ++i)
        gsl_vector_set(result, i, gsl_ran_ugaussian(r));

      gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, L, result);
      gsl_vector_add(result, mu);

      return GSL_SUCCESS;
    }
}
