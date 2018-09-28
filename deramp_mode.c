#include <stdio.h>
#include <stdlib.h>
 
#define ARRSIZE(arr) (sizeof(arr)/sizeof(*(arr)))

// compute moving average
void average(tr,windowlen,nx2)
double *tr;
int windowlen;
int nx2;
{
  int uk,bk,ck;
  double sum1,sum2;
  double *buffer;


  buffer = (double *)malloc((size_t)(windowlen*sizeof(double)));

  for(uk=nx2;uk<=windowlen-(nx2+1);uk++)
  {   
       sum1=0.0;
       sum2=0.0;
       for(bk=1;bk<=nx2;bk++)
       {
            sum1=sum1+(*(tr + uk-bk));
            sum2=sum2+(*(tr + uk+bk));
        }
        *(buffer+uk) = (sum1+sum2+(*(tr+uk)))/(float)(2*nx2+1);
  }
 
  for(uk=1;uk<=nx2-1;uk++)
  {
       sum1=0.0;
       sum2=0.0;
       for(bk=1;bk<=uk;bk++)
       {
           sum1=sum1+(*(tr + uk-bk));
           sum2=sum2+(*(tr + uk+bk));
           ck=bk;
        }
        *(buffer+uk) = (sum1+sum2+(*(tr+uk)))/(float)(2*ck+1);
  }
 
  for(uk=1;uk<=nx2-1;uk++)
  {
       sum1=0.0;
       sum2=0.0;
       for(bk=1;bk<=uk;bk++)
       {
            sum1=sum1+(*(tr + windowlen-1-uk-bk));
            sum2=sum2+(*(tr + windowlen-1-uk+bk));
            ck=bk;
       }
       *(buffer+windowlen-1-uk) =
            (sum1+sum2+(*(tr+windowlen-1-uk)))/(float)(2*ck+1);
  }
 
  *(buffer+windowlen-1) = *(tr+windowlen-1);
  *(buffer) = *(tr);
  for(uk=0;uk<windowlen;uk++)
         tr[uk] = buffer[uk];
 
  free((char *) buffer);
}



// this is the compare function for qsort 
int compare(const void *a, const void *b) {
    double x1 = *(const double*)a;
    double x2 = *(const double*)b;
    if (x1 > x2) return  1;
    if (x1 < x2) return -1;
    return 0;
}

// function to calculate mode
long calc_mode(int n, double a[]) {
    double maxValue = 0;
    int maxCount = 0, i, j;
    for (i = 0; i < n; ++i) {
        int count = 0;
              
        for (j = 0; j < n; ++j) {
            if (a[j] == a[i])
                ++count;
            }
                                          
        if (count > maxCount) {
            maxCount = count;
            maxValue = a[i];
            }
        }
   return maxValue;
}

// compute the average over all second row elements with same first row element
int sorted_avg(double matr[][2], double matr_sorted[][2], int N) {
	int i, j, k, l, i0, i1;
    double* arr;
	j = 0;
    i0 = 0;
    i1 = 0;
	for (i = 0; i < N-1; i++){
		if (matr[i][0] == matr[i+1][0]){
            i1 += 1;
		}
		if (matr[i][0] != matr[i+1][0]){
			matr_sorted[j][0] = matr[i][0];
            l = i1 - i0 + 1;
            arr = (double *)malloc((size_t)(l*sizeof(double)));
	        for (k = i0; k <= i1; k++){
                arr[k-i0] = matr[k][1];    
                }
			matr_sorted[j][1] = calc_mode(l, arr);
            free(arr);
            i1 += 1;
            i0 = i1;
			j += 1;
		}
	}
	return j;
}

// compute the average over all second row elements with same first row element
int sorted_avg_old(double matr[][2], double matr_sorted[][2], int N) {
	int i, j;
	double a, l;
	l = 1.;
	j = 0;
	a = matr[0][1];
	for (i = 0; i < N-1; i++){
		if (matr[i][0] == matr[i+1][0]){
			l += 1.;
			a += matr[i+1][1];
		}
		if (matr[i][0] != matr[i+1][0]){
			matr_sorted[j][0] = matr[i][0];
			matr_sorted[j][1] = a / l;
			l = 1.;
			a = matr[i+1][1];
			j += 1;
		}
	}
	return j;
}

// subtract sorted average from the original second row elements according to first row elements
int subtract_ramp(double matr[][2], double** matr_sorted, int N, int n){
	int i, j;
	for (i = 0; i < N; i++){
		for (j = 0; j < n; j++){
			if (matr[i][0] == matr_sorted[j][0]){
				matr[i][1] = matr[i][1] - matr_sorted[j][1];
			}
		}
	}
	return 0;
}



// final deramp function
int deramp_mode(double matr[][2], double matr_orig[][2], double matr_sorted[][2], double x[1], int N, int n_avg, double mean[1]){
	int n, i, j;
    double* arr;
    double** matr_sorted_short;

    arr = (double *)malloc((size_t)(N*sizeof(double)));
    // calculate overall mean value
	for (i = 0; i < N; i++){
        arr[i] = matr[i][1];
    }
    mean[0] = calc_mode(N, arr);

    qsort(matr, N, sizeof(*matr), compare);
 	n = sorted_avg(matr, matr_sorted, N);
    printf("\n (%d) \n", n);

    matr_sorted_short = (double **)malloc((size_t)(n*sizeof(double*)));
    for (i = 0; i < n; i++){
        matr_sorted_short[i] = (double *)malloc((size_t)(2*sizeof(double)));
    }
	for (i = 0; i < n; i++){
		for (j = 0; j < 2; j++){
			matr_sorted_short[i][j] = matr_sorted[i][j];
		}
	}
    
	for (i = 0; i < n; i++){
        x[i] = matr_sorted_short[i][1];
    }
    if (n_avg > 0){
        average(x, n, n_avg);
        }
	for (i = 0; i < n; i++){
	    matr_sorted_short[i][1] = x[i] - mean[0];
	}
	subtract_ramp(matr_orig, matr_sorted_short, N, n);

	return n;
}

