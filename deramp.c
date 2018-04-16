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
       for(bk=0;bk<=nx2;bk++)
       {
            sum1=sum1+(*(tr + uk-bk));
            sum2=sum2+(*(tr + uk+bk));
        }
        *(buffer+uk) = (sum1+sum2+(*(tr+uk)))/(float)(2*nx2+2);
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
        *(buffer+uk) = (sum1+sum2+(*(tr+uk)))/(float)(2*ck+2);
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
            (sum1+sum2+(*(tr+windowlen-1-uk)))/(float)(2*ck+2);
  }
 
  *(buffer+windowlen-1) = *(tr+windowlen-1);
  *(buffer) = *(tr);
  for(uk=0;uk<windowlen;uk++)
         tr[uk] = buffer[uk];
 
  free((char *) buffer);
}



// this is the compare function for qsort 
int compare(const void *a, const void *b) {
    long x1 = *(const long*)a;
    long x2 = *(const long*)b;
    if (x1 > x2) return  1;
    if (x1 < x2) return -1;
    return 0;
}

// compute the average over all second row elements with same first row element
int sorted_avg(long matr[][2], double matr_sorted[][2], int N) {
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
int deramp(long matr[][2], double matr_orig[][2], double matr_sorted[][2], int N, int n_avg){
	int n, i, j;
    double mean;
    double* x;
    double** matr_sorted_short;

    mean = 0.;
    // calculate overall mean value
	for (i = 0; i < N; i++){
        mean += matr[i][1];
    }
    mean = (mean / N);

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
    
    x = (double *)malloc((size_t)(n*sizeof(double)));
	for (i = 0; i < n; i++){
        x[i] = matr_sorted_short[i][1] - mean;
    }
    average(x, n, n_avg);
	for (i = 0; i < n; i++){
	    matr_sorted_short[i][1] = x[i];
	}
	subtract_ramp(matr_orig, matr_sorted_short, N, n);

	return 0;
}

 
int main(void) {
    long matrix[][2] = {{8,6}, {4,2}, {1,0}, {4,8}, {2,4},
                       {4,3}, {1,2}, {2,2}, {8,3}, {5,5}};
    double matrix_orig[][2] = {{8,6}, {4,2}, {1,0}, {4,8}, {2,4},
                       {4,3}, {1,2}, {2,2}, {8,3}, {5,5}};
	int N;
	size_t i;
 
	N = ARRSIZE(matrix);
	double m[N][2];

    printf("Original: ");
    for (i = 0; i < ARRSIZE(matrix_orig); i++)
        printf("(%f,%f) ", matrix_orig[i][0], matrix_orig[i][1]);
    putchar('\n');
 
	deramp(matrix, matrix_orig, m, N, 3);

    printf("Sorted  : ");
    for (i = 0; i < ARRSIZE(matrix); i++)
        printf("(%li,%li) ", matrix[i][0], matrix[i][1]);
    putchar('\n');
    printf("deramped: ");
    for (i = 0; i < ARRSIZE(matrix); i++)
        printf("(%f,%f) ", matrix_orig[i][0], matrix_orig[i][1]);
    putchar('\n');
 
    return 0;
}
