#include <stdio.h>
#include <string.h>
#include <malloc.h>


void average(tr,windowlen,nx2)
float *tr;
int windowlen;
int nx2;
{
  int uk,bk,ck;
  float sum1,sum2;
  float *buffer;


  buffer = (float *)malloc((size_t)(windowlen*sizeof(float)));

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


