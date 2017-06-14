/*
                             Developed by
                 Laboratory of Computational Astrophysics
                University of Illinois at Urbana-Champaign
*/
#ifdef UNICOS
/**********************************************************************/
BCDFLT(in,x,ioff,n)
int *in;
float *x;
int *ioff,*n;
/*
  Convert a character string into a floating point number.
----------------------------------------------------------------------*/
{
  char buf[80],*p,*q;
  int i;
  double atof();

  strncpy(buf,(char *)in + *ioff,*n);
  *x = atof(buf);
}
#else
/**********************************************************************/
void bcdflt_(in,x,ioff,n)
int *in;
float *x;
int *ioff,*n;
/*
  Convert a character string into a floating point number.
----------------------------------------------------------------------*/
{
  char buf[80],*p,*q;
  int i;
  double atof();

  strncpy(buf,(char *)in + *ioff,*n);
  *x = atof(buf);
}
#endif
