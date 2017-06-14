/*
                             Developed by
                 Laboratory of Computational Astrophysics
                University of Illinois at Urbana-Champaign
*/
#ifdef UNICOS
/***********************************************************************
    Check if input is pending on stdin. If so, read it and return
    it to the caller.
 
  Unicos version.
  On a Cray-2, compile with
    % cc -Dcray2 checkin.c
 
  Input:
    wait        Flag to indicate whether to wait for input.
 
  Output:
    string      Buffer containing the text read from stdin.
*/
 
#include <sys/time.h>
 
/* Define some macros which convert FORTRAN strings into something
   usable in C. */
 
typedef int STRING;
#ifdef cray2
#define len(string) ((string & 0x1FFFFFF800000000) >> 35 )
#define addr(string) ((char *)(string & 0xE0000007FFFFFFFF))
#else
#define len(string) ((string & 0x1FFFFFFF00000000) >> 32 )
#define addr(string) ((char *)(string & 0xE0000000FFFFFFFF))
#endif

/**********************************************************************/
int CHECKIN(string,wait)
STRING string;
int *wait;
{
  int mask,nfound,i,length;
  char *q,*msg;
  struct timeval timeout;

/* Check in the standard input is a terminal. */

  if(! isatty(0)) return(0);

/* Convert the FORTRAN Unicos calling sequence to something more
   usable. */

  msg = addr(string);
  length = len(string);

/* Something to read? */

  mask = 1;
  timeout.tv_sec = 0;
  timeout.tv_usec = 0;
  if(*wait)nfound = select(1,&mask,0,0,&timeout);
  else     nfound = 1;
  if(!nfound)return(0);

/* Do a read of the input, and pack it into an integer array. */

  nfound = read(0,msg,length) - 1;
  q = (char *)msg + nfound;
  for(i=nfound; i < length; i++) *q++ = ' ';
  return(nfound);
}
#else
/***********************************************************************
    Check if input is pending on stdin. If so, read it and return
    it to the caller.
  Output:
    msg		Buffer containing the text read from stdin.  */

#include <sys/time.h>

int checkin_(msg,wait,length)
int *msg,*wait,length;
{
/*  int mask,nfound,i; */
  int nfound,i;
  fd_set mask;
  char *q;
  struct timeval timeout;

/* Something to read? */

/*  mask = 1; */
  FD_SET(1,&mask);
  timeout.tv_sec = 0;
  timeout.tv_usec = 0;
  if(*wait)nfound = select(1,&mask,0,0,&timeout);
  else     nfound = 1;
  if(!nfound)return(0);

/* Do a read of the input, and pack it into an integer array. */

  nfound = read(0,msg,length) - 1;
  q = (char *)msg + nfound;
  for(i=nfound; i < length; i++) *q++ = ' ';
  return(nfound);
}
#endif
