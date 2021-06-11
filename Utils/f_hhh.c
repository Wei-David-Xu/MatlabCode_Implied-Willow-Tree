#include <stdio.h>
#include "nr.h"
#include "nrutil.c"
#include "mex.h"
#include <math.h>
#include "nrutil.h"

/**********************************************************************
**
**  C equivalent of the Fortran function "sign"
**
***********************************************************************/

double sign(double x, double y)
{
  if (y>=0.0)
     return fabs(x);
  else
     return -fabs(x);
}


/**********************************************************************
**
**  Subroutine Mom
**
***********************************************************************/


void mom(double g, double d, double *a, double *fault)
{

    /*declarations*/
      double *b,*c;
      double zz,vv,rttwo,rrtpi,w,e,r,h,t,u,y,x,v,f,z,s,p,q,aa,ab,expa,expb,zero,
             quart,half,p75,one,two,three,l,limit,k,m;
      int i;


    /*allocate space for vectors*/
      b=dvector(1,6); c=dvector(1,6);

    /*assign values*/
      zz=1.0e-5; vv=1.0e-8; limit=500.0;
      rttwo=1.414213562; rrtpi=0.5641895835; expa=80.0; expb=23.7;
      zero=0.0; quart=0.25; half=0.5; p75=0.75; one=1.0; two=2.0; three=3.0;
/*
//    main code starts here
*/
      *fault=0.0;
      for (i=1;i<=6;i++) {c[i]=zero;}
      w= g / d;
/*
//    Trial value of h
*/
	  if (w>expa) {goto g140;}
	  e=exp(w)+one;
	  r=rttwo / d;
	  h=p75;
	  if ( d<three) {h=quart * d;};
	  k=1.0;
	  goto g40;
/*
//    start of outer loop
*/
g20:  k=k+1.0;
	  if (k>limit) {goto g140;}
	  for (i=1;i<=6;i++) {c[i]=a[i];}
/*
//    no convergence yet - try smaller h
*/
      h=half*h;
g40:  t=w;
	  u=t;
	  y=h*h;
	  x=two*y;
	  a[1]=one/e;
	  for (i=2;i<=6;i++) {a[i]=a[i-1]/e;}
	  v=y;
	  f=r*h;
	  m=0.0;
/*
//    start inner loop to evaluate infinite series
*/
g60:  m=m+1.0;
	  if (m>limit) {goto g140;}
	  for (i=1;i<=6;i++) {b[i]=a[i];}
	  u=u-f;
	  z=one;
	  if (u>(-expb)) {z=exp(u)+z;}
	  t=t+f;
	  if (t>expb) {l=1.0;} else {l=0.0;}
	  if (l==0.0) {s=exp(t)+one;}
	  p=exp(-v);
	  q=p;
	  for (i=1;i<=6;i++) {
		  aa=a[i];
	      p=p/z;
		  ab=aa;
		  aa=aa+p;
		  if (aa==ab) {goto g100;}
		  if (l==1.0) {goto g80;}
		  q=q/s;
		  ab=aa;
		  aa=aa+q;
		  if (aa==ab) {l=1.0;} else {l=0.0;};
g80:      a[i]=aa;
	  }
g100: y=y+x;
	  v=v+y;
  	  for (i=1;i<=6;i++) {
		  if (a[i]==0.0) {goto g140;}
		  if (fabs((a[i]-b[i])/a[i])>vv) {goto g60;}
	  }
/*
//    end of inner loop
*/
      v=rrtpi*h;
	  for (i=1;i<=6;i++) {a[i]=v*a[i];}
	  for (i=1;i<=6;i++) {
		  if (a[i]==0.0) {goto g140;}
		  if (fabs((a[i]-c[i])/a[i])>zz) {goto g20;}
	  }
	  goto gmem;

g140: *fault=1.0;
      goto gmem;

gmem: free_dvector(b,1,6);
      free_dvector(c,1,6);

return;

}


/**********************************************************************
**
**  Subroutine SBfit	
**
***********************************************************************/


void sbfit(double xbar,double sigma, double rtb1,double b2,double tol,double *gamma, double *delta, double *xlam, double *xi, double *fault)
{

    /*declarations*/
      double *hmu,*deriv,*dd,tt,rb1,b1,e,u,x,y,w,f,d,g,s,h2,t,h2a,h2b,h3,h4,rbet,bet2,zero,one,two,three,four,six;
	  double half,quart,one5,neg,limit,m;
	  int j,k;
	  double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22;

      a1=0.0124;  a2=0.0623;  a3=0.4043;  a4=0.408;  a5=0.479;    a6=0.485;
      a7=0.5291;  a8=0.5955;  a9=0.626;   a10=0.64;  a11=0.7077;  a12=0.7466;
      a13=0.8;    a14=0.9281; a15=1.0614; a16=1.25;  a17=1.7973;  a18=1.8;
      a19=2.163;  a20=2.5;    a21=8.5245; a22=11.346;

    /*allocate space for vectors*/
      hmu=dvector(1,6); deriv=dvector(1,4); dd=dvector(1,4);

      zero=0.0; one=1.0; two=2.0; three=3.0; four=4.0; six=6.0; half=0.5; quart=0.25; one5=1.5; tt=1.0e-4; limit=50.0; 

/*
//Program begins here
*/
	  rb1=fabs(rtb1);
      b1=rb1*rb1;
	  if (rtb1<zero) 
	  {neg=1.0;}
	  else
	  {neg=0.0;}
/*
//Get d as a first estimate of delta
*/
	  e=b1+one;
	  x=half*b1+one;
	  y=fabs(rb1)*sqrt(quart*b1+one);
	  u=pow((x+y),one/three);
	  w=u+one/u-one;
	  f=w*w*(three+w*(two+w))-three;
      e=(b2-e)/(f-e);
	  if (fabs(rb1)>tol) {goto g5;}
      f=two;
	  goto g20;
g5:   d=one/sqrt(log(w));
	  if (d<a10) {goto g10;}
	  f=two-a21 / (d*(d*(d-a19)+a22));
	  goto g20;
g10:  f=a16*d;
g20:  f=e*f+one;
	  if (f<a18) {goto g25;}
	  d=(a9*f-a4)*pow((three-f),(-a5));
	  goto g30;
g25:  d=a13*(f-one);
/*
//    Get G as first estimate of gamma
*/
g30:  g=zero;
	  if (b1<tt) {goto g70;}
	  if (d>one) {goto g40;}
	  g=(a12*pow(d,a17)+a8)*pow(b1,a6);
	  goto g70;
g40:  if (d<=a20) {goto g50;}
	  u=a1;
	  y=a7;
	  goto g60;
g50:  u=a2;
	  y=a3;
g60:  g=pow(b1,(u*d+y))*(a14+d*(a15*d-a11));
g70:  m=0.0;
/*
//    main iteration starts here
*/
g80:  m=m+1.0;
	  if (m>limit) 
	  {*fault=1.0;}
	  else
  	  {*fault=0.0;}

	  if (*fault==1.0) {goto gmem;}
/*
//    Get first six moments for latest g and d values
*/
	  mom(g, d, hmu, fault);
	  if (*fault==1.0) {goto gmem;}
	  s=hmu[1]*hmu[1];
	  h2=hmu[2]-s;
	  if (h2<=0.0)
	  {*fault=1.0;}
	  else
	  {*fault=0.0;}
	  if (*fault==1.0) {goto gmem;}
	  t=sqrt(h2);
	  h2a=t*h2;
	  h2b=h2*h2;
	  h3=hmu[3]-hmu[1]*(three*hmu[2]-two*s);
	  rbet=h3/h2a;
	  h4=hmu[4]-hmu[1]*(four*hmu[3]-hmu[1]*(six*hmu[2]-three*s));
	  bet2=h4/h2b;
	  w=g*d;
	  u=d*d;
/*
//    Get derivatives
*/
   	  for (j=1;j<=2;j++) {
		  for (k=1;k<=4;k++) {
               t=(double)k;
			   if (j==1) {goto g90;}
               s=((w-t)*(hmu[k]-hmu[k+1])+(t+one)*(hmu[k+1]-hmu[k+2]))/u;
			   goto g100;
g90:           s=hmu[k+1]-hmu[k];
g100:          dd[k]=t*s/d;
		  }
		  t=two*hmu[1]*dd[1];
		  s=hmu[1]*dd[2];
		  y=dd[2]-t;
		  deriv[j]=(dd[3]-three*(s+hmu[2]*dd[1]-t*hmu[1])-one5*h3*y/h2)/h2a;
		  deriv[j + 2] = (dd[4] - four * (dd[3] * hmu[1] + dd[1] * hmu[3]) + six * (hmu[2] * t + hmu[1] * (s - t * hmu[1])) - two * h4 * y / h2) / h2b;
      }
	  t=one/(deriv[1]*deriv[4]-deriv[2]*deriv[3]);
	  u=(deriv[4]*(rbet-rb1)-deriv[2]*(bet2-b2))*t;
	  y=(deriv[1]*(bet2-b2)-deriv[3]*(rbet-rb1))*t;
/*
//    form new estimates of g and d
*/
	  g=g-u;
	  if ( (b1==zero) || (g<zero) ) {g=zero;}
	  d=d-y;
	  if ( (fabs(u)>tt) || (fabs(y)>tt) ) {goto g80;}
/*
//    end of iteration
*/
      *delta=d;
	  *xlam=sigma/sqrt(h2);
	  if (neg==1.0) {goto g130;}
	  *gamma=g;
	  goto g140;
g130: *gamma=-g;
      hmu[1]=one-hmu[1];
g140: *xi=xbar-*xlam*hmu[1];
	  goto gmem;
    /*
    //release memory
    */
gmem: free_dvector(hmu,1,6);
	  free_dvector(deriv,1,4);
	  free_dvector(dd,1,4);
	  return;
}




/**********************************************************************
**
**  Subroutine Sufit	
**
***********************************************************************/

void sufit(double xbar,double sd, double rb1,double b2,double tol,double *gamma, double *delta, double *xlam, double *xi)
{

    /*declarations*/
      double b1,b3,w,y,w1,wm1,z,v,a,b,x,zero,one,two,three,four,six,seven,eight,nine,ten,half,one5,two8;

      zero=0.0; one=1.0; two=2.0; three=3.0; four=4.0; six=6.0; seven=7.0; eight=8.0; nine=9.0; half=0.5; one5=1.5; two8=2.8; ten=10.0;

    /*program starts here*/
	  b1=rb1*rb1;
	  b3=b2-three;
/*
// W is first estimate of exp(delta ** (-2))
*/
      w=sqrt(sqrt(two*b2-two8*b1-two)-one);
	  if (fabs(rb1)>tol) {goto g10;}
/*
// Symmetrical case - results are known
*/
      y=zero;
	  goto g20;
/*
// Johnson iteration (using y for his M)
*/
g10:  w1=w+one;
	  wm1=w-one;
	  z=w1*b3;
	  v=w*(six+w*(three+w));
	  a=eight * (wm1*(three+w*(seven+v))-z);
	  b=16.0*(wm1*(six+v)-b3);
	  y=(sqrt(a*a-two*b*(wm1*(three+w*(nine+w*(ten+v)))-two*w1*z))-a)/b;
	  z=y*wm1*pow((four*(w+two)*y+three*w1*w1),2.0) / (two*pow(two*y+w1,3.0));
	  v=w*w;
	  w=sqrt(sqrt(one-two*(one5-b2+(b1*(b2-one5-v*(one+half*v)))/z))-one);
	  if (fabs(b1-z)>tol) {goto g10;}
/*
// End of iteration
*/
	  y=y/w;
	  y=log(sqrt(y)+sqrt(y+one));
	  if (rb1>zero) {y=-y;}
g20:  x=sqrt(one/log(w));
	  *delta=x;
	  *gamma=y*x;
	  y=exp(y);
	  z=y*y;
	  x=sd/sqrt(half*(w-one)*(half*w*(z+one/z)+one));
	  *xlam=x;
      *xi=(half*sqrt(w)*(y-one/y))*x+xbar;
	  return;
}



/***********************************************************************
**
**  Hill, Hill and Holder algorithm
**
***********************************************************************/

void hhh(double xbar,double sd,double rb1,double bb2,double *itype, double *gamma, double *delta, double *xlam, double *xi, double *ifault)
{

    /*declarations*/
      double tol,b1,b2,y,x,u,w,fault,zero,one,two,three,four,half,quart;

      tol=0.000001; zero=0.0; one=1.0; two=2.0; three=3.0; four=4.0; half=0.5; quart=0.25;

    /*main code*/
      *ifault=1.0;
      if (sd<zero) {return;}
      *ifault=zero;
      *xi=zero;
      *xlam=zero;
      *gamma=zero;
      *delta=zero;
      if (sd>zero) {goto g10;}
      *itype=5.0;
      *xi=xbar;
      return;
g10:  b1=rb1*rb1;
      b2=bb2;
      fault=zero;
/*
// Test whether lognormal (or normal) requested
*/
      if (b2>=zero) {goto g30;}
g20:  if (fabs(rb1)<=tol) {goto g70;}
	  goto g80;
/*
// Test for position relative to boundary line
*/
g30:  if (b2>(b1+tol+one)) {goto g60;}
	  if (b2<(b1+one)) {goto g50;}
/* 
// ST distribution
*/
g40:  *itype=5.0;
      y=half+half*sqrt(one-four/(b1+four));
	  if (rb1>zero) {y=one-y;}
	  x=sd/sqrt(y*(one-y));
	  *xi=xbar-y*x;
	  *xlam=*xi+x;
	  *delta=y;
	  return;
g50:  *ifault=2.0;
	  return;
g60:  if ((fabs(rb1)>tol) || (fabs(b2-three)>tol)) {goto g80;}
/*
// Normal distribution
*/
g70:  *itype=4.0;
	  *delta=one/sd;
	  *gamma=-xbar/sd;
      *xlam=1.0;
	  return;
/*
// Test for position relative to lognormal line
*/
g80:  x = half *b1 + one;
	  y = fabs(rb1) * sqrt(quart * b1 + one);
      u = pow( (x + y), (one / three) );
      w = u + one / u - one;
      u = w * w * (three + w * (two + w)) - three;
	  if ((b2<zero) || (fault==1.0)) {b2=u;}
	  x=u-b2;
	  if (fabs(x)>tol) {goto g90;}
/*
// Lognormal (SL) distribution
*/
	  *itype=1.0;
	  *xlam=sign(one,rb1);
	  u=*xlam*xbar;
	  x=one/sqrt(log(w));
	  *delta=x;
	  y=half*x*log(w*(w-one)/(sd*sd));
	  *gamma=y;
	  *xi=*xlam*(u-exp((half/x-y)/x));
	  return;
/*
// SB or SU distribution
*/
g90:  if (x > zero) {goto g100;}
	  *itype=2.0;
	  sufit(xbar,sd,rb1,b2,tol,gamma,delta,xlam,xi);
	  return;
g100: *itype=3.0;
	  sbfit(xbar,sd,rb1,b2,tol,gamma,delta,xlam,xi,&fault);
	  if (fault==0.0) {return;}
/*
// Failure - try to fit approximate result
*/
	  *ifault=3.0;
	  if (b2>b1+two) {goto g20;}
	  goto g40;

}






/***********************************************************************
**
**  Gateway function
**
***********************************************************************/

void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])

{
/*vectors passed to the Matlab function*/
  double xbar,sd,rb1,bb2;

/*vector returned to the Matlab function*/
  double *vitype,*vgamma,*vdelta,*vxlam,*vxi,*vifault;
  double   itype,  gamma,  delta,  xlam,  xi,  ifault;


/*Right hand side of the MEX instruction*/
  xbar     = mxGetScalar(prhs[0]);   
  sd       = mxGetScalar(prhs[1]);
  rb1      = mxGetScalar(prhs[2]);
  bb2      = mxGetScalar(prhs[3]);

/*Left hand side of the MEX instruction */
  plhs[0]=mxCreateDoubleMatrix((int) 1,(int) 1,mxREAL);
  plhs[1]=mxCreateDoubleMatrix((int) 1,(int) 1,mxREAL);
  plhs[2]=mxCreateDoubleMatrix((int) 1,(int) 1,mxREAL);
  plhs[3]=mxCreateDoubleMatrix((int) 1,(int) 1,mxREAL);
  plhs[4]=mxCreateDoubleMatrix((int) 1,(int) 1,mxREAL);
  plhs[5]=mxCreateDoubleMatrix((int) 1,(int) 1,mxREAL);


  vgamma  = mxGetPr(plhs[0]); 
  vdelta  = mxGetPr(plhs[1]); 
  vxlam   = mxGetPr(plhs[2]); 
  vxi     = mxGetPr(plhs[3]); 
  vitype  = mxGetPr(plhs[4]);
  vifault = mxGetPr(plhs[5]);


/*Call the main computational routine*/
  hhh(xbar,sd,rb1,bb2,&itype,&gamma,&delta,&xlam,&xi,&ifault);

  *vgamma=gamma; *vdelta=delta; *vxlam=xlam; *vxi=xi; *vifault=ifault; *vitype=itype;

}
