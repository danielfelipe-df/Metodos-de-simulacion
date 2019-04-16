#include<iostream>
#include<cmath>
#include"Bessel.h"
#include"rk4.h"
#include<fstream>
#include <cstdio>
const double tf=1000;
const double lf=15.0;
const double h =0.1;
const double rf=1.0;
const int size =2;
double Cero2(double a, double b,double r,double h,std::vector<double>&y);
double function(double x,double r,double h,std::vector<double>&y);
double Cero(double a,double b,int n);
string loop(int i);
int main()
{  std::vector<double> y(size);
  //condiciones iniciales
  y[0]=0.999;
  y[1]=0.001;
  int count =0;
  int n=0;
  double a,b,ce,ce2,r;
  for(double lambda=0.1;lambda<lf;lambda+=h)
    {
      a=lambda;
      b=lambda+h;
      ce=Cero2(a,b,r,h,y);
      ce2=Cero(a,b,n);
      
      /*std::ofstream fu("tabla.dat");
      fu<<ce<<'\t'<<ce2<<std::endl;
      fu.close();
      
      */
    
   for(double r=0.01;r<rf;r+=h/1.0){ 
     file<<r<<'\t'<<Bessel(n,r*lambda*lambda)<<endl;
	}
	count+=1;
	file.close();}

    std::cout << "set xrange[0.1:1.0]" << std::endl;
  std::cout << "set yrange[-1:1]" << std::endl;
  std::cout << "plot for[i=1:9] 'ans_'.i.'.dat'  using 1:2 w l title 'modo'.i " << std::endl;
  std::cout << "pause 10" << std::endl;
  return 0;
      return 0;     
}
double Cero(double a,double b,int n){
  double m,fa,fm;
  fa=Bessel(n,a);
  while(b-a>ERR){
    m=(a+b)/2; fm=Bessel(n,m); //calculo m;
    if(fa*fm<0)  //  if(f(a) y f(m) son de signo contrario)
      b=m;
    else
      {a=m; fa=fm;}
  }
  return (a+b)/2;
}
double Cero2(double a,double b,double r,double h,std::vector<double>&y){
  double m,fa,fm;
  fa=function(a,r,h,y);
  while(b-a>ERR){
    m=(a+b)/2; fm=function(m,r,h,y); //calculo m;
    if(fa*fm<0)  //  if(f(a) y f(m) son de signo contrario)
      b=m;
    else
      {a=m; fa=fm;}
  }
}

 
string loop(int i)
{ 
  char c=i;
  string s="ans_";
  s+=c;
  s+=".dat";
  return s;
}
 double function(double x,double r,double h,std::vector<double>&y)
{
  
  for(double lambda=0.1;lambda<x;lambda+=h)
    rk4(r,h,y,lambda);
  y[0]=1;//R(0)
  y[1]=0;//y(0), y=R'
  r=1;
  return y[0];
}
