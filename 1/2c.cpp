#include<cmath>
#include<iostream>
#include<vector>
#include <fstream>
#include"rk4.h"
using namespace std;
//constantes del problema//
const double h=0.1;//tamaño de paso
const double lf=15.0;//lambda final
const double rf=1.0;
const int size =2;//cuantas variables son independientes, y & R
const double ERR=1e-7;
//funciones
string loop(int i);
double Cero(double a, double b,double r,double h,std::vector<double>&y);
double function(double x,double r,double h,std::vector<double>&y);

double Cero(double a,double b,double r,double h,std::vector<double>&y){
  double m,fa,fm;
  fa=function(a,r,h,y);
  while(b-a>ERR){
    m=(a+b)/2; fm=function(m,r,h,y); //calculo m;
    if(fa*fm<0)  //  if(f(a) y f(m) son de signo contrario)
      b=m;
    else
      {a=m; fa=fm;}
  }
  return (a+b)/2;
}

int main(void)
{
 
  std::vector<double> y(size);
  //condiciones iniciales
  y[0]=1;//R(0)
  y[1]=0;//y(0), y=R'
  double r,ce,a,b;
  int count=1;
  for(double lambda=0.1;lambda<lf;lambda+=h)
    {
      a=lambda;
      b=lambda+h;
      ce=Cero(a,b,r,h,y);
      std::ofstream file(loop(count));
      for(double r=0.01;r<rf;r+=h/10){ file <<r << '\t'<<y[0]<<std::endl;
	  rk4(r,h,y,ce);
	    
      }
       file.close();
      
      count+=1;
   
      
    }

  std::cout << "set xrange[0.1:0.3]" << std::endl;
  std::cout << "set yrange[-0.1:0.1]" << std::endl;
  std::cout << "plot for[i=1:9] 'ans_'.i.'.dat'  using 1:2 w l title 'modo'.i " << std::endl;
  std::cout << "pause 10" << std::endl;
  return 0;
  
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
