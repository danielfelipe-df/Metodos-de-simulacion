//Mi primer programa en C++
#include<iostream>
#include<cmath>
using namespace std;

const double ERR=1e-7;
const double lambda = 2.4048;

double f(int n,double x,double t){
  return cos(n*t-x*sin(t));
}

double IntegralPorSimpson(int n,double x,double a,double b,int N){
  int DosN=2*N; double h=(b-a)/DosN; int i; double t,suma;
  for(suma=0,i=0;i<=DosN;i++){
    t=a+h*i;
    if(i==0 || i==DosN)
      suma+=f(n,x,t);
    else if(i%2==1)
      suma+=4*f(n,x,t);
    else
      suma+=2*f(n,x,t);
  }  
  return suma*h/3;
}

double Bessel(int n,double x){
  return 1/M_PI*IntegralPorSimpson(n,x,0,M_PI,50);
}

int main(void){
  int n=0; double x;
  
  for(x=0;x<10;x+=0.1)
    cout<<x<<" "<<Bessel(n,lambda*x)<<endl;
  
  return 0;
}
