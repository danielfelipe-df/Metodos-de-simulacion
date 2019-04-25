//Mi primer programa en C++
#include<iostream>
#include<cmath>
using namespace std;

const double ERR=1e-7;

double f(double x){
  return sin(x)/x;
}

double CeroPorBiseccion(double a,double b){
  double m,fa,fm;
  fa=f(a);
  while(b-a>ERR){
    m=(a+b)/2; fm=f(m); //calculo m;
    if(fa*fm<0)  //  if(f(a) y f(m) son de signo contrario)
      b=m;
    else
      {a=m; fa=fm;}
  }
  return (a+b)/2;
}

int main(void){
  double a,b;
  
  a=1; b=4;
  cout<<"f(x)=0 en x="<<CeroPorBiseccion(a,b)<<endl;

  return 0;
}
