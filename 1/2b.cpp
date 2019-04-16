#include<iostream>
#include"rk4.h"
#include <fstream>
//constantes del problema//
const double h=0.1;//tamaño de paso
const double lf=15.0;//lambda final
const int size =2;//cuantas variables son independientes, y & R
const double ERR=1e-7;
int main(void)
{
 
  std::vector<double> y(size);
  //condiciones iniciales
  y[0]=1;//R(0)
  y[1]=0;//y(0), y=R'
  double r=1;
  for(double lambda=0.1;lambda<lf;lambda+=h)
    { rk4(r,h,y,lambda);
      if (y[0]<ERR)
	std::cout<<lambda<<std::endl;
    }
  return 0;
  
}
