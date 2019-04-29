#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
#include"rk4.h"
//constantes del problema//
const double lambda =1;
const double h=0.001;//tamaño de paso
const double tf=10;//tiempo final
const int size =2;//cuantas variables son independientes, y & R

int main(void)
{
  
  std::vector<double> y(size);
  //condiciones iniciales
  y[0]=1;//R(0)
  y[1]=0;//y(0), y=R'
  std::ofstream file("2a.dat");
   std::cout<<"set xrange[0:10]"<<std::endl;
   std::cout<<"set yrange[-0.6:1]"<<std::endl;
  for(double r=0.01;r<tf;r+=h)
    {  file<<r<<" "<<y[0]<<std::endl;
      rk4(r,h,y,lambda);
      }
  std::cout<<"plot \"2a.dat\" w lines,'datos_2a.dat' w l "<<std::endl;
   std::cout << "pause 10" << std::endl;
  
  return 0;
  
}

