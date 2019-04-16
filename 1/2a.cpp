#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
#include"rk4.h"
//constantes del problema//
const lambda =1;
const double h=0.1;//tamaño de paso
const double tf=1000;//tiempo final
const int size =2;//cuantas variables son independientes, y & R

int main(void)
{
  std::ofstream file("data.dat");
  std::vector<double> y(size);
  //condiciones iniciales
  y[0]=1;//R(0)
  y[1]=0;//y(0), y=R'
  for(double r=0.01;r<tf;r+=h)
    { fout1 << r << " " << y[0]<<'\n';
      rk4(r,h,y,lambda);
      }
  file.close();
  std::cout << "set xrange[0.01:10]" << std::endl;
  std::cout << "set yrange[-1:1]" << std::endl;
  std::cout << "plot 'data.dat' using 1:2 w l title 'R(r)' " << std::endl;
  std::cout << "pause 10" << std::endl;
  
  return 0;
  
}

