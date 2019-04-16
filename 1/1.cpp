#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
//constantes del problema//
const double beta=0.35;
const double g=0.08;
const double h=0.1;
const double tf=1000;
const int size =2;
//funciones
void rk4(double t,double h,std::vector<double>&y);
double f(double t, const std::vector<double> & y, int id);

int main(void)
{

 std::vector<double> y(size);//se inicia el vector
  //condiciones iniciales
  y[0]=0.999;
  y[1]=0.001;
   std::ofstream file("1.dat");
  for(double t=0;t<tf;t+=h)
    {
      file << t << '\t' << y[0]<< '\t' << y[1]<< " " << 1-y[0]-y[1] << std::endl;      rk4(t,h,y);
      }
  file.close();
  //gnuplot
  std::cout << "set xrange[0:1000]" << std::endl;
  std::cout << "set yrange[0:1]" << std::endl;
  std::cout << "plot '1.dat' u 1:2 w l, '1.dat' u 1:3 w l, '1.dat' u 1:4 w l " << std::endl;
  std::cout << "pause 3" << std::endl;
  
  return 0;
  
}

//runge kutta
void rk4(double t,double h,std::vector<double>&y)
{
  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());
  // k1
  for(int i=0; i< y.size(); i++){
    k1[i]=h*f(t,y,i);
  }
  //k2 aux
  for(int i=0; i< y.size() ;i++){
    aux[i]=y[i]+k1[i]/2;
  }
  //k2
  for(int i=0; i< y.size();i++){
    k2[i]= h*f(t+h/2,aux,i);
  }
  //k3 aux
  for(int i=0; i< y.size(); i++){
    aux[i]=y[i]+k2[i]/2;
  }
  //k3
  for(int i=0; i< y.size(); i++){
    k3[i]=h*f(t+h/2,aux,i);
  }
  //k4 aux
  for(int i=0; i< y.size(); i++){
    aux[i]=y[i]+k3[i];
  }
  //k4
  for(int i=0; i< y.size(); i++){
    k4[i]=h*f(t+h,aux,i);
  }
  //runge kutta
  for(int i=0; i<y.size();i++)
    {
      
      y[i]+= (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    
    }
}
//ec dif.
double f(double t, const std::vector<double> & y,int id)
{
  if(0==id)
    {
      return -beta*y[0]*y[1];
    }
  else if (1==id)
    {
      return (beta*y[0]*y[1])-(g*y[1]);
    }
  else{
    std::cerr<<"Error!Id does not exist->"<<id<<std::endl;
    exit(1);
  }
}
