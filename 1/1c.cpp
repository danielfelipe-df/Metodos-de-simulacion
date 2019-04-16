#include<iostream>
#include<cmath>
#include<vector>
#include <fstream>
//constantes del problema//
const double g=0.08;
const double h=0.1;
const double tf=1000;
const int size =2;
//funciones
void rk4(double t,double h,std::vector<double>&y,double beta);
double f(double t, const std::vector<double> & y, int id,double beta);
int main(void)
{
 std::ofstream fout2("datos_c.dat");
  std::vector<double> y(size);
  //condiciones iniciales
  y[0]=0.999;
  y[1]=0.001;
  double beta=0.1;
  double beta_max = 0.4;
  for(beta=0.1; beta<beta_max; beta+=0.005){
    y[0]=0.999;    y[1]=0.001;
    std::ofstream fout2("datos_c.dat");
    for(double t=0; t<1000; t+=h){
      rk4(t,h,y,beta);
    }
    fout2 << beta/g << " " << y[0] << std::endl;
  }
  fout2.close();
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula.gif'" << std::endl;
  std::cout << "unset key" << std::endl;
  std::cout << "set xrange[0:1000]" << std::endl;
  std::cout << "set yrange[0:1]" << std::endl;
  std::cout << "plot 'datos_c.dat' using 1:2 w l " << std::endl;
  std::cout << "pause 3" << std::endl;
 
  return 0;
  
}
void rk4(double t,double h,std::vector<double>&y,double beta)
{
  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());
  // k1
  for(int i=0; i< y.size(); i++){
    k1[i]=h*f(t,y,i,beta);
  }
  //k2 aux
  for(int i=0; i< y.size() ;i++){
    aux[i]=y[i]+k1[i]/2;
  }
  //k2
  for(int i=0; i< y.size();i++){
    k2[i]= h*f(t+h/2,aux,i,beta);
  }
  //k3 aux
  for(int i=0; i< y.size(); i++){
    aux[i]=y[i]+k2[i]/2;
  }
  //k3
  for(int i=0; i< y.size(); i++){
    k3[i]=h*f(t+h/2,aux,i,beta);
  }
  //k4 aux
  for(int i=0; i< y.size(); i++){
    aux[i]=y[i]+k3[i];
  }
  //k4
  for(int i=0; i< y.size(); i++){
    k4[i]=h*f(t+h,aux,i,beta);
  }
  //runge kutta
  for(int i=0; i<y.size();i++)
    {
      
      y[i]+= (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    
    }
}
double f(double t, const std::vector<double> & y,int id,double beta)
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
