#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

const double ERR=1e-7;
const double omega=1.0,betha=0.35,gama=0.08;

double f1(double s ,double i ,double t){
  return -betha*s*i;}

double f2(double s ,double i ,double t){
  return betha*s*i-gama*i;}

void UnPasoDeRungeKutta4(double & s,double & i,double & t, double dt){
  double ds1,ds2,ds3,ds4; 
  double di1,di2,di3,di4;
  ds1=f1(s,i,t)*dt;                            di1=f2(s,i,t)*dt;
  ds2=f1(s+0.5*ds1, i+0.5*di1, t+0.5*dt)*dt;   di2=f2(s+0.5*ds1, i+0.5*di1, t+0.5*dt)*dt;
  ds3=f1(s+0.5*ds2, i+0.5*di2, t+0.5*dt)*dt;   di3=f2(s+0.5*ds2, i+0.5*di2, t+0.5*dt)*dt;
  ds4=f1(s+ds3, i+di3, t+dt)*dt;               di4=f2(s+ds3, i+di3, t+dt)*dt;

  t+=dt;  s+=(ds1+2*(ds2+ds3)+ds4)/6;          i+=(di1+2*(di2+di3)+di4)/6;
} 


int main(void){
  double sin,t,s,i,dt=0.1;

  ofstream f;
  f.open ("datos1a.dat");
  
  for(t=0,s=0.999,i=0.001;t<120; ){
    //cout<<t<<" "<<s<<" "<<i<<" "<<1.0-s-i<<endl;
    f << t<<" "<<s<<" "<<i<<" "<<1.0-s-i<<"\n";
    UnPasoDeRungeKutta4(s,i,t,dt);  
  }

  f.close();
  //---------------Configurando Imagen------------
  cout<<"set key at graph .92, .92 spacing 2 font 'Helvetica, 14'"<<endl;
  cout<<"set key box lt -1 lw 2"<<endl;
  cout<<"set xlabel 'Tiempo (s)' "<<endl;
  cout<<"set ylabel 'Proporcion' "<<endl;
  cout<<"set title 'MODELO SIR'"<<endl;
  cout<<"plot \'datos1a.dat\' u 1:2 title 'Suceptibles' w l,\'datos1a.dat\' u 1:3  w l title 'Infectados',\'datos1a.dat\' u 1:4 w l title 'Retirados'"<<endl;
  //---------------Creando Imagen------------
  cout<<"set term png"<<endl;
  cout<<"set output 'grafica1-a.png'"<<endl;
  cout<<"replot"<<endl;
  cout<<"set term x11"<<endl;
  return 0;
}
