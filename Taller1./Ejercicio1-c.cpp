#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

const double ERR=1e-7;
const double omega=1.0,gama=0.08;

double f1(double s ,double i ,double t, double betha){
  return -betha*s*i;}

double f2(double s ,double i ,double t, double betha){
  return betha*s*i-gama*i;}

void UnPasoDeRungeKutta4(double & s,double & i,double & t, double dt, double betha){
  double ds1,ds2,ds3,ds4; 
  double di1,di2,di3,di4;
  ds1=f1(s,i,t,betha)*dt;                            di1=f2(s,i,t,betha)*dt;
  ds2=f1(s+0.5*ds1, i+0.5*di1, t+0.5*dt,betha)*dt;   di2=f2(s+0.5*ds1, i+0.5*di1, t+0.5*dt, betha)*dt;
  ds3=f1(s+0.5*ds2, i+0.5*di2, t+0.5*dt,betha)*dt;   di3=f2(s+0.5*ds2, i+0.5*di2, t+0.5*dt, betha)*dt;
  ds4=f1(s+ds3, i+di3, t+dt,betha)*dt;               di4=f2(s+ds3, i+di3, t+dt, betha)*dt;

  t+=dt;  s+=(ds1+2*(ds2+ds3)+ds4)/6;          i+=(di1+2*(di2+di3)+di4)/6;
} 

  
int main(void){
  double betha,t,s,i,dt=0.1; ;
  ofstream f;
  f.open ("datos1c.dat");
  for(betha=0.1;betha<0.4;betha+=0.0005){
    for(t=0,s=0.999,i=0.001;t<120;) {
      UnPasoDeRungeKutta4(s, i, t, dt, betha);
    }
    f<<betha/gama<<" "<<s<<endl;
  }
  f.close();
  //-----------Configurando grafica------
  cout<<"set key at graph .95, .95 spacing 1 font 'Helvetica, 14'"<<endl;
  cout<<"set key box lt -1 lw 2"<<endl;
  cout<<"set xlabel '{/Symbol b}/{/Symbol g}' "<<endl;
  cout<<"set ylabel 's_{/Symbol \245' "<<endl;
  cout << "plot 'datos1c.dat' u 1:2  w l title ' s_{/Symbol \245} vs {/Symbol b}/{/Symbol g}' " << endl;
    //---------------Creando Imagen------------
  cout<<"set term png"<<endl;
  cout<<"set output 'grafica1-c.png'"<<endl;
  cout<<"replot"<<endl;
  cout<<"set term x11"<<endl;
  return 0;
}
