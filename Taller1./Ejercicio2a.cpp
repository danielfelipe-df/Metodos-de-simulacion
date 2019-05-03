#include<iostream>
#include<fstream>
#include<cmath>
#include<tuple>
using namespace std;

const double ERR=1e-7;
const double lamda=1.0, a=1.0;

double f1(double R1 ,double R2 ,double r){ return R2/r;} //f1 =dR1/dr=R2/r
double f2(double R1 ,double R2 ,double r){ return -lamda*lamda*r*R1;} //f2=dR2/dr=-gama**2*r*R1

void UnPasoDeRungeKutta4(double & R1,double & R2,double & r,double dr){
  double dR11,dR12,dR13,dR14;   double dR21,dR22,dR23,dR24;
  dR11=f1(R1,R2,r)*dr;                            dR21=f2(R1,R2,r)*dr;
  dR12=f1(R1+0.5*dR11, R2+0.5*dR21, r+0.5*dr)*dr;   dR22=f2(R1+0.5*dR11, R2+0.5*dR21, r+0.5*dr)*dr;
  dR13=f1(R1+0.5*dR12, R2+0.5*dR22, r+0.5*dr)*dr;   dR23=f2(R1+0.5*dR12, R2+0.5*dR22, r+0.5*dr)*dr;
  dR14=f1(R1+dR13, R2+dR23, r+dr)*dr;               dR24=f2(R1+dR13, R2+dR23, r+dr)*dr;

  r+=dr;  R1+=(dR11+2*(dR12+dR13)+dR14)/6;          R2+=(dR21+2*(dR22+dR23)+dR24)/6; 
} 


int main(void){
  double r,R1,R2,dr=0.1;
  ofstream f,g; 
  f.open("datos2a.dat"); 
    for(r=0.01,R1=1.0,R2=0;r<=10.0;){
      f<<r<<" "<<R1<<" "<<endl;
      UnPasoDeRungeKutta4(R1,R2,r,dr);
      }
  f.close();
  
  //---------------Configurando Imagen------------
  cout<<"set key at graph .92, .92 spacing 2 font 'Helvetica, 14'"<<endl;
  cout<<"set key box lt -1 lw 2"<<endl;
  cout<<"set xlabel 'radio' "<<endl;
  cout<<"set ylabel 'R(r)' "<<endl;
  cout<<"set title 'ECUACION DE BESSEL'"<<endl;
  cout<<"plot \'datos2a.dat\' w l title 'Bessel"<<endl;
  //cout<<"plot \'datos2a-F.dat\' w l title 'F'"<<endl;
  //cout<<"pause 15"<<endl;
  //---------------Creando Imagen------------
  cout<<"set term png"<<endl;
  cout<<"set output 'grafica2a.png'"<<endl;
  cout<<"replot"<<endl;
  cout<<"set term x11"<<endl;
  return 0;
}
