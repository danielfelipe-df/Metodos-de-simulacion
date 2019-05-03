#include<iostream>
#include<fstream>
#include<cmath>
#include<tuple>
using namespace std;

const double ERR=1e-7;
//const double lamda=1.0, a=1.0;

double f1(double R1 ,double R2 ,double r){ return R2/r;} //f1 =dR1/dr=R2/r
double f2(double R1 ,double R2 ,double r, double lamda){ return -lamda*lamda*r*R1;} //f2=dR2/dr=-lamda**2*r*R1
double f(double x){
  double r, R1, R2, dr=0.1, lamda;
  for(r=0.01,R1=1.0,R2=0;r<1.0;){
      UnPasoDeRungeKutta4(R1,R2,r,dr, lamda);
      if(r==x){
	break;
      }
  }
  return R1
}

void UnPasoDeRungeKutta4(double & R1,double & R2,double & r,double dr, double lamda){
  double dR11,dR12,dR13,dR14;   double dR21,dR22,dR23,dR24;
  dR11=f1(R1,R2,r)*dr;                            dR21=f2(R1,R2,r,lamda)*dr;
  dR12=f1(R1+0.5*dR11, R2+0.5*dR21, r+0.5*dr)*dr;   dR22=f2(R1+0.5*dR11, R2+0.5*dR21, r+0.5*dr,lamda)*dr;
  dR13=f1(R1+0.5*dR12, R2+0.5*dR22, r+0.5*dr)*dr;   dR23=f2(R1+0.5*dR12, R2+0.5*dR22, r+0.5*dr, lamda)*dr;
  dR14=f1(R1+dR13, R2+dR23, r+dr)*dr;               dR24=f2(R1+dR13, R2+dR23, r+dr,lamda)*dr;

  r+=dr;  R1+=(dR11+2*(dR12+dR13)+dR14)/6;          R2+=(dR21+2*(dR22+dR023)+dR24)/6; 
} 

//----------Método de Bisección----------
double CeroPorBiseccion(double a,double b, double lamda){
  double m,fa,fm, lamda;
  
  fa=f(a,lamda);
  while(b-a>ERR){
    m=(a+b)/2; fm=f(m,lamda); //calculo m;
    if(fa*fm<0)  //  if(f(a) y f(m) son de signo contrario)
      b=m;
    else
      {a=m; fa=fm;}
  }
  return (a+b)/2;
}



int main(void){
  double r,R1,R2,dr=0.1, lamda; int j;
  double a=1,b=4;
  ofstream f; //    float funcion[2][2981];
  f.open("datos2b.dat");
  for(lamda=0.1, j=0; lamda<=15.0; lamda=lamda+0.005, j++){
    for(r=0.01,R1=1.0,R2=0;r<1.0;){
      UnPasoDeRungeKutta4(R1,R2,r,dr, lamda);
      }
    f<<lamda<<" "<<R1<<" "<<endl;
    //funcion[0][j]=lamda;     funcion[1][j]=R1;     
  }
  f.close();
  cout<<"f(x)=0 en x="<<CeroPorBiseccion(a, b, lamda)<<endl;

  //---------------Configurando Imagen------------
  cout<<"set key at graph .92, .92 spacing 2 font 'Helvetica, 14'"<<endl;
  cout<<"set key box lt -1 lw 2"<<endl;
  cout<<"set xlabel '{/Symbol l}' "<<endl;
  cout<<"set ylabel 'f({/Symbol l})' "<<endl;
  //cout<<"set title ''"<<endl;
  cout<<"plot \'datos2b.dat\' w l title 'f({/Symbol l};r=1)"<<endl;
  //cout<<"pause 15"<<endl;
  //---------------Creando Imagen------------
  cout<<"set term png"<<endl;
  cout<<"set output 'grafica2b.png'"<<endl;
  cout<<"replot"<<endl;
  cout<<"set term x11"<<endl;
  return 0;
}
