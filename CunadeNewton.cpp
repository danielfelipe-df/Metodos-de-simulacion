#include<iostream>
#include<fstream> 
#include <cmath>
#include "Vector.h"
using namespace std; 

const double g=9.8;
const double K=1e10;

const double Zeta=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Xi=-0.06626458266981849;
const int N = 3;

class Cuerpo;
class Colisionador;


//___________________________Clase Cuerpo_____________________________
class Cuerpo{
private:
  double tau, omega, theta, m, R, L, I, xcorrido; 
public:
  void Inicie(double theta0, double omega0,double m0, double R0, double L0, double x0corrido);
  void InicieFuerza(void);
  void AgregueFuerza(double tau0);
  void Mueva_r(double dt, double Constante);
  void Mueva_V(double dt, double Constante);
  void Dibujese(void);
  double x(void){return xcorrido+L*sin(theta);};
  double y(void){return -L*cos(theta);};
  //double Gettau(void){return -L*cos(th;};
  
  friend class Colisionador; 
}; 

void Cuerpo::Inicie(double theta0, double omega0, double m0, double R0, double L0, double x0corrido){
  theta=theta0; omega=omega0; m=m0; R=R0; L=L0; xcorrido=x0corrido; I=m*L*L;
}


void Cuerpo::InicieFuerza(void){
  tau=-m*g*L*sin(theta);
}

void Cuerpo::Mueva_r(double dt, double Constante){
  theta+=omega*(Constante*dt);
}

void Cuerpo::Mueva_V(double dt, double Constante){
  omega+=tau*(Constante*dt/I);
}

void Cuerpo::Dibujese(void){
  cout<<", "<<xcorrido+L*sin(theta)<<"+"<<R<<"*cos(t),"<<-L*cos(theta)<<"+"<<R<<"*sin(t)," << xcorrido<< "+" << L*sin(theta)/7.0 << "*t" << "," << -L*cos(theta)/7.0 << "*t";
}

void Cuerpo::AgregueFuerza(double tau0){
  tau+=tau0;
}

//------------Clase Colisionador----------
class Colisionador{
private:
 
public:
  void CalculeTodasLasFuerzas(Cuerpo* Pendulo );
  void CalculeLaFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);
  
};

void Colisionador::CalculeTodasLasFuerzas(Cuerpo* Pendulo){
  int i,j;
  //Borrar todas las fuerzas
  for(i=0; i<N;i++) Pendulo[i].InicieFuerza();
  //Agregue fuerzas externas
  for(i=0; i<N; i++) Pendulo[i].AgregueFuerza(-Pendulo[i].m*g*sin(Pendulo[i].theta));
  //Calcular todas las fuerzas entre parejas de pendulos
  for(i=0; i<N;i++)
    for(j=i+1; j<N;j++)
      CalculeLaFuerzaEntre(Pendulo[i], Pendulo[j]);
  
}

void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  double s, d21, F2, tau1, tau2;
  d21= Pendulo2.x() -Pendulo1.x();
  s= (Pendulo2.R+Pendulo1.R) - d21;
  if(s>0){
    F2= K*pow(s,1.5)*d21/abs(d21); tau2 =F2*Pendulo2.L ; tau1= -F2*Pendulo1.L;
    Pendulo2.AgregueFuerza(tau2);
    Pendulo1.AgregueFuerza(tau1);
  }
  
}

//________________________Funciones Globales_________________________

void InicieAnimacion(void){
  // cout<<"set terminal gif animate"<<endl;
  // cout<<"set output 'MiPendulo.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-12:22]"<<endl;
  cout<<"set yrange[-12:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
}

void InicieCuadro(void){
  cout<<"plot 0,0 ";
}

void TermineCuadro(void){
  cout<<endl;
}

int main(void)
{
  double t, dt=0.01;
  double tdibujo; int Ndibujos;
  Cuerpo Pendulo[N];
  int i;
  Colisionador Newton;

  double m0=100, L0=12;
  double R0= 1.5;
  double x0corrido=0, theta0= -15*M_PI/180;

  double T=2*M_PI*sqrt(L0/g), tmax=1*T;
  
  InicieAnimacion(); Ndibujos=1000;
  
  //variables de la funcion Inicie(theta0, omega0, m0, R0, L0, x0corrido)
                 Pendulo[0].Inicie(theta0,      0, m0, R0, L0, x0corrido);
  for(i=1; i< N; i++ ) Pendulo[i].Inicie(0, 0,m0,R0,L0,x0corrido+ 2*i*R0);
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo>2*tmax/Ndibujos){
      InicieCuadro();
      for(i=0; i<N;i++)   Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    // for(i=0;i<N; i++)    cout<<Pendulo[i].x()<<" "<<Pendulo[i].y()<<endl;

    //Muevase con el algoritmo Omelyan PEFRL
    for(i=0; i<N;i++)  Pendulo[i].Mueva_r(dt,Zeta);
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0; i<N;i++)  Pendulo[i].Mueva_V(dt,(1-2*Lambda)/2); for(i=0; i<N;i++) Pendulo[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0; i<N;i++) Pendulo[i].Mueva_V(dt,Lambda); for(i=0; i<N;i++) Pendulo[i].Mueva_r(dt,1-2*(Xi+Zeta));
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0; i<N;i++) Pendulo[i].Mueva_V(dt,Lambda); for(i=0; i<N;i++) Pendulo[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0; i<N;i++) Pendulo[i].Mueva_V(dt,(1-2*Lambda)/2);    for(i=0; i<N;i++) Pendulo[i].Mueva_r(dt,Zeta);
  }

}
