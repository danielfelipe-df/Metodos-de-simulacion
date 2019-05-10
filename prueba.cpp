#include <fstream>
#include <cmath>
#include <iostream>
#include "Vector.h"
using namespace std;

const double g=9.8;
const int N=3;
const double K=0.1e10;

const double Zeta=0.1786178958448091;
const double Xi=-0.06626458266981849;
const double Lambda = -0.2123418310626054;
const double theta=1.0/(2-pow(2.0,1.0/3));

class Cuerpo ;
class Colisionador;

//---------------------------------------clase cuerpo
class Cuerpo{
private:
  double tau, omega, theta, m, R, L, I, xcorrido;
public:
  void Inicie(double theta0,double omega0, double m0,double R0,double L0,double x0corrido);
  void BorrarFuerza(void);
  void InicieFuerza(void);
  void AgregarFuerza(double tau0);
  void Mueva_theta(double dt,double Constante);
  void Mueva_omega(double dt,double Constante);
  void Dibujese();
  double x(void){return xcorrido+L*sin(theta);};
  double y(void){return -L*cos(theta);};
  
  friend class Colisionador;
};

void Cuerpo::Inicie(double theta0,double omega0, double m0,double R0,double L0,double x0corrido){
  omega=omega0;
  theta=theta0;
  m=m0;	R=R0; L=L0;
  xcorrido=x0corrido;
  I=m*L*L;
}

//void Cuerpo::BorrarFuerza(void){
//  tau=0;
//}
void Cuerpo::InicieFuerza(void){
  tau=-m*g*L*sin(theta);
}
void Cuerpo::AgregarFuerza(double tau0){
  tau+=tau0;
}

void Cuerpo::Mueva_theta(double dt,double Constante){
  theta+=omega*(Constante*dt);
}
void Cuerpo::Mueva_omega(double dt,double Constante){
  omega+=tau*(Constante*dt)/I;
}

void Cuerpo::Dibujese(void){
   cout<<", "<<xcorrido+L*sin(theta)<<"+"<<R<<"*cos(t),"<<-L*cos(theta)<<"+"<<R<<"*sin(t)," << xcorrido<< "+" << L*sin(theta)/7.0 << "*t" << "," << -L*cos(theta)/7.0 << "*t";   
}

//---------------------------------------clase colisionador

class Colisionador{
private:
  
public:
  void CalculeTodasLasFuerzas(Cuerpo* Pendulo);
  void CalculeLaFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);
};

void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Pendulo){
  int i,j;
  for(i=0; i<N; i++){
    Pendulo[i].InicieFuerza();
  }/*
  for(i=0; i<N; i++){
    Pendulo[i].AgregarFuerza(-Pendulo[i].m*g*Pendulo[i].L*sin(Pendulo[i].theta));
    }*/
  for(i=0;i<N;i++)
    for(j=i+1; j<N;j++){
      CalculeLaFuerzaEntre(Pendulo[i],Pendulo[j]);
    }
}

void Colisionador::CalculeLaFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  double s, d21, F2, tau1, tau2;
  d21= Pendulo2.x()-Pendulo1.x();
  s=(Pendulo2.R+Pendulo1.R) - d21;
  if(s>0){
    F2= K*pow(s,1.5); tau2 =F2*Pendulo2.L ; //tau1= -F2*Pendulo1.L;
    Pendulo2.AgregarFuerza(tau2);
    Pendulo1.AgregarFuerza(-tau2);
  }
}

//-----------------------------------Animación 
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'planeta.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-0.12:0.12]"<<endl;
  cout<<"set yrange [-0.15:0]"<<endl;
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

int main(void){
  double t, dt=0.001;
  Cuerpo Pendulo[N];
  Colisionador Newton;
  int i;
    
  double m0=0.1, L0=0.12, R0=0.015;
  double x0corrido=0, theta0=-15*M_PI/180, T= 2*M_PI*sqrt(L0/g), tmax=0.5*T;
  
  InicieAnimacion();
  
  Pendulo[0].Inicie(theta0, 0, m0, R0, L0, x0corrido);
  for(i=1;i<N;i++)
    Pendulo[i].Inicie(0, 0, m0, R0, L0,x0corrido+2*i*R0);
  
  
  double tdibujo; int Ndibujos=200;
  
  for (t=tdibujo=0; t<3*tmax; t+=dt, tdibujo+=dt){
    if (tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for (i=0; i<N; i++)
	Pendulo[i].Dibujese();
      TermineCuadro();
      std::cout << "pause 1" << std::endl;
      tdibujo=0;
    }
    
    //cout<<Pendulo[i].x()<<" "<<Pendulo[i].y()<<endl;
    
    //Muevase con Forest-Ruth
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,0.5*theta);//dt, Zeta
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,theta);//dt, (1-2*Lambda)/2
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(1-theta)/2);//dt, Xi 
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,1-2*theta);//dt, Lambda
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(1-theta)/2);//dt, 1-2*(Xi+Zeta)
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,theta);//dt, Lambda
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(theta)/2);//dt, Xi
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,theta);//dt, (1-2*Lambda)/2
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(theta)/2);//dt, Zeta
  }
  return 0;
}

