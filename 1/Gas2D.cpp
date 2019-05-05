//Mi Primer Programa en C++
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;
const double K=1.0e4;
const int Nx=5,Ny=5,N=Nx*Ny;
const double Lx=60;
const double Ly=60;
const double KT=10;
const double c1=0.1786178958448091;
const double c2=-0.2123418310626054;
const double c3=-0.06626458266981849;
const double umdc2=1-2*c2;
const double umdc3=(1-2*(c3+c1));

// Declaración de clases

class Cuerpo;
class Colisionador;
//clase cuerpo
class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AgregueFuerza(vector3D dF){F+=dF;}; 
  void Mueva_r(double dt,double Coeficiente);
  void Mueva_V(double dt,double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  double GetVx(void){return V.x();};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double Coeficiente){
  r+=V*(dt*Coeficiente);
}
void Cuerpo::Mueva_V(double dt,double Coeficiente){
  V+=F*(dt*Coeficiente)/m;
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//-------------------Clase colisionador-----------------
class Colisionador{
private:
  
public:
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void CalculeTodasLasFuerzas(Cuerpo * Molecula);
};
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2){
  vector3D dr= Molecula2.r-Molecula1.r;
  double norma_dr=norma(dr);
  double s=(Molecula2.R+Molecula1.R)-norma_dr;
  if(s>0){
    vector3D F2=dr*(K*pow(s,1.5)/norma_dr);
    Molecula2.AgregueFuerza(F2);
    Molecula1.AgregueFuerza(F2*(-1));
  }
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i ,j;
  for(i=0;i<(N+4);i++)Molecula[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=i+1;j<(N+4);j++)
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
}
//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";//pared de abajo
    cout<<" , "<<Ly/7<<"*t,"<<Ly;//pared de arriba 
    cout<<", 0,"<<Ly/7<<"*t";//pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t";//pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

int main(void){
  Cuerpo Molecula[N+4];
  Crandom ran64(1);
  Colisionador Newton;
  double t,tdibujo,dt=1.0e-3;
  double m0=1,R0=3;
  int i,j;
  InicieAnimacion();
  
  //Paredes
  double Mpared=100*m0,Rpared=10000;
  //------(x0,y0,z0,vx0,vy0,vz0,m0,R0)M
  //Pared Arriba
  Molecula[N].Inicie(Lx/2,Ly+Rpared,0,0,0,0,Mpared,Rpared);
  //Pared Abajo
  Molecula[N+1].Inicie(Lx/2,-Rpared,0,0,0,0,Mpared,Rpared);
  //pared derecha
  Molecula[N+2].Inicie(Lx+Rpared,Ly/2,0,0,0,0,Mpared,Rpared);
  //Pared Izquiera
  Molecula[N+3].Inicie(-Rpared,Ly/2,0,0,0,0,Mpared,Rpared);

  
  //Moleculas
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1),x0,y0, theta;
  double V0=sqrt(2*KT/m0);
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      { x0=dx*(i+1); y0=(j+1)*dy;theta=2*M_PI*ran64.r();
	Molecula[i*Ny+j].Inicie(x0,y0,0,V0*cos(theta),V0*sin(theta),0,m0,R0);
	 
      }
  //------------(x0,y0,z0,Vx0,Vy0  ,Vz0,  m0, R0)
  double T=100;
  for(t=tdibujo=0;t<T;t+=dt,tdibujo+=dt){
    //    cout<<t<<" "<<Molecula[0].Getx()<<" "<<Molecula[0].Gety()<<endl;
    
    if(tdibujo>T/1000){
      InicieCuadro();
      for(int i =0; i<N;i++) Molecula[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    for(int i=0;i<N;i++)Molecula[i].Mueva_r(dt,c1);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(int i=0;i<N;i++)Molecula[i].Mueva_V(dt,umdc2);
    
    for(int i=0;i<N;i++)Molecula[i].Mueva_r(dt,c3);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(int i=0;i<N;i++)Molecula[i].Mueva_V(dt,c2);
    
    for(int i=0;i<N;i++)Molecula[i].Mueva_r(dt,umdc3);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(int i=0;i<N;i++)Molecula[i].Mueva_V(dt,c2);
   

    for(int i=0;i<N;i++)Molecula[i].Mueva_r(dt,c3);

    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(int i=0;i<N;i++)Molecula[i].Mueva_V(dt,umdc2);

    for(int i=0;i<N;i++)Molecula[i].Mueva_r(dt,c1);
    
  }
  
  //----------IMPRIMIR RESULTADO-----------------------
  for(i=0;i<Nx;i++)
    {
      //    cout<<Molecula[N].GetVx()<<endl;
    }

  return 0;
}
