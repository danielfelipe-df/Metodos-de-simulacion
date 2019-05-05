#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
#include <cstdlib>
#include <time.h>
using namespace std;
const double k=1.0e4;
const int Nx=5,Ny=5,N=Nx*Ny;
const double Lx=60;
const double Ly=60;
const double KT=10;
const double c1=0.1786178958448091;
const double c2=-0.2123418310626054;
const double c3=-0.06626458266981849;
const double umdc2=1-2*c2;
const double umdc3=(1-2*(c3+c1));
const double E=1.0;
const double r0=10.0;
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
  double aux=12.0*E*pow(norma2(dr),-1)*((pow(r0,12)*pow(norma2(dr),-6))-(pow(r0,6)*pow(norma2(dr),-3)));
 
    Molecula2.AgregueFuerza( dr* aux);
    Molecula1.AgregueFuerza(dr*(-1)*aux);
 

}

void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i ,j;  double h;vector3D Punto1,Punto2,F,a,b;
  Punto1.cargue(0,120,0);Punto2.cargue(60,0,0);

  for(i=0;i<(N);i++)Molecula[i].BorreFuerza();
  for(i=0;i<N;i++){
    	a=Punto1-Molecula[i].r;
	if(abs(a.x())-Molecula[i*Ny+j].R<=0)
	  {
	    
	    h=Molecula[i].R-abs(a.x());
	    F.cargue(k*pow(h,1.5),0,0);
	    Molecula[i].AgregueFuerza(F);
	  }
	else if(abs(a.y())-Molecula[i].R<=0)
	  {
	   
	    h=Molecula[i].R-abs(a.y());
	    F.cargue(0,-k*pow(h,1.5),0);
	    	Molecula[i].AgregueFuerza(F);
	  }
	else
	  {
	    b=Punto2-Molecula[i].r;
	    if(abs(b.x())-Molecula[i*Ny+j].R<=0)
	  {
	    
	     h=Molecula[i].R-abs(b.x());
	    F.cargue(-k*pow(h,1.5),0,0);
	    	Molecula[i].AgregueFuerza(F);
	  }
	else if(abs(b.y())-Molecula[i].R<=0)
	  {
	     h=Molecula[i].R-abs(b.y());
	    F.cargue(0,k*pow(h,1.5),0);
	    	Molecula[i].AgregueFuerza(F);
	  }
        
	    
	  }
	
      }
    for(j=i+1;j<(N);j++)
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
}
//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<2*Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";//pared de abajo
    cout<<" , "<<Ly/7<<"*t,"<<2*Ly;//pared de arriba 
    cout<<", 0,"<<2*Ly/7<<"*t";//pared de la izquierda
    cout<<" , "<<Lx<<","<<2*Ly/7<<"*t";//pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

int main(void){
  clock_t tStart = clock();
  Cuerpo Molecula[N];
  Crandom ran64(10);
  Colisionador Newton;
  double t,tdibujo,dt=1.0e-3;
  double m0=1,R0=2.5;
  int i,j;
  
  //Paredes

  
  //Moleculas
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1),x0,y0, theta;
  double V0=sqrt(2*KT/m0);
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      { x0=dx*(i+1); y0=(j+1)*dy;theta=2*M_PI*ran64.r();
	Molecula[i*Ny+j].Inicie(x0,y0,0,V0*cos(theta),V0*sin(theta),0,m0,R0);
	 
      }
  //------------(x0,y0,z0,Vx0,Vy0  ,Vz0,  m0, R0)
  double T=400;
  double sum,sum_old;
  int c=0;
  ofstream file("dat.dat");
  for(t=tdibujo=0;t<T;t+=dt,tdibujo+=dt){
    sum_old=sum;
    for(i=0;i<N;i++)
      	sum+=Molecula[i].Gety()/N;

    file<<t<<" "<<sum<<endl;
    sum=0;
     
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
   
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    file.close();

  return 0;
}
