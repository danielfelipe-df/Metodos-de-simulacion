//Mi Primer Programa en C++
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;

const double k=1.0e4;
const double Lx=60,Ly=120;
const int Nx=5,Ny=5, N=Nx*Ny;
const double E = 1.0;
const double r0 = 10.0;
const double R1 = 2.5;


const double Zi=0.1786178958448091e0;
const double Lambda=0.2123418310626054*(-1);
const double Xi=0.06626458266981849*(-1);

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=(1-2*(Xi+Zi));

//------------ Declaracion de las clases --------
class Cuerpo;
class Colisionador;
//------------ Clase Cuerpo --------
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
  double GetVx(void){return V.x();}; //Inline
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
  V+=F*(dt*Coeficiente/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------------------ Clase Colisionador -----------------
class Colisionador{
private:
public:
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void CalculeTodasLasFuerzas(Cuerpo * Molecula);
};
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2){
  vector3D dr=Molecula1.r - Molecula2.r;
  double dr2=norma2(dr);
  double aux=12.0*E*pow(dr2,-1)*((pow(r0,12)*pow(dr2,-6))-(pow(r0,6)*pow(dr2,-3)));
  vector3D F2=aux*dr;
  Molecula2.AgregueFuerza((-1)*F2);  Molecula1.AgregueFuerza(F2);
  
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i,j;
 
  for(i=0;i<N;i++) Molecula[i].BorreFuerza();
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
     CalculeFuerzaEntre(Molecula[i],Molecula[j]);
    }
  }
  
  //Calcular Fuerza Con las Paredes

  vector3D S1, S2, S3, S4, F3, P13, P24;
  double h, FP;
  
  P13.cargue(0, Ly, 0);
  P24.cargue(Lx, 0, 0);
  
  for (i=0; i<N; i++){
       
  S1 = P24 - Molecula[i].r;
  S2 = (-1)*S1;
  
  S3 = P13 - Molecula[i].r;
  S4 = (-1)*S3;
  
  
    
  if(S1.x() <= R1){
    h = R1-S1.x();
    FP = -k*pow(h,1.5);
    
    F3.cargue(FP,0,0);}
  
  else if (S2.y() <= R1 ){
    h = R1-S2.y();
    FP = k*pow(h,1.5);
    
    F3.cargue(0,FP,0); }
  
  else  if (S3.y() <= R1 ){
    h = R1-S3.y();
    FP = -k*pow(h,1.5);
    
    F3.cargue(0,FP,0); }
  
  else   if  (S4.x() <= R1 ){
    h =R1-S4.x();
    FP = k*pow(h,1.5);
    
    F3.cargue(FP,0,0); }
  
  else {F3.cargue(0,0,0);}
  
  Molecula[i].AgregueFuerza(F3);
  }
    
  }
  //------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+40<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+40<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

int main(void){
  Cuerpo Molecula[N];
  Colisionador Newton;
  Crandom ran64(1);
  double t,tdibujo,dt=1.0e-3;
  double m0=1.0,R0=2.5;
  int i,j,n;
  double kT=10.0;
  double VelocidadParticulaenteq[N];

  // InicieAnimacion();

  //--------------INICIALIZACION----------------------------

  

  //MOLECULAS
  double dx=Lx/(Nx+1), dy= 0.5*Ly/(Ny+1), x0,y0;
  double V0=sqrt(2*kT/m0), theta,Vx0,Vy0;
  
  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++){
      x0=(i+1)*dx; y0=(j+1)*dy; 
      theta=2*M_PI*ran64.r(); Vx0=V0*cos(theta); Vy0=V0*sin(theta);
      //---------------------(x0,y0,z0,Vx0,Vy0,Vz0, m0, R0)
      Molecula[i*Ny+j].Inicie(x0,y0, 0,Vx0,Vy0,  0, m0, R0);
    }

  double teq=55, T=200;
  for(t=0;t<=teq;t+=dt){
   
    //Moverlo Segun PEFRL Orden 4
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Zi);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Coeficiente2);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Zi);
    
    if(t<=(teq+dt) && t>=(teq-dt)){
      for(i=0; i<N; i++){
	
	VelocidadParticulaenteq[i]= Molecula[i].GetVx();
      }
    }
  } 

  double Total;
  double Promedio;
  for (i=0; i<N; i++){
    Total +=  VelocidadParticulaenteq[i];
  }
  
  Promedio = Total/N;
  
  //Desviación estandar

  double sigma;
  double dif2, suma2;
  //Sumas al cuadrado

  for(i=0; i<N; i++){
  dif2 =pow(VelocidadParticulaenteq[i]-Promedio, 2);
  suma2 +=dif2;
}

  sigma = pow((suma2/N),0.5);
  
  cout<<sigma<<endl;

  cout<<"Error porcentual"<<endl;
  cout<< (sigma-3.16228)*100.0/3.16228; 
  
  //------------IMPRIMIR RESULTADOS-------------------------
  return 0;
}
