#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
using namespace std;

const int N=3;//cantidad de planetas
const double G=1.0;
//Para PEFRL
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
  vector3D r,V,F; double m,R,omega,l;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AgregueFuerza(vector3D dF){F+=dF;};
  void Mueva_r(double dt,double Coeficiente);
  void Mueva_V(double dt,double Coeficiente);
  void DibujeseRotado(void);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
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
void Cuerpo::DibujeseRotado(void){
  cout<<" , "<<r.x()<<"*cos(omega*l)-"<<r.y()<<"*sin(omega*l)+"<<R<<"*cos(t), "<<r.x()<<"*sin(omega*l)+"<<r.y()<<"*cos(omega*l)+"<<R<<"*sin(t) , "<<1000<<"*cos(t),"<<1000<<"*sin(t)";
}
//------------------ Clase Colisionador -----------------
class Colisionador{
private: 
public:
  void CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);
  void CalculeTodasLasFuerzas(Cuerpo * Planeta,double omega);
};
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  vector3D dr=Planeta2.r-Planeta1.r;
  double aux=G*Planeta1.m*Planeta2.m*pow(norma2(dr),-1.5);
  vector3D F1=dr*aux;
  Planeta1.AgregueFuerza(F1);  Planeta2.AgregueFuerza(F1*(-1));
 
}
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Planeta,double omega){
  int i,j;
  for(i=0;i<N;i++) Planeta[i].BorreFuerza();
    
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]); 
  
}


int main(void){
  Cuerpo Planeta[N];
  Colisionador Newton;
  int i,c=0;
  double t,tdibujo,dt=50.0;
  double r=1000,m0=1047,m1=1,m2=0.005;
  double M,omega,T,x0,x1,V0,V1,l,x2,y2,Vx2, Vy2,x,x_old,x_old_old,l_old,sum;
  double R01=5,R02=5;
  double a=-M_PI/3,perturbacion=3e-3;
  M=m0+m1; omega=sqrt(G*M/(r*r*r)); T=2*M_PI/omega;
  x0=m1*r/M; x1=x0-r; V0=-omega*x0; V1=-omega*x1;
  x2=x0-r*cos(a);y2=-r*sin(a);
  Vx2 = omega*y2+perturbacion;   Vy2 = -omega*x2;
  //---------------(x0,y0,z0,Vx0,Vy0,Vz0, m0, R0)
  Planeta[0].Inicie(x0, 0, 0, 0, V0, 0, m0, R01);
  Planeta[1].Inicie(x1, 0, 0, 0, V1,  0, m1,  R02);
  Planeta[2].Inicie(x2, y2,0,Vx2, Vy2,   0, m2,  R02);
  std::ofstream fout("data.dat");
  std::ofstream fout2("puntos.dat");
   x= Planeta[2].Getx()*cos(omega*l)-Planeta[2].Gety()*sin(omega*l);
   x_old=x;
   x_old_old=x_old;
  for(l=0;l<30*T;l+=dt){
    x= Planeta[2].Getx()*cos(omega*l)-Planeta[2].Gety()*sin(omega*l);
    if(x_old<x && x_old<x_old_old){
      fout2 <<x<<"  "<< l <<'\n';
      if(c==0){l_old=l;c+=1;}
      else
	{ 
	  sum+=l-l_old;
	  l_old=l;
	  c+=1;
	}
    }
    else
       fout <<x<<"  "<< l <<'\n';
     
    //Moverlo Segun PEFRL Orden 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zi);
    Newton.CalculeTodasLasFuerzas(Planeta,omega);
    
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Planeta,omega);
    
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Coeficiente2);
    Newton.CalculeTodasLasFuerzas(Planeta,omega);
    
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Xi);
    Newton.CalculeTodasLasFuerzas(Planeta,omega);
   
    for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Coeficiente1);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zi);
    x_old_old=x_old;
    x_old=x;
  }
  fout.close();
  fout2.close();
  cout<<sum/c;
  return 0;
}
