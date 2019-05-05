//Mi Primer Programa en C++
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
using namespace std;
const int N=2;
const double l=12;
const double K=1.0e4;
const double c1=0.1786178958448091;
const double c2=-0.2123418310626054;
const double c3=-0.06626458266981849;
const double umdc2=1-2*c2;
const double umdc3=(1-2*(c3+c1));
const double g=980;
const double x0=10;
const double Y0=2.0;

// Declaración de clases

class Cuerpo;
class Colisionador;
//clase cuerpo
class Cuerpo{
private:
  vector3D r,V,F; double m,R,theta;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0,double R0,double t);
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
		    double m0,double R0,double t){
  r.cargue(x0,y0,z0); V.cargue(Vx0,Vy0,Vz0); m=m0; R=R0; theta=t;
							  
}
void Cuerpo::Mueva_r(double dt,double Coeficiente){
  r+=V*(dt*Coeficiente);
}
void Cuerpo::Mueva_V(double dt,double Coeficiente){
  V+=F*(dt*Coeficiente)/m;
}

void Cuerpo::Dibujese(void){
  std::cout << " , " << R << "*cos(2*pi*t) + " << x0 << " + " << l << "*sin(" << theta << ")," << R << "*sin(2*pi*t) + " << l << "*(1-cos(" << theta << "))";
  std::cout << " , " << l << "*sin(" << theta << ")*(1-t) + " << x0 << "," << l+R << "-" << l << "*cos(" << theta << ")*(1-t)";
}
//-------------------Clase colisionador-----------------
class Colisionador{
private:
  
public:
  void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2);
  void CalculeTodasLasFuerzas(Cuerpo * Molecula);
  void Calcule_T(Cuerpo * Molecula);
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
void Colisionador::Calcule_T(Cuerpo * Molecula){
  for(int i=0; i<N;i++)
    {
      Molecula[i].theta=Molecula[i].Getx()/l;
    }
   }
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Molecula){
  int i ,j;vector3D Gravedad, Tension;double aux;
  for(i=0;i<(N);i++)Molecula[i].BorreFuerza();
  for(i=0;i<N;i++){
    aux=(-1)*Molecula[i].m*g;
    Tension.cargue(sin(Molecula[i].theta)*aux,cos(Molecula[i].theta)*aux,0);
    Gravedad.cargue(0,aux,0);
    Molecula[i].AgregueFuerza(Gravedad);
    Molecula[i].AgregueFuerza(Tension);
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
  }
}

//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<10<<"]"<<endl;
  cout<<"set yrange[-10:"<<10<<"]"<<endl;
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
  Cuerpo Molecula[N];
  Crandom ran64(1);
  Colisionador Newton;
  double t,tdibujo,dt=1.0e-3;
  double m0=100,R0=1.5;
  int i,j;
  InicieAnimacion();
  double Theta_1=0,Theta_2=M_PI/3;
 
  
  //Moleculas
   
  Molecula[1].Inicie(x0*sin(Theta_1),Y0*cos(Theta_1),0,0,0,0,m0,R0,Theta_1);
  Molecula[1].Inicie(x0*sin(Theta_2)+R0,Y0*cos(Theta_2),0,0,0,0,m0,R0,Theta_2);

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
    
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c1);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,umdc2);
    
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c3);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,c2);
    
    for( i=0;i<N;i++)Molecula[i].Mueva_r(dt,umdc3);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,c2);
   

    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c3);

    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for( i=0;i<N;i++)Molecula[i].Mueva_V(dt,umdc2);

    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c1);

    Newton.Calcule_T(Molecula);
  }
  


  return 0;
}
