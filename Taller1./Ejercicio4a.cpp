#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

//const int N = 2;  //al creal los granos, envez de crear 2 se crean 6, n4. Son las paredes y los ultimos 4 son las paredes
const double K = 10e9;
const double Lx=100,Ly=100;   //lonitud de las paredes
const int Nx=2,Ny=1,N=Nx*Ny;  //Numero de granos en fila, columna y total

const double G = 1.0;
const double zetha = 0.1786178958448091;
const double lambda = -0.2123418310626054;
const double xhi = -0.06626458266981849;

//-------- Declaración de las clases -------
class Cuerpo;
class Colisionador;
//-------- Declaración de las clases -------

class Cuerpo{
private:
  vector3D r, V, F;
  double m, R;

public:
  void Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AgregueFuerza(vector3D df){F+=df;};
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  void Arranque(double dt);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  double GetVx(void){return V.x();}; //Inline
  friend class Colisionador; //Esto permite que la clase Colisionador pueda usar los datos privados de Cuerpo
};


class Colisionador{
private:

public:
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
  void CalculeTodasLasFuerzas(Cuerpo * Grano);

};

//-------Funciones Globales------
void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);

int main(void)
{
  Cuerpo Grano[N+4];
  Colisionador Newton;
  Crandom ran64(1);
  double t, tdibujo, dt=1.0;
  int i,j;
  double m0=100, R0=1.5;
  double kT=10;
  
  InicieAnimacion();

  //INICIALIZACION-----------------------------------
  
  //std::cout << T << std::endl;
  //PAREDE
  double Mpared=100*m0, Rpared=10000;
  //Pared Arriba
  Grano[N].Inicie(Lx/2,Ly+Rpared, 0, 0, 0, 0, Mpared, Rpared);
  //Pared Abajo
  Grano[N+1].Inicie(Lx/2,-Rpared, 0, 0, 0, 0, Mpared, Rpared);
  //Pared Derecha
  Grano[N+2].Inicie(Lx+Rpared,Ly/2, 0, 0, 0, 0, Mpared, Rpared);
  //Pared Izquierda
  Grano[N+3].Inicie(-Rpared,Ly/2, 0, 0, 0, 0, Mpared, Rpared);

  
  //----------------GRANOS-----------------
  double dx=Lx/(Nx+1), dy=Ly/(Ny+1), x0, y0, Vx0, Vy0;   //espacio uniforme entre moleculas
  double V0=sqrt(2*kT/m0), theta;


  for(i=0;i<Nx;i++)
    for(j=0;j<Ly;j++){
      //-----(x0, y0, z0, Vx0, Vy0, Vz0, m0, R0)
      x0=(1+i)*dx ; y0=(j+1)*dy;      theta=2*M_PI*ran64.r(); Vx0= V0*cos(theta); Vy0=V0*sin(theta);
      Grano[i*Ny+j].Inicie(x0, 4.0, 0, Vx0, Vy0,0, m0, R0);
    }
      
  //CORRER LA SIMULACION----------------
  double T=50,Teq=20;
  for(t=tdibujo=0;t<Teq+T;t+=dt,tdibujo+=dt){
    
    if(tdibujo>T/1000){
      // if(t>Teq) for(int n=0;n<N;n++) cout<<Grano[n].GetVx()<<endl;
      
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }
    //Moverlo según Forest Ruth Orden 4
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt, zetha);
    Newton.CalculeTodasLasFuerzas(Grano);
    for(i=0;i<N;i++) Grano[i].Mueva_V(dt, (1-2*lambda)/2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt, xhi);
    Newton.CalculeTodasLasFuerzas(Grano);
    for(i=0;i<N;i++) Grano[i].Mueva_V(dt, lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt, 1-2*(xhi+zetha));
    Newton.CalculeTodasLasFuerzas(Grano);
    for(i=0;i<N;i++) Grano[i].Mueva_V(dt, lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt, xhi);
    Newton.CalculeTodasLasFuerzas(Grano);
    for(i=0;i<N;i++) Grano[i].Mueva_V(dt, (1-2*lambda)/2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt, zetha);
  }

  return 0;
}


//--------- Funciones Cuerpo------
void Cuerpo::Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0)
{
  r.cargue(x0,y0,z0);  V.cargue(Vx0,Vy0,Vz0);  m=m0;  R=R0;
}

void Cuerpo::Mueva_r(double dt, double coeficiente)
{
  r += V*(dt*coeficiente);
}

void Cuerpo::Mueva_V(double dt, double coeficiente)
{
  V += F*(dt*coeficiente/m);
}

void Cuerpo::Dibujese(void)
{
  std::cout << " , " << r.x() << "+" << R << "*cos(t)," << r.y() << "+" << R << "*sin(t)";
}


//-------Funciones Colisionador------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano)
{
  int i,j;
  for(i=0;i<(N+4);i++) Grano[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=i+1;j<(N+4);j++)
      CalculeFuerzaEntre(Grano[i], Grano[j]);
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2)
{
  vector3D dr = Grano2.r - Grano1.r;
  double Norma_dr=norma(dr);
  double s=(Grano2.R+Grano1.R)-Norma_dr;
  if(s>0){
    vector3D F2=dr*(K*pow(s,1.5)/Norma_dr);
    Grano2.AgregueFuerza(F2);    Grano1.AgregueFuerza(F2*(-1));
  }
}


//-------Funciones Globales------
void InicieAnimacion(void)
{
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula.gif'" << std::endl;
  std::cout << "unset key" << std::endl;
  std::cout << "set xrange[0:60]" << std::endl;
  std::cout << "set yrange[0:60]" << std::endl;
  std::cout << "set size ratio -1" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange[0:7]" << std::endl;
  std::cout << "set isosamples 12" << std::endl;
}

void InicieCuadro(void)
{
  std::cout << "plot 0,0 ";
  cout<<" , "<<Lx/7<<"*t,0";
  cout<<" , "<<Lx/7<<"*t,"<<Ly;
  cout<<" , 0,"<<Ly/7<<"*t";
  cout<<" , "<<Lx<<","<<Ly/7<<"*t";
}

void TermineCuadro(void)
{
  std::cout << std::endl;
}
