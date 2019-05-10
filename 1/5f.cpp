#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"
#include <cstdlib>
#include <time.h>
using namespace std;
const double ka=1.0e4;
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
  double GetVy(void){return V.y();};
  double GetVz(void){return V.z();};
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
  int k ,l;  double h;vector3D Punto1,Punto2,F,a,b;
  Punto1.cargue(0,120,0);Punto2.cargue(60,0,0);

  for(k=0;k<(N);k++)Molecula[k].BorreFuerza();
  for(k=0;k<N;k++){
    a=Punto1-Molecula[k].r;
    if(abs(a.x())-Molecula[k].R<=0)
      {	
	h=Molecula[k].R-abs(a.x());
	F.cargue(ka*pow(h,1.5),0,0);
	Molecula[k].AgregueFuerza(F);
      }
    else if(abs(a.y())-Molecula[k].R<=0)
      {       
	h=Molecula[k].R-abs(a.y());
	F.cargue(0,-ka*pow(h,1.5),0);
	Molecula[k].AgregueFuerza(F);
      }
    else
      {
	b=Punto2-Molecula[k].r;
	if(abs(b.x())-Molecula[k].R<=0)
	  {
	    
	    h=Molecula[k].R-abs(b.x());
	    F.cargue(-ka*pow(h,1.5),0,0);
	    Molecula[k].AgregueFuerza(F);
	  }
	else if(abs(b.y())-Molecula[k].R<=0)
	  {
	    h=Molecula[k].R-abs(b.y());
	    F.cargue(0,ka*pow(h,1.5),0);
	    Molecula[k].AgregueFuerza(F);
	  }         
      }
    for(l=k+1;l<N;l++)
      CalculeFuerzaEntre(Molecula[k],Molecula[l]);
  }
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
  Crandom ran64(1);
  Colisionador Newton;
  double t,dt=1.0e-3;
  double m0=1,R0=2.5;
  int i,j;
  double teq = 70;
  //Moleculas
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1),x0,y0, theta;
  double V0=sqrt(2*KT/m0);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      x0=dx*(i+1); y0=(j+1)*dy;theta=2*M_PI*ran64.r();
      Molecula[i*Ny+j].Inicie(x0,y0,0,V0*cos(theta),V0*sin(theta),0,m0,R0); 
    }
  }

  //------------(x0,y0,z0,Vx0,Vy0  ,Vz0,  m0, R0)
  double T=200;
  int M = int(T/dt);
  double array_Vx[M];
  double sum1;
  
  //ofstream file("data_5f.dat");
  for(t=0, j=0;t<T;t+=dt, j++){
    
    if(t<=(teq+dt) && t>=(teq-dt)){
      ofstream file("data_5f.dat");
      for(i=0;i<N;i++){
	file << i << " " << Molecula[i].GetVx() << "\n";
	array_Vx[i] = Molecula[i].GetVx();
      }
      file.close();
      break;
    }
    
    /*
    if(t>teq){
      sum1=0;
      for(i=0; i<N; i++) sum1 += Molecula[i].GetVx();
      file << t << " " << sum1/N << "\n";
      array_Vx[j] = sum1/N;
    }
    */
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c1);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,umdc2);
    
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c3);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,c2);
    
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,umdc3);
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,c2);
    
    
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c3);
    
    
    Newton.CalculeTodasLasFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,umdc2);
    
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,c1);    
  }
  //file.close();
  //printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC)

  std::cout << "set terminal png" << std::endl;
  std::cout << "set output 'vx.png'" << std::endl;
  std::cout << "plot 'data_5f.dat' w l" << std::endl;
  //std::cout << "pause 10" << std::endl;

  std::cout << "set terminal png" << std::endl;
  std::cout << "set output 'histogram_vx.png'" << std::endl;
  std::cout << "binwidth=1" << std::endl;
  std::cout << "bin(x,width)=width*floor(x/width)" << std::endl;
  std::cout << "plot 'data_5f.dat' using (bin($2,binwidth)):(1.0) smooth freq with boxes" << std::endl;
  //std::cout << "pause 10" << std::endl;

  double sum=0, average, s;
  double sigma = sqrt(KT/m0);
  for(i=0; i<N; i++) sum += array_Vx[i];
  average = sum/N;

  sum=0;
  for(i=0; i<N;i++) sum += pow(array_Vx[i] - average, 2);
  s = sqrt(sum/(N-1));

  std::cout << sigma << "\t" << s << "\t" << sigma-s << "\t" << abs(sigma-s)/sigma << std::endl;
  
  return 0;
}
