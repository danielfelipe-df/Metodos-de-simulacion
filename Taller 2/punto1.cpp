#include <iostream>
#include <fstream>
#include <cmath>
#include "Random64.h"
#include "Vector.h"
using namespace std;
  
const int Lx=256; 
const int Ly=256;
const double p0=0.25;
const double p=0.25;
const double p1=1-(2*p0)-p;
const double p2=p;
class LatticeGas{
private:
  int V[4][1];
  int n[Lx][4], nnew[Lx][4];//4 direcciones arriba abajo derecha izquierda
public:
  void Inicie(int B, double mu, double sigma, Crandom & ran64);
  void Show(bool ImprimirNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double GetSigma2(void);
};
double LatticeGas::GetSigma2(void)
{
  double N, xprom, sigma2;  int i,k;
  //cuantas bolitas 
  for(N=0, i=0; i<Lx; i++){
    for(k=0; k<4; k++){
      N+=n[i][k];
    }
  }
  //Distancia promedio
  for(xprom=0, i=0; i<Lx; i++){
    for(k=0; k<4; k++){
      xprom+=i*n[i][k];
    }
  }
  xprom /= N;
  //distancia al cuadrado
  for(sigma2=0, i=0; i<Lx; i++){
    for(k=0; k<4; k++){
      sigma2+=i*i*n[i][k];
    }
  }
  //distancia al cuadrado menos N por el promedio al cuadrado sobre n-1. asi es la definicion.
  sigma2 = (sigma2-N*xprom*xprom)/(N-1);

  return sigma2;
}
void LatticeGas::Adveccione(void){
  for(int i=0;i<Lx;i++)//para cada celda
    for(int k=0;k<4;k++)
      n[(i+Lx+V[k][0])%Lx][k]=nnew[i][k]; // el primer parentesis con modulo genera fronteras periodicas
}  
 
void LatticeGas::Colisione(Crandom & ran64){
  for(int i=0;i<Lx;i++)//ir celda por celda
    if(ran64.r()<p0){
      for(int k=0;k<4;k++)
	nnew[i][k]=n[i][k];//se deja igual   
    }
    else if(ran64.r()<(p0+p)){
      for(int k=0;k<4;k++)
	nnew[(5+i)%4][k]=n[i][k];//90 grados
    }
    else if(ran64.r()<(p0+p+p1)){
      for(int k=0;k<4;k++)
	nnew[(6+i)%4][k]=n[i][k];//180 grados
    }
    else{
      for(int k=0;k<4;k++)
	nnew[(7+i)%4][k]=n[i][k];//270 grados
    }
 
}

void LatticeGas::Show(bool ImprimirNew){
  for(int k=0;k<4;k++){
    for(int i=0;i<Lx;i++)
      if(ImprimirNew) cout<<nnew[i][k]; else cout<<n[i][k];
    cout<<endl;
      }
}
void LatticeGas::Inicie(int B, double mu, double sigma, Crandom & ran64){
  V[0][0]=1;   V[1][0]= 1; V[2][0]=-1;V[3][0]=-1;
  //iniciar los contenidos
  int i,k,b;
  for(i=0;i<Lx;i++)
    for(k=0;k<4;k++)
      n[i][k]=nnew[i][k]=0;
  while(B>0){ // colocar B bolitas
    i=(int)ran64.gauss(mu,sigma); if(i<0) i=0; if(i>Lx-1) i=Lx-1; // genera un numero que no se sale de la celda
    k=(int)4*ran64.r(); // genera un mumero entre 0 y 3

    if(n[i][k]==0){
      n[i][k]=1;B--; // reduce a B
    }
  }
}

int main(void){
  LatticeGas Difusion;
  Crandom ran64(1);
  int B=400; double mu=Lx/2.0,sigma=16;
  int t, tmax=350;
  Difusion.Inicie(B,mu,sigma,ran64);
  // Difusion.Show(false);
  ofstream file("dat.dat");
  for(t=0;t<tmax;t++){
    file << t << " " << Difusion.GetSigma2() << endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  file.close();
  cout<<"plot \"dat.dat\" w l "<<endl;
  cout<<"pause 10"<<endl;
  //cout<<endl;
  //Difusion.Show(true);
  return 0;
}
