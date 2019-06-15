#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include "Random64.h"
#include "Vector.h"
using namespace std;
  
const int Lx=256; 
const int Ly=256;
const double p0=0.25;
const double p=0.25;
const double pp=1-(2*p)-p0;
const double pp1=p;
const double p1=p+p0;
const double p2=p+p0+pp;
const int Q=4;
class LatticeGas{
private:
  int V[4][2];
  int n[Lx][Ly][Q], nnew[Lx][Ly][Q];//4 direcciones derecha 0 , arriba 1, izquierda 2, abajo 3
public:
  void Inicie(int B, double mu, double sigma, Crandom & ran64,double mu2);
  void Show(bool ImprimirNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double GetSigma2(void);
};
double LatticeGas::GetSigma2(void)
{
  double N, xprom, sigma2;  int i,k, j;
  for(xprom=0,sigma2=0,N=0, i=0; i<Lx; i++)
    for(j=0; j<Ly; j++)
      for(k=0; k<Q; k++){
	N+=n[i][j][k];
	xprom+=pow((i*i)+(j*j),0.5)*n[i][j][k];
	sigma2+=((i*i)+(j*j))*n[i][j][k];
      }
  xprom /= N;
  sigma2 = (sigma2-N*xprom*xprom)/(N-1);
  return sigma2;
}
  void LatticeGas::Adveccione(void){
  for(int i=0;i<Lx;i++)//para cada celda
    for(int j=0; j<Ly; j++)
      for(int k=0;k<Q;k++)
        n[(i+Lx+V[k][0])%Lx][(j+Ly+V[k][1])%Ly][k]=nnew[i][j][k]; // el primer parentesis con modulo genera fronteras periodicas
}  
 
void LatticeGas::Colisione(Crandom & ran64){
  int i,j, k;
  for(i=0;i<Lx;i++)
    for(j=0; j<Ly; j++){//ir celda por celda
      if(ran64.r()<p0){
        for( k=0;k<Q;k++)
          nnew[i][j][k]=n[i][j][k];//se deja igual
      }
      else if(ran64.r()<p1 ){
        for( k=0;k<Q;k++)
          nnew[i][j][(1+k)%Q]=n[i][j][k];//90 grados
      }
      else if(ran64.r()<p2){
        for( k=0;k<Q;k++)
          nnew[i][j][(2+k)%Q]=n[i][j][k];//180 grados
      }
      else{
        for( k=0;k<Q;k++)
          nnew[i][j][(3+k)%Q]=n[i][j][k];//270 grados
      }
    }
}

void LatticeGas::Show(bool ImprimirNew){
  for(int k=0;k<Q;k++){
    for(int i=0;i<Lx;i++)
      for(int j=0;j<Ly;j++)
        if(ImprimirNew) cout<<nnew[i][j][k]; else cout<<n[i][j][k];
      cout<<endl;
  }
}
void LatticeGas::Inicie(int B, double mu, double sigma, Crandom & ran64,double mu2){
  V[0][0]=1;   V[1][0]=0; V[2][0]=-1;  V[3][0]=0;
  V[0][1]=0;   V[1][1]=1; V[2][1]=0;  V[3][1]=-1;
  //iniciar los contenidos
  int i,k,j;
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<Q;k++)
        n[i][j][k]=nnew[i][j][k]=0;
  
  while(B>0){ // colocar B bolitas
    i=(int)ran64.gauss(mu,sigma); if(i<0) i=0; if(i>Lx-1) i=Lx-1; // genera un numero que no se sale de la celda
    j=(int)ran64.gauss(mu2,sigma); if(j<0) j=0; if(j>Ly-1) j=Ly-1;
    k=(int)Q*ran64.r(); // genera un mumero entre 0 y Q-1
    
    if(n[i][j][k]==0){
      n[i][j][k]=1;B--; // reduce a B
    }
  }
}

int main(void){
  LatticeGas Difusion;
  Crandom ran64(1);
  int B=2400; double mu=Lx/2.0,sigma=32,mu2=Ly/2.0;
  int t, tmax=350;
  Difusion.Inicie(B,mu,sigma,ran64,mu2);
  ofstream file("dat.dat");
  for(t=0;t<tmax;t++){
    file << t << " " << Difusion.GetSigma2() << endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  file.close();
  cout<<"plot \"dat.dat\" w l  "<<endl;
  cout<<"pause 100"<<endl;
  return 0;
}
