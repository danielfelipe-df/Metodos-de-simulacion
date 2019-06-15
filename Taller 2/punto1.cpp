#include <iostream>
#include <fstream>
#include <cmath>
#include "Random64.h"
#include "Vector.h"
using namespace std;
  
const int Lx=256; 
const int Ly=256;
const double p0=0.5;
const double p=0.25;
const double p1=1-(2*p0)-p;
const double p2=p;
class LatticeGas{
private:
  int V[4][2];
  int n[Lx][Ly][4], nnew[Lx][Ly][4];//4 direcciones derecha 0 , arriba 1, izquierda 2, abajo 3
public:
  void Inicie(int B, double mu, double sigma, Crandom & ran64);
  void Show(bool ImprimirNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double GetSigma2(void);
};
double LatticeGas::GetSigma2(void)
{
  double N, xprom, sigma2;  int i,k, j;
  //cuantas bolitas 
  for(N=0, i=0; i<Lx; i++){
    for(j=0; j<Ly; j++){
      for(k=0; k<4; k++){
        N+=n[i][j][k];
      }
    }
  }
  //Distancia promedio
  for(xprom=0, i=0; i<Lx; i++){
    for(j=0; j<Ly; j++){
      for(k=0; k<4; k++){
        xprom+=pow((i*i)+(j*j),0.5)*n[i][j][k];
      }
    }
  }
  xprom /= N;
  //distancia al cuadrado
  for(sigma2=0, i=0; i<Lx; i++){
    for(j=0; j<Ly; j++){
      for(k=0; k<4; k++){
        sigma2+=((i*i)+(j*j))*n[i][j][k];
      }
    }
  }
  //distancia al cuadrado menos N por el promedio al cuadrado sobre n-1. asi es la definicion.
  sigma2 = (sigma2-N*xprom*xprom)/(N-1);

  return sigma2;
}
  void LatticeGas::Adveccione(void){
  for(int i=0;i<Lx;i++)//para cada celda
    for(int j=0; j<Ly; j++)
      for(int k=0;k<4;k++)
        n[(i+Lx+V[k][0])%Lx][(i+Ly+V[k][1])%Ly][k]=nnew[i][j][k]; // el primer parentesis con modulo genera fronteras periodicas
}  
 
void LatticeGas::Colisione(Crandom & ran64){
  for(int i=0;i<Lx;i++){
    for(int j=0; j<Ly; j++){//ir celda por celda
      if(ran64.r()<p0){
        for(int k=0;k<4;k++)
          nnew[i][j][k]=n[i][j][k];//se deja igual
      }
      else if(ran64.r()<(p0+p) && ran64.r()>p0 ){
        for(int k=0;k<4;k++)
          nnew[i][j][(1+k)%4]=n[i][j][k];//90 grados
      }
      else if(ran64.r()<(p0+p+p1)&&ran64.r()>(p0+p) ){
        for(int k=0;k<4;k++)
          nnew[i][j][(2+k)%4]=n[i][j][k];//180 grados
      }
      else{
        for(int k=0;k<4;k++)
          nnew[i][j][(3+k)%4]=n[i][j][k];//270 grados
      }
    }
  }
}

void LatticeGas::Show(bool ImprimirNew){
  for(int k=0;k<4;k++){
    for(int i=0;i<Lx;i++)
      for(int j=0;j<Ly;j++)
        if(ImprimirNew) cout<<nnew[i][j][k]; else cout<<n[i][j][k];
      cout<<endl;
  }
}
void LatticeGas::Inicie(int B, double mu, double sigma, Crandom & ran64){
  V[0][0]=1;   V[1][0]=0; V[2][0]=-1;  V[3][0]=0;
  V[0][1]=0;   V[1][1]=1; V[2][1]=0;  V[3][1]=-1;
  //iniciar los contenidos
  int i,k,b,j;
  for(i=0;i<Lx;i++)
    for(j=0;j<Ly;j++)
      for(k=0;k<4;k++)
        n[i][j][k]=nnew[i][j][k]=0;

  while(B>0){ // colocar B bolitas
    i=(int)ran64.gauss(mu,sigma); if(i<0) i=0; if(i>Lx-1) i=Lx-1; // genera un numero que no se sale de la celda
    k=(int)4*ran64.r(); // genera un mumero entre 0 y 3
    
    if(n[i][j][k]==0){
      n[i][j][k]=1;B--; // reduce a B
    }
  }
}

int main(void){
  LatticeGas Difusion;
  Crandom ran64(1);
  int B=600; double mu=Lx/2.0,sigma=16;
  int t, tmax=10;
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
