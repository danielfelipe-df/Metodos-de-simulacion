#include <iostream>
#include <fstream>
#include "Random64.h"

const int Lx = 256;
const int Ly = 256;
const double p0=0.25;
const double p=0.25;
const double pp=1-(2*p)-p0;
const double p1=p+p0;
const double p2=p+p0+pp;
const int Q=4;

class LatticeGas
{
private:
  int V[Q][2];
  double n[Lx][Ly][4], nnew[Lx][Ly][Q];
public:
  void Inicie(int B, double mu1, double sigma, Crandom & ran64, double mu2);
  void Show(bool ImprimirNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double GetSigma2(void);
};

void LatticeGas::Inicie(int B, double mu1, double sigma, Crandom & ran64, double mu2)
{
  int i,j,k,b;
  V[0][0]=1;  V[1][0]=0;  V[2][0]=-1;   V[3][0]=0;
  V[0][1]=0;  V[1][1]=1;  V[2][1]=0;    V[3][1]=-1;

  for(i=0; i<Lx; i++)
    for(j=0; j<Ly; j++)
      for(k=0; k<Q; k++)
        n[i][j][k]=nnew[i][j][k]=0;

  while (B>0){
    i = (int) ran64.gauss(mu1,sigma);
    j = (int) ran64.gauss(mu2,sigma);
    if(i<0) i=0;    if(j<0) j=0;   if(i>Lx-1) i=Lx-1;   if(j>Ly-1) j=Ly-1;
    k = (int) Q*ran64.r();
    if(n[i][j][k]==0){n[i][j][k]=1;  B--;}
  }
}


void LatticeGas::Show(bool ImprimirNew)
{
  for(int k=0; k<Q; k++)
    for(int i=0; i<Lx; i++)
      for(int j=0; j<Ly; j++){
        if(ImprimirNew){std::cout << nnew[i][j][k];}
        else{std::cout << n[i][j][k];}
      }
      std::cout << std::endl;
}


void LatticeGas::Colisione(Crandom &ran64)
{
  int i,j,k;
  for(i=0; i<Lx; i++)
    for(j=0; j<Ly; j++){//Ir celda por celda
      nnew[i][j][0]=p0*n[i][j][0] + p*n[i][j][1] + p*n[i][j][3] + pp*n[i][j][2];
      nnew[i][j][1]=p0*n[i][j][1] + p*n[i][j][2] + p*n[i][j][0] + pp*n[i][j][3];
      nnew[i][j][2]=p0*n[i][j][2] + p*n[i][j][3] + p*n[i][j][1] + pp*n[i][j][0];
      nnew[i][j][3]=p0*n[i][j][3] + p*n[i][j][0] + p*n[i][j][2] + pp*n[i][j][1];
      //std::cout << nnew[i][j][0] << "\t" << nnew[i][j][1] << "\t" << nnew[i][j][2] << "\t" << nnew[i][j][3] << std::endl;

    }
}


void LatticeGas::Adveccione(void)
{
  for(int i=0; i<Lx;i++)
    for(int j=0; j<Ly; j++)
      for(int k=0; k<Q; k++)
        n[(i+Lx+V[k][0])%Lx][(j+Ly+V[k][1])%Ly][k]=nnew[i][j][k];
}

double LatticeGas::GetSigma2(void)
{
  int i,k,j;
  double xprom=0, yprom =0, N=0;
  double varx= 0, vary=0, covarxy=0, sigma2=0;

  for(i=0; i<Lx; i++){
    for(j=0; j<Ly; j++){
      for(k=0; k<Q; k++){
        N += n[i][j][k];
        xprom += i*n[i][j][k];	yprom += j*n[i][j][k];
        varx += i*i*n[i][j][k];	vary += j*j*n[i][j][k];
        covarxy += j*i*n[i][j][k];
        //std::cout << xprom << std::endl;
      }
    }
  }

  xprom /= N;
  yprom /= N;
  varx = (varx-N*xprom*xprom)/(N-1);
  vary = (vary-N*yprom*yprom)/(N-1);
  covarxy = (covarxy-N*xprom*yprom)/(N-1);
  sigma2 = varx + vary + 2*covarxy;
  /*
  std::cout << N << std::endl;
  std::cout << xprom << std::endl;
  std::cout << yprom << std::endl;
  std::cout << varx << std::endl;
  std::cout << vary << std::endl;
  std::cout << covarxy << std::endl;
  std::cout << sigma2 << std::endl;
  */
  return sigma2;
}

int main(void)
{
  LatticeGas Difusion;
  Crandom ran64(1);
  double D=p1/(2*(1-p1));
  double Pendiente=4*D;
  int B=2400;   double mu1=Lx/2.0, mu2=Ly/2.0, sigma=16;
  int t, tmax=350;
  Difusion.Inicie(B,mu1,sigma,ran64, mu2);
  //Difusion.Show(true);
  std::ofstream file("datos_2.csv");
  for(t=0; t<tmax;t++){
    //file << t << " " << Difusion.GetSigma2() << '\n';
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
    //Difusion.Show(true);
    file << t << " " << Difusion.GetSigma2() << '\n';
  }
  file.close();
  std::cout<<"f(x)=m*x+b "<<std::endl;
  std::cout<<"a="<<Pendiente<<std::endl;
  std::cout<<"fit f(x) 'datos_2.csv' u 1:2 via m,b"<<std::endl;
  std::cout<<"title_f(m,b)=sprintf('f(x)=%.2fx+ %.2f',m,b)"<<std::endl;
  std::cout<<"g(x)=a*x+b"<<std::endl;
  std::cout<<"title_g(a,b)=sprintf('g(x)=%.2gx+ %.2g',a,b)"<<std::endl;
  std::cout<<"set terminal png"<<std::endl;
  std::cout<<"set output 'P0="<<p0<<"P="<<p<<"_2.png'" <<std::endl;
  std::cout<<"set xlabel 't' "<<std::endl;
  std::cout<<"set ylabel 'sigma cuadrado'"<<std::endl;
  std::cout<<"plot f(x) t title_f(m,b), g(x) t title_g(a,b)"<<std::endl;
}
