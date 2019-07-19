#ifndef LATTICEBOLTZMANN_H
#define LATTICEBOLTZMANN_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "omp.h"
/*
  Para proportion=12 se manejaron valores de lambda=20 y A=150.
  Además de que se comienza a graficar en x=5, para no ver la fuente inicial.
*/
const int proportion=2;
const int Lx=40*proportion;
const int Ly=21*proportion;
const int Lz=11*proportion;

const int Q=7;
const double W0=1.0/4;
//Constante de reflexión
//const double k=1;
const double kx = 0; //Constante de reflexión para las zonas donde x es constante. Es decir, en la pared del tablero y el fondo.
const double ky = 0; //Constante de reflexión para las zonas donde y es constante. Es decir, en los otros costados del salón.
const double kz = 0; //Constante de reflexión para las zonas donde z es constante. Es decir, en el piso y el techo.
const double k_1 = 0; //Constante de reflexión para la columna.
const double k_2 = 0; //Constante de reflexión para la puerta.

const double C=0.5; // C<0.707 celdas/click
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann
{
private:
  double w[Q];
  int V[3][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy, V[2][i]=V_iz
  double f[Lx][Ly][Lz][Q], fnew[Lx][Ly][Lz][Q]; // f[ix][iy][iz][iz][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,int iz,bool UseNew);
  double Jx(int ix,int iy,int iz,bool UseNew);
  double Jy(int ix,int iy,int iz,bool UseNew);
  double Jz(int ix,int iy,int iz,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,double Jz0,int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0,double Jx0,double Jy0, double Jz0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo);
  void Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo);
  bool Columna(int x1, int x2, int x);
};

#endif // LATTICEBOLTZMANN_H
