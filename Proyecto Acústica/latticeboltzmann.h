#ifndef LATTICEBOLTZMANN_H
#define LATTICEBOLTZMANN_H

#include <iostream>
#include <fstream>
#include <cmath>

/*
  Para proportion=12 se manejaron valores de lambda=20 y A=150.
  Además de que se comienza a graficar en x=5, para no ver la fuente inicial.
*/


const int proportion=3;
const int Lx=40*proportion;
const int Ly=21*proportion;
const int Lz=11*proportion;

const int Q=7;
const double W0=1.0/4;
//Constante de reflexión
//const double k=1;
const double kx_1 = 0; //Constante de reflexión para la zona del tablero.
const double kx_2 = 0; //Constante de reflexión para la zona del fondo.
const double ky_1 = 0; //Constante de reflexión para la zona de las puertas.
const double ky_2 = 0; //Constante de reflexión para la zona de las ventanas.
const double kz_1 = 0; //Constante de reflexión para la zona del piso.
const double kz_2 = 0; //Constante de reflexión para la zona del techo.
const double k_1 = 0; //Constante de reflexión para la columna.
const double k_2 = 0; //Constante de reflexión para la puerta.
const double k_3 = 0; //Constante de reflexión para el tablero.
const double k_4 = 0; //Constante de reflexión para la ventanita del fondo.

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
  bool Rectangulo(int x1, int x2, int x, int y1, int y2, int y);
};

#endif // LATTICEBOLTZMANN_H
