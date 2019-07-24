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
const int Lx=(40*proportion)+2;
const int Ly=(21*proportion)+2;
const int Lz=(11*proportion)+2;

const int Q=7;
const double W0=1.0/4;
//Constante de reflexión
//const double k=1;

const double k_ladrillo = 0.0; //Constante de reflexión para el ladrillo normal
const double k_zigzag = 0.0; //Constante de reflexión para el ladrillo en zig-zag
const double k_pared = 0.0; //Constante de reflexión para la pared blanca
const double k_baldosa = 0.0; //Constante de reflexión para la baldosa
const double k_techo = 0.0; //Constante de reflexión para el techo
const double k_ventanita = 0.0; //Constante de reflexión para la ventanita del fondo
const double k_tablero = 0.0; //Constante de reflexión para el tablero
const double k_ventanas = 0.0; //Constante de reflexión para las ventanas grandes
const double k_puerta = 0.0; //Constante de reflexión para las puertas

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
  double GetRho(int ix, int iy, int iz, bool algo){double rho0 = rho(ix,iy,iz,algo);  return rho0;};
};

#endif // LATTICEBOLTZMANN_H
