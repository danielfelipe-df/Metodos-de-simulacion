#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=400;
const int Ly=200;

const int Q=5;
const double W0=1.0/3;

//const double C=0.5; // C<0.707 celdas/click
//const double TresC2=3*C*C;
//const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,int i, double Ccelda); //Ahora se coloca Ccelda como argumento para calcular siempre C, para cada casilla
  double Ccelda(int ix, int iy);  //Se define la función que calcula C para cada celda
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0,double Jx0,double Jy0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1-W0)/4.0;
  //Cargar los vectores
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
  return suma;
}
double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i, double Ccelda){ //Se agrega Ccelda como argumento de la función
  double TresC2 = 3*Ccelda*Ccelda, Aux = 1-TresC2*(1-W0);
  if(i==0)
    return rho0*Aux; //AUX en forma explicita
  else
    return w[i]*(TresC2*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0)); //Se colocan las constantes explicitamente, TresC2 y AUX0
}
void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Jx0,Jy0,c;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscópicas
      rho0=rho(ix,iy,false);  Jx0=Jx(ix,iy,false);  Jy0=Jy(ix,iy,false);  c = Ccelda(ix,iy);
      for(i=0;i<Q;i++)
        fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Jx0,Jy0,i,c); //Se agregua Ccelda como argumento
    }
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
        f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  double c; int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      c = Ccelda(ix,iy);
      for(i=0;i<Q;i++)
        f[ix][iy][i]=feq(rho0,Jx0,Jy0,i,c); //Se agregua Ccelda y se calcula para todas las casillas siguiendo el suavizado definida, tanh
    }
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy; double lambda,omega,rho0,Jx0,Jy0,c;
  lambda=10; ix=0; //Primera columna, todos los valores de iy (filas)
  for(iy=0;iy<Ly;iy++){
    c = Ccelda(ix,iy); omega=2*M_PI*c/lambda; rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
    for(i=0;i<Q;i++){
      fnew[ix][iy][i]=feq(rho0,Jx0,Jy0,i,c);
    } //Se agregua Ccelda como argumento y se calcula para cada celda de la condición inicial
  }
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  ofstream MiArchivo(NombreArchivo); double rho0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
      //cout<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

//Se define la función para calcular C en cada celda
double LatticeBoltzmann::Ccelda(int ix, int iy){
 float angulo=M_PI/4; int ix0=iy*tan(angulo); double n = 0.5*tanh(ix-ix0)+1.5;
 return 0.5/n; //Suavizado realizado, antes de ix0 n=1, después n=2
}

int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=400;
/*
  // Estos comandos se descomentan si se quiere guardar el gif
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'OndasTest.gif'" << std::endl;
  
  //Estos comandos se descomentan para hacer el gif
  cout << "set pm3d" << endl;
  cout << "set size ratio 1" << endl;
  cout << "set palette defined (-5 \"red\", 0 \"white\", 5 \"blue\")" << endl;
  cout << "set cbrange[-5:5]" << endl;
  cout << "set xrange[0:400]; set yrange[0:200]; set zrange[-10:10]" << endl;
  cout << "set view 90,0" << endl;
*/

  Ondas.Inicie(0,0,0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
    //Ondas.Imprimase("ondas.dat");
    //cout << "splot 'ondas.dat'" << endl;
  }
  Ondas.Imprimase("Ondas.dat");
  //std::cout << "splot 'Ondas.dat'" << std::endl;

  return 0;
}

