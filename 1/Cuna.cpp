#include <fstream>
#include <cmath>
#include <iostream>
#include "Vector.h"
using namespace std;

//-------------Constantes----------------
const double g=9.8;
const int N=3;
const double K=0.1e10;
const double L=0.12; 
//----------Constantes del metodo--------
const double Zeta=0.1786178958448091;
const double Xi=-0.06626458266981849;
const double Lambda = -0.2123418310626054;
const double theta=1.0/(2-pow(2.0,1.0/3));

class Cuerpo ;
class Colisionador;

//---------------------------------------clase cuerpo
class Cuerpo{
private:
  double tau, omega, theta, m, R, I, xinicial;
public:
  void Inicie(double theta0,double omega0, double m0,double R0,double x0inicial);
  void BorrarTorque(void){tau=0;};
  void AgregarTorque(double dtau){tau+=dtau;};
  void Mueva_theta(double dt,double Constante);
  void Mueva_omega(double dt,double Constante);
  void Dibujese(void);
  double Gettheta(void){return theta;};
  double Getx(void){return xinicial+L*sin(theta);};
  double Gety(void){return -L*cos(theta);};
  double Gettau(void){return tau;};
  friend class Colisionador;
  friend int main ();
};

void Cuerpo::Inicie(double theta0,double omega0, double m0,double R0,double x0inicial){
  omega=omega0;
  theta=theta0;
  m=m0;	R=R0; 
  xinicial=x0inicial;
  I=m*L*L;

}
void Cuerpo::Mueva_theta(double dt,double Constante){
theta+=omega*(Constante*dt);
}
void Cuerpo::Mueva_omega(double dt,double Constante){
omega+=tau*(Constante*dt)/I;
}

void Cuerpo::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<xinicial<<"+"<<sin(theta)*L/7<<"*t,"<<-cos(theta)*L/7<<"*t";
   }

//---------------------------------------clase colisionador

class Colisionador{
private:

public:
	void CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);
        void CalculeTodosLosTorques(Cuerpo* Pendulo);

};

void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  double s, d21, F2, tau1, tau2;
  d21= Pendulo2.Getx()-Pendulo1.Getx();
  s=(Pendulo2.R+Pendulo1.R) - fabs(d21);
  if(s>=0){ // si hay choque 
    F2=d21*K*pow(s,1.5)/fabs(d21); tau2 =F2*L ;
    Pendulo2.AgregarTorque(tau2);
    Pendulo1.AgregarTorque(-tau2);
  }
} 
  void Colisionador::CalculeTodosLosTorques(Cuerpo* Pendulo){
  int i,j;
  
  for(i=0;i<N;i++){
   Pendulo[i].BorrarTorque();
  }
  //fuerza gravitacional 
   for(i=0; i<N; i++){
   Pendulo[i].AgregarTorque(-Pendulo[i].m*g*L*sin(Pendulo[i].theta));
    }
   //choque entre las bolas
  for(i=0;i<N;i++)
  for(j=i+1; j<N;j++){
    CalculeTorqueEntre(Pendulo[i],Pendulo[j]);
  }
  }
  


//-----------------------------------Animación 
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl;
  //cout<<"set output 'planeta.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange [-0.05:0.17]"<<endl;
  cout<<"set yrange [-0.15:0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
}

void InicieCuadro(void){
  cout<<"plot 0,0 ";
}

void TermineCuadro(void){
  cout<<endl;
}

int main(void){
  double t, dt=1e-6;
  Cuerpo Pendulo[N];
  Colisionador Newton;
  int i;
    
  double m0=1, R0=0.015;
  double x0inicial=0, theta0=-15*M_PI/180, T= 2*M_PI*sqrt(L/g), tmax=5*T;
  
  //InicieAnimacion();

//-----------Inicie(theta0,omega0, m0, R0, x0inicial);
  Pendulo[0].Inicie(theta0,   0,   m0, R0, x0inicial);
  for(i=1;i<N;i++)
  Pendulo[i].Inicie(     0,   0,   m0, R0,x0inicial+2*i*R0);
  
  
  double tdibujo; int Ndibujos=200;
   
  ofstream fout("data_cuna.dat");
  for (t=tdibujo=0; t<0.2; t+=dt, tdibujo+=dt){
    //cout<<t<<" "<<Pendulo[0].tau<<endl;
    // cout<<Pendulo[0].Getx()<<'\t'<<Pendulo[0].Gety()<<endl;
    //cout<<t<<" "<<Pendulo[0].Getx()<<endl;
    // cout<<t<<" "<<Pendulo[1].Getx()<<endl;
    //cout<<t<<" "<<Pendulo[2].Getx()<<endl;
    /*
     if (tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for (i=0; i<N; i++)
	Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;

      }
    */

    fout << t << "\t";
    for(i=0;i<N;i++){
       fout << Pendulo[i].Gettau() << "\t";
    }
    fout << "\n";

    //Muevase con Forest-Ruth
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,0.5*theta);  //dt, Zeta
    Newton.CalculeTodosLosTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,theta);      //dt, (1-2*Lambda)/2
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(1-theta)/2);//dt, Xi 
    Newton.CalculeTodosLosTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,1-2*theta);  //dt, Lambda
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(1-theta)/2);//dt, 1-2*(Xi+Zeta)
    Newton.CalculeTodosLosTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,theta);      //dt, Lambda
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(theta)/2);  //dt, Xi
    Newton.CalculeTodosLosTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,theta);      //dt, (1-2*Lambda)/2
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,(theta)/2);  //dt, Zeta
  }
  fout.close();
  cout << "plot 'data_cuna.dat' u 1:3 w l" << endl;
  cout << "pause 10" << endl;
  return 0;
}

