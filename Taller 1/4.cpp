#include <iostream>
#include <cmath>
#include "Vector.h"

void Iniciar();


class Pendulo
{
private:
  double r, l, x0, theta, omega, tau, I, m;

public:
  void First_ubication(double x0, double r0, double theta0, double l0, double I0, double m0);
  void Dibujarse();
  void AgregueTorque(double dtau){tau+=dtau;};
  void InicieTorque(void);
  void Mueva_theta(double dt,double Coeficiente);
  void Mueva_omega(double dt,double Coeficiente);
  //double Getx(void){return r.x();}; //Inline
  //double Gety(void){return r.y();}; //Inline
  //double Getz(void){return r.z();}; //Inline
  friend class Colisionador;
};


class Colisionador{
private:

public:
  void CalculeTodasLasFuerzas(Pendulo* Cuerpo );
  void CalculeLaFuerzaEntre(Pendulo & Pendulo1,Pendulo & Pendulo2);

};


int main(void)
{
  Pendulo objeto[4];
  int i;
  double x0 = 3, theta0 = 0, r0 = 1, l0 = 2; m0 = 1; I0 = m0*r0*r0/2;
  Iniciar();

  objeto[0].First_ubication(x0,r0,-M_PI/12,l0,I0,m0);
  for(i=1, x0=2*r0; i<4; x0+=2*r0, i++)  objeto[i].First_ubication(x0,r0,theta0,l0,I0,m0);
  Ndibujos=1000;

  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo>2*tmax/Ndibujos){
      InicieCuadro();
      for(i=0; i<N;i++)   objeto[i].Dibujarse();
      TermineCuadro();
      tdibujo=0;
    }
    // for(i=0;i<N; i++)    cout<<Pendulo[i].x()<<" "<<Pendulo[i].y()<<endl;

    //Muevase con el algoritmo Omelyan PEFRL
    for(i=0; i<N;i++)  objeto[i].Mueva_theta(dt,Zeta);
    Newton.CalculeTodasLasFuerzas(objeto);
    for(i=0; i<N;i++)  objeto[i].Mueva_omega(dt,(1-2*Lambda)/2);
    for(i=0; i<N;i++)  objeto[i].Mueva_theta(dt,Xi);
    Newton.CalculeTodasLasFuerzas(objeto);
    for(i=0; i<N;i++)  objeto[i].Mueva_omega(dt,Lambda);
    for(i=0; i<N;i++)  objeto[i].Mueva_theta(dt,1-2*(Xi+Zeta));
    Newton.CalculeTodasLasFuerzas(objeto);
    for(i=0; i<N;i++)  objeto[i].Mueva_omega(dt,Lambda);
    for(i=0; i<N;i++)  objeto[i].Mueva_theta(dt,Xi);
    Newton.CalculeTodasLasFuerzas(objeto);
    for(i=0; i<N;i++)  objeto[i].Mueva_omega(dt,(1-2*Lambda)/2);
    for(i=0; i<N;i++)  objeto[i].Mueva_theta(dt,Zeta);
  }
  return 0;
}

//------- Funciones de la clase Pendulo --------
void Pendulo::First_ubication(double x, double r0, double theta0, double l0, double I0, double m0)
{
  r = r0;  theta = theta0;  x0 = x;   l = l0;   I = I0;   m = m0;
}

void Pendulo::Dibujarse()
{
  std::cout << " , " << r << "*cos(2*pi*t) + " << x0 << " + " << l << "*sin(" << theta << ")," << r << "*sin(2*pi*t) + " << l << "*(1-cos(" << theta << "))";
  std::cout << " , " << l << "*sin(" << theta << ")*(1-t) + " << x0 << "," << l+r << "-" << l << "*cos(" << theta << ")*(1-t)";
}

void Pendulo::Mueva_theta(double dt, double Constante){
  theta+=omega*(Constante*dt);
}

void Pendulo::Mueva_omega(double dt, double Constante){
  omega+=tau*(Constante*dt/I);
}

void Pendulo::InicieTorque(void){
  tau=-m*g*l*sin(theta);
}

//Funciones de la clase Colisionador
void Colisionador::CalculeTodasLasFuerzas(Pendulo* Cuerpo){
  int i,j;
  //Borrar todas las fuerzas
  for(i=0; i<N;i++) Cuerpo[i].InicieFuerza();
  //Agregue fuerzas externas
  for(i=0; i<N; i++) Cuerpo[i].AgregueFuerza(-Cuerpo[i].m*g*sin(Cuerpo[i].theta));
  //Calcular todas las fuerzas entre parejas de pendulos
  for(i=0; i<N;i++)
    for(j=i+1; j<N;j++)
      CalculeLaFuerzaEntre(Cuerpo[i], Cuerpo[j]);

}

void Colisionador::CalculeLaFuerzaEntre(Pendulo & Pendulo1,Pendulo & Pendulo2){
  double s, d21, F2, tau1, tau2;
  d21= Pendulo2.x() -Pendulo1.x();
  s= (Pendulo2.R+Pendulo1.R) - d21;
  if(s>0){
    F2= K*pow(s,1.5)*d21/abs(d21); tau2 =F2*Pendulo2.L ; tau1= -F2*Pendulo1.L;
    Pendulo2.AgregueFuerza(tau2);
    Pendulo1.AgregueFuerza(tau1);
  }

}


//------ Funciones normales ------
void Iniciar()
{
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula2.gif'" << std::endl;
  std::cout << "set xrange[-10:10]" << std::endl;
  std::cout << "set yrange[-10:10]" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange[0:1]" << std::endl;
  //std::cout << "plot 0,0";
}

