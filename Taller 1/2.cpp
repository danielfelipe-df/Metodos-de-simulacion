#include <iostream>
#include <fstream>
#include <cmath>

const double ERR =1e-7;

double f1(double r1, double r2, double r, double lambda);
double f2(double r1, double r2, double r, double lambda);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda);
double Biseccion(double a, double b);
double f(double a);

int main(void)
{
  double r1 = 1, r2 = 0, r;
  double dr = 0.001;
  double lambda = 2.4048, dlambda;

  //-----Punto a-----

  std::ofstream fout1("datos_2a.dat");
  for(r=dr/10; r<10; r+=dr){
    fout1 << r << " " << r1 << '\n';
    UnPasoDeRungeKutta4(r1,r2,r,dr,lambda);
  }
  fout1.close();

  std::cout << "plot 'datos_2a.dat' w l, 0" << std::endl;
  std::cout << "pause 10" << std::endl;


  //-----Punto b-----

  dlambda = 0.01;
  /*
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula1.gif'" << std::endl;
  std::cout << "unset key" << std::endl;
  std::cout << "set xrange[0:2]" << std::endl;
  std::cout << "set yrange[-1:1]" << std::endl;
  */
  std::ofstream fout2("datos_2b.dat");
  for(lambda=dlambda; lambda<15; lambda+=dlambda){
    r1 = 1;   r2 = 0;
    for(r=dr/10; r<1; r+=dr){
      UnPasoDeRungeKutta4(r1,r2,r,dr,lambda);
    }
    fout2 << lambda << " " << r1 << '\n';
  }
  fout2.close();

  std::cout << "plot 'datos_2b.dat' w l, 0" << std::endl;
  std::cout << "pause 10" << std::endl;

  //----- Punto c------



  //----- Punto d------


  return 0;
}


double f1(double r1, double r2, double r, double lambda){return r2;}

double f2(double r1, double r2, double r, double lambda){return -(lambda*lambda*r1) - (r2/r);}

void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda)
{
  double dr11, dr12, dr13, dr14;
  double dr21, dr22, dr23, dr24;

  dr11=f1(r1,r2,r,lambda);                                            dr21=f2(r1,r2,r,lambda);
  dr12=f1(r1+0.5*dr11*dr, r2+0.5*dr21*dr, r+0.5*dr, lambda);          dr22=f2(r1+0.5*dr11*dr, r2+0.5*dr21*dr, r+0.5*dr, lambda);
  dr13=f1(r1+0.5*dr12*dr, r2+0.5*dr22*dr, r+0.5*dr, lambda);          dr23=f2(r1+0.5*dr12*dr, r2+0.5*dr22*dr, r+0.5*dr, lambda);
  dr14=f1(r1+dr13*dr, r2+dr23*dr, r+dr*dr, lambda);                   dr24=f2(r1+dr13*dr, r2+dr23*dr, r+dr, lambda);

  r += dr;  r1 += dr*(dr11+2*(dr12+dr13)+dr14)/6;          r2 += dr*(dr21+2*(dr22+dr23)+dr24)/6;
}

double Biseccion(double a, double b)
{
  double m,fa,fm;
  fa=f(a);
  while(b-a>ERR){
    m=(a+b)/2; fm=f(m); //calculo m;
    if(fa*fm<0)  //  if(f(a) y f(m) son de signo contrario)
      b=m;
    else
      {a=m; fa=fm;}
  }
  return (a+b)/2;
}

double f(double a){return a;}
