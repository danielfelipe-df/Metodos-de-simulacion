#include <iostream>
#include <fstream>
#include <cmath>

double f1(double r1, double r2, double r, double lambda);
double f2(double r1, double r2, double r, double lambda);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda);

int main(void)
{
  double r1 = 1, r2 = 0, r;
  double dr = 0.001;
  double lambda = 1.0, dlambda;

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

  r = 1;    dlambda = 0.01;
  std::ofstream fout2("datos_2b.dat");
  for(lambda=dlambda*10; lambda<15; lambda+=dlambda){
    fout2 << lambda << " " << r1 << '\n';
    UnPasoDeRungeKutta4(r1,r2,r,dlambda,lambda);
  }
  fout2.close();
  
  std::cout << "plot 'datos_2b.dat' w l, 0" << std::endl;
  std::cout << "pause 10" << std::endl;

  return 0;
}


double f1(double r1, double r2, double r, double lambda){return r2;}

double f2(double r1, double r2, double r, double lambda){return -lambda*lambda*r1 - r2/r;}

void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda)
{
  double dr11, dr12, dr13, dr14;
  double dr21, dr22, dr23, dr24;
  
  dr11=f1(r1,r2,r,lambda)*dr;                                   dr21=f2(r1,r2,r,lambda)*dr;
  dr12=f1(r1+0.5*dr11, r2+0.5*dr21, r+0.5*dr, lambda)*dr;       dr22=f2(r1+0.5*dr11, r2+0.5*dr21, r+0.5*dr, lambda)*dr;
  dr13=f1(r1+0.5*dr12, r2+0.5*dr22, r+0.5*dr, lambda)*dr;       dr23=f2(r1+0.5*dr12, r2+0.5*dr22, r+0.5*dr, lambda)*dr;
  dr14=f1(r1+dr13, r2+dr23, r+dr, lambda)*dr;                   dr24=f2(r1+dr13, r2+dr23, r+dr, lambda)*dr;

  r += dr;  r1 += (dr11+2*(dr12+dr13)+dr14)/6;          r2 += (dr21+2*(dr22+dr23)+dr24)/6;
}
