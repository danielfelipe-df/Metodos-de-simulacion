#include <iostream>
#include <fstream>
#include <cmath>

const double lambda = 1.0;

double f1(double r1, double r2, double r);
double f2(double r1, double r2, double r);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr);

int main(void)
{
  //Definimos r1 = R(r) y r2 = dR/dr
  //Así obtenemos que dr1/dr = r2 y dr2/dr1 = -(lambda*lambda*r1) - (r2/r)
  //Planteamos las condiciones iniciales
  double r1 = 1, r2 = 0, r;
  double dr = 0.01;

  //Hacemos el Runge-Kutta para el intervalo pedido e imprimos los datos para después poder gaficar
  std::ofstream fout("datos_2a.dat");
  for(r=dr; r<(10+dr); r+=dr){
    fout << r << " " << r1 << '\n';
    UnPasoDeRungeKutta4(r1,r2,r,dr);
  }
  fout.close();

  //Generamos la gráfica y la guardamos
  std::cout << "set terminal pdf" << std::endl;
  std::cout << "set output 'figure_2a.pdf'" << std::endl;
  std::cout << "set title 'Imagen del ejercio 2a'" << std::endl;
  std::cout << "set xlabel 'r'" << std::endl;
  std::cout << "set ylabel 'R(r)'" << std::endl;
  std::cout << "plot 'datos_2a.dat' w l title 'R(r)', 0 title ''" << std::endl;

  return 0;
}

//Definimos f1(r) = dr1/dr
double f1(double r1, double r2, double r){return r2;}

//Definimos f2(r) = dr2/dr
double f2(double r1, double r2, double r){return -(lambda*lambda*r1) - (r2/r);}

//Planteamos el Runge-Kutta respectivo
//Nótese que como hacemos el for sobre r con r+=dr el paso de loop, no es necesario definir en esta función r+=dr
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr)
{
  double dr11, dr12, dr13, dr14;
  double dr21, dr22, dr23, dr24;

  dr11=f1(r1,r2,r);                                            dr21=f2(r1,r2,r);
  dr12=f1(r1+0.5*dr11*dr, r2+0.5*dr21*dr, r+0.5*dr);          dr22=f2(r1+0.5*dr11*dr, r2+0.5*dr21*dr, r+0.5*dr);
  dr13=f1(r1+0.5*dr12*dr, r2+0.5*dr22*dr, r+0.5*dr);          dr23=f2(r1+0.5*dr12*dr, r2+0.5*dr22*dr, r+0.5*dr);
  dr14=f1(r1+dr13*dr, r2+dr23*dr, r+dr*dr);                   dr24=f2(r1+dr13*dr, r2+dr23*dr, r+dr);

  r1 += dr*(dr11+2*(dr12+dr13)+dr14)/6;          r2 += dr*(dr21+2*(dr22+dr23)+dr24)/6;
}
