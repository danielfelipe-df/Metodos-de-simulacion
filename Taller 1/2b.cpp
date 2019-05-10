#include <iostream>
#include <fstream>
#include <cmath>

double f1(double r1, double r2, double r, double lambda);
double f2(double r1, double r2, double r, double lambda);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda);

int main(void)
{
  double r1, r2, r, r_max = 1.0;
  double dr = 0.01;
  double lambda, dlambda = 0.01, lambda_max = 15;

  //Vamos a resolver la misma ecuación diferencial hasta r_max = 1.0 para distintos lambda
  //Así, el último valor de r será la función f(r=rmax=1.0; lambda)
  std::ofstream fout("datos_2b.dat");
  for(lambda=dlambda*10; lambda<(lambda_max+dlambda); lambda+=dlambda){
    r1 = 1;   r2 = 0;
    for(r=dr; r<(r_max+dr); r+=dr){
      UnPasoDeRungeKutta4(r1,r2,r,dr,lambda);
    }
    fout << lambda << " " << r1 << '\n';
  }
  fout.close();

  //Generamos la gráfica y la guardamos
  std::cout << "set terminal pdf" << std::endl;
  std::cout << "set output 'figure_2b.pdf'" << std::endl;
  std::cout << "set title 'Imagen del ejercio 2b'" << std::endl;
  std::cout << "set xlabel 'lambda'" << std::endl;
  std::cout << "set ylabel 'f(lambda)'" << std::endl;
  std::cout << "plot 'datos_2b.dat' w l title 'f(lambda)', 0 title ''" << std::endl;

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

  r1 += dr*(dr11+2*(dr12+dr13)+dr14)/6;          r2 += dr*(dr21+2*(dr22+dr23)+dr24)/6;
}
