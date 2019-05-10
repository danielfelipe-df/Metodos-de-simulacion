#include <iostream>
#include <fstream>
#include <cmath>

double f1(double r1, double r2, double r, double lambda);
double f2(double r1, double r2, double r, double lambda);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda);
double Biseccion(int a, int b, double *f);
double f(int n,double x,double t);
double IntegralPorSimpson(int n,double x,double a,double b,int N);
double Bessel(int n,double x);


int main(void)
{
  double r1 = 1, r2 = 0, r_max = 10, r;
  double dr = 0.001;
  double lambda = 1, dlambda = 0.001;

  int N = int(15/dr);
  double data_bessel[N];
  int i;

  //Se hallan los valores de la función de Bessel de orden cero.
  //Se utiliza el array data_bessel para utilizar la misma función Bisección. (Se puede hacer una nueva función Bisección para este caso).
  std::ofstream fout("datos_2d.dat");
  for(r=dr, i=0; r<(r_max+dr); r+=dr, i++){
    fout << r << " " << Bessel(0,r) << '\n';
    data_bessel[i] = Bessel(0,r);
  }
  fout.close();

  //Se hace la gŕafica de la función Bessel
  std::cout << "set terminal pdf" << std::endl;
  std::cout << "set output 'figure_2d.pdf'" << std::endl;
  std::cout << "set title 'Función de Bessel de orden 0'" << std::endl;
  std::cout << "set xlabel 'r'" << std::endl;
  std::cout << "set ylabel 'J0(r)'" << std::endl;
  std::cout << "set xrange[0:10]" << std::endl;
  std::cout << "set yrange[-1:1]" << std::endl;
  std::cout << "plot 'datos_2d.dat' w l title '', 0 title ''" << std::endl;

  //Se hallan los ceros de los datos data_bessel
  int intervalos[6] = {1,3,6,9,12,15};
  double pre_lambda, a, b;
  for(i=0; i<(6-1); i++){
    a = intervalos[i]/dlambda;   b = intervalos[i+1]/dlambda;
    pre_lambda = Biseccion(a, b, data_bessel);
    std::cout << "Cero " << i+1 << " de J0(r): " << pre_lambda*dlambda << std::endl;
  }

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

//Aquí le paso como argumentos las posiciones inciales en el array para hacer el método.
//Después le paso la dirección (tipo puntero) de la primera casilla del array.
double Biseccion(int a, int b, double *f)
{
  //Defino la variable 'm', la cual será mi variable de intercambio en el método. (Nótese que representará el número de una casilla).
  //Los punteros fa y fm serán las direcciones de la posición 'a' y 'm'
  int m;     double *fa, *fm;
  //Aquí se hace álgebra de punteros para asignarle la dirección a fa
  fa = f + a;
  //Aquí el while del método. Y lo que hago es definir que mientras la diferencia de casillas sea 5 o menor termine de buscar la raíz.
  while(b-a > 2){
    //Asignación de las variables
    m=(int)(a+b)/2;  fm = f + m;
    //Aplicación del método de bisección
    if((*fa)*(*fm)<0){
      b=m;
    }
    else{
      a=m;  fa=fm;
    }
  }
  //Devuelvo la ubicación de la raíz.
  return (a+b)/2.0;
}


double f(int n,double x,double t){return cos(n*t-x*sin(t));}

double IntegralPorSimpson(int n,double x,double a,double b,int N)
{
  int DosN=2*N; double h=(b-a)/DosN; int i; double t,suma;
  for(suma=0,i=0;i<=DosN;i++){
    t=a+h*i;
    if(i==0 || i==DosN)
      suma+=f(n,x,t);
    else if(i%2==1)
      suma+=4*f(n,x,t);
    else
      suma+=2*f(n,x,t);
  }
  return suma*h/3;
}

double Bessel(int n,double x){return 1/M_PI*IntegralPorSimpson(n,x,0,M_PI,50);}
