#include <iostream>
#include <fstream>
#include <cmath>

double f1(double r1, double r2, double r, double lambda);
double f2(double r1, double r2, double r, double lambda);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda);
double Biseccion(int a, int b, double *f);

int main(void)
{
  double r1 = 1, r2 = 0, r;
  double dr = 0.001;
  double lambda = 1, dlambda=0.001, lambda_max=15;

  //----- Punto c------

  //Este N se define de esta forma porque lo que hago es generar un array que tenga el tamaño necesario para guardar los datos
  int N = int(lambda_max/dlambda);
  double data[N];    int i;
  for(lambda=dlambda, i=0; lambda<lambda_max; lambda+=dlambda, i++){
    r1 = 1;   r2 = 0;
    for(r=dr/10; r<1; r+=dr){
      UnPasoDeRungeKutta4(r1,r2,r,dr,lambda);
    }
    //Aquí se guarda cada dato de la función en un array
    data[i] = r1;
  }

  //Definimos los intervalos donde se hallarán los ceros. Además, se define el array donde se van a guardar los lambdas a los que corresponde cada cero.
  int intervalos[6] = {1,3,6,9,12,15};
  double cero_lambda[5];
  double pre_lambda, a, b;

  //Hallamos cada cero y lo guardamos
  for(int i=0; i<(6-1); i++){
    a = intervalos[i]/dlambda;   b = intervalos[i+1]/dlambda;
    pre_lambda = Biseccion(a, b, data);
    cero_lambda[i] = pre_lambda*dlambda;
  }

  //Definimos las condiciones iniciales para cada lambda y resolvemos la ecuación diferencial para cada uno.
  double rr1[5] = {1,1,1,1,1}, rr2[5] = {0,0,0,0,0};
  double r_max = 1.2;

  std::ofstream fout("datos_2c.dat");
  for(r=dr; r<(r_max+dr); r+=dr){
    fout << r << " " << rr1[0] << " " << rr1[1] << " " << rr1[2] << " " << rr1[3] << " " << rr1[4] << '\n';
    UnPasoDeRungeKutta4(rr1[0],rr2[0],r,dr,cero_lambda[0]);
    UnPasoDeRungeKutta4(rr1[1],rr2[1],r,dr,cero_lambda[1]);
    UnPasoDeRungeKutta4(rr1[2],rr2[2],r,dr,cero_lambda[2]);
    UnPasoDeRungeKutta4(rr1[3],rr2[3],r,dr,cero_lambda[3]);
    UnPasoDeRungeKutta4(rr1[4],rr2[4],r,dr,cero_lambda[4]);
  }
  fout.close();

  std::cout << "set terminal pdf" << std::endl;
  std::cout << "set output 'figure_2c.pdf'" << std::endl;
  std::cout << "set title 'Imagen del ejercio 2c'" << std::endl;
  std::cout << "set xlabel 'r'" << std::endl;
  std::cout << "set ylabel 'R(r)'" << std::endl;
  std::cout << "plot 0 title '' ,'datos_2c.dat' using 1:2 w l title 'lambda = " << cero_lambda[0] << "', ";
  std::cout << "'datos_2c.dat' using 1:3 w l title 'lambda = " << cero_lambda[1] << "', ";
  std::cout << "'datos_2c.dat' using 1:4 w l title 'lambda = " << cero_lambda[2] << "', ";
  std::cout << "'datos_2c.dat' using 1:5 w l title 'lambda = " << cero_lambda[3] << "', ";
  std::cout << "'datos_2c.dat' using 1:6 w l title 'lambda = " << cero_lambda[4] << std::endl;

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
