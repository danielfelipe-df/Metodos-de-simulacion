#include <iostream>
#include <fstream>
#include <cmath>

const double ERR =1e-7;

double f1(double r1, double r2, double r, double lambda);
double f2(double r1, double r2, double r, double lambda);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda);
double Biseccion(int &a, int &b, double *f);


int main(void)
{
  double r1 = 1, r2 = 0, r;
  double dr = 0.001;
  double lambda = 2.4048, dlambda, lambda_max;

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

  dlambda = 0.001;
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
  dlambda = 0.001;  lambda_max = 15;
  //Este N se define de esta forma porque lo que hago es generar un array que tenga el tamaño necesario para guardar los datos
  int N = (int)lambda_max/dlambda;
  double data[N];
  int i;
  for(lambda=dlambda, i=0; lambda<lambda_max; lambda+=dlambda, i++){
    r1 = 1;   r2 = 0;
    for(r=dr/10; r<1; r+=dr){
      UnPasoDeRungeKutta4(r1,r2,r,dr,lambda);
    }
    //Aquí se guarda cada dato de la función en un array
    data[i] = r1;
  }

  std::cout << N << std::endl;
  //Las variables int son las que me indican la posición que voy a mirar en el array
  //En el caso de 'a' estoy interesado que casilla ocupa el valor lambda=3 en el array. Y en el caso de 'b' es lambda=6
  double algo;   int a = 3/dlambda, b=6/dlambda;
  algo = Biseccion(a,b,data);
  //Aquí imprimo, inicialmente, el número de la casilla, y después el valor de lambda al que corresponde
  std::cout << algo << " " << algo*dlambda << std::endl;
  int cosa = (int)algo;
  //Aquí imprimo el valor del cero
  std::cout << data[cosa] << std::endl;
  //Aquí grafico en donde quedaría ubicado el cero de la función
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange[-0.6:1]" << std::endl;
  std::cout << "plot " << algo*dlambda << ",t , 'datos_2b.dat' w l , t*10,0" << std::endl;
  std::cout << "pause 5" << std::endl;

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

  r1 += dr*(dr11+2*(dr12+dr13)+dr14)/6;          r2 += dr*(dr21+2*(dr22+dr23)+dr24)/6;
}

//Aquí le paso como argumentos las posiciones inciales en el array para hacer el método.
//Después le paso la dirección (tipo puntero) de la primera casilla del array.
double Biseccion(int &a, int &b, double *f)
{
  //Defino la variable 'm', la cual será mi variable de intercambio en el método. (Nótese que representará el número de una casilla).
  //Los punteros fa y fm serán las direcciones de la posición 'a' y 'm'
  int m;     double *fa, *fm;
  //Aquí se hace álgebra de punteros para asignarle la dirección a fa
  fa = f + a;
  //Aquí el while del método. Y lo que hago es definir que mientras la diferencia de casillas sea 5 o menor termine de buscar la raíz.
  while(b-a > 5){
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
