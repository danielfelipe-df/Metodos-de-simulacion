#include <iostream>
#include <fstream>
#include <cmath>

double f1(double r1, double r2, double r, double lambda);
double f2(double r1, double r2, double r, double lambda);
void UnPasoDeRungeKutta4(double &r1, double &r2, double &r, double &dr, double &lambda);
double Biseccion(int &a, int &b, double *f);
double f(int n,double x,double t);
double IntegralPorSimpson(int n,double x,double a,double b,int N);
double Bessel(int n,double x);


int main(void)
{
  double r1 = 1, r2 = 0, r;
  double dr = 0.001;
  double lambda = 1, dlambda, lambda_max;

  //-----Punto a-----

  std::ofstream fout("datos_2a.dat");
  for(r=dr/10; r<10; r+=dr){
    fout << r << " " << r1 << '\n';
    UnPasoDeRungeKutta4(r1,r2,r,dr,lambda);
  }
  fout.close();

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
  fout.open("datos_2b.dat");
  for(lambda=dlambda; lambda<15; lambda+=dlambda){
    r1 = 1;   r2 = 0;
    for(r=dr/10; r<1; r+=dr){
      UnPasoDeRungeKutta4(r1,r2,r,dr,lambda);
    }
    fout << lambda << " " << r1 << '\n';
  }
  fout.close();

  std::cout << "plot 'datos_2b.dat' w l, 0" << std::endl;
  std::cout << "pause 10" << std::endl;


  //----- Punto c------

  dlambda = 0.001;  lambda_max = 15;   dr=0.001;
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

  //Las variables int son las que me indican la posición que voy a mirar en el array
  //En el caso de 'a' estoy interesado que casilla ocupa el valor lambda=1 en el array. Y en el caso de 'b' es lambda=3
  double algo;   int a = 1/dlambda, b=3/dlambda;
  algo = Biseccion(a,b,data);

  //Aquí grafico en donde quedaría ubicado el cero de la función
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange[-0.6:1]" << std::endl;
  std::cout << "plot " << algo*dlambda << ",t , 'datos_2b.dat' w l , t*10,0" << std::endl;
  std::cout << "pause 5" << std::endl;

  //Definimos los intervalos donde se hallarán los ceros. Además, se define el array donde se van a guardar los lambdas a los que corresponde cada cero.
  int intervalos[6] = {1,3,6,9,12,15};
  double cero_lambda[5];

  //Hallamos cada cero y lo guardamos
  for(int i=0; i<(6-1); i++){
    a = intervalos[i]/dlambda;   b = intervalos[i+1]/dlambda;
    algo = Biseccion(a, b, data);
    cero_lambda[i] = algo*dlambda;
  }

  //Definimos las condiciones iniciales para cada lambda y resolvemos la ecuación diferencial para cada uno.
  double rr1[5] = {1,1,1,1,1}, rr2[5] = {0,0,0,0,0};
  fout.open("datos_2c.dat");
  for(r=dr/10; r<10; r+=dr){
    fout << r << " " << rr1[0] << " " << rr1[1] << " " << rr1[2] << " " << rr1[3] << " " << rr1[4] << '\n';
    UnPasoDeRungeKutta4(rr1[0],rr2[0],r,dr,cero_lambda[0]);
    UnPasoDeRungeKutta4(rr1[1],rr2[1],r,dr,cero_lambda[1]);
    UnPasoDeRungeKutta4(rr1[2],rr2[2],r,dr,cero_lambda[2]);
    UnPasoDeRungeKutta4(rr1[3],rr2[3],r,dr,cero_lambda[3]);
    UnPasoDeRungeKutta4(rr1[4],rr2[4],r,dr,cero_lambda[4]);
  }
  fout.close();

  std::cout << "plot 'datos_2c.dat' using 1:2 w l, 'datos_2c.dat' using 1:3 w l, 'datos_2c.dat' using 1:4 w l, 'datos_2c.dat' using 1:5 w l, 'datos_2c.dat' using 1:6 w l" << std::endl;
  std::cout << "pause 10" << std::endl;

  //Se hace un zoom en el intervalo x[0:2]. Aquí se ve que efectivamente todos se cruzan en r=1.
  std::cout << "set xrange[0:2]" << std::endl;
  std::cout << "set yrange[-0.6:1]" << std::endl;
  std::cout << "plot 'datos_2c.dat' using 1:2 w l, 'datos_2c.dat' using 1:3 w l, 'datos_2c.dat' using 1:4 w l, 'datos_2c.dat' using 1:5 w l, 'datos_2c.dat' using 1:6 w l" << std::endl;
  std::cout << "pause 10" << std::endl;


  //----- Punto d------

  N = int(15/dr);
  double data_bessel[N];

  //Se hallan los valores de la función de Bessel de orden cero.
  //Se utiliza el array data_bessel para utilizar la misma función Bisección. (Se puede hacer una nueva función Bisección para este caso).
  fout.open("datos_2d.dat");
  for(r=dr, i=0; r<10; r+=dr, i++){
    fout << r << " " << Bessel(0,r) << '\n';
    data_bessel[i] = Bessel(0,r);
  }
  fout.close();

  //Se hace la gŕafica de la función Bessel
  std::cout << "set xrange[0:10]" << std::endl;
  std::cout << "set yrange[-1:1]" << std::endl;
  std::cout << "plot 'datos_2d.dat' w l, 0" << std::endl;
  std::cout << "pause 10" << std::endl;

  //Se hallan los ceros de los datos data_bessel
  double cero_bessel[5];
  for(i=0; i<(6-1); i++){
    a = intervalos[i]/dlambda;   b = intervalos[i+1]/dlambda;
    algo = Biseccion(a, b, data_bessel);
    cero_bessel[i] = algo*dlambda;
  }

  //Se halla el error relativo entre los ceros
  for(i=0; i<5; i++){
    std::cout << "Cero de Bessel: " << cero_bessel[i] << "\tCero de Derivada: " << cero_lambda[i] << "\tError relativo: " << 100*std::abs(cero_bessel[i]-cero_lambda[i])/cero_bessel[i] << "%" << std::endl;
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
double Biseccion(int &a, int &b, double *f)
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
