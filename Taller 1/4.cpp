#include <iostream>
#include <cmath>


void Iniciar();

class Pendulo
{
private:
  double r, l, x0;

public:
  double theta;
  void First_ubication(double x0, double r0, double theta0, double l0);
  void Dibujarse();
};


int main(void)
{
  Pendulo objeto[4];
  int i;
  double x0 = 3, theta0 = -M_PI/2, r0 = 1, l0 = 2;
  Iniciar();

  for(i=0, x0=0; i<4; x0+=2*r0, i++){
    objeto[i].First_ubication(x0,r0,theta0,l0);
  }

  for(theta0=-M_PI/2; theta0<=M_PI/2; theta0+=0.01){
    std::cout << "plot 0,0";
    for(i=0; i<4; i++){
      objeto[i].theta = theta0;
      objeto[i].Dibujarse();
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  //std::cout << "pause 10" << std::endl;
  return 0;
}

//------- Funciones de la clase Pendulo --------
void Pendulo::First_ubication(double x, double r0, double theta0, double l0)
{
  r = r0;  theta = theta0;  x0 = x;   l = l0;
}

void Pendulo::Dibujarse()
{
  std::cout << " , " << r << "*cos(2*pi*t) + " << x0 << " + " << l << "*sin(" << theta << ")," << r << "*sin(2*pi*t) - " << l << "*(cos(" << theta << ")-1)";
  std::cout << " , " << l << "*sin(" << theta << ")*(1-t) + " << x0 << "," << l+r << "-" << l << "*cos(" << theta << ")*(1-t)";
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
