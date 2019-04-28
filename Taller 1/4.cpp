#include <iostream>
#include <cmath>

void Iniciar();

class Pendulo
{
private:
  double r, theta, l, x0;

public:
  void First_ubication(double x0, double r0, double theta0, double l0);
  void Dibujarse();
};


int main(void)
{
  Pendulo objeto[3];
  double x = 3, r0 = 1, theta0 = M_PI/2, l0 = 2;
  Iniciar();
  for(int i=0, x=0; i<4; x+=2*r0, i++){
    objeto[i].First_ubication(x, r0, theta0, l0);
    objeto[i].Dibujarse();
  }
  std::cout << std::endl;
  std::cout << "pause 10" << std::endl;
  return 0;
}

//------- Funciones de la clase Pendulo --------
void Pendulo::First_ubication(double x, double r0, double theta0, double l0)
{
  r = r0;  theta = theta0;  x0 = x;   l = l0;
}

void Pendulo::Dibujarse()
{
  std::cout << " , " << r << "*cos(t) + " << x0 << " + " << l << "*sin(" << theta << ")," << r << "*sin(t) + " << l << "*(1+cos(" << theta << "))";
  std::cout << " , " << x0 << "," << l/10 << "*t + " << r;
}

//------ Funciones normales ------
void Iniciar()
{
  std::cout << "set xrange[-10:10]" << std::endl;
  std::cout << "set yrange[-10:10]" << std::endl;
  std::cout << "set parametric" << std::endl;
  std::cout << "set trange[0:10]" << std::endl;
  std::cout << "plot 0,0";
}
