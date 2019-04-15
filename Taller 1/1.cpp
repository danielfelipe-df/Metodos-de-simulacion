#include <iostream>
#include <fstream>
#include <cmath>

const double my_gamma = 0.08;

double f1(double s, double i, double t, double beta);
double f2(double s, double i, double t, double beta);
void UnPasoDeRungeKutta4(double &s, double &i, double &t, double &dt, double &beta);

int main(void)
{
  double s = 0.999, i = 0.001;
  double dt = 0.1;
  double beta;

  //-----Punto a-----

  std::ofstream fout1("datos_1a.dat");
  beta = 0.35;
  for(double t=0; t<1000; t+=dt){
    fout1 << t << " " << s << " " << i << " " << 1-s-i << '\n';
    UnPasoDeRungeKutta4(s,i,t,dt,beta);
  }
  fout1.close();
  
  std::cout << "plot 'datos_1a.dat' using 1:2 w l, 'datos_1a.dat' using 1:3 w l, 'datos_1a.dat' using 1:4 w l " << std::endl;
  std::cout << "pause 3" << std::endl;

  //-----Punto c-----
  
  /*
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula.gif'" << std::endl;
  std::cout << "unset key" << std::endl;
  std::cout << "set xrange[0:1000]" << std::endl;
  std::cout << "set yrange[0:1]" << std::endl;
  */

  double beta_max = 0.4;
  std::ofstream fout3("datos_1c.dat");
  for(beta=0.1; beta<beta_max; beta+=0.0005){
    s=0.999;    i=0.001;
    //std::ofstream fout2("datos_1c.dat");
    for(double t=0; t<1000; t+=dt){
      //fout2 << t << " " << s << " " << i << " " << 1-s-i << '\n';
      UnPasoDeRungeKutta4(s,i,t,dt,beta);
    }
    //fout2.close();
    fout3 << beta/my_gamma << " " << s << std::endl;
    //std::cout << "plot 'datos_1c.dat' using 1:2 w l, 'datos_1c.dat' using 1:3 w l, 'datos_1c.dat' using 1:4 w l" << std::endl;
  }
  fout3.close();
  
  std::cout << "plot 'datos_1c.dat' w l" << std::endl;
  std::cout << "pause 3" << std::endl;

  return 0;
}


double f1(double s, double i, double t, double beta){return -beta*s*i;}

double f2(double s, double i, double t, double beta){return beta*s*i - my_gamma*i;}

void UnPasoDeRungeKutta4(double &s, double &i, double &t, double &dt, double &beta)
{
  double ds1, ds2, ds3, ds4;
  double di1, di2, di3, di4;
  
  ds1=f1(s,i,t,beta)*dt;                                 di1=f2(s,i,t,beta)*dt;
  ds2=f1(s+0.5*ds1, i+0.5*di1, t+0.5*dt, beta)*dt;       di2=f2(s+0.5*ds1, i+0.5*di1, t+0.5*dt, beta)*dt;
  ds3=f1(s+0.5*ds2, i+0.5*di2, t+0.5*dt, beta)*dt;       di3=f2(s+0.5*ds2, i+0.5*di2, t+0.5*dt, beta)*dt;
  ds4=f1(s+ds3, i+di3, t+dt, beta)*dt;                   di4=f2(s+ds3, i+di3, t+dt, beta)*dt;

  t += dt;  s += (ds1+2*(ds2+ds3)+ds4)/6;          i += (di1+2*(di2+di3)+di4)/6;
}
