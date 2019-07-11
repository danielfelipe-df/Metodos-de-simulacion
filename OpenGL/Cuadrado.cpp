#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/gl.h>
#include <GL/glut.h>

const int Lx=400;
const int Ly=400;

const int Q=5;
const double W0=1.0/3;
const double k=0.0;

const double C=0.5; // C<0.707 celdas/click
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;


bool Cuadrado(int x, int y, int ix, int iy, int L);
bool Circunferencia(int x, int y, int ix, int iy, int R);
bool Elipse(int x, int y, double ix, double iy, double Ra, double Rb);

class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0,double Jx0,double Jy0);
  void ImponerCampos(int t, int ix, int iy);
  void Imprimase(const char * NombreArchivo);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1-W0)/4.0;
  //Cargar los vectores
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[0][i]; else suma+=f[ix][iy][i]*V[0][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]*V[1][i]; else suma+=f[ix][iy][i]*V[1][i];
  return suma;
}
double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i){
  if(i==0)
    return rho0*AUX0;
  else
    return w[i]*(TresC2*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
}
void LatticeBoltzmann::Colisione(){
  int ix,iy,i; double rho0,Jx0,Jy0;
  int cuadrado_x=50, cuadrado_y=50, lado=25;
  int centro_x=100, centro_y=100, radio=25;
  double elipse_x=50, elipse_y=150, radio_a=25, radio_b=30;

  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly; iy++){
      rho0=rho(ix,iy,false);  Jx0=Jx(ix,iy,false);  Jy0=Jy(ix,iy,false);
      if(Cuadrado(ix,iy,cuadrado_x,cuadrado_y,lado) || Circunferencia(ix,iy,centro_x,centro_y,radio) || Elipse(ix,iy,elipse_x,elipse_y,radio_a,radio_b)){
        fnew[ix][iy][0]=UmUtau*f[ix][iy][0]+Utau*feq(rho0,Jx0,Jy0,0);
        fnew[ix][iy][3]=k*f[ix][iy][1];  fnew[ix][iy][1]=k*f[ix][iy][3];
        fnew[ix][iy][4]=k*f[ix][iy][2];  fnew[ix][iy][2]=k*f[ix][iy][4];
      }
      else
        {for(i=0; i<Q; i++){fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Jx0,Jy0,i);}}
    }
  }
  /*
  for(ix=cuadrado_x; ix<=cuadrado_x+lado; ix++){
    for(iy=cuadrado_y; iy<=cuadrado_y+lado; iy++){
      rho0=rho(ix,iy,false);  Jx0=Jx(ix,iy,false);  Jy0=Jy(ix,iy,false);
      if(Cuadrado(ix,iy,cuadrado_x,cuadrado_y,lado)){
        fnew[ix][iy][0]=UmUtau*f[ix][iy][0]+Utau*feq(rho0,Jx0,Jy0,0);
        fnew[ix][iy][3]=k*f[ix][iy][1];  fnew[ix][iy][1]=k*f[ix][iy][3];
        fnew[ix][iy][4]=k*f[ix][iy][2];  fnew[ix][iy][2]=k*f[ix][iy][4];
      }
    }
  }*/
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
        if(ix+V[0][i]<Lx && ix+V[0][i]>=0 && iy+V[1][i]<Ly && iy+V[1][i]>=0){
          f[ix+V[0][i]][iy+V[1][i]][i]=fnew[ix][iy][i];
        }
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
        f[ix][iy][i]=feq(rho0,Jx0,Jy0,i);
}
void LatticeBoltzmann::ImponerCampos(int t, int ix, int iy){
  int i; double lambda,omega,rho0,Jx0,Jy0;
  lambda=10.0; omega=2*M_PI*C/lambda;
  rho0=30*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][i]=feq(rho0,Jx0,Jy0,i);
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double rho0, Jx0, Jy0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true); //Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      MiArchivo<<ix<<'\t'<<iy<<'\t'<<rho0<<'\n';
      //MiArchivo<<ix<<" "<<iy<<" "<<Jx0<<" "<<Jy0<<endl;
    }
    MiArchivo<<'\n';
  }
  MiArchivo.close();
}

bool Cuadrado(int x, int y, int ix, int iy, int L)
{
  if(x==ix && y>=iy && y<=iy+L){return true;}//Costado de la izquierda
  else if(x==ix+L && y>=iy && y<=iy+L){return true;}//Costado de la derecha
  else if(y==iy && x>=ix && x<=ix+L){return true;}//Costado de abajo
  else if(y==iy+L && x>=ix && x<=ix+L){return true;}//Costado de arriba
  else{return false;}
}

bool Circunferencia(int x, int y, int ix, int iy, int R)
{
  if((x-ix)*(x-ix) + (y-iy)*(y-iy) <= R*R){return true;}
  else{return false;}
}

bool Elipse(int x, int y, double ix, double iy, double Ra, double Rb)
{
  if(((x-ix)*(x-ix)/(Ra*Ra)) + ((y-iy)*(y-iy)/(Rb*Rb)) <= 1){return true;}
  else{return false;}
}

//----------- Funciones Globales ----------

double rho1[Lx][Ly];
LatticeBoltzmann Ondas;
int t,t_otro;

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);
  //glColor3f(0.5, 1.0, 1.0);

  Ondas.Colisione();
  Ondas.ImponerCampos(t,5,67);
  Ondas.Adveccione();
  double rho0, normalizacion=0.5;
  glPointSize(3.0);
  glBegin(GL_POINTS);
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
      rho0 = rho1[ix][iy];
      if(rho0>normalizacion){glColor3f(0.0, 0.0, 1.0);}
      else if(rho0<-normalizacion){glColor3f(1.0, 0.0, 0.0);}
      else if(rho0>0){glColor3f(1.0-rho0/normalizacion, 1.0-rho0/normalizacion, 1.0);}
      else if(rho0<0){glColor3f(1.0, 1.0+rho0/normalizacion, 1.0+rho0/normalizacion);}
      else{glColor3f(1.0, 1.0, 1.0);}
      glVertex3f(ix*0.0025,iy*0.0025,0);
    }
  }

  //glVertex3d(0.25, 0.25, 1.0);
  //  glVertex3d(0.75, 0.25, 0.0);
  //glVertex3d(0.75, 0.75, 1.0);
  //glVertex3d(0.25, 0.75, 0.0);
  glEnd();
  glFlush();
}

void init(void)
{
  /* select clearing (background) color
  glClearColor(0.0, 0.0, 0.0, 0.0);
  /*
  /*initialize viewing values */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
}

void AmplitudDisplay(void)
{
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
      rho1[ix][iy] = Ondas.rho(ix,iy,true);
    }
  }
  t++;  t_otro++;
  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
  switch (button) {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN)
      glutIdleFunc(AmplitudDisplay);
      break;
    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN)
      glutIdleFunc(NULL);
      break;
    default:
      break;
  }
}


int main(int argc, char** argv)
{
  //OpenGL

  Ondas.Inicie(0,0,0);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(600, 600);
  glutInitWindowPosition(500, 100);
  glutCreateWindow("Animation");
  init();
  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  glutMainLoop();

  //Gnuplot
  /*
  //std::cout <<"set terminal gif animate" << std::endl;
  //std::cout << "set output 'pelicula_3.gif'" << std::endl;
  std::cout << "set pm3d map" << std::endl;
  std::cout << "set palette defined (-1 \"red\", 0 \"white\", 1 \"blue\")" << std::endl;
  //std::cout << "set cbrange[-0.05:0.05]" << std::endl;
  std::cout << "set xrange[0:200]; set yrange[0:80]; set zrange[-1:1]" << std::endl;

  int tmax=600;
  int pared_x[4]={0,50,100,200};
  int huecos_y[5]={0,50,100,150,200};
  Ondas.Inicie(0,0,0);

  for(t=0;t<tmax;t++){
    Ondas.Colisione(pared_x,4,huecos_y,5);
    //if(t<20)
    Ondas.ImponerCampos(t,5,40);
    //Ondas.ImponerCampos(t,5,37);
    Ondas.Adveccione();
    Ondas.Imprimase("Ondas.dat");
    std::cout << "splot 'Ondas.dat' title ''" << std::endl;
  }
  //Ondas.Imprimase("Ondas.dat");
  */
  return 0;
}
