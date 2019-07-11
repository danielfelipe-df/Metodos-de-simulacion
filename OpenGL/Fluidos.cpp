#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/gl.h>
#include <GL/glut.h>
using namespace std;

/*
  La animación que se ve bien pero es lenta tiene:
  tau =0.505;  Lx=900; Ly=200; iyc=100;  R=5;
  Rhoinicial=1.0; Vventilador=0.1;

*/
const int Lx=900;
const int Ly=200;

const int Q=9;

const double tau=0.6;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

double Rhoinicial=1.0, Vventilador=0.1, Rhocontinua=1.0;

bool Cuadrado(int x, int y, int ix, int iy, int L);
bool Circunferencia(int x, int y, int ix, int iy, int R);
bool Elipse(int x, int y, double ix, double iy, double Ra, double Rb);
bool Ala(int x, int y, double ix, double iy,double lx, double ly);
bool Corazon(int x, int y, double ix, double iy, double lx, double ly);

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
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Colisione(double Vventilador);
  void Adveccione(void);
  void Inicie(double rho0,double Ux0,double Uy0);
  //void ImponerCampos(double Vventilador);
  void Imprimase(const char * NombreArchivo,double Vventilador);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=4/9.0; 
  w[1]=w[2]=w[3]=w[4]=1/9.0;
  w[5]=w[6]=w[7]=w[8]=1/36.0;
  //Cargar los vectores
  V[0][0]=0;
  V[1][0]=0;

  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;

  V[0][5]=1;  V[0][6]=-1; V[0][7]=-1; V[0][8]=1;
  V[1][5]=1;  V[1][6]=1;  V[1][7]=-1; V[1][8]=-1;
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
double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
void LatticeBoltzmann::Colisione(double Vventilador){
  int ix,iy,i; double rho0,Ux0,Uy0;

  for(ix=0; ix<3; ix++){
    for(iy=0; iy<Ly; iy++){
      for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(Rhocontinua,Vventilador,0,i);
    }
  }
  //Para cada celda
  for(ix=3;ix<250;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscópicas
      rho0=rho(ix,iy,false);  Ux0=Jx(ix,iy,false)/rho0;  Uy0=Jy(ix,iy,false)/rho0;
      //if(Elipse(ix,iy,30.0,100.0,20.0,50.0))
      //if(Cuadrado(ix,iy,10,50,100))
      //if(Circunferencia(ix,iy,60.0,100.0,50.0))
      if(Ala(ix,iy,10.0,100.0,200.0,50.0))
      //if(Corazon(ix,iy,10,100,100,50))
        for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);
      else if(ix==50 && iy==100+50+1)
        for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);
      else
      for(i=0;i<Q;i++)  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i);
    }

  for(ix=250;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      //Calcular las cantidades macroscópicas
      rho0=rho(ix,iy,false);  Ux0=Jx(ix,iy,false)/rho0;  Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++)
        fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i);
    }
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
        if(ix+V[0][i]<Lx && ix+V[0][i]>=0 && iy+V[1][i]<Ly && iy+V[1][i]>=0){
          f[ix+V[0][i]][iy+V[1][i]][i]=fnew[ix][iy][i];
        }
}
void LatticeBoltzmann::Inicie(double rho0,double Ux0,double Uy0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
        f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
}
/*
void LatticeBoltzmann::ImponerCampos(double Vventilador){
  int i,ix,iy; double rho0; int ixc=60,iyc=100; int R=50, R2=R*R;
  
  for(ix=0;ix<250;ix++)
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //ventilador
      if(ix==1 || ix==2)
        for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(Rhocontinua,Vventilador,0,i);
      else if(Elipse(ix,iy,30.0,100.0,20.0,50.0))
      //else if(Cuadrado(ix,iy,10,50,100))
      //else if(Circunferencia(ix,iy,ixc,iyc,R))
      //else if(Ala(ix,iy,10.0,100.0,200.0,50.0))
      //else if(Corazon(ix,iy,10,100,100,50))
        for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);
      /*else if(ix==ixc && iy==iyc+150+3)
        for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);*//*
      else if(iy == 0 || iy ==Ly-1)
        for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,0,0,i);
    }
}
*/
void LatticeBoltzmann::Imprimase(const char * NombreArchivo,double Vventilador){
  ofstream MiArchivo(NombreArchivo); double rho0,Ux0,Uy0, vx, vy;
  for(int ix=0;ix<Lx;ix+=4){
    for(int iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true);  Ux0=Jx(ix,iy,true)/rho0;  Uy0=Jy(ix,iy,true)/rho0;
      vx = 4*(Ux0+Vventilador)/Vventilador; vy = 4*Uy0/Vventilador;
      //MiArchivo<<ix<<" "<<iy<<" "<<4*(Ux0-Vventilador)/Vventilador<<" "<<4*Uy0/Vventilador<<endl;
      //MiArchivo << ix << " " << iy << " " << rho0 << endl;
      MiArchivo << ix << " " << iy << " " << sqrt(vx*vx + vy*vy) << endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

bool Cuadrado(int x, int y, int ix, int iy, int L)
{
  if(x>=ix && x<=ix+L && y>=iy && y<=iy+L){return true;}
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

bool Ala(int x, int y, double ix, double iy,double lx, double ly)
{
  if(x>=ix && x<=lx+ix){
    double sup = ly*sin(M_PI*x/lx);
    if(y <= sup+iy && y>= iy-sup){return true;}
    else{return false;}
  }
  else{return false;}
}

bool Corazon(int x, int y, double ix, double iy, double lx, double ly)
{
  if(x>=ix && x<=lx+ix){
    double sup = ly*abs(sin(2*M_PI*(x-ix)/lx));
    double inf = ly*0.03*(abs((x-ix)-0.5*lx)-0.5*lx);
    if(y<=sup+iy && y>=inf+iy){return true;}
    else{return false;}
  }
  else{return false;}
}

//---------- Funciones Generales --------

double rho1[Lx][Ly], Jx1[Lx][Ly], Jy1[Lx][Ly];
LatticeBoltzmann Aire;
int t=0;

void RenderString(float x, float y, void *font, const unsigned char* string)
{
  char *c;

  //glColor3f(rgb.r, rgb.g, rgb.b);
  glRasterPos2f(x, y);

  glutBitmapLength(font, string);
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);

  double rho0, Jx0, Jy0, Ux0, Uy0, normalization=3.0;
  double x,y,r;
  double red0 = 160.0/256.0, green0 = 32.0/256.0, blue0 = 240.0/256.0;
  double red1 = 256.0/256.0, green1 = 256.0/256.0, blue1 = 0.0;
  Aire.Colisione(Vventilador);
  //Aire.ImponerCampos(Vventilador);
  Aire.Adveccione();
  glPointSize(1.0);
  glBegin(GL_POINTS);
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
      rho0 = rho1[ix][iy];  Jx0 = Jx1[ix][iy];  Jy0 = Jy1[ix][iy];
      Uy0 = Jx0/rho0; Ux0/rho0;
      x = (Ux0-Vventilador)/Vventilador;  y=Uy0/Vventilador;
      r=sqrt(x*x+y*y)/normalization;
      if(r>1.0){glColor3f(red0,green0,blue0);}
      else{glColor3f(red0*r,green0*r,blue0*r);}//red0*r,green0*r,blue0*r para una combinación linda
      glVertex3f(ix*0.0011,iy*0.005,0);
    }
  }
  glEnd();
  //std::cout << t << std::endl;
  glFlush();
}

void init(void)
{
  /* select clearing (background) color
  glClearColor(0.0, 0.0, 0.0, 0.0);
  */
  /*initialize viewing values */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
}

void AmplitudDisplay(void)
{
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
      rho1[ix][iy] = Aire.rho(ix,iy,true);
      Jx1[ix][iy] = Aire.Jx(ix,iy,true);
      Jy1[ix][iy] = Aire.Jy(ix,iy,true);
    }
  }
  t++;
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

  Aire.Inicie(Rhoinicial,Vventilador,0);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(900, 200);
  glutInitWindowPosition(500, 100);
  glutCreateWindow("Ala");
  init();
  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  //const unsigned char * strings = reinterpret_cast<const unsigned char *>( "123" );
  //RenderString(600.0f, 100.0f, GLUT_BITMAP_TIMES_ROMAN_24, strings);
  glutMainLoop();
  //std::cout << "Elipse: " << t << std::endl;

  //Gnuplot
  /*
  int t,tmax=200;
  double RHOinicial=1.0, Vventilador=0.1;
  
  Aire.Inicie(RHOinicial,Vventilador,0);
  std::cout << "set pm3d map" << std::endl;
  //std::cout << "set term gif animate" << std::endl;
  //std::cout << "set output 'fluido_2.gif'" << std::endl;
  std::cout << "set parametric; set urange[0:2*pi]; set vrange[0:2*pi]" << std::endl;
  std::cout << "set xrange[0:300]; set yrange[0:70]" << std::endl;

  for(t=0;t<tmax;t++){
    Aire.Colisione();
    Aire.ImponerCampos(Vventilador);
    Aire.Adveccione();
    Aire.Imprimase("Aire.dat",Vventilador);
    //std::cout << "plot 'Aire.dat' w vec, 12*cos(t)+32,12*sin(t)+32" << std::endl;
    std::cout << "splot 'Aire.dat', 12*cos(u)+32,12*sin(v)+32,0" << std::endl;
  }
  
  Aire.Imprimase("Aire.dat",Vventilador);
  */
  return 0;
}
