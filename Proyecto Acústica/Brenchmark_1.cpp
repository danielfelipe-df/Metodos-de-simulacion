#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/gl.h>
#include <GL/glut.h>

const int Lx=200;
const int Ly=1;
const int Lz=1;

const int Q=7;
const double W0=1.0/4;

const double C=0.5; // C<0.707 celdas/click
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double w[Q];
  int V[3][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy, V[2][i]=V_iz
  double f[Lx][Ly][Lz][Q], fnew[Lx][Ly][Lz][Q]; // f[ix][iy][iz][iz][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,int iz,bool UseNew);
  double Jx(int ix,int iy,int iz,bool UseNew);
  double Jy(int ix,int iy,int iz,bool UseNew);
  double Jz(int ix,int iy,int iz,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,double Jz0,int i);
  void Colisione(double k);
  void Adveccione(void);
  void Inicie(double rho0,double Jx0,double Jy0, double Jz0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo);
  void Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=w[5]=w[6]=W0/2.0;
  //Cargar los vectores
  V[0][0]=0;  V[1][0]=0;  V[2][0]=0;

  V[0][1]=1;  V[0][2]=-1;  V[0][3]=0;  V[0][4]=0;  V[0][5]=0;  V[0][6]=0;
  V[1][1]=0;  V[1][2]=0;  V[1][3]=1;  V[1][4]=-1; V[1][5]=0;  V[1][6]=0;
  V[2][1]=0;  V[2][2]=0;  V[2][3]=0;  V[2][4]=0;  V[2][5]=1;  V[2][6]=-1;
}
double LatticeBoltzmann::rho(int ix,int iy, int iz,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]; else suma+=f[ix][iy][iz][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy, int iz,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]*V[0][i]; else suma+=f[ix][iy][iz][i]*V[0][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy,int iz,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]*V[1][i]; else suma+=f[ix][iy][iz][i]*V[1][i];
  return suma;
}
double LatticeBoltzmann::Jz(int ix,int iy,int iz,bool UseNew){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]*V[2][i]; else suma+=f[ix][iy][iz][i]*V[2][i];
  return suma;
}
double LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,double Jz0,int i){
  if(i==0)
    return rho0*(1-TresC2);
  else
    return 4*w[i]*(C*C*rho0+(V[0][i]*Jx0+V[1][i]*Jy0+V[2][i]*Jz0));
}
void LatticeBoltzmann::Colisione(double k){
  int ix,iy=0,iz=0,i; double rho0,Jx0,Jy0,Jz0;
  //Para cada celda

  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Calcular las cantidades macroscópicas
        rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);

        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
        if(ix==Lx-1){fnew[ix][iy][iz][2]=k*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k*f[ix][iy][iz][2];}
        else{fnew[ix][iy][iz][1]=UmUtau*f[ix][iy][iz][1]+Utau*feq(rho0,Jx0,Jy0,Jz0,1);
            fnew[ix][iy][iz][2]=UmUtau*f[ix][iy][iz][2]+Utau*feq(rho0,Jx0,Jy0,Jz0,2);}
        for(i=3; i<Q; i++){
          fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);
        }
        /*
        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
        if(ix==0 || ix==Lx-1){fnew[ix][iy][iz][2]=k*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k*f[ix][iy][iz][2];}
        else{fnew[ix][iy][iz][1]=UmUtau*f[ix][iy][iz][1]+Utau*feq(rho0,Jx0,Jy0,Jz0,1);
            fnew[ix][iy][iz][2]=UmUtau*f[ix][iy][iz][2]+Utau*feq(rho0,Jx0,Jy0,Jz0,2);}
        /*
        if(iy==0 || iy==Ly-1){fnew[ix][iy][iz][4]=0; fnew[ix][iy][iz][3]=0;}
        else{fnew[ix][iy][iz][3]=UmUtau*f[ix][iy][iz][3]+Utau*feq(rho0,Jx0,Jy0,Jz0,3);
            fnew[ix][iy][iz][4]=UmUtau*f[ix][iy][iz][4]+Utau*feq(rho0,Jx0,Jy0,Jz0,4);}

        if(iz==0 || iz==Lz-1){fnew[ix][iy][iz][5]=0; fnew[ix][iy][iz][6]=0;}
        else{fnew[ix][iy][iz][5]=UmUtau*f[ix][iy][iz][5]+Utau*feq(rho0,Jx0,Jy0,Jz0,5);
            fnew[ix][iy][iz][6]=UmUtau*f[ix][iy][iz][6]+Utau*feq(rho0,Jx0,Jy0,Jz0,6);}
        *//*
        fnew[ix][iy][iz][3]=UmUtau*f[ix][iy][iz][3]+Utau*feq(rho0,Jx0,Jy0,Jz0,3);
        fnew[ix][iy][iz][4]=UmUtau*f[ix][iy][iz][4]+Utau*feq(rho0,Jx0,Jy0,Jz0,4);
        fnew[ix][iy][iz][5]=UmUtau*f[ix][iy][iz][5]+Utau*feq(rho0,Jx0,Jy0,Jz0,5);
        fnew[ix][iy][iz][6]=UmUtau*f[ix][iy][iz][6]+Utau*feq(rho0,Jx0,Jy0,Jz0,6);*/
      }

    /*
    for(ix=0;ix<Lx;ix++){
      rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
      fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
      if(ix==0 || ix==Lx-1){fnew[ix][iy][iz][2]=k*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k*f[ix][iy][iz][2];}
      else{fnew[ix][iy][iz][1]=UmUtau*f[ix][iy][iz][1]+Utau*feq(rho0,Jx0,Jy0,Jz0,1);
          fnew[ix][iy][iz][2]=UmUtau*f[ix][iy][iz][2]+Utau*feq(rho0,Jx0,Jy0,Jz0,2);}
    }
    */

}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
        for(int i=0;i<Q;i++)
          if(ix+V[0][i]>=0){
            f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][(iz+V[2][i]+Lz)%Lz][i]=fnew[ix][iy][iz][i];
          }
  /*
  for(int ix=0; ix<Lx; ix++){
    for(int i=0; i<Q; i++){
      if(ix+V[0][i]<Lx && ix+V[0][i]>=0){
        f[ix+V[0][i]][0][0][i]=fnew[ix][0][0][i];
      }
    }
  }*/
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0, double Jz0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
        for(int i=0;i<Q;i++)
          f[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix=10,iy=0,iz=0; double lambda,omega,rho0,Jx0,Jy0,Jz0;
  lambda=50; omega=2*M_PI*C/lambda; //ix=Lx/2; //iy=Ly/2; iz=Lz/2;
  rho0=sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0,Jz0;
  for(int ix=0;ix<Lx;ix++){
    //for(int iy=0;iy<Ly;iy++){
      //for(int iz=0;iz<Lz;iz++){
        int iz=0, iy=0;
        rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        MiArchivo<<ix<<'\t'<<rho0<<'\n';
      //}
      //MiArchivo<<'\n';
    //}
    //MiArchivo<<'\n';
  }
  MiArchivo.close();
}
void LatticeBoltzmann::Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo){
  double rho0 = rho(ix,iy,iz,false);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << rho0 << '\n';
  ofs.close();
}

//------ Funciones Generales-----

const int number=5;
double distance = 1.0/(number+1.0);
double rho1[number][Lx][Ly][Lz];
LatticeBoltzmann Ondas[number];
int t,t_otro;
double k[number]={1.0,0.75,0.5,0.25,0.0};

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);

  int iz=0; double rho0;
  for(int i=0; i<number; i++){
    Ondas[i].Colisione(k[i]);
    Ondas[i].ImponerCampos(t);
    Ondas[i].Adveccione();
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for(int ix=0; ix<Lx; ix++){
      for(int iy=0; iy<Ly; iy++){
        rho0 = rho1[i][ix][iy][iz];
        glColor3f(1.0,1.0,1.0);
        glVertex3f(ix*0.005,rho0*0.05+(distance*(i+1)),0);
      }
    }
    glEnd();
  }

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
      for(int iz=0; iz<Lz; iz++){
        for(int i=0; i<number; i++){
          rho1[i][ix][iy][iz] = Ondas[i].rho(ix,iy,iz,true);
        }
      }
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

  for(int i=0; i<number; i++){
    Ondas[i].Inicie(0,0,0,0);
  }

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(600, 50);
  glutInitWindowPosition(500, 100);
  glutCreateWindow("Animation");
  init();
  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  glutMainLoop();


  //Gnuplot
  /*
  int t,tmax=1000;
  double k=0, dk=0.05;

  //std::cout << "set terminal gif animate" << std::endl;
  //std::cout << "set output 'pelicula2.gif'" << std::endl;
  //std::cout << "set pm3d explicit" << std::endl;
  //std::cout << "set palette defined (-0.000002 \"red\", 0 \"white\", 0.000002 \"blue\")" << std::endl;
  //std::cout << "set cbrange[-0.000002:0.000002]" << std::endl;
  //std::cout << "set xrange[-1:100]; set yrange[-1:1]; set zrange[-0.000005:0.000005]" << std::endl;
  //std::cout << "set view 90,0" << std::endl;
  std::cout << "set xrange[0:160]; set yrange[-0.25:0.25]" <<std::endl;

  //std::cout << "set xrange[18:22]; set yrange[18:22]; set zrange[18:22]" << std::endl;

  //for(k=0;k<(0+dk); k+=dk){
  Ondas.Inicie(0,0,0,0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione(k);
    //if(t<21)
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
    //Ondas.Imprimir(t,25,25,25,"datos.dat");
    Ondas.Imprimase("Ondas.dat");
    //std::cout << "set pm3d explicit" << std::endl;
    std::cout << "plot 'Ondas.dat'" << std::endl;
  }
//}
  */
  return 0;
}

