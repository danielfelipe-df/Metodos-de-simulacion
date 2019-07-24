#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/gl.h>
#include <GL/glut.h>
#include "latticeboltzmann.h"

double rho1[Lx][Ly][Lz];
LatticeBoltzmann Ondas;
int t,t_otro;
/*
void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT);

  double rho0;  int iz=Lz/2;
  Ondas.Colisione();
  Ondas.ImponerCampos(t);
  Ondas.Adveccione();
  glPointSize(3.0);
  glBegin(GL_POINTS);
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
      rho0 = rho1[ix][iy][iz];
      /*
      if(rho0>0){glColor3f(1.0-rho0*40.0,1.0-rho0*40.0,1.0);}
      else if(rho0<0){glColor3f(1.0,1.0+rho0*40.0,1.0+rho0*40.0);}
      else{glColor3f(1.0,1.0,1.0);}
      glVertex3d(ix*0.0025-0.15,iy*0.0075-0.25,rho0*2.0);
      *//*
      if(rho0>0){glColor3f(1.0-rho0*40.0,1.0-rho0*40.0,1.0);}
      else if(rho0<0){glColor3f(1.0,1.0+rho0*40.0,1.0+rho0*40.0);}
      else{glColor3f(1.0,1.0,1.0);}
      glVertex3d(ix*0.0075+0.025,iy*0.0075+0.25,0.0);
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
  /*initialize viewing values *//*
  //glMatrixMode(GL_PROJECTION);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
  //glOrtho(0.0, 1.0, -2.0, 2.0, -10.0, 1.0);
  //gluLookAt(0.25, -0.5, 2.0, 10.5, 10.5, 70.0, 1.0, 1.0, -1.0);
}

void AmplitudDisplay(void)
{
  for(int ix=0; ix<Lx; ix++){
    for(int iy=0; iy<Ly; iy++){
      for(int iz=0; iz<Lz; iz++){
        rho1[ix][iy][iz] = Ondas.rho(ix,iy,iz,true);
      }
    }
  }
  t++;
  std::cout << t << std::endl;
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
*/

int main(int argc, char** argv)
{
  //OpenGL
  /*
  Ondas.Inicie(0,0,0,0);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(400, 400);
  glutInitWindowPosition(500, 100);
  glutCreateWindow("Animation");
  init();
  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  glutMainLoop();
  */
  //Gnuplot

  int t,tmax=1000;

  // Estos comandos se descomentan si se quiere guardar el gif
  /*
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set output 'pelicula0.gif'" << std::endl;

  //Estos comandos se descomentan para hacer el gif

  std::cout << "set pm3d" << std::endl;
  std::cout << "set palette defined (-1 \"red\", 0 \"white\", 1 \"blue\")" << std::endl;
  std::cout << "set cbrange[-1:1]" << std::endl;
  std::cout << "set xrange[-1:41]; set yrange[-1:41]; set zrange[-1:5]" << std::endl;
  */

  Ondas.Inicie(0,0,0,0);
  std::ofstream fout("Datos.csv");
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
    fout << t << '\t' ;
    fout << Ondas.GetRho(13*proportion,(3+5)*proportion,2*proportion,true) << '\t'; //Primera posición
    fout << Ondas.GetRho(13*proportion,(21-4)*proportion,2*proportion,true) << '\t'; //Segunda posición
    fout << Ondas.GetRho(28*proportion,(3+5)*proportion,2*proportion,true) << '\t'; //Tercera posición
    fout << Ondas.GetRho(28*proportion,(21-4)*proportion,2*proportion,true) << '\t'; //Cuarta posición
    fout << Ondas.GetRho(39*proportion,(3+5)*proportion,2*proportion,true) << '\t'; //Quinta posición
    fout << Ondas.GetRho(39*proportion,(21-4)*proportion,2*proportion,true) << '\t'; //Sexta posición
    fout << Ondas.GetRho(39*proportion,10*proportion,2*proportion,true) << '\t'; //Séptima posición
    fout << Ondas.GetRho(39*proportion,2*proportion,2*proportion,true) << '\t'; //Octava posición
    fout << Ondas.GetRho(39*proportion,(21-1)*proportion,2*proportion,true) << '\t'; //Novena posición
    fout << Ondas.GetRho(20*proportion,10*proportion,2*proportion,true) << '\t'; //Décima posición
    fout << Ondas.GetRho(2*proportion,2*proportion,2*proportion,true) << '\t'; //Onceava posición
    fout << Ondas.GetRho(2*proportion,(21-2)*proportion,2*proportion,true) << '\t'; //Novena posición
    fout << '\n';
  }
  fout.close();
  //Ondas.Imprimase("Ondas.dat");

  return 0;
}

