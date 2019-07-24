#include <omp.h>
#include "latticeboltzmann.h"

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

void LatticeBoltzmann::Colisione(void){
  int ix,iy,iz,i; double rho0,Jx0,Jy0,Jz0;

#pragma omp parallel for private(rho0, Jx0, Jy0, Jz0, ix, iy, iz, i)
  //Región sin frontera ni obstáculos dentro del auditorio
  for(iy=3*proportion;iy<Ly-(2*proportion);iy++){
    for(ix=1;ix<(Lx-1);ix++){
      for(iz=1;iz<(Lz-1);iz++){
        //Calcular las cantidades macroscópicas
        rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        for(i=0; i<Q; i++){
          fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);
        }
      }
    }
  }


#pragma omp parallel for private(rho0, Jx0, Jy0, Jz0, ix, iy, iz, i)
  //Región de las columnas
  for(ix=0; ix<Lx; ix++){
    //Caso de las columnas
    if(Columna(1,1*proportion,ix) || Columna(13*proportion,16*proportion,ix) || Columna(27*proportion,30*proportion,ix)){//Defino las regiones de la columna
      for(iz=0; iz<Lz; iz++){
        for(iy=0; iy<(3*proportion); iy++){
          rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
          fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
          fnew[ix][iy][iz][2]=k_ladrillo*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_ladrillo*f[ix][iy][iz][2];
          fnew[ix][iy][iz][4]=k_ladrillo*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_ladrillo*f[ix][iy][iz][4];
          fnew[ix][iy][iz][6]=k_ladrillo*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_ladrillo*f[ix][iy][iz][6];
        }
      }
    }
    else{
      for(iz=0; iz<Lz; iz++){
        //Caso de las ventanas
        if(Columna(3*proportion,9*proportion,iz)){
          iy=0; rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
          fnew[ix][iy][iz][3]=k_ventanas*f[ix][iy][iz][4]; fnew[ix][iy][iz][4]=k_ventanas*f[ix][iy][iz][3];
          for(i=0; i<3; i++){fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);}
          for(i=5; i<Q; i++){fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);}

          for(iy=1; iy<(3*proportion);iy++){
            rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
            for(i=0; i<Q; i++){fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);}
          }
        }
        else{
          //Otros casos
          for(iy=0; iy<(3*proportion); iy++){
            rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
            for(i=0; i<3; i++){fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);}

            if(iz==0){fnew[ix][iy][iz][6]=k_baldosa*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_baldosa*f[ix][iy][iz][6];}
            else if(iz==Lz-1){fnew[ix][iy][iz][6]=k_techo*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_techo*f[ix][iy][iz][6];}
            else{fnew[ix][iy][iz][5]=UmUtau*f[ix][iy][iz][5]+Utau*feq(rho0,Jx0,Jy0,Jz0,5);
                fnew[ix][iy][iz][6]=UmUtau*f[ix][iy][iz][6]+Utau*feq(rho0,Jx0,Jy0,Jz0,6);}

            if(iy==0){fnew[ix][iy][iz][3]=k_ladrillo*f[ix][iy][iz][4]; fnew[ix][iy][iz][4]=k_ladrillo*f[ix][iy][iz][3];}
            else{fnew[ix][iy][iz][3]=UmUtau*f[ix][iy][iz][3]+Utau*feq(rho0,Jx0,Jy0,Jz0,3);
                fnew[ix][iy][iz][4]=UmUtau*f[ix][iy][iz][4]+Utau*feq(rho0,Jx0,Jy0,Jz0,4);}
          }
        }
      }
    }
  }

#pragma omp parallel for private(rho0, Jx0, Jy0, Jz0, ix, iy, iz, i)
  //Región de las puertas
  for(ix=0; ix<Lx; ix++){
    //Defino las regiones de las columnas
    if(Columna(6*proportion,13*proportion,ix)){
      for(iz=0; iz<Lz; iz++){
        //Defino la región de la puerta en ese columna
        if(Rectangulo(7*proportion,12*proportion,ix,1,6*proportion,iz)){
          for(iy=Ly-(2*proportion); iy<Ly; iy++){
            rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
            fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
            fnew[ix][iy][iz][2]=k_puerta*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_puerta*f[ix][iy][iz][2];
            fnew[ix][iy][iz][4]=k_puerta*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_puerta*f[ix][iy][iz][4];
            fnew[ix][iy][iz][6]=k_puerta*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_puerta*f[ix][iy][iz][6];
          }
        }
        //Defino el resto de la columna
        else{
          for(iy=Ly-(2*proportion); iy<Ly; iy++){
            rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
            fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
            fnew[ix][iy][iz][2]=k_pared*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_pared*f[ix][iy][iz][2];
            fnew[ix][iy][iz][4]=k_pared*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_pared*f[ix][iy][iz][4];
            fnew[ix][iy][iz][6]=k_pared*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_pared*f[ix][iy][iz][6];
          }
        }
      }
    }
    else if(Columna(28*proportion,35*proportion,ix)){
      for(iz=0; iz<Lz; iz++){
        //Defino la región de la puerta en ese columna
        if(Rectangulo(29*proportion,34*proportion,ix,1,6*proportion,iz)){
          for(iy=Ly-(2*proportion); iy<Ly; iy++){
            rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
            fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
            fnew[ix][iy][iz][2]=k_puerta*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_puerta*f[ix][iy][iz][2];
            fnew[ix][iy][iz][4]=k_puerta*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_puerta*f[ix][iy][iz][4];
            fnew[ix][iy][iz][6]=k_puerta*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_puerta*f[ix][iy][iz][6];
          }
        }
        //Defino el resto de la columna
        else{
          for(iy=Ly-(2*proportion); iy<Ly; iy++){
            rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
            fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
            fnew[ix][iy][iz][2]=k_pared*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_pared*f[ix][iy][iz][2];
            fnew[ix][iy][iz][4]=k_pared*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_pared*f[ix][iy][iz][4];
            fnew[ix][iy][iz][6]=k_pared*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_pared*f[ix][iy][iz][6];
          }
        }
      }
    }
    else{//Defino las otras regiones
      for(iz=0; iz<Lz; iz++){
        for(iy=Ly-(2*proportion); iy<Ly; iy++){
          rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
          for(i=0; i<3; i++){fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);}

          if(iz==0){fnew[ix][iy][iz][6]=k_baldosa*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_baldosa*f[ix][iy][iz][6];}
          else if(iz==Lz-1){fnew[ix][iy][iz][6]=k_techo*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_techo*f[ix][iy][iz][6];}
          else{fnew[ix][iy][iz][5]=UmUtau*f[ix][iy][iz][5]+Utau*feq(rho0,Jx0,Jy0,Jz0,5);
                fnew[ix][iy][iz][6]=UmUtau*f[ix][iy][iz][6]+Utau*feq(rho0,Jx0,Jy0,Jz0,6);}

          if(iy==Ly-1){fnew[ix][iy][iz][3]=k_zigzag*f[ix][iy][iz][4]; fnew[ix][iy][iz][4]=k_zigzag*f[ix][iy][iz][3];}
          else{fnew[ix][iy][iz][3]=UmUtau*f[ix][iy][iz][3]+Utau*feq(rho0,Jx0,Jy0,Jz0,3);
                fnew[ix][iy][iz][4]=UmUtau*f[ix][iy][iz][4]+Utau*feq(rho0,Jx0,Jy0,Jz0,4);}
        }
      }
    }
  }

  //Región del tablero
  ix=0;
#pragma omp parallel for private(rho0, Jx0, Jy0, Jz0, iy, iz, i)
  for(iz=0; iz<Lz; iz++){
    if(Columna(3*proportion,7*proportion,iz)){
      for(iy=3*proportion;iy<Ly-(2*proportion);iy++){
        rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
        fnew[ix][iy][iz][2]=k_tablero*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_tablero*f[ix][iy][iz][2];
        fnew[ix][iy][iz][4]=k_tablero*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_tablero*f[ix][iy][iz][4];
        fnew[ix][iy][iz][6]=k_tablero*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_tablero*f[ix][iy][iz][6];
      }
    }
    else{
      for(iy=3*proportion;iy<Ly-(2*proportion);iy++){
        rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
        fnew[ix][iy][iz][2]=k_pared*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_pared*f[ix][iy][iz][2];
        fnew[ix][iy][iz][4]=k_pared*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_pared*f[ix][iy][iz][4];
        fnew[ix][iy][iz][6]=k_pared*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_pared*f[ix][iy][iz][6];
      }
    }
  }

  //Región opuesta del tablero
  ix=Lx-1;
#pragma omp parallel for private(rho0, Jx0, Jy0, Jz0, iy, iz, i)
  for(iy=3*proportion;iy<Ly-(2*proportion);iy++){
    for(iz=0; iz<Lz; iz++){
      rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
      if(Rectangulo(9*proportion,12*proportion,iy,3*proportion,6*proportion,iz)){
        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
        fnew[ix][iy][iz][2]=k_ventanita*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_ventanita*f[ix][iy][iz][2];
        fnew[ix][iy][iz][4]=k_ventanita*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_ventanita*f[ix][iy][iz][4];
        fnew[ix][iy][iz][6]=k_ventanita*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_ventanita*f[ix][iy][iz][6];
      }
      else{
        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
        fnew[ix][iy][iz][2]=k_pared*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k_pared*f[ix][iy][iz][2];
        fnew[ix][iy][iz][4]=k_pared*f[ix][iy][iz][3]; fnew[ix][iy][iz][3]=k_pared*f[ix][iy][iz][4];
        fnew[ix][iy][iz][6]=k_pared*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_pared*f[ix][iy][iz][6];
      }
    }
  }

  //Para el piso
  iz=0;
#pragma omp parallel for private(rho0, Jx0, Jy0, Jz0, ix, iy, i)
  for(iy=3*proportion;iy<Ly-(2*proportion);iy++){
    for(ix=1; ix<(Lx-1); ix++){
      rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
      for(i=0; i<5; i++){fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);}
      fnew[ix][iy][iz][6]=k_baldosa*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_baldosa*f[ix][iy][iz][6];
    }
  }

  //Para el techo
  iz=Lz-1;
#pragma omp parallel for private(rho0, Jx0, Jy0, Jz0, ix, iy, i)
  for(iy=3*proportion;iy<Ly-(2*proportion);iy++){
    for(ix=1; ix<(Lx-1); ix++){
      rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
      for(i=0; i<5; i++){fnew[ix][iy][iz][i]=UmUtau*f[ix][iy][iz][i]+Utau*feq(rho0,Jx0,Jy0,Jz0,i);}
      fnew[ix][iy][iz][6]=k_techo*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k_techo*f[ix][iy][iz][6];
    }
  }
}

void LatticeBoltzmann::Adveccione(void){
  int ix,iy,iz,i;
#pragma omp parallel for private(ix, iy, iz, i)
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
        for(i=0;i<Q;i++)
          if(ix+V[0][i]<Lx && iy+V[1][i]<Ly && iz+V[2][i]<Lz && ix+V[0][i]>=0 && iy+V[1][i]>=0 && iz+V[2][i]>=0)
            f[ix+V[0][i]][iy+V[1][i]][iz+V[2][i]][i]=fnew[ix][iy][iz][i];
}

void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0, double Jz0){
  int ix,iy,iz,i;
#pragma omp parallel for private(ix, iy, iz, i)
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++)
        for(i=0;i<Q;i++)
          f[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}

void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy,iz,A; double lambda,omega,rho0,Jx0,Jy0,Jz0;
  lambda=7; omega=2*M_PI*C/lambda; A=2;
  //Posiciones fuente normal
  ix=1*proportion; iy=9*proportion; iz=3*proportion;
  //Posicionnes fuente rosa
  /*ix=1*proportion; iy=9*proportion; iz=3*proportion;*/
  rho0=A*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}

void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0,Jz0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      //for(int iz=0;iz<Lz;iz++){
        int iz=20;
        rho0=rho(ix,iy,iz,true);  //Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        MiArchivo<<ix<<" "<<iy<<" "<<rho0<<'\n';
      //}
      //MiArchivo<<'\n';
    }
    MiArchivo<<'\n';
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo){
  double rho0 = rho(ix,iy,iz,true);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << rho0 << '\n';
  ofs.close();
}

bool LatticeBoltzmann::Columna(int x1, int x2, int x)
{
  if(x>=x1 && x<=x2){return true;}
  else{return false;}
}

bool LatticeBoltzmann::Rectangulo(int x1, int x2, int x, int y1, int y2, int y)
{
  if(x>=x1 && x<=x2){
    if(y>=y1 && y<=y2){return true;}
    else{return false;}
  }
  else{return false;}
}
