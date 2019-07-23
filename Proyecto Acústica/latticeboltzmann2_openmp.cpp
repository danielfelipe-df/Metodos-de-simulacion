#include "latticeboltzmann2.h"
#include "omp.h"

void LatticeBoltzmann::Psquared(void)
{
  int ix=5;
  #pragma omp paralel for
  {
    for(int iy=0; iy<Ly;iy++){
      for(int iz=0; iz<Lz;iz++){
      Psq[ix][iy][iz]= Paverage[ix][iy][iz]*Paverage[ix][iy][iz];
      }
    }
  }
}



void LatticeBoltzmann::Pminfunction(bool init)
{
  double rho0,min; int ix=5;
  if(init)
    // for(int ix=0; ix<Lx;ix++)
#pragma omp paralel for private (rho0,min)
    {    
      for(int iy=0; iy<Ly;iy++){
	for(int iz=0; iz<Lz;iz++)
	  {
	  min=Pmin[ix][iy][iz];
	  rho0=rho(ix,iy,iz,true);
	    if(rho0<min)
	      Pmin[ix][iy][iz]=rho0;
	}
      }
    }
  else
    // for(int ix=0; ix<Lx;ix++)
#pragma omp paralel for
    {
      for(int iy=0; iy<Ly;iy++){
	for(int iz=0; iz<Lz;iz++){
	  Pmin[ix][iy][iz]=rho(ix,iy,iz,false);}
      }
    }
}


void LatticeBoltzmann::Paveragefunction(bool init)
{
  int ix=5;
  if(init)
#pragma omp paralel for
    {
    //for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++){
	for(int iz=0; iz<Lz;iz++){
	  if(ix==10)Paverage[ix][iy][iz]=Paverage[ix-1][iy][iz];
	  else Paverage[ix][iy][iz]=(Pmax[ix][iy][iz]-Pmin[ix][iy][iz])/2;
	}
      }
    }
  else
#pragma omp paralel for
    {
    //for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++){
	for(int iz=0; iz<Lz;iz++){
	  Paverage[ix][iy][iz]=0;
	}
      }
    }
}
void LatticeBoltzmann::Pmaxfunction(bool init)
{ int ix=5;
  double rho0,max;
  if(init)
#pragma omp paralel
    {
      //for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++){
	for(int iz=0; iz<Lz;iz++)
	  {
	    max=Pmax[ix][iy][iz];
	    rho0=rho(ix,iy,iz,true);
	    if(rho0>max)    Pmax[ix][iy][iz]=rho0;
	  }
      }
    }
  else
#pragma omp paralel for
    {
      //for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++){
	for(int iz=0; iz<Lz;iz++){
	  Pmax[ix][iy][iz]=rho(ix,iy,iz,false);
	}
      }
    }
}

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
  //Para cada celda
  #pragma omp paralel for
  {
  for(iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
      for(iz=0;iz<Lz;iz++){
        //Calcular las cantidades macroscópicas
        rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
	if(ix==Lx-1 || ix==0){fnew[ix][iy][iz][2]=k*f[ix][iy][iz][1]; fnew[ix][iy][iz][1]=k*f[ix][iy][iz][2];}
        else{fnew[ix][iy][iz][1]=UmUtau*f[ix][iy][iz][1]+Utau*feq(rho0,Jx0,Jy0,Jz0,1);
            fnew[ix][iy][iz][2]=UmUtau*f[ix][iy][iz][2]+Utau*feq(rho0,Jx0,Jy0,Jz0,2);}

        if(iz==Lz-1){fnew[ix][iy][iz][6]=k*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k*f[ix][iy][iz][6];}
       	else if(iz==0){fnew[ix][iy][iz][6]=k2*f[ix][iy][iz][5]; fnew[ix][iy][iz][5]=k2*f[ix][iy][iz][6];}
	else{fnew[ix][iy][iz][5]=UmUtau*f[ix][iy][iz][5]+Utau*feq(rho0,Jx0,Jy0,Jz0,5);
	  fnew[ix][iy][iz][6]=UmUtau*f[ix][iy][iz][6]+Utau*feq(rho0,Jx0,Jy0,Jz0,6);}
      
	if(iy==Ly-1 || iy==0){fnew[ix][iy][iz][3]=k*f[ix][iy][iz][4]; fnew[ix][iy][iz][4]=k*f[ix][iy][iz][3];}
        else{fnew[ix][iy][iz][3]=UmUtau*f[ix][iy][iz][3]+Utau*feq(rho0,Jx0,Jy0,Jz0,3);
	  fnew[ix][iy][iz][4]=UmUtau*f[ix][iy][iz][4]+Utau*feq(rho0,Jx0,Jy0,Jz0,4);}
      }
    }
  }
  }
    
}
void LatticeBoltzmann::Adveccione(void){
  #pragma omp paralel for
  {
    for(int ix=0;ix<Lx;ix++){
      for(int iy=0;iy<Ly;iy++){
	for(int iz=0;iz<Lz;iz++){
        for(int i=0;i<Q;i++)
          if(ix+V[0][i]<Lx && iy+V[1][i]<Ly && iz+V[2][i]<Lz && ix+V[0][i]>=0 && iy+V[1][i]>=0 && iz+V[2][i]>=0)
            f[ix+V[0][i]][iy+V[1][i]][iz+V[2][i]][i]=fnew[ix][iy][iz][i];
	}
      }
    }
  }
}

void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0, double Jz0){
  #pragma omp paralel for
  {
    for(int ix=0;ix<Lx;ix++){
      for(int iy=0;iy<Ly;iy++){
	for(int iz=0;iz<Lz;iz++){
        for(int i=0;i<Q;i++)
          f[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
	}
      }
    }
  }
}

void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy,iz,A; double lambda,omega,rho0,Jy0,Jx0,Jz0;
  lambda=7; omega=2*M_PI*C/lambda;A=1;
  ix=1; iy=3; iz=11; 
  rho0=A*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
  
  ix=36; iy=37; iz=45; 
  rho0=A*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
  
  ix=47; iy=11; iz=16; 
  rho0=A*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);

  ix=12; iy=30; iz=25; 
  rho0=A*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);

  ix=23; iy=25; iz=30; 
  rho0=A*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);

  ix=18; iy=15; iz=32; 
  rho0=A*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);

  
}

void LatticeBoltzmann::Imprimase(int t,const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double squared;
  int ix=5,iy=5,iz=5;
  squared=Psq[ix][iy][iz];
  MiArchivo<<t<<" "<<squared<<std::endl;
  MiArchivo.close();
}
void LatticeBoltzmann::Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo){
  double squaredd = Psq[ix][iy][iz];
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << squaredd << '\n';
  ofs.close();
}
bool LatticeBoltzmann::Columna(int x1, int x2, int x)
{
  if(x>=x1 && x<=x2){return true;}
  else{return false;}
}
