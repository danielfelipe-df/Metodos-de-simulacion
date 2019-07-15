#include <iostream>
#include <fstream>
#include <cmath>

const int Lx=80;
const int Ly=20;
const int Lz=20;

const int Q=7;
const float W0=1.0/4;
//Constante de reflexión
const float k=0.1;

const float C=0.5; // C<0.707 celdas/click
const float TresC2=3*C*C;
const float AUX0=1-TresC2*(1-W0);

const float tau=0.5;
const float Utau=1.0/tau;
const float UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  float w[Q];
  int V[3][Q]; // V[0][i]=V_ix,  V[1][i]=V_iy, V[2][i]=V_iz
  float f[Lx][Ly][Lz][Q], fnew[Lx][Ly][Lz][Q], Pmax[Lx][Ly][Lz]; // f[ix][iy][iz][iz][i]
public:
  LatticeBoltzmann(void);
  float rho(int ix,int iy,int iz,bool UseNew);
  float Jx(int ix,int iy,int iz,bool UseNew);
  float Jy(int ix,int iy,int iz,bool UseNew);
  float Jz(int ix,int iy,int iz,bool UseNew);
  float feq(float rho0,float Jx0,float Jy0,float Jz0,int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(float rho0,float Jx0,float Jy0, float Jz0);
  void ImponerCampos(int t);
  void Imprimase(const char * NombreArchivo);
  void Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo);
  void Pmaxfunction(bool init);
  
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
float LatticeBoltzmann::rho(int ix,int iy, int iz,bool UseNew){
  int i; float suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]; else suma+=f[ix][iy][iz][i];
  return suma;
}
float LatticeBoltzmann::Jx(int ix,int iy, int iz,bool UseNew){
  int i; float suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]*V[0][i]; else suma+=f[ix][iy][iz][i]*V[0][i];
  return suma;
}
float LatticeBoltzmann::Jy(int ix,int iy,int iz,bool UseNew){
  int i; float suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]*V[1][i]; else suma+=f[ix][iy][iz][i]*V[1][i];
  return suma;
}
float LatticeBoltzmann::Jz(int ix,int iy,int iz,bool UseNew){
  int i; float suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][iz][i]*V[2][i]; else suma+=f[ix][iy][iz][i]*V[2][i];
  return suma;
}
float LatticeBoltzmann::feq(float rho0,float Jx0,float Jy0,float Jz0,int i){
  if(i==0)
    return rho0*(1-TresC2);
  else
    return 4*w[i]*(C*C*rho0+(V[0][i]*Jx0+V[1][i]*Jy0+V[2][i]*Jz0));
}
void LatticeBoltzmann::Colisione(void){
  int ix,iy,iz,i; float rho0,Jx0,Jy0,Jz0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Calcular las cantidades macroscópicas
        rho0=rho(ix,iy,iz,false);  Jx0=Jx(ix,iy,iz,false);  Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
	
        fnew[ix][iy][iz][0]=UmUtau*f[ix][iy][iz][0]+Utau*feq(rho0,Jx0,Jy0,Jz0,0);
	
        if(ix==Lx-1 || ix==0){
	  fnew[ix][iy][iz][2]=k*f[ix][iy][iz][1];
	  fnew[ix][iy][iz][1]=k*f[ix][iy][iz][2];}
	else{
	  fnew[ix][iy][iz][1]=UmUtau*f[ix][iy][iz][1]+Utau*feq(rho0,Jx0,Jy0,Jz0,1);
	  fnew[ix][iy][iz][2]=UmUtau*f[ix][iy][iz][2]+Utau*feq(rho0,Jx0,Jy0,Jz0,2);}

        if(iy==Ly-1 || iy==0){
	  fnew[ix][iy][iz][4]=k*f[ix][iy][iz][3];
	  fnew[ix][iy][iz][3]=k*f[ix][iy][iz][4];}
	
        else{
	  fnew[ix][iy][iz][3]=UmUtau*f[ix][iy][iz][3]+Utau*feq(rho0,Jx0,Jy0,Jz0,3);
	  fnew[ix][iy][iz][4]=UmUtau*f[ix][iy][iz][4]+Utau*feq(rho0,Jx0,Jy0,Jz0,4);}
	
        if(iz==Lz-1 || iz==0){
	  fnew[ix][iy][iz][6]=k*f[ix][iy][iz][5];
	  fnew[ix][iy][iz][5]=k*f[ix][iy][iz][6];}
        else{
	  fnew[ix][iy][iz][5]=UmUtau*f[ix][iy][iz][5]+Utau*feq(rho0,Jx0,Jy0,Jz0,5);
	  fnew[ix][iy][iz][6]=UmUtau*f[ix][iy][iz][6]+Utau*feq(rho0,Jx0,Jy0,Jz0,6);}
      }
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
        for(int i=0;i<Q;i++)
          if(ix+V[0][i]<Lx && iy+V[1][i]<Ly && iz+V[2][i]<Lz && ix+V[0][i]>=0 && iy+V[1][i]>=0 && iz+V[2][i]>=0)
            f[ix+V[0][i]][iy+V[1][i]][iz+V[2][i]][i]=fnew[ix][iy][iz][i];
}
void LatticeBoltzmann::Inicie(float rho0,float Jx0,float Jy0, float Jz0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
        for(int i=0;i<Q;i++)
          f[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy,iz; float lambda,omega,rho0,Jx0,Jy0,Jz0;
  lambda=10; omega=2*M_PI/lambda; ix=1; iy=Ly/2; iz=Lz/2;
  rho0=10*sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}

void LatticeBoltzmann::Pmaxfunction(bool init)
{
  double rho0,max;
  if(init)
    for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++)
	for(int iz=0; iz<Lz;iz++)
	  {
	    max=Pmax[ix][iy][iz];
	    rho0=rho(ix,iy,iz,true);
	    if(rho0>max)
	      Pmax[ix][iy][iz]=rho0;
	  }
  else
    for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++)
	for(int iz=0; iz<Lz;iz++)
	  Pmax[ix][iy][iz]=rho(ix,iy,iz,false);
}


void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double rho0, max;
  int iy = Ly/2;
  int iz = Lz/2;

  for(int ix=1;ix<Lx;ix++){
	max=Pmax[ix][iy][iz];
	MiArchivo<<ix<<" "<<max<<std::endl;
  }
    MiArchivo.close();
}


  void LatticeBoltzmann::Imprimir(int t, int ix, int iy, int iz, const char * NombreArchivo){
  float rho0 = rho(ix,iy,iz,false);
  std::ofstream ofs;
  ofs.open(NombreArchivo, std::ofstream::out | std::ofstream::app);
  ofs << t << '\t' << rho0 << '\n';
  ofs.close();
}


int main(void)
{
  LatticeBoltzmann Ondas;
  int t,tmax=1000;
 
  Ondas.Inicie(0,0,0,0);
  Ondas.Pmaxfunction(false);
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
    Ondas.Pmaxfunction(true);
   }
  Ondas.Imprimase("cuadrado.dat");

  return 0;
}
