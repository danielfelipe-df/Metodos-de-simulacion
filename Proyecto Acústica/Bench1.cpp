#include <iostream>
#include <fstream>
#include <cmath>

const int Lx=2000;
const int Ly=1;
const int Lz=1;

const int Q=7;
const double W0=1.0/4;
//Constante de reflexión

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
  double f[Lx][Ly][Lz][Q],Pmax[Lx][Ly][Lz],Paverage[Lx][Ly][Lz],Pmin[Lx][Ly][Lz], fnew[Lx][Ly][Lz][Q]; // f[ix][iy][iz][iz][i]
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
  void Pmaxfunction(bool init);
 void Pminfunction(bool init);
  void Paveragefunction(bool init);
};
void LatticeBoltzmann::Paveragefunction(bool init)
{
  if(init)
    for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++)
	for(int iz=0; iz<Lz;iz++){
	  if(ix==10)Paverage[ix][iy][iz]=Paverage[ix-1][iy][iz];
	  else Paverage[ix][iy][iz]=(Pmax[ix][iy][iz]-Pmin[ix][iy][iz])/2;
      }
  else
   for(int ix=0; ix<Lx;ix++)
     for(int iy=0; iy<Ly;iy++)
       for(int iz=0; iz<Lz;iz++) 
	 Paverage[ix][iy][iz]=0;
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
void LatticeBoltzmann::Pminfunction(bool init)
{
  double rho0,min;
  if(init)
    for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++)
	for(int iz=0; iz<Lz;iz++)
	  {
	    min=Pmin[ix][iy][iz];
	    rho0=rho(ix,iy,iz,true);
	    if(rho0<min)
	      Pmin[ix][iy][iz]=rho0;
	  }
  else
    for(int ix=0; ix<Lx;ix++)
      for(int iy=0; iy<Ly;iy++)
	for(int iz=0; iz<Lz;iz++)
	  Pmin[ix][iy][iz]=rho(ix,iy,iz,false);
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
void LatticeBoltzmann::Colisione(double k){
  int ix,iy,iz,i; double rho0,Jx0,Jy0,Jz0;
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
      }
}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int iz=0;iz<Lz;iz++)
        for(int i=0;i<Q;i++)
          if(ix+V[0][i]>=0){
            f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][(iz+V[2][i]+Lz)%Lz][i]=fnew[ix][iy][iz][i];
          }
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
  lambda=1000; omega=2*M_PI*C/lambda; //ix=Lx/2; //iy=Ly/2; iz=Lz/2;
  rho0=sin(omega*t); Jx0=Jx(ix,iy,iz,false); Jy0=Jy(ix,iy,iz,false); Jz0=Jz(ix,iy,iz,false);
  for(i=0;i<Q;i++)
    fnew[ix][iy][iz][i]=feq(rho0,Jx0,Jy0,Jz0,i);
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  std::ofstream MiArchivo(NombreArchivo); double rho0,average;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
	average=Paverage[ix][iy][iz];
	MiArchivo<<ix<<" "<<average<<std::endl;
      }
     }
  }
    MiArchivo.close();
}
int main(void)
{
  LatticeBoltzmann Ondas;
  int t,tmax=100000,tmax2=10000;
  double k,kmax=1;
  Ondas.Inicie(0,0,0,0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione(k);
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
 
  }
  std::cout<<"set print \"wo.dat\" "<<std::endl;
  for(k=0;k<kmax;k+=0.05){
    Ondas.Pmaxfunction(false);
    Ondas.Pminfunction(false);
    Ondas.Paveragefunction(false);
    for(t=tmax;t<(tmax2+tmax);t++){
      Ondas.Colisione(k);
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
    Ondas.Pmaxfunction(true);
    Ondas.Pminfunction(true);
    }
  
    Ondas.Paveragefunction(true);
    Ondas.Imprimase("datos.dat");
    std::cout<<"stats \"datos.dat\" using 2 nooutput"<<std::endl;
    std::cout<<"y_min=STATS_min"<<std::endl;
    std::cout<<"y_max=STATS_max"<<std::endl;
    std::cout<<"SWR =y_max/y_min"<<std::endl;
    std::cout<<"SWRMU=SWR-1.0"<<std::endl;
    std::cout<<"SWRMS=SWR+1.0"<<std::endl;
    std::cout<<"cociente=SWRMU/SWRMS"<<std::endl;
    std::cout<<"cuadrado=cociente**2"<<std::endl;
    std::cout<<"T"<<k*100<<"=1.0-cuadrado"<<std::endl;
    std::cout<<"print T"<<k*100<<","<<k<<std::endl;
     }
  std::cout<<"set terminal jpeg enhanced"<<std::endl;
  std::cout<<"set output \"Grafica1.jpg\""<<std::endl;
  std::cout<<"f(x)=m*x+b "<<std::endl;
  std::cout<<"fit f(x) 'wo.dat' u 1:2 via m,b"<<std::endl;
  std::cout<<"title_f(m,b)=sprintf('f(x)=%.2fx+ %.2f',m,b)"<<std::endl;
  std::cout<<"set xlabel 'T' "<<std::endl;
  std::cout<<"set ylabel 'K'"<<std::endl;
  std::cout<<"plot f(x) t title_f(m,b), 'wo.dat' w p notitle"<<std::endl;
  return 0;
}

