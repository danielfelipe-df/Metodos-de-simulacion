#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=512;
const int Ly=64;

const int Q=9;

const double tau=1.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double visc = (tau - 1/2.0)/3.0;


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
  double sigmaxx(int ix, int iy, double rho0);
  double sigmayy(int ix, int iy, double rho0);
  double sigmaxy(int ix, int iy, double rho0);
  double Fuerza_x(double x, double y, double dAx,double dAy);
  double Fuerza_y(double x, double y, double dAx,double dAy);
  
  void Colisione(void);
  void Adveccione(void);
  void Inicie(double rho0,double Ux0,double Uy0);
  void ImponerCampos(double Vventilador);
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
double LatticeBoltzmann::sigmaxx(int ix, int iy, double rho0){
  double Ux0 = 0, Uxx = 0;
  Ux0 = Jx(ix,iy,true)/rho0;
  
  for(int i=0;i<Q;i++){  Uxx += w[i]*V[0][i]*Ux0; }
  
  return -rho0/3.0 + 6*visc*Uxx;
}

double LatticeBoltzmann::sigmayy(int ix, int iy, double rho0){
  double Uy0, Uyy = 0;
  Uy0 = Jy(ix,iy,true)/rho0;
  
  for(int i=0;i<Q;i++){  Uyy += w[i]*V[1][i]*Uy0; }
  
  return -rho0/3.0 + 6*visc*Uyy;
}


double LatticeBoltzmann::sigmaxy(int ix, int iy, double rho0){
  double Uy0, Ux0, Uxy = 0, Uyx = 0;
    Uy0 = Jy(ix,iy,true)/rho0;
    Ux0 = Jx(ix,iy,true)/rho0;

    for(int i=0;i<Q;i++){
    
     Uyx += w[i]*V[0][i]*Uy0;
     Uxy += w[i]*V[1][i]*Ux0;
   }
   
   return 3*visc*Uyx + 3*visc*Uxy;
}


double LatticeBoltzmann::Fuerza_x(double x, double y, double dAx,double dAy){

  double sigmaxxP=0, sigmaxyP=0, sigmayyP=0;
  double dFx=0;
  double rho0=0, rho1=0,rho2=0,rho3=0;

  int ix = floor(x);
  int iy = floor(y);
  
  double u =x-ix, v =y-iy ;
  
     rho0=rho(   ix,    iy,false);
     rho1=rho( ix+1,    iy,false);
     rho2=rho(   ix,  iy+1,false);
     rho3=rho( ix+1,  iy+1,false);

    
    sigmaxxP = sigmaxx(ix,iy,rho0)*(1-v)*(1-u)+sigmaxx(ix+1,iy,rho1)*u*(1-v)+sigmaxx(ix,iy+1,rho2)*v*(1-u)+sigmaxx(ix+1,iy+1,rho3)*u*v;
    
    sigmayyP = sigmayy(ix,iy,rho0)*(1-v)*(1-u)+sigmayy(ix+1,iy,rho1)*u*(1-v)+sigmayy(ix,iy+1,rho2)*v*(1-u)+sigmayy(ix+1,iy+1,rho3)*u*v;
    
    sigmaxyP = sigmaxy(ix,iy,rho0)*(1-v)*(1-u)+sigmaxy(ix+1,iy,rho1)*u*(1-v)+sigmaxy(ix,iy+1,rho2)*v*(1-u)+sigmaxy(ix+1,iy+1,rho3)*u*v;

  return sigmaxxP*dAx + sigmaxyP*dAy; 
}

double LatticeBoltzmann::Fuerza_y(double x, double y, double dAx,double dAy){

  double sigmaxxP=0, sigmaxyP=0, sigmayyP=0;
  double dFx=0;
  double rho0=0, rho1=0,rho2=0,rho3=0;

  int ix = floor(x);
  int iy = floor(y);
  
  double u =x-ix, v =y-iy ;
  
     rho0=rho(   ix,    iy,false);
     rho1=rho( ix+1,    iy,false);
     rho2=rho(   ix,  iy+1,false);
     rho3=rho( ix+1,  iy+1,false);

    
    sigmaxxP = sigmaxx(ix,iy,rho0)*(1-v)*(1-u)+sigmaxx(ix+1,iy,rho1)*u*(1-v)+sigmaxx(ix,iy+1,rho2)*v*(1-u)+sigmaxx(ix+1,iy+1,rho3)*u*v;
    
    sigmayyP = sigmayy(ix,iy,rho0)*(1-v)*(1-u)+sigmayy(ix+1,iy,rho1)*u*(1-v)+sigmayy(ix,iy+1,rho2)*v*(1-u)+sigmayy(ix+1,iy+1,rho3)*u*v;
    
    sigmaxyP = sigmaxy(ix,iy,rho0)*(1-v)*(1-u)+sigmaxy(ix+1,iy,rho1)*u*(1-v)+sigmaxy(ix,iy+1,rho2)*v*(1-u)+sigmaxy(ix+1,iy+1,rho3)*u*v;

  return  sigmaxyP*dAx + sigmayyP*dAy; 
}






void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Ux0,Uy0;
  //Para cada celda
  for(ix=0;ix<Lx;ix++)
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
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Inicie(double rho0,double Ux0,double Uy0){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[ix][iy][i]=feq(rho0,Ux0,Uy0,i);
}
void LatticeBoltzmann::ImponerCampos(double Vventilador){
  int i,ix,iy; double rho0; int ixc=Lx/2,iyc=Ly/2; int R=Ly/8, R2=R*R;
  double omega=2*M_PI/1000;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      //ventilador
      if(ix==0) 
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,Vventilador,0,i);
      //Obstáculo cilíndrico
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2)
	Ux0=-omega*(iy-iyc); Uy0=-omega*(ix-iyx);
	for(i=0;i<Q;i++)  fnew[ix][iy][i]=feq(rho0,Ux0,Uy0,i);

    }

}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo,double Vventilador){
  ofstream MiArchivo(NombreArchivo); double rho0,Ux0,Uy0;
  for(int ix=0;ix<Lx;ix+=4){
    for(int iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true);  Ux0=Jx(ix,iy,true)/rho0;  Uy0=Jy(ix,iy,true)/rho0;
      MiArchivo<<ix<<" "<<iy<<" "<<4*(Ux0-Vventilador)/Vventilador<<" "<<4*Uy0/Vventilador<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}


int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=600;
  double RHOinicial=1.0, Vventilador=0.1;
  int ixc=Lx/2,iyc=Ly/2; int R=Ly/8;
  double theta;
  double thetamax = 2*M_PI;
  double dtheta = thetamax/12;
  double x, y, dAx, dAy, Ffx, Ffy;  
  Aire.Inicie(RHOinicial,Vventilador,0);
  
  for(t=0;t<tmax;t++){
    Aire.Colisione();
    Aire.ImponerCampos(Vventilador);
    Aire.Adveccione();
    //    Aire.Imprimase("Aire.dat",Vventilador);
  
    //  cout<<"plot 'Aire.dat' w vec"<<endl;
  }
  
  for(theta=0;theta<=thetamax; theta+=dtheta ){
    
    x=R*cos(theta+dtheta/2.0) + ixc;
    y=R*sin(theta+dtheta/2.0) + iyc;
    
    dAx=sin(theta+dtheta/2.0)*M_PI*R/30;
    dAy=cos(theta+dtheta/2.0)*M_PI*R/30;

    Ffx +=Aire.Fuerza_x(x,y,dAx,dAy);
    Ffy +=Aire.Fuerza_y(x,y,dAx,dAy);

  }
  cout <<Ffx<<' ' <<Ffy<<endl;
  return 0;
}
