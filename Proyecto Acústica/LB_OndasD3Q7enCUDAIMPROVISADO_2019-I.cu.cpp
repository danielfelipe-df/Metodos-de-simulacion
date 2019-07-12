// Programa que hace c=a+b para a,b,c vectores en CUDA
#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 128
#define Ly 128
#define Lz 128

#define Nx 8
#define Ny 8
#define Nz 8


const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Ny;
const int Mz=(Lz+Nz-1)/Nz;


#define Q 7
#define C (0.5) // C<0.707 celdas/click
#define TresC2 (0.75)
#define AUX0 (0.4375)

#define tau (0.5)
#define Utau (2.0)
#define UmUtau (-1.0)

//--------------------KERNELS----------------
__constant__ float d_w[Q];
__constant__ int d_Vx[Q];
__constant__ int d_Vy[Q];
__constant__ int d_Vz[Q];


//Calcular Campos Macroscopicos y la funcion de equilibrio en Device
__device__ float d_rho(float f0,float f1,float f2,float f3,float f4, float f5, float f6){
  return f0+f1+f2+f3+f4+f5+f6;
}
__device__ float d_Jx(float f0,float f1,float f2,float f3,float f4, float f5, float f6){
  return f0*d_Vx[0]+f1*d_Vx[1]+f2*d_Vx[2]+f3*d_Vx[3]+f4*d_Vx[4]+f5*d_Vx[5]+f6*d_Vx[6];
}
__device__ float d_Jy(float f0,float f1,float f2,float f3,float f4, float f5, float f6){
  return f0*d_Vy[0]+f1*d_Vy[1]+f2*d_Vy[2]+f3*d_Vy[3]+f4*d_Vy[4]+f5*d_Vy[5]+f6*d_Vy[6];
}
  __device__ float d_Jz(float f0,float f1,float f2,float f3,float f4, float f5, float f6){
  return f0*d_Vz[0]+f1*d_Vz[1]+f2*d_Vz[2]+f3*d_Vz[3]+f4*d_Vz[4]+f5*d_Vz[5]+f6*d_Vz[6];
}
__device__ float d_feq0(float rho0,float Jx0,float Jy0, float Jz0){
  return rho0*AUX0;
}
__device__ float d_feqi(float rho0,float Jx0,float Jy0,float Jz0, int i){
  return 4*d_w[i]*(C*C*rho0+(d_Vx[i]*Jx0 + d_Vy[i]*Jy0 + d_Vz[i]*Jz0));
}
//Kernels de Evolucion
__global__ void d_Colisione(float *d_f0,size_t pitchf0,
			    float *d_f1,size_t pitchf1,
			    float *d_f2,size_t pitchf2,
			    float *d_f3,size_t pitchf3,
			    float *d_f4,size_t pitchf4,
			    float *d_f5,size_t pitchf5,
			    float *d_f6,size_t pitchf6,
			    float *d_f0new,size_t pitchf0new,
			    float *d_f1new,size_t pitchf1new,
			    float *d_f2new,size_t pitchf2new,
			    float *d_f3new,size_t pitchf3new,
			    float *d_f4new,size_t pitchf4new,
			    float *d_f5new,size_t pitchf5new,
			    float *d_f6new,size_t pitchf6new){
  //Definir variables
  int ix,iy, iz; float *f0,*f1,*f2,*f3,*f4,*f5,*f6,*f0new,*f1new,*f2new,*f3new,*f4new,*f5new,*f6new;

  //Determinar la posicion que atiende el hilo (thread) = ¿Quién soy yo?
  ix=blockIdx.x*blockDim.x+threadIdx.x;   iy=blockIdx.y*blockDim.y+threadIdx.y;     iz=blockIdx.z*blockDim.z+threadIdx.z;

  //Identificar las posiciones de lectura y de escritura
  f0=d_f0+(ix*pitchf0)/sizeof(float)+iy+iz;  f0new=d_f0new+(ix*pitchf0new)/sizeof(float)+iy+iz;
  f1=d_f1+(ix*pitchf1)/sizeof(float)+iy+iz;  f1new=d_f1new+(ix*pitchf1new)/sizeof(float)+iy+iz;
  f2=d_f2+(ix*pitchf2)/sizeof(float)+iy+iz;  f2new=d_f2new+(ix*pitchf2new)/sizeof(float)+iy+iz;
  f3=d_f3+(ix*pitchf3)/sizeof(float)+iy+iz;  f3new=d_f3new+(ix*pitchf3new)/sizeof(float)+iy+iz;
  f4=d_f4+(ix*pitchf4)/sizeof(float)+iy+iz;  f4new=d_f4new+(ix*pitchf4new)/sizeof(float)+iy+iz;  
  f5=d_f5+(ix*pitchf5)/sizeof(float)+iy+iz;  f5new=d_f5new+(ix*pitchf5new)/sizeof(float)+iy+iz;
  f6=d_f6+(ix*pitchf6)/sizeof(float)+iy+iz;  f6new=d_f6new+(ix*pitchf6new)/sizeof(float)+iy+iz; //Ni idea de como va esta wea
  
  //Procesar los datos

  //Calcular las cantidades macroscópicas
  float rho0,Jx0,Jy0, Jz0;
  rho0=d_rho(*f0,*f1,*f2,*f3,*f4,*f5,*f6); Jx0=d_Jx(*f0,*f1,*f2,*f3,*f4,*f5,*f6); Jy0=d_Jy(*f0,*f1,*f2,*f3,*f4,*f5,*f6);

  //Evolucionar
  (*f0new)=UmUtau*(*f0)+Utau*d_feq0(rho0,Jx0,Jy0,Jz0);
  (*f1new)=UmUtau*(*f1)+Utau*d_feqi(rho0,Jx0,Jy0,Jz0,1);
  (*f2new)=UmUtau*(*f2)+Utau*d_feqi(rho0,Jx0,Jy0,Jz0,2);
  (*f3new)=UmUtau*(*f3)+Utau*d_feqi(rho0,Jx0,Jy0,Jz0,3);
  (*f4new)=UmUtau*(*f4)+Utau*d_feqi(rho0,Jx0,Jy0,Jz0,4);
  (*f5new)=UmUtau*(*f5)+Utau*d_feqi(rho0,Jx0,Jy0,Jz0,5);
  (*f6new)=UmUtau*(*f6)+Utau*d_feqi(rho0,Jx0,Jy0,Jz0,6);
}
__global__ void d_ImponerCampos(float *d_f0,size_t pitchf0,
				float *d_f1,size_t pitchf1,
				float *d_f2,size_t pitchf2,
				float *d_f3,size_t pitchf3,
				float *d_f4,size_t pitchf4,
				float *d_f5,size_t pitchf5,
				float *d_f6,size_t pitchf6,
				float *d_f0new,size_t pitchf0new,
				float *d_f1new,size_t pitchf1new,
				float *d_f2new,size_t pitchf2new,
				float *d_f3new,size_t pitchf3new,
				float *d_f4new,size_t pitchf4new,
				float *d_f5new,size_t pitchf5new,
				float *d_f6new,size_t pitchf6new,
				float rhoFuente){
  //Definir variables
  int ix,iy,iz; float *f0,*f1,*f2,*f3,*f4,*f5,*f6,*f0new,*f1new,*f2new,*f3new,*f4new,*f5new,*f6new;
  float rho0,Jx0,Jy0,Jz0;

  //Determinar la posicion que atiende el hilo (thread) = ¿Quién soy yo?
  ix=blockIdx.x*blockDim.x+threadIdx.x;   iy=blockIdx.y*blockDim.y+threadIdx.y;    iz=blockIdx.z*blockDim.z+threadIdx.z;
  if(ix==Lx/2 && iy==Ly/2 && iz==Lz/2){  //Esto tambien hay que modificarlo, i guess

  //Identificar las posiciones de lectura y de escritura
    f0=d_f0+(ix*pitchf0)/sizeof(float)+iy+iz;  f0new=d_f0new+(ix*pitchf0new)/sizeof(float)+iy+iz;
    f1=d_f1+(ix*pitchf1)/sizeof(float)+iy+iz;  f1new=d_f1new+(ix*pitchf1new)/sizeof(float)+iy+iz;
    f2=d_f2+(ix*pitchf2)/sizeof(float)+iy+iz;  f2new=d_f2new+(ix*pitchf2new)/sizeof(float)+iy+iz;
    f3=d_f3+(ix*pitchf3)/sizeof(float)+iy+iz;  f3new=d_f3new+(ix*pitchf3new)/sizeof(float)+iy+iz;
    f4=d_f4+(ix*pitchf4)/sizeof(float)+iy+iz;  f4new=d_f4new+(ix*pitchf4new)/sizeof(float)+iy+iz; //Ni idea again

    //Procesar los datos

    //Calcular las cantidades macroscópicas
    rho0=rhoFuente; Jx0=d_Jx(*f0,*f1,*f2,*f3,*f4,*f5,*f6); Jy0=d_Jy(*f0,*f1,*f2,*f3,*f4,*f5,*f6); Jz0=d_Jz(*f0,*f1,*f2,*f3,*f4,*f5,*f6);

    //Imponer las condiciones de frontera
    (*f0new)=d_feq0(rho0,Jx0,Jy0,Jz0);
    (*f1new)=d_feqi(rho0,Jx0,Jy0,Jz0,1);
    (*f2new)=d_feqi(rho0,Jx0,Jy0,Jz0,2);
    (*f3new)=d_feqi(rho0,Jx0,Jy0,Jz0,3);
    (*f4new)=d_feqi(rho0,Jx0,Jy0,Jz0,4);
    (*f5new)=d_feqi(rho0,Jx0,Jy0,Jz0,5);
    (*f6new)=d_feqi(rho0,Jx0,Jy0,Jz0,6);
  }
}
__global__ void d_Adveccione(float *d_f0,size_t pitchf0,
			     float *d_f1,size_t pitchf1,
			     float *d_f2,size_t pitchf2,
			     float *d_f3,size_t pitchf3,
			     float *d_f4,size_t pitchf4,
			     float *d_f5,size_t pitchf5,
			     float *d_f6,size_t pitchf6,
			     float *d_f0new,size_t pitchf0new,
			     float *d_f1new,size_t pitchf1new,
			     float *d_f2new,size_t pitchf2new,
			     float *d_f3new,size_t pitchf3new,
			     float *d_f4new,size_t pitchf4new,
			     float *d_f5new,size_t pitchf5new,
			     float *d_f6new,size_t pitchf6new){
  //Definir variables
  int ix,iy,iz; float *f0,*f1,*f2,*f3,*f4,*f5,*f6; float *f0new,*f1new,*f2new,*f3new,*f4new,*f5new,*f6new;

  //Determinar la posicion que atiende el hilo (thread) = ¿Quién soy yo?
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;  iz=blockIdx.z*blockDim.z+threadIdx.z;

  //LEE EN fnew Y ESCRIBE EN f
  //Leer los datos
  f0new=d_f0new+(ix*pitchf0new)/sizeof(float)+iy+iz;
  f1new=d_f1new+(ix*pitchf0new)/sizeof(float)+iy+iz;
  f2new=d_f2new+(ix*pitchf0new)/sizeof(float)+iy+iz;
  f3new=d_f3new+(ix*pitchf0new)/sizeof(float)+iy+iz;
  f4new=d_f4new+(ix*pitchf0new)/sizeof(float)+iy+iz;  //Revisar
  f5new=d_f5new+(ix*pitchf0new)/sizeof(float)+iy+iz;  //Revisar
  f6new=d_f6new+(ix*pitchf0new)/sizeof(float)+iy+iz;  //Revisar
 
  //Determinar las posiciones de escritura
  f0=d_f0+(((ix+d_Vx[0]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[0]+Ly)%Ly)+((iz+d_Vz[0]+Lz)%Lz); 
  f1=d_f1+(((ix+d_Vx[1]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[1]+Ly)%Ly)+((iz+d_Vz[1]+Lz)%Lz); 
  f2=d_f2+(((ix+d_Vx[2]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[2]+Ly)%Ly)+((iz+d_Vz[2]+Lz)%Lz); 
  f3=d_f3+(((ix+d_Vx[3]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[3]+Ly)%Ly)+((iz+d_Vz[3]+Lz)%Lz); 
  f4=d_f4+(((ix+d_Vx[4]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[4]+Ly)%Ly)+((iz+d_Vz[4]+Lz)%Lz);   //Revisar
  f5=d_f5+(((ix+d_Vx[5]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[5]+Ly)%Ly)+((iz+d_Vz[5]+Lz)%Lz);   //Revisar
  f6=d_f6+(((ix+d_Vx[6]+Lx)%Lx)*pitchf0)/sizeof(float)+((iy+d_Vy[6]+Ly)%Ly)+((iz+d_Vz[6]+Lz)%Lz);   //Revisar
 
  //Escribirlos en las nuevas posiciones
  (*f0)=(*f0new);
  (*f1)=(*f1new);
  (*f2)=(*f2new);
  (*f3)=(*f3new);
  (*f4)=(*f4new);
  (*f5)=(*f5new);
  (*f6)=(*f6new);
}


//------------------- CLASES ----------------
class LatticeBoltzmann{
private:
  float h_w[Q]; int h_Vx[Q], h_Vy[Q], h_Vz[Q];

  float h_f0[Lx][Ly][Lz]; float*d_f0; size_t pitchf0; 
  float h_f1[Lx][Ly][Lz]; float*d_f1; size_t pitchf1; 
  float h_f2[Lx][Ly][Lz]; float*d_f2; size_t pitchf2; 
  float h_f3[Lx][Ly][Lz]; float*d_f3; size_t pitchf3; 
  float h_f4[Lx][Ly][Lz]; float*d_f4; size_t pitchf4; 
  float h_f5[Lx][Ly][Lz]; float*d_f5; size_t pitchf5; 
  float h_f6[Lx][Ly][Lz]; float*d_f6; size_t pitchf6; 

  
  float h_f0new[Lx][Ly][Lz]; float*d_f0new; size_t pitchf0new; 
  float h_f1new[Lx][Ly][Lz]; float*d_f1new; size_t pitchf1new; 
  float h_f2new[Lx][Ly][Lz]; float*d_f2new; size_t pitchf2new; 
  float h_f3new[Lx][Ly][Lz]; float*d_f3new; size_t pitchf3new; 
  float h_f4new[Lx][Ly][Lz]; float*d_f4new; size_t pitchf4new; 
  float h_f5new[Lx][Ly][Lz]; float*d_f5new; size_t pitchf5new; 
  float h_f6new[Lx][Ly][Lz]; float*d_f6new; size_t pitchf6new;
  
public:
  //Constructor y Destructor
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  //Cantidades macroscópicas
  double h_rho(int ix,int iy,int iz,bool UseNew);
  double h_Jx(int ix,int iy,int iz,bool UseNew);
  double h_Jy(int ix,int iy,int iz,bool UseNew);
  double h_Jz(int ix,int iy,int iz,bool UseNew);
  double h_feq(double rho0,double Jx0,double Jy0, double Jz0, int i);
  //Funciones de evolución
  void Inicie(double rho0,double Jx0,double Jy0, double Jz0);
  void Colisione(void);
  void ImponerCampos(int t);
  void Adveccione(void);
  void Imprimase(const char * NombreArchivo);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Cargar los pesos
  h_w[0]=1.0/4; h_w[1]=h_w[2]=h_w[3]=h_w[4]=h_w[5]=h_w[6]=1.0/8; //cargarlos en el Host
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice); //enviarlos al Device

  //Cargar los vectores
  h_Vx[0]=0;  //cargarlos en el Host
  h_Vy[0]=0;
  h_Vz[0]=0;

  h_Vx[1]=1;  h_Vx[2]=-1;  h_Vx[3]=0; h_Vx[4]=0;  h_Vx[5]=0;  h_Vx[6]=0;
  h_Vy[1]=0;  h_Vy[2]=0;  h_Vy[3]=1;  h_Vy[4]=-1;  h_Vy[5]=0;  h_Vy[6]=0;
  h_Vz[1]=0;  h_Vz[2]=0;  h_Vz[3]=1;  h_Vz[4]=-1;  h_Vz[5]=1;  h_Vz[6]=-1;

  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice); //enviarlos al Device
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice); //enviarlos al Device
  cudaMemcpyToSymbol(d_Vz,h_Vz,Q*sizeof(int),0,cudaMemcpyHostToDevice); //enviarlos al Device

  //Construir las matrices f en el Device
  cudaMallocPitch((void**) &d_f0,&pitchf0,Ly*sizeof(float), Lz*sizeof(float), Lx);
  cudaMallocPitch((void**) &d_f1,&pitchf1,Ly*sizeof(float), Lz*sizeof(float), Lx);
  cudaMallocPitch((void**) &d_f2,&pitchf2,Ly*sizeof(float), Lz*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f3,&pitchf3,Ly*sizeof(float), Lz*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f4,&pitchf4,Ly*sizeof(float), Lz*sizeof(float),Lx); //Revisar
  cudaMallocPitch((void**) &d_f5,&pitchf5,Ly*sizeof(float), Lz*sizeof(float),Lx); //Revisar
  cudaMallocPitch((void**) &d_f6,&pitchf6,Ly*sizeof(float), Lz*sizeof(float),Lx); //Revisar

  cudaMallocPitch((void**) &d_f0new,&pitchf0new,Ly*sizeof(float),Lz*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f1new,&pitchf1new,Ly*sizeof(float),Lz*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f2new,&pitchf2new,Ly*sizeof(float),Lz*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f3new,&pitchf3new,Ly*sizeof(float),Lz*sizeof(float),Lx);
  cudaMallocPitch((void**) &d_f4new,&pitchf4new,Ly*sizeof(float),Lz*sizeof(float),Lx);//Revisar
  cudaMallocPitch((void**) &d_f5new,&pitchf5new,Ly*sizeof(float),Lz*sizeof(float),Lx);//Revisar
  cudaMallocPitch((void**) &d_f6new,&pitchf6new,Ly*sizeof(float),Lz*sizeof(float),Lx);//Revisar
}

LatticeBoltzmann::~LatticeBoltzmann(void){
  cudaFree(d_f0);  cudaFree(d_f1);  cudaFree(d_f2);  cudaFree(d_f3);  cudaFree(d_f4);cudaFree(d_f5);cudaFree(d_f6);
  cudaFree(d_f0new);cudaFree(d_f1new);cudaFree(d_f2new);cudaFree(d_f3new);cudaFree(d_f4new);cudaFree(d_f5new);cudaFree(d_f6new);
}

//Calcular Campos Macroscopicos y la funcion de equilibrio en Host
double LatticeBoltzmann::h_rho(int ix,int iy,int iz,bool UseNew){
  if(UseNew) 
    return h_f0new[ix][iy][iz]+h_f1new[ix][iy][iz]+h_f2new[ix][iy][iz]+h_f3new[ix][iy][iz]+h_f4new[ix][iy][iz]+h_f5new[ix][iy][iz]+h_f6new[ix][iy][iz];
  else
    return h_f0[ix][iy][iz]+ h_f1[ix][iy][iz]+h_f2[ix][iy][iz]+h_f3[ix][iy][iz]+h_f4[ix][iy][iz]+h_f5[ix][iy][iz]+h_f6[ix][iy][iz];
}
double LatticeBoltzmann::h_Jx(int ix,int iy, int iz,bool UseNew){
  if(UseNew) 
    return h_Vx[0]*h_f0new[ix][iy][iz]+h_Vx[1]*h_f1new[ix][iy][iz]+h_Vx[2]*h_f2new[ix][iy][iz]
      +h_Vx[3]*h_f3new[ix][iy][iz]+h_Vx[4]*h_f4new[ix][iy][iz]+h_Vx[5]*h_f5new[ix][iy][iz]+ h_Vx[6]*h_f6new[ix][iy][iz];
  else
    return h_Vx[0]*h_f0[ix][iy][iz]+h_Vx[1]*h_f1[ix][iy][iz]+h_Vx[2]*h_f2[ix][iy][iz]
      +h_Vx[3]*h_f3[ix][iy][iz]+h_Vx[4]*h_f4[ix][iy][iz]+h_Vx[5]*h_f5[ix][iy][iz]+h_Vx[6]*h_f6[ix][iy][iz];
}



double LatticeBoltzmann::h_Jy(int ix,int iy, int iz,bool UseNew){
  if(UseNew) 
    return h_Vy[0]*h_f0new[ix][iy][iz]+h_Vy[1]*h_f1new[ix][iy][iz]+h_Vy[2]*h_f2new[ix][iy][iz]
      +h_Vy[3]*h_f3new[ix][iy][iz]+h_Vy[4]*h_f4new[ix][iy][iz]+h_Vy[5]*h_f5new[ix][iy][iz]+ h_Vy[6]*h_f6new[ix][iy][iz];
  else
    return h_Vy[0]*h_f0[ix][iy][iz]+h_Vy[1]*h_f1[ix][iy][iz]+h_Vy[2]*h_f2[ix][iy][iz]
      +h_Vy[3]*h_f3[ix][iy][iz]+h_Vy[4]*h_f4[ix][iy][iz]+h_Vy[5]*h_f5[ix][iy][iz]+h_Vy[6]*h_f6[ix][iy][iz];
}



double LatticeBoltzmann::h_Jz(int ix,int iy, int iz,bool UseNew){
  if(UseNew) 
    return h_Vz[0]*h_f0new[ix][iy][iz]+h_Vz[1]*h_f1new[ix][iy][iz]+h_Vz[2]*h_f2new[ix][iy][iz]
      +h_Vz[3]*h_f3new[ix][iy][iz]+h_Vz[4]*h_f4new[ix][iy][iz]+h_Vz[5]*h_f5new[ix][iy][iz]+ h_Vz[6]*h_f6new[ix][iy][iz];
  else
    return h_Vz[0]*h_f0[ix][iy][iz]+h_Vz[1]*h_f1[ix][iy][iz]+h_Vz[2]*h_f2[ix][iy][iz]
      +h_Vz[3]*h_f3[ix][iy][iz]+h_Vz[4]*h_f4[ix][iy][iz]+h_Vz[5]*h_f5[ix][iy][iz]+h_Vz[6]*h_f6[ix][iy][iz];
}


double LatticeBoltzmann::h_feq(double rho0,double Jx0,double Jy0,double Jz0,int i){
  if(i==0)
    return rho0*AUX0;
  else
    return 4*h_w[i]*(C*C*rho0+(h_Vx[i]*Jx0 + h_Vy[i]*Jy0 + h_Vz[i]*Jz0));
}
//Funciones para la Evolucion
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0, double Jz0){
  //Cargarlas en el Host
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++){
      for(int iz=0;iz<Lz;iz++){
      
	h_f0[ix][iy][iz]=h_feq(rho0,Jx0,Jy0,Jz0,0);
	h_f1[ix][iy][iz]=h_feq(rho0,Jx0,Jy0,Jz0,1);
	h_f2[ix][iy][iz]=h_feq(rho0,Jx0,Jy0,Jz0,2);
	h_f3[ix][iy][iz]=h_feq(rho0,Jx0,Jy0,Jz0,3);
	h_f4[ix][iy][iz]=h_feq(rho0,Jx0,Jy0,Jz0,4);
	h_f5[ix][iy][iz]=h_feq(rho0,Jx0,Jy0,Jz0,5);
	h_f6[ix][iy][iz]=h_feq(rho0,Jx0,Jy0,Jz0,6);
      }
    }

  
  //Enviarlas al Device
  cudaMemcpy3D(d_f0,pitchf0,h_f0,Ly*sizeof(float),Lz*sizeof(float),Ly*sizeof(float),Lz*sizeof(float), Lx,cudaMemcpyHostToDevice);
  cudaMemcpy3D(d_f1,pitchf1,h_f1,Ly*sizeof(float),Lz*sizeof(float),Ly*sizeof(float),Lz*sizeof(float), Lx,cudaMemcpyHostToDevice);
  cudaMemcpy3D(d_f2,pitchf2,h_f2,Ly*sizeof(float),Lz*sizeof(float),Ly*sizeof(float),Lz*sizeof(float), Lx,cudaMemcpyHostToDevice);
  cudaMemcpy3D(d_f3,pitchf3,h_f3,Ly*sizeof(float),Lz*sizeof(float),Ly*sizeof(float),Lz*sizeof(float), Lx,cudaMemcpyHostToDevice);
  cudaMemcpy3D(d_f4,pitchf4,h_f4,Ly*sizeof(float),Lz*sizeof(float),Ly*sizeof(float),Lz*sizeof(float), Lx,cudaMemcpyHostToDevice);
  cudaMemcpy3D(d_f5,pitchf5,h_f5,Ly*sizeof(float),Lz*sizeof(float),Ly*sizeof(float),Lz*sizeof(float), Lx,cudaMemcpyHostToDevice);
  cudaMemcpy3D(d_f6,pitchf6,h_f6,Ly*sizeof(float),Lz*sizeof(float),Ly*sizeof(float),Lz*sizeof(float), Lx,cudaMemcpyHostToDevice);
}
void LatticeBoltzmann::Colisione(void){
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,Nz);
  dim3 BlocksPerGrid(Mx,My,Mz);
  d_Colisione<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
						 d_f1,pitchf1,
						 d_f2,pitchf2,
						 d_f3,pitchf3,
						 d_f4,pitchf4,
						 d_f5,pitchf5,
						 d_f6,pitchf6,
						 d_f0new,pitchf0new,
						 d_f1new,pitchf1new,
						 d_f2new,pitchf2new,
						 d_f3new,pitchf3new,
						 d_f4new,pitchf4new,
						 d_f5new,pitchf5new,
						 d_f6new,pitchf6new);
}
void LatticeBoltzmann::ImponerCampos(int t){
  double A=10, lambda=10, omega=2*M_PI/lambda;
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,Nz);
  dim3 BlocksPerGrid(Mx,My,Mz);
  d_ImponerCampos<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
						 d_f1,pitchf1,
						 d_f2,pitchf2,
						 d_f3,pitchf3,
						 d_f4,pitchf4,
						 d_f5,pitchf5,
						 d_f6,pitchf6,
						 d_f0new,pitchf0new,
						 d_f1new,pitchf1new,
						 d_f2new,pitchf2new,
						 d_f3new,pitchf3new,
						 d_f4new,pitchf4new,
						 d_f5new,pitchf5new,
						 d_f6new,pitchf6new,
						 A*sin(omega*t));
}
void LatticeBoltzmann::Adveccione(void){
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,Nz);
  dim3 BlocksPerGrid(Mx,My,Mz);
  d_Adveccione<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f0,pitchf0,
						 d_f1,pitchf1,
						 d_f2,pitchf2,
						 d_f3,pitchf3,
						 d_f4,pitchf4,
						 d_f5,pitchf5,
						 d_f6,pitchf6,
						 d_f0new,pitchf0new,
						 d_f1new,pitchf1new,
						 d_f2new,pitchf2new,
						 d_f3new,pitchf3new,
						 d_f4new,pitchf4new,
						 d_f5new,pitchf5new,
						 d_f6new,pitchf6new);
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  //traer los datos al Host
  cudaMemcpy3D(h_f0new,Ly*sizeof(float),Lz*sizeof(float),d_f0new,pitchf0new,Ly*sizeof(float),Lz*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy3D(h_f1new,Ly*sizeof(float),Lz*sizeof(float),d_f1new,pitchf1new,Ly*sizeof(float),Lz*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy3D(h_f2new,Ly*sizeof(float),Lz*sizeof(float),d_f2new,pitchf2new,Ly*sizeof(float),Lz*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy3D(h_f3new,Ly*sizeof(float),Lz*sizeof(float),d_f3new,pitchf3new,Ly*sizeof(float),Lz*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy3D(h_f4new,Ly*sizeof(float),Lz*sizeof(float),d_f4new,pitchf4new,Ly*sizeof(float),Lz*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy3D(h_f5new,Ly*sizeof(float),Lz*sizeof(float),d_f5new,pitchf5new,Ly*sizeof(float),Lz*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  cudaMemcpy3D(h_f6new,Ly*sizeof(float),Lz*sizeof(float),d_f6new,pitchf6new,Ly*sizeof(float),Lz*sizeof(float),Lx,cudaMemcpyDeviceToHost);

  //Imprimirlos
  ofstream MiArchivo(NombreArchivo); double rho0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      rho0=h_rho(ix,iy,true);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}
// ----------------- FUNCIONES GLOBALES ----------
int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=100;

  Ondas.Inicie(0,0,0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
  }
  Ondas.Imprimase("Ondas.dat");

  return 0;
}
