#include<vector>
#include<cmath>
double f(double r, const std::vector<double> & y,int id,double lambda);
void rk4(double r,double h,std::vector<double>&y,double lambda);
void rk4(double r,double h,std::vector<double>&y,double lambda)
{
  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());
  // k1
  for(int i=0; i< y.size(); i++){
    k1[i]=h*f(r,y,i,lambda);
  }
  //k2 aux
  for(int i=0; i< y.size() ;i++){
    aux[i]=y[i]+k1[i]/2;
  }
  //k2
  for(int i=0; i< y.size();i++){
    k2[i]= h*f(r+h/2,aux,i,lambda);
  }
  //k3 aux
  for(int i=0; i< y.size(); i++){
    aux[i]=y[i]+k2[i]/2;
  }
  //k3
  for(int i=0; i< y.size(); i++){
    k3[i]=h*f(r+h/2,aux,i,lambda);
  }
  //k4 aux
  for(int i=0; i< y.size(); i++){
    aux[i]=y[i]+k3[i];
  }
  //k4
  for(int i=0; i< y.size(); i++){
    k4[i]=h*f(r+h,aux,i,lambda);
  }
  //runge kutta
  for(int i=0; i<y.size();i++)
    {
      
      y[i]+= (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    
    }
}
double f(double r, const std::vector<double> & y,int id,double lambda)
{
  if(0==id)
    {
      return y[1];
    }
  else if (1==id)
    {
      return -(y[1]/r)-(y[0]*lambda*lambda);
    }
  else{
    std::cerr<<"Error!Id does not exist->"<<id<<std::endl;
    exit(1);
 
}
}
