#pragma GCC diagnostic error "-std=c++11"
#include <iostream>
#include "nqs.hh"

using namespace std;
using namespace nqs;

int main(){

  int nspins=20;

  typedef Ising1d Hamiltonian;
  double h_field=1;
  double Jz=1;
  Hamiltonian hamiltonian(nspins,h_field,Jz);

  //Defining the Rbm State
  typedef Rbm RbmState;
  int nhidden=20;
  RbmState rbm(nspins,nhidden);   

  int seed=12345;
  rbm.InitRandomPars(seed,0.01); 

  //The Gibbs sampler
  typedef Gibbs<RbmState> Sampler;
  Sampler sampler(rbm);

  //Using a simple Stochastic Gradient Descent optimizer
  typedef Sgd Optimizer;
  double eta=0.2;   
  Sgd opt(eta);

  Variational<Hamiltonian,RbmState,Sampler,Optimizer> var(hamiltonian,sampler,opt);  

  int batch_size=100;
  int max_iter=10000000;  
  var.Run(batch_size,max_iter);

}
