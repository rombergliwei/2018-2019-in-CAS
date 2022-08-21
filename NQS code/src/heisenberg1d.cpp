//
//  heisenberg1d.cpp
//  NQS
//
//  Created by PI on 2018/11/16.
//  Copyright © 2018 PI. All rights reserved.
//

#include <iostream>
#include <vector>
#include <complex>
#include "nqs_paper.h"

//Anti-ferromagnetic Heisenberg model in 1d
class Heisenberg1d{
    
    //number of spins  旋转次数
    const int nspins_;
    
    //option to use periodic boundary conditions
    const bool pbc_;
    
    //coupling constant  耦合常数
    const double jz_;
    
public:
    
    Heisenberg1d(int nspins,double jz,bool pbc=true):nspins_(nspins),pbc_(pbc),jz_(jz){
        std::cout<<"# Using the 1d Heisenberg model with J_z = "<<jz_<<std::endl;
    }
    
    
    //Finds the non-zero matrix elements of the hamiltonian
    //on the given state
    //i.e. all the state' such that <state'|H|state> = mel(state') \neq 0
    //state' is encoded as the sequence of spin flips to be performed on state
    void FindConn(const std::vector<int> & state,std::vector<std::vector<int> > & flipsh,std::vector<std::complex<double> > & mel){
        
        mel.resize(1);
        flipsh.resize(1);
        
        //computing interaction part Sz*Sz
        mel[0]=0.;
        
        for(int i=0;i<(nspins_-1);i++){
            mel[0]+=double(state[i]*state[i+1]);
        }
        
        if(pbc_){
            mel[0]+=double(state[nspins_-1]*state[0]);
        }
        
        mel[0]*=jz_;
        
        //Looks for possible spin flips
        for(int i=0;i<(nspins_-1);i++){
            if(state[i]!=state[i+1]){
                mel.push_back(-2);
                flipsh.push_back(std::vector<int>({i,i+1}));
            }
        }
        
        if(pbc_){
            if(state[nspins_-1]!=state[0]){
                mel.push_back(-2);
                flipsh.push_back(std::vector<int>({nspins_-1,0}));
            }
        }
        
    }
    
    int MinFlips()const{
        return 2;
    }
    
    
};
