//
//  ising1d.cpp
//  NQS
//
//  Created by PI on 2018/11/16.
//  Copyright © 2018 PI. All rights reserved.
//

#include <iostream>
#include <vector>
#include <complex>
#include "nqs_paper.h"

//Transverse-field Ising model in 1d  横场Ising模型
class Ising1d{
    
    //number of spins
    const int nspins_;
    
    //value of the transverse field
    const double hfield_;
    
    //option to use periodic boundary conditions  PBCs
    const bool pbc_;
    
    //pre-computed quantities
    std::vector<std::complex<double> > mel_;  //复数向量
    std::vector<std::vector<int> > flipsh_;   //整型向量
    
public:
    
    Ising1d(int nspins,double hfield,bool pbc=true):nspins_(nspins),hfield_(hfield),pbc_(pbc){   //默认有PBC  ising1d模型参数有 自旋粒子，hfield的强度，pbc边界
        Init();
    }
    
    void Init(){
        mel_.resize(nspins_+1,0.);
        flipsh_.resize(nspins_+1);
        
        for(int i=0;i<nspins_;i++){
            mel_[i+1]=-hfield_;
            flipsh_[i+1]=std::vector<int>(1,i);
        }
        std::cout<<"# Using the 1d Transverse-field Ising model with h = "<<hfield_<<std::endl;
    }
    
    //Finds the non-zero matrix elements of the hamiltonian
    //on the given state
    //i.e. all the state' such that <state'|H|state> = mel(state') \neq 0  不等于0
    //state' is encoded as the sequence of spin flips to be performed on state
    void FindConn(const std::vector<int> & state,std::vector<std::vector<int> > & flipsh,std::vector<std::complex<double> > & mel){
        //函数的参数为整型一维向量state，二维整型向量flipsh，一维复数double型mel
        mel.resize(nspins_+1);
        flipsh.resize(nspins_+1);
        
        //assigning pre-computed matrix elements and spin flips 分配预先计算的矩阵元素和旋转翻转
        mel=mel_;  //矩阵元素
        flipsh=flipsh_;  //旋转翻转
        
        //computing interaction part Sz*Sz
        mel[0]=0.;
        
        for(int i=0;i<(nspins_-1);i++){
            mel[0]-=double(state[i]*state[i+1]);
        }
        
        if(pbc_){
            mel[0]-=double(state[nspins_-1]*state[0]);
        }
        
    }
    
    int MinFlips()const{
        return 1;
    }
    
};
