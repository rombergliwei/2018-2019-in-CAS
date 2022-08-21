//
//  main.cpp
//  NQS
//
//  Created by PI on 2018/11/16.
//  Copyright © 2018 PI. All rights reserved.
//

#include "nqs_paper.h"

int main(int argc, char *argv[]){
    
    auto opts=ReadOptions(argc,argv);  //ReadOptions是一个定义的函数
    
    //Definining the neural-network wave-function
    Nqs wavef(opts["filename"]);   //Nqs为新定义的一个class wavef是Nqs类的一个对象
    
    int nsweeps=std::stod(opts["nsweeps"]);   //nsweep = 扫描的次数
    int nspins=wavef.Nspins();   //nspins = 可见层的元素个数
    
    //Problem hamiltonian inferred from file name  选择的模型
    std::string model=opts["model"];
    
    bool printastes=opts.count("filestates");
    
    int seed=std::stoi(opts["seed"]);  //随机数种子
    
    if(model=="Ising1d"){
        double hfield=std::stod(opts["hfield"]);  //定义hfield的大小
        Ising1d hamiltonian(nspins,hfield);  //Ising1d是新定义的一个class ,hamiltonian 为Ising1d的一个对象，参数是napins和hfield
        
        //Defining and running the sampler   选择模型后运行sampler
        Sampler<Nqs,Ising1d> sampler(wavef,hamiltonian,seed);   //采样函数的参数为选择的波函数，给定的哈密顿量，随机数种子
        if(printastes){
            sampler.SetFileStates(opts["filestates"]);
        }
        sampler.Run(nsweeps);  //只有一个参数nsweeps
    }
    else if(model=="Heisenberg1d"){
        double jz=std::stod(opts["jz"]);
        Heisenberg1d hamiltonian(nspins,jz);   //Heisenberg1d是新定义的一个class
        
        //Defining and running the sampler
        Sampler<Nqs,Heisenberg1d> sampler(wavef,hamiltonian,seed);
        if(printastes){
            sampler.SetFileStates(opts["filestates"]);
        }
        sampler.Run(nsweeps);
    }
    else if(model=="Heisenberg2d"){
        double jz=std::stod(opts["jz"]);
        Heisenberg2d hamiltonian(nspins,jz);   //Heisenberg2d是新定义的一个class
        
        //Defining and running the sampler
        Sampler<Nqs,Heisenberg2d> sampler(wavef,hamiltonian,seed);
        if(printastes){
            sampler.SetFileStates(opts["filestates"]);
        }
        sampler.Run(nsweeps);
    }
    else{
        std::cerr<<"#The given input file does not correspond to one of the implemented problem hamiltonians";
        std::abort();
    }
    
}
