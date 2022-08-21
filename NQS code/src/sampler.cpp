//
//  sampler.cpp
//  NQS
//
//  Created by PI on 2018/11/16.
//  Copyright © 2018 PI. All rights reserved.
//

#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ctime>
#include "nqs_paper.h"

//Simple Monte Carlo sampling of a spin  蒙特卡罗采样
//Wave-Function
template<class Wf,class Hamiltonian> class Sampler{
    
    //wave-function
    Wf & wf_;
    
    //Hamiltonian
    Hamiltonian & hamiltonian_;
    
    //number of spins  自旋粒子
    const int nspins_;
    
    //current state in the sampling
    std::vector<int> state_;
    
    //random number generators and distributions  随机数生成器和发布器
    std::mt19937 gen_;
    std::uniform_real_distribution<> distu_;
    std::uniform_int_distribution<> distn_;
    
    //sampling statistics     抽样统计
    double accept_;
    double nmoves_;
    
    //container for indices of randomly chosen spins to be flipped
    std::vector<int> flips_;
    
    //option to write the sampled configuration on a file
    bool writestates_;
    std::ofstream filestates_;
    
    //quantities needed by the hamiltonian
    //non-zero matrix elements  非零矩阵元素
    std::vector<std::complex<double> > mel_;
    
    //flip connectors for the hamiltonian(see below for details)
    std::vector<std::vector<int> > flipsh_;
    
    //storage for measured values of the energy   存储能量的测量值
    std::vector<std::complex<double> > energy_;
    
public:
    
    Sampler(Wf & wf,Hamiltonian & hamiltonian,int seed):
    wf_(wf),hamiltonian_(hamiltonian),distu_(0,1),nspins_(wf.Nspins()),distn_(0,nspins_-1)
    {
        
        writestates_=false;
        Seed(seed);
        ResetAv();
    }
    
    ~Sampler(){     //析构函数
        if(writestates_){
            filestates_.close();
        }
    }
    
    //Uniform random number in [0,1) uniform函数返回的是0和1之间的随机数
    inline double Uniform(){
        return distu_(gen_);
    }
    
    inline void Seed(int seed){
        if(seed<0){
            gen_.seed(std::time(nullptr));
        }
        else{
            gen_.seed(seed);
        }
    }
    
    //Random spin flips (max 2 spin flips in this implementation)  这部分并不是很懂
    //if mag0=true, when doing 2 spin flips the total magnetization is kept equal to 0
    //if mag0=true,当进行2次旋转翻转时，总磁化强度保持等于0
    inline bool RandSpin(std::vector<int> & flips,int nflips,bool mag0=true){
        flips.resize(nflips);
        
        flips[0]=distn_(gen_);
        if(nflips==2){
            flips[1]=distn_(gen_);
            if(!mag0){
                return flips[1]!=flips[0];
            }
            else{
                return state_[flips[1]]!=state_[flips[0]];
            }
        }
        
        return true;
    }
    
    //Initializes a random spin state
    //if mag0=true, the initial state is prepared with zero total magnetization
    //if mag0=true, 则初始状态准备为零总磁化
    void InitRandomState(bool mag0=true){
        state_.resize(nspins_);
        for(int i=0;i<nspins_;i++){
            state_[i]=(Uniform()<0.5)?(-1):(1);
        }
        
        if(mag0){
            int magt=1;
            if(nspins_%2){
                std::cerr<<"# Error : Cannot initializate a random state with zero magnetization for odd number of spins"<<std::endl;
            std:abort();
            }
            while(magt!=0){
                magt=0;
                for(int i=0;i<nspins_;i++){
                    magt+=state_[i];
                }
                if(magt>0){
                    int rs=distn_(gen_);
                    while(state_[rs]<0){
                        rs=distn_(gen_);
                    }
                    state_[rs]=-1;
                    magt-=1;
                }
                else if(magt<0){
                    int rs=distn_(gen_);
                    while(state_[rs]>0){
                        rs=distn_(gen_);
                    }
                    state_[rs]=1;
                    magt+=1;
                }
            }
        }
    }
    
    void ResetAv(){  //重置Av
        accept_=0;
        nmoves_=0;
    }
    
    inline double Acceptance()const{
        return accept_/nmoves_;
    }
    
    void Move(int nflips){
        
        //Picking "nflips" random spins to be flipped
        if(RandSpin(flips_,nflips)){
            
            //Computing acceptance probability
            double acceptance=std::norm(wf_.PoP(state_,flips_));
            
            //Metropolis-Hastings test  测试MH算法  SM--s11附近
            if(acceptance>Uniform()){
                
                //Updating look-up tables in the wave-function  更新查找表
                wf_.UpdateLt(state_,flips_);
                
                //Moving to the new configuration  转到新的configuration
                for(const auto& flip : flips_){
                    state_[flip]*=-1;
                }
                
                accept_+=1;
            }
        }
        
        nmoves_+=1;
    }
    
    void SetFileStates(std::string filename){
        writestates_=true;
        filestates_.open(filename.c_str());
        if(!filestates_.is_open()){
            std::cerr<<"# Error : Cannot open file "<<filename<<" for writing"<<std::endl;
            std::abort();
        }
        else{
            std::cout<<"# Saving sampled configuration to file "<<filename<<std::endl;
        }
    }
    
    void WriteState(){
        for(const auto & spin_value : state_){
            filestates_<<std::setw(2)<<spin_value<<" ";
        }
        filestates_<<std::endl;
    }
    
    //Measuring the value of the local energy  在当前状态测量能量的值
    //on the current state
    void MeasureEnergy(){
        std::complex<double> en=0.;
        
        //Finds the non-zero matrix elements of the hamiltonian
        //on the given state
        //i.e. all the state' such that <state'|H|state> = mel(state') \neq 0   不等于0
        //state' is encoded as the sequence of spin flips to be performed on state
        hamiltonian_.FindConn(state_,flipsh_,mel_);
        
        for(int i=0;i<flipsh_.size();i++){
            en+=wf_.PoP(state_,flipsh_[i])*mel_[i];
        }
        
        energy_.push_back(en);
    }
    
    
    //Run the Monte Carlo sampling
    //nsweeps is the total number of sweeps to be done
    //thermfactor is the fraction of nsweeps to be discarded during the initial equilibration
    //sweepfactor set the number of single spin flips per sweeps to nspins*sweepfactor
    //nflipss is the number of random spin flips to be done, it is automatically set to 1 or 2 depending on the hamiltonian
    //运行蒙特卡罗采样
       // nsweeps是要完成的扫描总数
       // thermfactor是在初始平衡期间要丢弃的nsweeps的分数
       // sweepfactor将每次扫描的单个自旋翻转次数设置为nspins * sweepfactor
       // nflipss是要完成的随机旋转翻转次数，根据汉密尔顿函数自动设置为1或2
    void Run(double nsweeps,double thermfactor=0.1,int sweepfactor=1,int nflipss=-1){  //后面三个参数都是默认设置
        
        int nflips=nflipss;
        
        if(nflips==-1){
            nflips=hamiltonian_.MinFlips();  //设置最小翻转次数为 1
        }
        
        //checking input consistency  检查输入的一致性
        if(nflips>2 || nflips<1){  //旋转反转次数应该等于1或2
            std::cerr<<"# Error : The number of spin flips should be equal to 1 or 2.";
            std::cerr<<std::endl;
            std::abort();
        }
        if(thermfactor>1 || thermfactor<0){  //热化因子应该在0和1之间
            std::cerr<<"# Error : The thermalization factor should be a real number between 0 and 1";
            std::cerr<<std::endl;
            std::abort();
        }
        if(nsweeps<50){
            std::cerr<<"# Error : Please enter a number of sweeps sufficiently large (>50)";
            std::cerr<<std::endl;
            std::abort();
        }
        
        std::cout<<"# Starting Monte Carlo sampling"<<std::endl;
        std::cout<<"# Number of sweeps to be performed is "<<nsweeps<<std::endl;
        
        InitRandomState();
        
        flips_.resize(nflips);
        
        //initializing look-up tables in the wave-function
        wf_.InitLt(state_);   //state最开始的入口
        
        ResetAv();
        
        std::cout<<"# Thermalization... ";
        std::flush(std::cout);
        
        //thermalization
        for(double n=0;n<nsweeps*thermfactor;n+=1){
            for(int i=0;i<nspins_*sweepfactor;i++){
                Move(nflips);
            }
        }
        std::cout<<" DONE "<<std::endl;
        std::flush(std::cout);
        
        ResetAv();
        
        std::cout<<"# Sweeping... ";
        std::flush(std::cout);
        
        //sequence of sweeps
        for(double n=0;n<nsweeps;n+=1){
            for(int i=0;i<nspins_*sweepfactor;i++){
                Move(nflips);
            }
            if(writestates_){
                WriteState();
            }
            MeasureEnergy();
        }
        std::cout<<" DONE "<<std::endl;
        std::flush(std::cout);
        
        OutputEnergy();
        
    }
    
    void OutputEnergy(){
        int nblocks=50;
        
        int blocksize=std::floor(double(energy_.size())/double(nblocks));
        
        double enmean=0;
        double enmeansq=0;
        
        double enmean_unblocked=0;
        double enmeansq_unblocked=0;
        
        for(int i=0;i<nblocks;i++){
            double eblock=0;
            for(int j=i*blocksize;j<(i+1)*blocksize;j++){
                eblock+=energy_[j].real();
                assert(j<energy_.size());
                
                double delta=energy_[j].real()-enmean_unblocked;
                enmean_unblocked+=delta/double(j+1);
                double delta2=energy_[j].real()-enmean_unblocked;
                enmeansq_unblocked+=delta*delta2;
            }
            eblock/=double(blocksize);
            double delta=eblock-enmean;
            enmean+=delta/double(i+1);
            double delta2=eblock-enmean;
            enmeansq+=delta*delta2;
        }
        
        enmeansq/=(double(nblocks-1));
        enmeansq_unblocked/=(double((nblocks*blocksize-1)));
        
        double estav=enmean/double(nspins_);
        double esterror=std::sqrt(enmeansq/double(nblocks))/double(nspins_);
        
        int ndigits=std::log10(esterror);
        if(ndigits<0){
            ndigits=-ndigits+2;
        }
        else{
            ndigits=0;
        }
        
        std::cout<<"# Estimated average energy per spin : "<<std::endl;
        std::cout<<"# "<<std::scientific<<std::setprecision(ndigits)<<estav;
        std::cout<<" +/-  "<<std::setprecision(0)<<esterror<<std::endl;
        std::cout<<"# Error estimated with binning analysis consisting of ";
        std::cout<<nblocks<<" bins "<<std::endl;
        std::cout<<"# Block size is "<<blocksize<<std::endl;
        std::cout<<"# Estimated autocorrelation time is ";
        std::cout<<std::setprecision(0);
        std::cout<<0.5*double(blocksize)*enmeansq/enmeansq_unblocked<<std::endl;
    }
    
};
