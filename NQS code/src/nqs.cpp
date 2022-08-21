//
//  nqs.cpp
//  NQS
//
//  Created by PI on 2018/11/16.
//  Copyright © 2018 PI. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <cassert>
#include "nqs_paper.h"

class Nqs{
    
    //Neural-network weights   网络权重一般取复数，能完整描述波函数的振幅和相位
    std::vector<std::vector<std::complex<double> > > W_;    //二维复数向量
    
    //Neural-network visible bias  可见层的偏置值
    std::vector<std::complex<double> > a_;
    
    //Neural-network hidden bias   隐含层的偏置值
    std::vector<std::complex<double> > b_;
    
    //Number of hidden units   隐含层的元素数量
    int nh_;
    
    //Number of visible units    可见层的元素数量
    int nv_;
    
    //look-up tables  查找表
    std::vector<std::complex<double> > Lt_;
    
    //Useful quantities for safe computation of ln(cosh(x))
    const double log2_;
    
public:
    
    Nqs(std::string filename):log2_(std::log(2.)){
        LoadParameters(filename);
    }
    
    //computes the logarithm of the wave-function  计算波函数的对数
    inline std::complex<double> LogVal(const std::vector<int> & state)const{
        
        std::complex<double> rbm(0.,0.);   //rbm 一个复数
        
        for(int v=0;v<nv_;v++){
            rbm+=a_[v]*double(state[v]);
        }
        
        for(int h=0;h<nh_;h++){
            
            std::complex<double> thetah=b_[h];
            
            for(int v=0;v<nv_;v++){
                thetah+=double(state[v])*(W_[v][h]);
            }
            
            rbm+=Nqs::lncosh(thetah);
        }
        
        return rbm;
    }
    
    //computes the logarithm of Psi(state')/Psi(state)  Psi就是wave-function
    //where state' is a state with a certain number of flipped spins
    //the vector "flips" contains the sites to be flipped
    //look-up tables are used to speed-up the calculation   查找表用于提高计算速度
    inline std::complex<double> LogPoP(const std::vector<int> & state,const std::vector<int> & flips)const{
        //函数的返回一个double的复数，参数是state，flips
        if(flips.size()==0){
            return 0.;
        }
        
        std::complex<double> logpop(0.,0.);
        
        //Change due to the visible bias
        for(const auto & flip : flips){   //就是循环，定义了变量begin是flip=0，end是flips flips类似下标，指的是有多少个偏置值，state指的是可见层的configurations
            logpop-=a_[flip]*2.*double(state[flip]);
        }
        
        //Change due to the interaction weights
        for(int h=0;h<nh_;h++){   //循环次数是隐含层元素的个数
            const std::complex<double> thetah=Lt_[h];  //定义新变量thetah，用于保存查找表，表的大小为隐含层元素的个数
            std::complex<double> thetahp=thetah;    //定义新变量thetahp
            
            for(const auto & flip : flips){
                thetahp-=2.*double(state[flip])*(W_[flip][h]);
            }
            logpop+= ( Nqs::lncosh(thetahp)-Nqs::lncosh(thetah) );
        }
        
        return logpop;
    }
    
    //
    inline std::complex<double> PoP(const std::vector<int> & state,const std::vector<int> & flips)const{
        return std::exp(LogPoP(state,flips));  //exp是计算e的x次方的函数
    }
    
    //initialization of the look-up tables  查找表的初始化，函数的参数是state一维整型向量
    void InitLt(const std::vector<int> & state){  //
        Lt_.resize(nh_);    //Lt_是一个一维数组，大小为隐含层的元素个数
        
        for(int h=0;h<nh_;h++){
            Lt_[h]=b_[h];
            for(int v=0;v<nv_;v++){
                Lt_[h]+=double(state[v])*(W_[v][h]);
            }
        }
        
    }
    
    //updates the look-up tables after spin flips  查找表的更新
    //the vector "flips" contains the indices of sites to be flipped  函数的参数为state一维向量，还有整型向量flips (SM-s12)
    void UpdateLt(const std::vector<int> & state,const std::vector<int> & flips){
        if(flips.size()==0){
            return;
        }
        
        for(int h=0;h<nh_;h++){
            for(const auto & flip : flips){
                Lt_[h]-=2.*double(state[flip])*W_[flip][h];  //就是SM-s12的公式
            }
        }
    }
    
    //loads the parameters of the wave-function from a given file  加载wf的参数
    void LoadParameters(std::string filename){  //将文件名作为参数，读取内容到相应的变量里
        
        std::ifstream fin(filename.c_str());
        
        if(!fin.good()){
            std::cerr<<"# Error : Cannot load from file "<<filename<<" : file not found."<<std::endl;
            std::abort();
        }
        
        fin>>nv_;
        fin>>nh_;
        
        if(!fin.good() || nv_<0 || nh_<0){
            std::cerr<<"# Trying to load from an invalid file.";
            std::cerr<<std::endl;
            std::abort();
        }
        
        a_.resize(nv_);   //将可见层的偏置值数量设为可见层的数量
        b_.resize(nh_);   //将隐含层的偏置值数量设为隐含层的数量
        W_.resize(nv_,std::vector<std::complex<double> > (nh_));  //将权值矩阵的行X列设为：可见层X隐含层
        
        for(int i=0;i<nv_;i++){  //将可见层的偏置值放到a_[]中
            fin>>a_[i];
        }
        for(int j=0;j<nh_;j++){   //将隐含层的偏置值放到b_[]中
            fin>>b_[j];
        }
        for(int i=0;i<nv_;i++){    //将权值放到W_[]中
            for(int j=0;j<nh_;j++){
                fin>>W_[i][j];
            }
        }
        
        if(!fin.good()){
            std::cerr<<"# Trying to load from an invalid file.";
            std::cerr<<std::endl;
            std::abort();
        }
        
        std::cout<<"# NQS loaded from file "<<filename<<std::endl;
        std::cout<<"# N_visible = "<<nv_<<"  N_hidden = "<<nh_<<std::endl;
    }
    
    //ln(cos(x)) for real argument
    //for large values of x we use the asymptotic expansion  求双曲余弦函数
    inline double lncosh(double x)const{
        const double xp=std::abs(x);
        if(xp<=12.){
            return std::log(std::cosh(xp));
        }
        else{
            return xp-log2_;
        }
    }
    
    //ln(cos(x)) for complex argument
    //the modulus is computed by means of the previously defined function
    //for real argument   复数的求双曲余弦函数
    inline std::complex<double> lncosh(std::complex<double> x)const{
        const double xr=x.real();
        const double xi=x.imag();
        
        std::complex<double> res=Nqs::lncosh(xr);
        res +=std::log( std::complex<double>(std::cos(xi),std::tanh(xr)*std::sin(xi)) );
        
        return res;
    }
    
    //total number of spins
    //equal to the number of visible units   所有的自旋粒子等于可见层的元素个数
    inline int Nspins()const{
        return nv_;
    }
    
};
