//
//  readoptions.cpp
//  NQS
//
//  Created by PI on 2018/11/16.
//  Copyright © 2018 PI. All rights reserved.
//
#include <getopt.h>
#include <map>
#include <iostream>
#include <string>
#include "nqs_paper.h"

//Various utilities to read the command line options

std::string FindModel(std::string strarg){   //找模型
    std::size_t found = strarg.find("Ising");  //ising模型
    if (found!=std::string::npos){             //npos表示最大值
        return "Ising1d";
    }
    found = strarg.find("Heisenberg1d");    //海森堡1维模型
    if (found!=std::string::npos){
        return "Heisenberg1d";
    }
    
    found = strarg.find("Heisenberg2d");   //海森堡2维模型
    if (found!=std::string::npos){
        return "Heisenberg2d";
    }
    return "None";
}

std::string FindCoupling(std::string strarg){   //找不同模型对应的一个参数hfield or Jz
    size_t found = strarg.find("_");
    size_t found1= std::string::npos;
    size_t found2= std::string::npos;
    if (found!=std::string::npos){
        found1=strarg.find("_",found+1);
        if (found1!=std::string::npos){
            found2=strarg.find("_",found1+1);
        }
    }
    if(found1!=std::string::npos && found2!=std::string::npos){
        return(strarg.substr(found1+1,found2-found1-1));
    }
    else{
        std::cerr<<"# Error : the filename is not in the format specified for the Ising/Heisenberg model"<<std::endl;
        std::abort();  //abort()中断程序的执行
    }
    std::exit(0);
    return "error";
}

void PrintHeader(){   //打印介绍信息
    std::cout<<std::endl;
    std::cout<<"\t|   Neural-network quantum states sampler   |"<<std::endl;
    std::cout<<"\t| written by Giuseppe Carleo, December 2016 |"<<std::endl<<std::endl;
}

void PrintInfoMessage(){   //打印介绍信息
    std::cout<<"Usage : ./nqs_run OPTIONS"<<std::endl<<std::endl;
    
    std::cout<<"Allowed OPTIONS are : "<<std::endl<<std::endl;
    
    std::cout<<"--filename=...  "<<std::endl;
    std::cout<<"\tname of the file containing neural-network weights"<<std::endl;
    std::cout<<"\t(chose one in directories Ground/ or Unitary/)"<<std::endl<<std::endl;
    
    std::cout<<"--nsweeps=... "<<std::endl;
    std::cout<<"\tnumber of Monte Carlo sweeps"<<std::endl;
    std::cout<<"\t(default value is 1.0e4)"<<std::endl<<std::endl;   //1.0e4 表示0.0001
    
    std::cout<<"--seed=... "<<std::endl;
    std::cout<<"\tinteger seed for pseudo-random numbers"<<std::endl;
    std::cout<<"\tseed<0 sets it to the internal clock value"<<std::endl;
    std::cout<<"\t(default value is -1)"<<std::endl<<std::endl;
    
    std::cout<<"--filestates=... "<<std::endl;
    std::cout<<"\tname of the file to print sampled configurations"<<std::endl;
    std::cout<<"\t(by default it is not set)"<<std::endl<<std::endl;
}

std::map<std::string,std::string> ReadOptions(int argc,char *argv[]){  //ReadOptions 函数的定义 map模板类-红黑树
    
    PrintHeader();
    
    std::map<std::string,std::string> options;
    
    if(argc==1){
        PrintInfoMessage();
        std::exit(0);
    }
    
    while (1)
    {
        static struct option long_options[] =      //定义一个option结构体
        {
            /* These options don’t set a flag.
             We distinguish them by their indices. */
            {"filename",  required_argument, 0, 'a'},
            {"nsweeps",  required_argument, 0, 'b'},
            {"seed",    required_argument, 0, 'c'},
            {"filestates",    required_argument, 0, 'd'},
            {0, 0, 0, 0}
        };
        
        /* getopt_long stores the option index here. */
        int option_index = 0;
        
        int c = getopt_long (argc, argv, "a:b:c:d:",
                             long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
        
        switch (c)
        {
            case 'a':
                options["filename"]=optarg;
                break;
                
            case 'b':
                options["nsweeps"]=optarg;
                break;
                
            case 'c':
                options["seed"]=optarg;
                break;
                
            case 'd':
                options["filestates"]=optarg;
                break;
                
            case '?':
                PrintInfoMessage();
                break;
                
            default:
                std::abort ();
        }
    }
    
    if(options.count("filename")==0){     //count函数是STL里面的 统计容器中等于value元素的个数
        std::cerr<<"# Error: Option filename must be specified with the option --filename=FILENAME"<<std::endl;
        std::abort();
    }
    
    if(options.count("nsweeps")==0){
        options["nsweeps"]="1.0e4";
    }
    
    if(options.count("seed")==0){
        options["seed"]="-1";
    }
    
    options["model"]=FindModel(options["filename"]);  //将读取到的模型放到options中
    
    if(options["model"]=="Ising1d"){
        options["hfield"]=FindCoupling(options["filename"]);
    }
    else{
        options["jz"]=FindCoupling(options["filename"]);
    }
    
    return options;
}
