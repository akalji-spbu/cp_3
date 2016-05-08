//
//  main.cpp
//  cp_3
//
//  Created by Nikolay Tikhonov on 13.04.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "matrix.hpp"
#include "SNLE.hpp"

void MN();
void MMN();
void MMNW();
void MMNH();

int main(int argc, const char * argv[]){

    MN();
    std::cout<<"==========================================="<<std::endl;
    
    MMN();
    std::cout<<"==========================================="<<std::endl;
    
    MMNW();
    std::cout<<"==========================================="<<std::endl;
    
    MMNH();
    std::cout<<"==========================================="<<std::endl;
    
    return 0;
}

void MN(){
    Matrix  X0(10,1);
    //INITIAL APPROXIMATION
    X0.Add(0, 0, 0.5);
    X0.Add(1, 0, 0.5);
    X0.Add(2, 0, 1.5);
    X0.Add(3, 0, -1);
    X0.Add(4, 0, -0.5);
    X0.Add(5, 0, 1.5);
    X0.Add(6, 0, 0.5);
    X0.Add(7, 0, -0.5);
    X0.Add(8, 0, 1.5);
    X0.Add(9, 0, -1.5);
    //END INITIAL APPROXIMATION
    
    Matrix X_Newton(Newton(X0));
    std::cout<<"Newton X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    X_Newton.Show();
    std::cout<<"Check Newton X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    Matrix Check_Newton(Func(X_Newton));
    Check_Newton.Show();
    
}
void MMN(){
    Matrix  X0(10,1);
    //INITIAL APPROXIMATION
    X0.Add(0, 0, 0.5);
    X0.Add(1, 0, 0.5);
    X0.Add(2, 0, 1.5);
    X0.Add(3, 0, -1);
    X0.Add(4, 0, -0.5);
    X0.Add(5, 0, 1.5);
    X0.Add(6, 0, 0.5);
    X0.Add(7, 0, -0.5);
    X0.Add(8, 0, 1.5);
    X0.Add(9, 0, -1.5);
    //END INITIAL APPROXIMATION
    
    Matrix X_Modified_Newton(ModifiedNewton(X0));
    std::cout<<"Modified Newton X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    X_Modified_Newton.Show();
    std::cout<<"Check Modified Newton X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    Matrix Check_Modified_Newton(Func(X_Modified_Newton));
    Check_Modified_Newton.Show();
    
}
void MMNW(){
    Matrix  X0(10,1);
    //INITIAL APPROXIMATION
    X0.Add(0, 0, 0.5);
    X0.Add(1, 0, 0.5);
    X0.Add(2, 0, 1.5);
    X0.Add(3, 0, -1);
    X0.Add(4, 0, -0.5);
    X0.Add(5, 0, 1.5);
    X0.Add(6, 0, 0.5);
    X0.Add(7, 0, -0.5);
    X0.Add(8, 0, 1.5);
    X0.Add(9, 0, -1.5);
    //END INITIAL APPROXIMATION
    
    unsigned k = 10;
    std::cout<<"To "<<k<<std::endl;
    Matrix X_Modified_Newton_Which(Modified_Newton_Which(X0, k));
    std::cout<<"Modified Newton which X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    X_Modified_Newton_Which.Show();
    std::cout<<"Check Modified Newton which X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    Matrix Check_Modified_Which_Newton(Func(X_Modified_Newton_Which));
    Check_Modified_Which_Newton.Show();
    std::cout<<"////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
    
}
void MMNH(){
    Matrix  X0(10,1);
    //INITIAL APPROXIMATION
    X0.Add(0, 0, 0.5);
    X0.Add(1, 0, 0.5);
    X0.Add(2, 0, 1.5);
    X0.Add(3, 0, -1);
    X0.Add(4, 0, -0.5);
    X0.Add(5, 0, 1.5);
    X0.Add(6, 0, 0.5);
    X0.Add(7, 0, -0.5);
    X0.Add(8, 0, 1.5);
    X0.Add(9, 0, -1.5);
    //END INITIAL APPROXIMATION
    
    unsigned k = 10;
    std::cout<<"To "<<k<<std::endl;
    Matrix X_Modified_Newton_Hibrid(Modified_Newton_Hibrid(X0, k));
    std::cout<<"Modified Newton hibrid X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    X_Modified_Newton_Hibrid.Show();
    std::cout<<"Check Modified Newton hibrid X:"<<std::endl;
    std::cout<<"___________________________________________"<<std::endl;
    Matrix Check_Modified_Hibrid_Newton(Func(X_Modified_Newton_Hibrid));
    Check_Modified_Hibrid_Newton.Show();
    std::cout<<"////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
}
