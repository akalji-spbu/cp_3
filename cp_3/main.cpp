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

int main(int argc, const char * argv[]){
    Matrix X0(10,1);
    //END INITIAL APPROXIMATION
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

    
    //Matrix M(10,10);
    //M=MJacobi(X);
    //M.Show();
    Matrix P(10,10), L(10,10), U(10,10), Ja(10,10),
    
    return 0;
}
