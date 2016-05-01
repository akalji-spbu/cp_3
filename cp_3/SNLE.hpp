//
//  SNLE.hpp
//  cp_3
//
//  Created by Nikolay Tikhonov on 21.04.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#ifndef SNLE_hpp
#define SNLE_hpp
#include "matrix.hpp"

Matrix MJacobi(Matrix &X);

//THE SYSTEM
/////////////////////////////////////////////////////////////////////////////////////////////////////////
double  f1(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f2(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f3(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f4(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f5(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f6(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f7(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f8(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f9(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
double  f10(double  x1, double  x2, double  x3, double  x4, double  x5, double  x6, double  x7, double  x8, double  x9, double  x10);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//END THE SYSTEM
#endif /* SNLE_hpp */
