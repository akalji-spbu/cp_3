//
//  SNLE.hpp
//  cp_3
//
//  Created by Nikolay Tikhonov on 21.04.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#ifndef SNLE_hpp
#define SNLE_hpp
#include <ctime>
#include "matrix.hpp"
#include "SLE.h"

Matrix MJacobi(const Matrix &X);
Matrix Func(const Matrix &X);
Matrix Newton(const Matrix &X0);
Matrix ModifiedNewton(const Matrix &X0);
Matrix Modified_Newton_Which(const Matrix &X0, unsigned k);
Matrix Modified_Newton_Hibrid(const Matrix &X0, unsigned k);

#endif /* SNLE_hpp */
