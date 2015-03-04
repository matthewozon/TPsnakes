#ifndef C_CURVE_H
#define C_CURVE_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <C_point.h>
#include <C_imgMatrix.h>

#define PI 3.14159265359

class C_curve
{
public:
    //default ctor: do nothing
    C_curve();
    //ctor: create a point vector of length N
    C_curve(unsigned short _N);
    //ctor: copy the vector _C into C
    C_curve(std::vector<C_point> _C);
    //copy ctors
    C_curve(const C_curve& pC);
    C_curve(C_curve& pC);

    C_curve operator= (C_curve const& otherC);


    //store the collection of points that constitute the curve
    std::vector<C_point> C;

    //init curves (change the number of point in the curve)
    void circle(double radius, double cx, double cy, unsigned short N);
    void square(double edge, double cx, double cy, unsigned short N);

    //save x and y coordinates of C into two separated files
    void saveToFile(std::string fileNameVectX, std::string fileNameVectY);

};

#endif // C_CURVE_H
