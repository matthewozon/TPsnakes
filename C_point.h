#ifndef C_POINT_H
#define C_POINT_H

#include <iostream>

class C_point //or a struct
{
public:
    C_point();
    C_point(double _x, double _y);
    C_point(const C_point& p);
    C_point(C_point& p);

    //position: x are the columns and y the rows
    double x,y;

    //speed
    double vx, vy;

    //external energy
    double Px, Py;

    C_point operator= (C_point const& c);

    void show();
};

#endif // C_POINT_H
