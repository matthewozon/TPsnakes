#include <iostream>
#include <stdio.h>
#include <sstream>

#include <C_imgMatrix.h>

#include <C_point.h>
#include <C_curve.h>
#include <C_snakes.h>
#include <vector>





int main(int argc, char *argv[])
{
    if(argc!=6)
    {
        std::cout << "commande " << argv[0] << " dt<double> lambda1<double> lambda2<double> lambda3<double> nb iter<int> fileName<char string>" << std::endl;
        return -2;
    }

    //creation d'une courbe
    C_curve my_curve;

    //chargement de l'image
    C_imgMatrix<double> imgDoub(argv[5]);


    return 1;
}
