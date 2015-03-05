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
    if(argc!=7)
    {
        std::cout << "commande " << argv[0] << " dt<double> lambda1<double> lambda2<double> lambda3<double> nb iter<int> fileName<char string>" << std::endl;
        return -2;
    }

    //chargement de l'image
    C_imgMatrix<double> imgDoub(argv[6]);

    //creation d'une courbe
    C_curve my_curve;

    //initialisation de la courbe a une forme circulaire (attention la courbe oit etre contenu dans l'image)
    int nb_point = 80;
    double r=50.0, cx=25.0, cy=25.0;
    my_curve.circle(r,cx,cy,nb_point);

    //instantiation d'un objet snake
    C_snakes A(my_curve,imgDoub);

    //initialisation de l'agorithme

    //evolution temporelle de la courbe
    for(int i=0 ; i<atoi(argv[5]) ; i++)
    {
        //faire un pas temporel
    }

    return 1;
}
