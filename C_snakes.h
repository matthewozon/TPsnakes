#ifndef C_SNAKES_H
#define C_SNAKES_H

#include <C_imgMatrix.h>
#include <C_curve.h>

class C_snakes
{
public:
    C_snakes(); //constructeur par defaut : initialise les attributs a une valeur par defaut
    C_snakes(C_curve snakeInit, C_imgMatrix<double> _img); //initialise le sanke avec la courbe snakeInit et initialise l'image avec _img
    C_snakes(const C_snakes& pS);//recopie
    C_snakes(C_snakes& pS);//recopie


    //snake evolution parameter
    double dt;
    double lambda1, lambda2, lambda3;

    //methode pour le calcul de l'energie du snake
    double getIntEnergy();
    double getExtEnergy();

    //methodes pour l'evolution du snake
    void initAlgorithm(void);
    void timeIteration(void);

    //le snake a faire evoluer
    C_curve S;

    //les attribut qui peuvent etre utile
    C_matrix<double> Ainv; //inverse du system
    C_matrix<double> D1; //matrice de derivation
    C_matrix<double> D2; //matrice de derivation seconde
    C_matrix<double> D4; //matrice de derivation 4eme
    void buildAInvMatrix(void);
    void buildD1(void);
    void buildD2(void);
    void buildD4(void);

    //methodes pour les seconds menbres
    C_matrix<double> getBxk(void);
    C_matrix<double> getByk(void);

    //attributs concernant l'image et son energie
    C_imgMatrix<double> img;
    C_imgMatrix<double> PX;
    C_imgMatrix<double> PY;
    C_imgMatrix<double> imgEnergy;

    //methode qui met a jours l'energie image en chaque point de la courbe
    void updateP(void); //to be implemented by sutdent

    //form column vectors containing x or y coordinates
    C_matrix<double> getXcoordinates(void);
    C_matrix<double> getYcoordinates(void);

};

#endif // C_SNAKES_H
