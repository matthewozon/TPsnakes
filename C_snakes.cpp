#include <C_snakes.h>

C_snakes::C_snakes()
{
    //set default settings
    dt = 1.0;
    lambda1 = 1.0;
    lambda2 = 1.0;
    lambda3 = 1.0;
}


C_snakes::C_snakes(C_curve snakeInit, C_imgMatrix<double> _img)
{
    dt = 1.0;
    lambda1 = 1.0;
    lambda2 = 1.0;
    lambda3 = 1.0;
    img.resize(_img.getNbRow(),_img.getNbColumn()); //met l'image a la bonne taille
    img = _img;
    S = snakeInit;
}


C_snakes::C_snakes(const C_snakes& pS)
{
    //
}

C_snakes::C_snakes(C_snakes& pS)
{
    //
}

void C_snakes::initAlgorithm(void)
{
    return;
}

void C_snakes::timeIteration(void)
{
    return;
}

C_matrix<double> C_snakes::getBxk(void)
{
    C_matrix<double> Bxk;
    return Bxk;
}

C_matrix<double> C_snakes::getByk(void)
{
    C_matrix<double> Byk;
    return Byk;
}

void C_snakes::buildAInvMatrix(void)
{
    return;
}


void C_snakes::buildD1(void)
{
    return;
}
void C_snakes::buildD2(void)
{
    return;
}
void C_snakes::buildD4(void)
{
    return;
}


void C_snakes::updateP(void)
{
    return;
}

C_matrix<double> C_snakes::getXcoordinates(void)
{
    C_matrix<double> Xpos;
    return Xpos;
}

C_matrix<double> C_snakes::getYcoordinates(void)
{
    C_matrix<double> Ypos;
    return Ypos;
}

double C_snakes::getIntEnergy()
{
    double E = 0.0;
    return E;
}

double C_snakes::getExtEnergy()
{
    double E = 0.0;
    return -E;
}
