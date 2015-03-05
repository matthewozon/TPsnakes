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
    //met des parametre par defaut qu'il faut estimer... manuellement
    dt = 1.0;
    lambda1 = 1.0;
    lambda2 = 1.0;
    lambda3 = 1.0;
    //
    img.resize(_img.getNbRow(),_img.getNbColumn()); //met l'image a la bonne taille
    //operation surchargee
    img = _img;
    //operation surchargee
    S = snakeInit;
}

//cptor
C_snakes::C_snakes(const C_snakes& pS)
{
    //il faudrait le faire si ca sert a qqch...
}
//cptor
C_snakes::C_snakes(C_snakes& pS)
{
    //il faudrait le faire si ca sert a qqch...
}

//initialisation de l'algorithme (on suppose que l'on a deja la forme initiale de la courbe)
void C_snakes::initAlgorithm(void)
{
    //construit le systeme inverxe pour faire evoluer la courbe
    buildAInvMatrix();

    //construit aussi la matrice de derivation premiere (pour le calcul de l'energie interne)
    buildD1();

    //lisse l'image par convolution avec un noyau gaussien 13*13
    C_matrix<double> G(13,13);
    double tx, ty, sum=0.0;
    for(int i=0 ; i<13 ; i++)
    {
        tx = -3.0+ ((double)i)*(6.0/12.0);
        for(int j=0 ; j<13 ; j++)
        {
            ty = -3.0+ ((double)j)*(6.0/12.0);
            G(i,j) = exp(-SQR(tx)-SQR(ty));
            sum+=G(i,j);
        }
    }
    G = G*(1.0/sum);
    //redimensionnement de imgEnergy (qui contiendra l'energie image, mais pour l'instant sert a stocker l'image lissee)
    imgEnergy.resize(img.getNbRow(),img.getNbColumn());
    imgEnergy = img.conv2(G);

    //definition des variable contenant le gradient de l'image
    C_imgMatrix<double> gradx(img.getNbRow(),img.getNbColumn());
    C_imgMatrix<double> grady(img.getNbRow(),img.getNbColumn());

    //derive selon : column/x (cf C_matrix.h)

    //derive selon : row/y (cf C_matrix.h)

    //calcule de la norme du gradient

    //calcule le gradient de l'energie image (Eext) et stock les composante dans Px et PY

    return;
}

//un pas temporelle de l'algorithme
void C_snakes::timeIteration(void)
{
    //met a jour les composante du gradient de l'energie image pour chaque point de la courbe

    //resolution des systeme equation 6 et 7 du TP (un pas temporel)

    //mettre a jour la courbe
    for(unsigned int i=0 ; i<S.C.size() ; i++)
    {
        S.C.at(i);//do something
    }
    return;
}

//construction des seconds membres (cf equation 6 et 7 du sujet de TP)
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

//cf equation 6 et 7 du TP
void C_snakes::buildAInvMatrix(void)
{
    //construit la matrice du systeme
    //inversion de cette matrice puis stockage dans Ainv (cf C_matrix.h)
    return;
}

//construction des matrice de derivation par rapport au parametre de la courbe (ce n'est pas par rapport a x ou y)
void C_snakes::buildD1(void)
{
    D1.resize(S.C.size(),S.C.size());
    D1 = 0.0;
    for(unsigned short i=0 ; i<S.C.size() ; i++)
    {
        //diag sup
        if((i+1)<S.C.size()) D2(i,(unsigned short) (i+1)) = 1.0;
        //diag inf
        if(i>0) D2(i,(unsigned short) (i-1)) = -1.0;
    }
    D1((unsigned short) 0,D1.endC) = -1.0;
    D1(D1.endL,(unsigned short) 0) = 1.0;
    return;
}
void C_snakes::buildD2(void)
{
    D2.resize(S.C.size(),S.C.size());
    D2 = 0.0;
    for(unsigned short i=0 ; i<S.C.size() ; i++)
    {
        D2(i,i) = -2.0;
        if((i+1)<S.C.size()) D2(i,(unsigned short) (i+1)) = 1.0;
        if(i>0) D2(i,(unsigned short) (i-1)) = 1.0;
    }
    D2((unsigned short) 0,D2.endC) = 1.0;
    D2(D2.endL,(unsigned short) 0) = 1.0;
    return;
}
void C_snakes::buildD4(void)
{
    D4.resize(S.C.size(),S.C.size());
    D4 = 0.0;
    for(unsigned short i=0 ; i<S.C.size() ; i++)
    {
        D4(i,i) = 6.0;
        if((i+1)<S.C.size()) D4(i,(unsigned short) (i+1)) = -4.0;
        if(i>0) D4(i,(unsigned short) (i-1)) = -4.0;
        if((i+2)<S.C.size()) D4(i,(unsigned short) (i+2)) = 1.0;
        if(i>1) D4(i,(unsigned short) (i-2)) = 1.0;
    }
    D4((unsigned short) 0,D4.endC) = -4.0;
    D4(D4.endL,(unsigned short) 0) = -4.0;
    D4((unsigned short) 0,(unsigned short) (D4.endC-1)) = 1.0;
    D4((unsigned short) (D4.endL-1),(unsigned short) 0) = 1.0;
    D4((unsigned short) 1,D4.endC) = 1.0;
    D4(D4.endL,(unsigned short) 1) = 1.0;
    return;
}

//mis a jour de Px et Py pour chaque point de la courbe (gradient de l'energie image)
void C_snakes::updateP(void)
{
    //attention a ne pas melanger les repere ligne/colone et x/y
    return;
}

//extrait tout les coordonnees x de la courbe et retourne une matrice
C_matrix<double> C_snakes::getXcoordinates(void)
{
    C_matrix<double> Xpos(S.C.size(),1);
    for(unsigned short i=0 ; i<S.C.size() ; i++)
    {
        Xpos(i,(unsigned short)1) = S.C.at(i).x;
    }
    return Xpos;
}

//extrait tout les coordonnees y de la courbe et retourne une matrice
C_matrix<double> C_snakes::getYcoordinates(void)
{
    C_matrix<double> Ypos(S.C.size(),1);
    for(unsigned short i=0 ; i<S.C.size() ; i++)
    {
        Ypos(i,(unsigned short)1) = S.C.at(i).y;
    }
    return Ypos;
}

//calcul de l'energie interieure de la courbe (cf page 2 de l'enonce de TP)
double C_snakes::getIntEnergy()
{
    double E = 0.0;
    return E;
}

//calcul de l'energie image de la courbe (contribution de l'energie image le long de la courbe) (cf page 2 de l'enonce de TP)
double C_snakes::getExtEnergy()
{
    double E = 0.0;
    return -E;
}
