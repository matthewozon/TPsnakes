#include <C_curve.h>

C_curve::C_curve()
{
}


C_curve::C_curve(unsigned short _N)
{
    for(unsigned short i=0 ; i<_N ; i++)
    {
        C_point p;
        C.push_back(p);
    }
}

C_curve::C_curve(std::vector<C_point> _C)
{
    C = _C;
}

C_curve::C_curve(const C_curve& pC)
{
    C = pC.C;
}

C_curve::C_curve(C_curve& pC)
{
    C = pC.C;
}


C_curve C_curve::operator= (C_curve const& otherC)
{
    if(this==&otherC) return *this;

    C = otherC.C;

    return *this;
}




void C_curve::saveToFile(std::string fileNameVectX, std::string fileNameVectY)
{
    std::ofstream myfileX, myfileY;
    myfileX.open (fileNameVectX.data(), std::ios::out | std::ios::app);
    myfileY.open (fileNameVectY.data(), std::ios::out | std::ios::app);
    if(myfileX.is_open() && myfileY.is_open())
    {
        for(unsigned short i=0 ; i<C.size() ; i++)
        {
            myfileX << C.at(i).x << ", ";
            myfileY << C.at(i).y << ", ";
        }
        myfileX << std::endl;
        myfileY << std::endl;
        myfileX.close();
        myfileY.close();
    }
}


void C_curve::circle(double radius, double cx, double cy, unsigned short N)
{
    C.clear();
    C_point p;
    double theta;
    for(unsigned short i=0 ; i<N ; i++)
    {
        theta = ((double) i)*2.0*PI/((double)N);
        p.x = cx+radius*cos(theta);
        p.y = cy+radius*sin(theta);
        C.push_back(p);
    }
    return;
}

void C_curve::square(double edge, double cx, double cy, unsigned short N)
{
    C.clear();
    C_point p;
    double s;
    for(unsigned short i=0 ; i<N ; i++)
    {
        s = ((double) i)/((double)N);
        if(s<0.25)
        {
            p.x = cx + 0.5*edge;
            p.y = cy + (((s-0.0)/0.25)-0.5)*edge;
        }
        else if(s<0.5 && s>=0.25)
        {
            p.x = cx + (0.5-((s-0.25)/0.25))*edge;
            p.y = cy + 0.5*edge;
        }
        else if(s<0.75 && s>=0.5)
        {
            p.x = cx - 0.5*edge;
            p.y = cy + (0.5-((s-0.5)/0.25))*edge;
        }
        else
        {
            p.x = cx + (((s-0.75)/0.25)-0.5)*edge;
            p.y = cy - 0.5*edge;
        }
        C.push_back(p);
    }
    return;
}
