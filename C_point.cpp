#include <C_point.h>

C_point::C_point()
{
    x = 0.0;
    y = 0.0;
    vx = 0.0;
    vy = 0.0;
    Px = 0.0;
    Py = 0.0;
}

C_point::C_point(double _x, double _y)
{
    x = _x;
    y = _y;
    vx = 0.0;
    vy = 0.0;
    Px = 0.0;
    Py = 0.0;
}

C_point::C_point(const C_point& p)
{
    x = p.x;
    y = p.y;
    vx = p.vx;
    vy = p.vy;
    Px = p.Px;
    Py = p.Py;
}

C_point::C_point(C_point& p)
{
    x = p.x;
    y = p.y;
    vx = p.vx;
    vy = p.vy;
    Px = p.Px;
    Py = p.Py;
}

C_point C_point::operator= (C_point const& p)
{
    if(this==&p) return *this;

    x = p.x;
    y = p.y;
    vx = p.vx;
    vy = p.vy;
    Px = p.Px;
    Py = p.Py;

    return *this;
}

void C_point::show()
{
    std::cout << "position " << x << " " << y << std::endl;
    std::cout << "speed " << vx << " " << vy << std::endl;
    std::cout << "image energy " << Px << " " << Py << std::endl;
}
