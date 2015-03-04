//good to go!

#ifndef C_MATRIX_H
#define C_MATRIX_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <C_vector.h>
#include <cmath>
#include <cstring>


const double epsilon=0.0000000001;


#define ABS(a) ((a) >= 0 ? (a) : -(a))
#define SQR(a) ((a)*(a))
#define SWAP(x,y) do \
   { unsigned char swap_temp[sizeof(x) == sizeof(y) ? (signed)sizeof(x) : -1]; \
     memcpy(swap_temp,&y,sizeof(x)); \
     memcpy(&y,&x,       sizeof(x)); \
     memcpy(&x,swap_temp,sizeof(x)); \
    } while(0)


#define NEAREST 0
#define LINEAR  1
#define SMALL_NUM_F 1e-37
#define SMALL_NUM_D 1e-307

//this is a template, meaning that you need to specify the type when you instanciate a object.
template<class dataType> class C_matrix
{
public:
    //constructors and destructors
    //default ctor: create empty matrix (use resize to set a new size and allocate the image container)
    C_matrix();
    //ctor: create mtrix of size _L rows, _C columns
    C_matrix(unsigned short _L, unsigned short _C);
    //copy ctor
    C_matrix(const C_matrix &M);
    //ctor: create a matrix copying _M
    C_matrix(dataType **_M, unsigned short _L, unsigned short _C);
    //dtor: free memory
    virtual ~C_matrix();


    //IO methods
    dataType& operator()(const int l, const int c);
    const dataType& operator()(const int l, const int c) const;
    dataType& operator()(const long l, const long c);
    const dataType& operator()(const long l, const long c) const;
    dataType& operator()(const short l, const short c);
    const dataType& operator()(const short l, const short c) const;
    dataType& operator()(const unsigned int l, const unsigned int c);
    const dataType& operator()(const unsigned int l, const unsigned int c) const;
    dataType& operator()(const unsigned short l, const unsigned short c);
    const dataType& operator()(const unsigned short l, const unsigned short c) const;
    dataType& operator()(const unsigned long l, const unsigned long c);
    const dataType& operator()(const unsigned long l, const unsigned long c) const;
    dataType& operator()(const float l, const float c);
    const dataType& operator()(const float l, const float c, unsigned short INTERPOLATION=NEAREST) const;
    dataType& operator()(const double l, const double c);
    const dataType& operator()(const double l, const double c, unsigned short INTERPOLATION=NEAREST) const;

    //overload operators
    C_matrix<dataType> operator= (C_matrix const& c);
    C_matrix<dataType> operator= (dataType const& x);
    C_matrix<dataType> operator+ (C_matrix const& c);
    C_matrix<dataType> operator+ (dataType const& x);
    C_matrix<dataType> operator- (C_matrix const& c);
    C_matrix<dataType> operator- (dataType const& x);
    C_matrix<dataType> operator* (C_matrix const& c);
    //C_vector<dataType>& operator* (C_vector<dataType> const& c); //obselete
    C_matrix operator* (const dataType& x);

    C_matrix<dataType> operator== (C_matrix const& c);
    C_matrix<dataType> operator== (dataType const& x);

    C_matrix<dataType> operator> (C_matrix const& c);
    C_matrix<dataType> operator> (dataType const& x);
    C_matrix<dataType> operator>= (C_matrix const& c);
    C_matrix<dataType> operator>= (dataType const& x);

    C_matrix<dataType> operator< (C_matrix const& c);
    C_matrix<dataType> operator< (dataType const& x);
    C_matrix<dataType> operator<= (C_matrix const& c);
    C_matrix<dataType> operator<= (dataType const& x);

    //cast operator (later)
    //operator int();

    //usual operations
    C_matrix<dataType> m_abs(void); //must be multithreaded
    dataType maxVal(void);
    dataType minVal(void);
    double sum(void);
    double mean(void);
    double var(void);
    C_matrix<dataType> SQRT(void); //must be multithreaded
    C_matrix<dataType> dotPower(double alpha);//octave .^ //must be multithreaded
    C_matrix<dataType> dotExp(double alpha=1.0); //must be multithreaded
    C_matrix<dataType> dotLog(double alpha=1.0); //must be multithreaded
    C_matrix<dataType> dotProduct(C_matrix const& B);//octave .* //must be multithreaded
    C_matrix<dataType> dotDiv(C_matrix const& B);//octave ./ //must be multithreaded
    C_matrix<dataType> dotProduct(C_matrix& B); //must be multithreaded
    C_matrix<dataType> dotDiv(C_matrix& B); //must be multithreaded

    //distance map
    C_matrix<double> bwdistEuclidean(void); //Algorithme de Danielson

    //convolution //must be multithreaded
    C_matrix<dataType> conv2(C_matrix<dataType> const& h);
    C_matrix<dataType> conv2(C_matrix<dataType>& h);
    C_matrix<dataType> gradX(void);
    C_matrix<dataType> gradY(void);


    //mesh (later)
    void meshRow(dataType minValue, dataType maxValue);
    void meshColumn(dataType minValue, dataType maxValue);
    //linespace (later)

    //matrix inversion
    C_matrix<dataType> inv(void);

    //diagonalisation (later)
    //svd (later)
    //random generation (later)
    //fft (later)

    //tools for linear system solving
    C_matrix<dataType> Transpose(void);
    C_matrix<dataType> LU(void);
    C_matrix<dataType> LUP(C_vector<int> &Indice);
    C_matrix<dataType> LMU(void);
    C_vector<dataType> LineAlgEq_LU(C_vector<dataType> &B);



    //other tools
    unsigned short getNbRow(void)const {return m_L;}
    unsigned short getNbColumn(void)const {return m_C;}
    void show(void);
    void save(std::string fileName);

    //to use carefully
    bool resize(unsigned short newL, unsigned short newC);
    unsigned short endL;
    unsigned short endC;

protected:
    unsigned short m_L;
    unsigned short m_C;
    dataType** m_A;

    //privates tools to allocate and free memory
    dataType** allocate(unsigned short M, unsigned short N);
    void deallocation(dataType** B, unsigned short M, unsigned short N);

};
template<class dataType> C_matrix<dataType>::C_matrix()
{
    m_L=0;
    m_C=0;
    m_A=NULL;
    endL = m_L-1;
    endC = m_C-1;
}
template<class dataType> C_matrix<dataType>::C_matrix(unsigned short _L, unsigned short _C) : m_L(_L),m_C(_C)
{
    m_A = allocate(m_L,m_C);
    if(m_A==NULL)
    {
        m_L=0;
        m_C=0;
    }
    endL = m_L-1;
    endC = m_C-1;
}

template<class dataType> C_matrix<dataType>::C_matrix(dataType **_M, unsigned short _L, unsigned short _C) : m_L(_L),m_C(_C)
{
    m_A = allocate(m_L, m_C);
    if(m_A==NULL)
    {
        m_L=0;
        m_C=0;
    }
    else
    {
        for(unsigned short i=0 ; i<m_L ; i++)
        {
            for(unsigned short j=0 ; j<m_C ; j++)
                m_A[i][j] = _M[i][j];
        }
    }
    endL = m_L-1;
    endC = m_C-1;
}

template<class dataType> C_matrix<dataType>::C_matrix(const C_matrix &X)
{
    m_L = X.getNbRow();
    m_C = X.getNbColumn();
    m_A = allocate(m_L, m_C);
    if(m_A==NULL)
    {
        m_L=0;
        m_C=0;
    }
    else
    {
        for(unsigned short i=0 ; i<m_L ; i++)
        {
            for(unsigned short j=0 ; j<m_C ; j++)
            {
                m_A[i][j] = X(i,j);
            }
        }
    }
    endL = m_L-1;
    endC = m_C-1;

}

template<class dataType> C_matrix<dataType>::~C_matrix()
{
    deallocation(m_A,m_L,m_C);
}

//to use carefully
template<class dataType> bool C_matrix<dataType>::resize(unsigned short newL, unsigned short newC)
{
    deallocation(m_A,m_L,m_C);
    m_A = NULL;
    m_L = 0;
    m_C = 0;
    m_A = allocate(newL, newC);
    if(m_A==NULL) return false;
    m_L = newL;
    m_C = newC;
    endL = m_L-1;
    endC = m_C-1;
    return true;
}

template<class dataType> dataType** C_matrix<dataType>::allocate(unsigned short M, unsigned short N)
{
    dataType** A = new dataType*[M];
    if(A==NULL) return A;

    for(unsigned short i=0 ; i<M ; i++)
        A[i] = new dataType[N];

    //should check if all allocation is OK
    bool cond=true;
    for(unsigned short i=0 ; i<M ; i++)
    {
        if(A[i]==NULL)
        {
            cond=false;
            i=M;
        }
    }
    if(!cond)
    {
        deallocation(A,M,N);
    }
    return A;
}

template<class dataType> void C_matrix<dataType>::deallocation(dataType** B, unsigned short M, unsigned short N)
{
    if(B!=NULL)
    {
        for(unsigned short i=0 ; i<M ; i++)
        {
            if(B[i]!=NULL)
            {
                delete B[i];
            }
        }
        delete B;
        B=NULL;
    }
    return;
}

template<class dataType> void C_matrix<dataType>::show(void)
{
    std::cout << "--------------------------------------------------------" << std::endl;
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            std::cout << m_A[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

template<class dataType> void C_matrix<dataType>::save(std::string fileName)
{
    std::ofstream myfileX;
    myfileX.open (fileName.data(), std::ios::out );

    if(myfileX.is_open())
    {
        for(unsigned short i=0 ; i<m_L ; i++)
        {
            for(unsigned short j=0 ; j<m_C ; j++)
            {
                myfileX << m_A[i][j] << "\t";
            }
            myfileX << std::endl;
        }
        myfileX.close();
    }
}



template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator= (C_matrix const& c)
{
    if(this==&c)return *this;

    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";

    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            this->m_A[i][j] = c(i,j);
        }
    }
    return *this;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator= (dataType const& x)
{
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            m_A[i][j] = x;
        }
    }

    return *this;
}




template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator+ (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = m_A[i][j] + c(i,j);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator+ (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = m_A[i][j] + x;
        }
    }

    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator- (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);// = new C_matrix(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = m_A[i][j] - c(i,j);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator- (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);// = new C_matrix
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = m_A[i][j] - x;
        }
    }

    return B;
}


template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator* (C_matrix const& c)
{
    if(m_C!=c.getNbRow()) throw "mismatch dimension matrix";

    C_matrix B(m_L,c.getNbColumn());
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<c.getNbColumn() ; j++)
        {
            dataType S = (dataType) 0;
            for(unsigned short k=0 ; k<m_C ; k++)
            {
                S += m_A[i][k]*c(k,j);
            }
            B(i,j) = S;
        }
    }
    return B;
}


template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator* (const dataType& x)
{
    C_matrix B(m_L,m_C);
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = m_A[i][j]*x;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator== (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]==c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator== (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]==x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator> (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]>c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator> (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]>x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator>= (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]>=c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator>= (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]>=x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}




template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator< (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]<c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator< (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]<x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator<= (C_matrix const& c)
{
    if(c.getNbColumn()!=this->getNbColumn() || c.getNbRow()!=this->getNbRow()) throw "dimension matrix must agree";
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]<=c(i,j)) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::operator<= (dataType const& x)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(m_A[i][j]<=x) B(i,j)=1.0;
            else B(i,j) = 0.0;
        }
    }
    return B;
}






template<class dataType> dataType& C_matrix<dataType>::operator()(const int l, const int c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const int l, const int c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const long l, const long c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const long l, const long c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const short l, const short c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const short l, const short c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const unsigned int l, const unsigned int c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const unsigned int l, const unsigned int c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}
template<class dataType> dataType& C_matrix<dataType>::operator()(const unsigned short l, const unsigned short c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const unsigned short l, const unsigned short c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const unsigned long l, const unsigned long c)
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const unsigned long l, const unsigned long c) const
{
    if( l>= this->m_L || c>= this->m_C )
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    return m_A[l][c];
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const float l, const float c)
{
    if( (l<0.0) || (l>=(float) (this->m_L)) || (c<0.0) || (c>= (float) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    unsigned int ll = (unsigned int) floor((double) l);
    unsigned int cc = (unsigned int) floor((double) c);
    return m_A[ll][cc];

}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const float l, const float c, unsigned short INTERPOLATION) const
{
    if( (l<0.0) || (l>=(float) (this->m_L)) || (c<0.0) || (c>= (float) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    if(INTERPOLATION==NEAREST)
    {
        unsigned int ll = (unsigned int) floor((double) l);
        unsigned int cc = (unsigned int) floor((double) c);
        return m_A[ll][cc];
    }
    else
    {
        //not yet available
        unsigned int ll = (unsigned int) floor((double) l);
        unsigned int cc = (unsigned int) floor((double) c);
        //return m_A[ll][cc];
        if(l<SMALL_NUM_F && c<SMALL_NUM_F)
        {
            return m_A[ll][cc];
        }
        if(l<SMALL_NUM_F)//column interpolation
        {
            unsigned int cc1 = cc+1;
            return (dataType) (((double)m_A[ll][cc])*( ((double)cc1) - ((double)c) ) + ((double)m_A[ll][cc1])*( ((double)c) - ((double)cc) ));
        }
        if(c<SMALL_NUM_F)//row interpolation
        {
            unsigned int ll1 = ll+1;
            return (dataType) (((double)m_A[ll][cc])*( ((double)ll1) - ((double)l) ) + ((double)m_A[ll1][cc])*( ((double)l) - ((double)ll) ));
        }
        //main case: bilinear interpolation
        unsigned int ll1 = ll+1;
        unsigned int cc1 = cc+1;
        double valuell = (((double)m_A[ll][cc])*( ((double)cc1) - ((double)c) ) + ((double)m_A[ll][cc1])*( ((double)c) - ((double)cc) ));
        double valuell1 = (((double)m_A[ll1][cc])*( ((double)cc1) - ((double)c) ) + ((double)m_A[ll1][cc1])*( ((double)c) - ((double)cc) ));
        return (dataType) (valuell*( ((double)ll1) - ((double)l) ) + valuell1*( ((double)l) - ((double)ll) ));
    }
}

template<class dataType> dataType& C_matrix<dataType>::operator()(const double l, const double c)
{
    if( (l<0.0) || (l>=(double) (this->m_L)) || (c<0.0) || (c>= (double) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    unsigned int ll = (unsigned int) floor(l);
    unsigned int cc = (unsigned int) floor(c);
    return m_A[ll][cc];

}
template<class dataType> const dataType& C_matrix<dataType>::operator()(const double l, const double c, unsigned short INTERPOLATION) const
{
    if( (l<0.0) || (l>=(double) (this->m_L)) || (c<0.0) || (c>= (double) (this->m_C)))
    {
        throw "index must be larger than 0 and smaller than the matrix dimensions";
    }
    if(INTERPOLATION==NEAREST)
    {
        unsigned int ll = (unsigned int) floor(l);
        unsigned int cc = (unsigned int) floor(c);
        return m_A[ll][cc];
    }
    else
    {
        //not yet available
        unsigned int ll = (unsigned int) floor(l);
        unsigned int cc = (unsigned int) floor(c);
        //return m_A[ll][cc];
        if(l<SMALL_NUM_F && c<SMALL_NUM_F)
        {
            return m_A[ll][cc];
        }
        if(l<SMALL_NUM_F)//column interpolation
        {
            unsigned int cc1 = cc+1;
            return (dataType) (((double)m_A[ll][cc])*( ((double)cc1) - c ) + ((double)m_A[ll][cc1])*( c - ((double)cc) ));
        }
        if(c<SMALL_NUM_F)//row interpolation
        {
            unsigned int ll1 = ll+1;
            return (dataType) (((double)m_A[ll][cc])*( ((double)ll1) - l ) + ((double)m_A[ll1][cc])*( l - ((double)ll) ));
        }
        //main case: bilinear interpolation
        unsigned int ll1 = ll+1;
        unsigned int cc1 = cc+1;
        double valuell = (((double)m_A[ll][cc])*( ((double)cc1) - c ) + ((double)m_A[ll][cc1])*( c - ((double)cc) ));
        double valuell1 = (((double)m_A[ll1][cc])*( ((double)cc1) - c ) + ((double)m_A[ll1][cc1])*( c - ((double)cc) ));
        return (dataType) (valuell*( ((double)ll1) - l ) + valuell1*( l - ((double)ll) ));
    }
}



template<class dataType> C_matrix<dataType> C_matrix<dataType>::m_abs(void)
{
    C_matrix<dataType> B(this->m_L,this->m_C);// = new C_matrix
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(i,j) = ABS(m_A[i][j]);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::SQRT(void)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = sqrt(m_A[i][j]);
        }
    }
    return B;
}

template<class dataType> C_matrix<double> C_matrix<dataType>::bwdistEuclidean(void)
{
    C_matrix<double> B(m_L,m_C), B2(m_L,m_C), B3(m_L,m_C), B4(m_L,m_C);
    C_matrix<double> dL(m_L,m_C), dC(m_L,m_C), dL2(m_L,m_C), dC2(m_L,m_C), dL3(m_L,m_C), dC3(m_L,m_C), dL4(m_L,m_C), dC4(m_L,m_C);
    double BIGNUM = ((double)m_L)*((double)m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(ABS(m_A[i][j])<SMALL_NUM_F)//if not in the ocean
            {
                B(i,j) = 0.0; B2(i,j) = 0.0; B3(i,j) = 0.0; B4(i,j) = 0.0;
                dL(i,j) = 0.0; dL2(i,j) = 0.0; dL3(i,j) = 0.0; dL4(i,j) = 0.0;
                dC(i,j) = 0.0; dC2(i,j) = 0.0; dC3(i,j) = 0.0; dC4(i,j) = 0.0;
            }
            else
            {
                B(i,j) = BIGNUM; B2(i,j) = BIGNUM; B3(i,j) = BIGNUM; B4(i,j) = BIGNUM;
                dL(i,j) = BIGNUM; dL2(i,j) = BIGNUM; dL3(i,j) = BIGNUM; dL4(i,j) = BIGNUM;
                dC(i,j) = BIGNUM; dC2(i,j) = BIGNUM; dC3(i,j) = BIGNUM; dC4(i,j) = BIGNUM;
            }
        }
    }

    double v, vi, vj;
    bool firstIslandMet=false;
    //first sweep
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(!firstIslandMet && B(i,j)<BIGNUM)
            {
                firstIslandMet = true;
            }
            if(B(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                v = SQR(dL(i,j))+SQR(dC(i,j));
                if(i>0) vi = SQR(dL((unsigned short)(i-1),j)-1.0)+SQR(dC((unsigned short)(i-1),j));
                else vi = BIGNUM;

                if(j>0) vj = SQR(dL(i,(unsigned short)(j-1)))+SQR(dC(i,(unsigned short)(j-1))-1.0);
                else vj = BIGNUM;

                if(v<vj && v<vi)
                {
                    //no change
                }
                else if(vi<v && vi<vj)
                {
                    if(i>0)
                    {
                        dL(i,j) = dL((unsigned short)(i-1),j)-1.0;
                        dC(i,j) = dC((unsigned short)(i-1),j);
                        if(B(i,j)>sqrt(vi)) B(i,j) = sqrt(vi);
                    }
                }else
                {
                    if(j>0)
                    {
                        dL(i,j) = dL(i,(unsigned short)(j-1));
                        dC(i,j) = dC(i,(unsigned short)(j-1))-1.0;
                        if(B(i,j)>sqrt(vj)) B(i,j) = sqrt(vj);
                    }
                }
            }
        }
    }

    //second sweep
    firstIslandMet=false;
    for(long i=(long)(m_L-1) ; i>=0 ; i--)
    {
        for(long j=0 ; j<m_C ; j++)
        {
            if(!firstIslandMet && ABS(m_A[i][j])<SMALL_NUM_F)
            {
                firstIslandMet = true;
            }
            if(B2(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                try
                {
                    v = SQR(dL2(i,j))+SQR(dC2(i,j));
                    if(i<(long)(m_L-1)) vi = SQR(dL2((long)(i+1),j)+1.0)+SQR(dC2((long)(i+1),j));
                    else vi = BIGNUM;

                    if(j>0) vj = SQR(dL2(i,(long)(j-1)))+SQR(dC2(i,(long)(j-1))-1.0);
                    else vj = BIGNUM;

                    if(v<vj && v<vi)
                    {
                        //no change
                    }
                    else if(vi<v && vi<vj)
                    {
                        if(i<(long)(m_L-1))
                        {
                            dL2(i,j) = dL2((long)(i+1),j)+1.0;
                            dC2(i,j) = dC2((long)(i+1),j);
                            B2(i,j) = sqrt(vi);
                        }
                    }
                    else
                    {
                        if(j>0)
                        {
                            dL2(i,j) = dL2(i,(long)(j-1));
                            dC2(i,j) = dC2(i,(long)(j-1))-1.0;
                            B2(i,j) = sqrt(vj);
                        }
                    }
                }
                catch(const char* a)
                {
                    std::cout << a << std::endl;
                }
            }
            if(B(i,j)>B2(i,j)) B(i,j)=B2(i,j);
        }
    }

    //third sweep
    firstIslandMet=false;
    for(long i=(long)(m_L-1) ; i>=0 ; i--)
    {
        for(long j=(m_C-1) ; j>=0 ; j-- )
        {
            if(!firstIslandMet && ABS(m_A[i][j])<SMALL_NUM_F)
            {
                firstIslandMet = true;
            }
            if(B3(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                try
                {
                    v = SQR(dL3(i,j))+SQR(dC3(i,j));
                    if(i<(long)(m_L-1)) vi = SQR(dL3((long)(i+1),j)+1.0)+SQR(dC3((long)(i+1),j));
                    else vi = BIGNUM;

                    if(j<(long)(m_C-1)) vj = SQR(dL3(i,(long)(j+1)))+SQR(dC3(i,(long)(j+1))+1.0);
                    else vj = BIGNUM;

                    if(v<vj && v<vi)
                    {
                        //no change
                    }
                    else if(vi<v && vi<vj)
                    {
                        if(i<(long)(m_L-1))
                        {
                            dL3(i,j) = dL3((long)(i+1),j)+1.0;
                            dC3(i,j) = dC3((long)(i+1),j);
                            B3(i,j) = sqrt(vi);
                        }
                    }
                    else
                    {
                        if(j<(long)(m_C-1))
                        {
                            dL3(i,j) = dL3(i,(long)(j+1));
                            dC3(i,j) = dC3(i,(long)(j+1))+1.0;
                            B3(i,j) = sqrt(vj);
                        }
                    }
                }
                catch(const char* a)
                {
                    std::cout << a << std::endl;
                }
            }
            if(B(i,j)>B3(i,j)) B(i,j)=B3(i,j);
        }
    }


    //fourth sweep
    firstIslandMet=false;
    std::cout << (unsigned short)(m_L-1) << " " << (unsigned short)(m_C-1) << std::endl;
    for(long i=0 ; i<m_L ; i++)
    {
        for(long j=(long)(m_C-1) ; j>=0 ; j--)
        {
            if(!firstIslandMet && ABS(m_A[i][j])<SMALL_NUM_F)
            {
                firstIslandMet = true;
            }
            if(B4(i,j)>SMALL_NUM_F && firstIslandMet)
            {
                try
                {
                    v = SQR(dL4(i,j))+SQR(dC4(i,j));
                    if(i>0) vi = SQR(dL4((long)(i-1),j)-1.0)+SQR(dC4((long)(i-1),j));
                    else vi = BIGNUM;

                    if(j<(long)(m_C-1)) vj = SQR(dL4(i,(long)(j+1)))+SQR(dC4(i,(long)(j+1))+1.0);
                    else vj = BIGNUM;

                    if(v<vj && v<vi)
                    {
                        //no change
                    }
                    else if(vi<v && vi<vj)
                    {
                        if(i>0)
                        {
                            dL4(i,j) = dL4((long)(i-1),j)-1.0;
                            dC4(i,j) = dC4((long)(i-1),j);
                            B4(i,j) = sqrt(vi);
                        }
                    }
                    else
                    {
                        if(j<(long)(m_C-1))
                        {
                            dL4(i,j) = dL4(i,(long)(j+1));
                            dC4(i,j) = dC4(i,(long)(j+1))+1.0;
                            B4(i,j) = sqrt(vj);
                        }
                    }
                }
                catch(const char* a)
                {
                    std::cout << a << std::endl;
                }
            }
            if(B(i,j)>B4(i,j)) B(i,j)=B4(i,j);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotPower(double alpha)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = (dataType) pow((double) m_A[i][j], alpha);
        }
    }
    return B;
}
template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotExp(double alpha)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = (dataType) exp(alpha*(double) m_A[i][j]);
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotLog(double alpha)
{
    C_matrix<dataType> B(m_L,m_C);
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B(i,j) = (dataType) log(alpha*(double) m_A[i][j]);
        }
    }
    return B;
}

template<class dataType> dataType C_matrix<dataType>::maxVal(void)
{
    dataType B = m_A[0][0];
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(B<m_A[i][j]) B = m_A[i][j];
        }
    }
    return B;
}

template<class dataType> dataType C_matrix<dataType>::minVal(void)
{
    dataType B = m_A[0][0];
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            if(B>m_A[i][j]) B = m_A[i][j];
        }
    }
    return B;
}

template<class dataType>double C_matrix<dataType>::sum(void)
{
    double B = 0.0;
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            B += (double) m_A[i][j];
        }
    }
    return B;
}

template<class dataType>double C_matrix<dataType>::mean(void)
{
    return sum()/(((double)m_L)*((double)m_C));
}

template<class dataType>double C_matrix<dataType>::var(void)
{
    double B = mean(), V=0.0;
    for(unsigned short i=0 ; i<m_L ; i++)
    {
        for(unsigned short j=0 ; j<m_C ; j++)
        {
            V += SQR(((double) m_A[i][j]) - B);
        }
    }
    return V/(((double)m_L)*((double)m_C) - 1.0);
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::Transpose(void)
{
    C_matrix<dataType> B(this->m_C,this->m_L);
    for(unsigned short i=0 ; i<this->m_L ; i++)
    {
        for(unsigned short j=0 ; j<this->m_C ; j++)
        {
            B(j,i) = m_A[i][j];
        }
    }
    return B;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotProduct(C_matrix const& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = m_A[i][j]*B(i,j);
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotProduct(C_matrix& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = (m_A[i][j])*(B(i,j));
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotDiv(C_matrix const& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix<dataType> C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = m_A[i][j]/(B(i,j));
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::dotDiv(C_matrix& B)
{
    if(this->getNbRow()!=B.getNbRow() || this->getNbColumn()!=B.getNbColumn()) throw "mismatch dimension matrix";
    C_matrix C(m_L,m_C);
    for(unsigned short i=0 ; i<C.getNbRow() ; i++)
    {
        for(unsigned short j=0 ; j<C.getNbColumn() ; j++)
        {
            C(i,j) = (m_A[i][j])/(B(i,j));
        }
    }
    return C;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::conv2(C_matrix<dataType> const& h)
{
    if((h.getNbRow() % 2 == 0) || (h.getNbColumn() % 2 == 0))
    {
        std::cout << "even kernel size are not handled. Try to change the size of the convolution kernel to an odd size" << std::endl;
    }

    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    tmp = 0.0;


    //compute convolution
    unsigned long K = h.getNbRow()/2;
    unsigned long L = h.getNbColumn()/2;
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            for(unsigned long ih=0 ; ih<2*K+1 ; ih++)
            {
                for(unsigned long jh=0 ; jh<2*L+1 ; jh++)
                {
                    if((i+K)>=ih && (j+L)>=jh && (i-ih+K)<m_L && (j-jh+L)<m_C) //check image boundaries
                    {
                        tmp(i,j) += m_A[i-ih+K][j-jh+L]*h(ih,jh);
                    }
                }
            }
        }
    }
    return tmp;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::conv2(C_matrix<dataType>& h)
{
    if((h.getNbRow() % 2 == 0) || (h.getNbColumn() % 2 == 0))
    {
        std::cout << "even kernel size are not handled. Try to change the size of the convolution kernel to an odd size" << std::endl;
    }
    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    //std::cout << m_L << " " << m_C << std::endl;
    tmp = 0.0;


    //compute convolution
    unsigned long K = h.getNbRow()/2;
    unsigned long L = h.getNbColumn()/2;
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            for(unsigned long ih=0 ; ih<2*K+1 ; ih++)
            {
                for(unsigned long jh=0 ; jh<2*L+1 ; jh++)
                {
                    if((i+K)>=ih && (j+L)>=jh && (i-ih+K)<m_L && (j-jh+L)<m_C) //check image boundaries
                    {
                        //std::cout << "doing the job " << i << " " << j << " " << ih << " " << jh << std::endl;
                        tmp(i,j) += m_A[i-ih+K][j-jh+L]*h(ih,jh);
                    }
                }
            }
        }
    }
    return tmp;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::gradX(void)
{
    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    tmp = 0.0;


    //compute convolution
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=1 ; j<m_C-1 ; j++)
        {
            tmp(i,j) = 0.5*(m_A[i][j+1]-m_A[i][j-1]);
        }
    }
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        tmp(i,(unsigned long)0) = m_A[i][1]-m_A[i][0];
        tmp(i,(unsigned long)m_C-1) = m_A[i][m_C-1]-m_A[i][m_C-2];
    }
    return tmp;
}

template<class dataType> C_matrix<dataType> C_matrix<dataType>::gradY(void)
{
    //returned matrix
    C_matrix<dataType> tmp(m_L,m_C);
    tmp = 0.0;


    //compute convolution
    for(unsigned long i=1 ; i<m_L-1 ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            tmp(i,j) = 0.5*(m_A[i+1][j]-m_A[i-1][j]);
        }
    }
    for(unsigned long i=0 ; i<m_C ; i++)
    {
        tmp((unsigned long) 0,i) = m_A[1][i]-m_A[0][i];
        tmp((unsigned long) m_L-1,i) = m_A[m_L-1][i]-m_A[m_L-2][i];
    }
    return tmp;
}
template<class dataType> void C_matrix<dataType>::random(void)
{
    srand(time(NULL));
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            m_A[i][j] = (dataType) rand();
        }
    }
    return;
}

template<class dataType> void C_matrix<dataType>::randomf(void)
{
    srand(time(NULL));
    for(unsigned long i=0 ; i<m_L ; i++)
    {
        for(unsigned long j=0 ; j<m_C ; j++)
        {
            m_A[i][j] = (dataType) (((double) rand())/((double) RAND_MAX));
        }
    }
    return;
}
template<class dataType> void C_matrix<dataType>::meshRow(dataType minValue, dataType maxValue)
{
    if(minValue>=maxValue)
    {
        std::cout << "min value is larger than max value" << std::endl;
        return;
    }
    double delta = (((double) maxValue) - ((double) minValue))/((double) (m_L-1));
    for(unsigned i=0 ; i<m_C ; i++)
    {
        for(unsigned j=0 ; j<m_L ; j++)
        {
            m_A[j][i] = (dataType) (((double) minValue) + ((double) j)*delta);
        }
    }
    return;
}

template<class dataType> void C_matrix<dataType>::meshColumn(dataType minValue, dataType maxValue)
{
    if(minValue>=maxValue)
    {
        std::cout << "min value is larger than max value" << std::endl;
        return;
    }
    double delta = (((double) maxValue) - ((double) minValue))/((double) (m_C-1));
    for(unsigned i=0 ; i<m_C ; i++)
    {
        for(unsigned j=0 ; j<m_L ; j++)
        {
            m_A[j][i] = (dataType) (((double) minValue) + ((double) i)*delta);
        }
    }
    return;
}



















//**********************
//Décomposition LU
//**********************
template<class dataType> C_matrix<dataType> C_matrix<dataType>::LU(void)
{
    C_matrix<dataType> MLU(*this);// = new C_matrix<dataType>(*this);
    C_matrix<dataType> tmp(m_L,m_C);


    std::cout<<"\nLU Decomposition Steps\n";
    std::cout<<"------------------------\n";

    if(m_L != m_C)
    {
        std::cout<<"\nColumn and Rows must be equal"<<std::endl;
    }
    else
    {
        for(int j=0; j<m_C; j++)
        {
            if(ABS(MLU(j,j))<epsilon)
            {
                std::cout<<"Pivot is equal to zero..."<<std::endl;
                break; //pivot nul
            }
            else
            {

                MLU.show();
                //Division colonne j par pivot
                for(int i=j+1; i<m_L; i++)
                {
                    MLU(i,j) = MLU(i,j)/MLU(j,j);
                }

                //Complément de Shur
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        tmp(ii,jj) = MLU(j,jj)*MLU(ii,j);
                    }
                }
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        MLU(ii,jj) = MLU(ii,jj)-tmp(ii,jj);
                    }
                }

            }
        }
    }
    std::cout<<"\t\tEnd LU Decomposition\n";
    //MLU.Affiche();
    return MLU;
}












//******************************
//Décomposition LU + Permutation
//******************************
template<class dataType> C_matrix<dataType> C_matrix<dataType>::LUP(C_vector<int> &Indice)
{
    C_matrix<dataType> MLUP(*this);
    C_matrix<dataType> tmp(m_L,m_C);
    dataType Tmp;

    std::cout<<"\nLUP Decomposition Steps\n";
    std::cout<<"------------------------\n";

    if(m_L != m_C){
        std::cout<<"\nColumn and Rows must be equal"<<std::endl;
    }
    else
    {

        //Indice contiendra les permutations
        for(int j=0; j<m_C; j++)
            Indice.set(j,(dataType) j);

        //Pour tous les pivots
        for(int j=0; j<m_C; j++)
        {
            if(ABS(MLUP(j,j))<epsilon)
            {
                std::cout<<"Pivot is equal to zero..."<<std::endl;
                break; //pivot nul
            }
            else
            {
                MLUP.show();

                //Recherche du pivot le plus grand
                dataType max=MLUP(j,j);
                int indmax=j;
                for(int i=j+1; i<m_L; i++)
                {
                    if(MLUP(i,j)>max)
                    {
                        max=MLUP(i,j);
                        indmax=i;
                    }
                }

                //Echange des lignes
                for(int jj=0;jj<m_C;jj++)
                {
                    Tmp = MLUP(j,jj);
                    MLUP(j,jj) = MLUP(indmax,jj);
                    MLUP(indmax,jj) = Tmp;
                }
                //Mise à jour du tableau Indice (=les perumtations)
                Indice.set(j,indmax);
                Indice.set(indmax,j);

                //Division colonne j par pivot
                for(int i=j+1; i<m_L; i++)
                {
                    MLUP(i,j) = MLUP(i,j)/MLUP(j,j);
                }

                //Complément de Shur
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        tmp(ii,jj) = MLUP(j,jj)*MLUP(ii,j);
                    }
                }
                for(int jj=j+1; jj<m_C; jj++)
                {
                    for(int ii=j+1; ii<m_L; ii++)
                    {
                        MLUP(ii,jj) = MLUP(ii,jj)-tmp(ii,jj);
                    }
                }

            }
        }
    }

    std::cout<<"\t\tEnd LUP Decomposition\n";
    //MLUP.Affiche();
    return MLUP;
}







//**********************
//Produit L*U
//**********************
template<class dataType> C_matrix<dataType> C_matrix<dataType>::LMU(void)
{
    C_matrix<dataType> L(m_L,m_C);
    C_matrix<dataType> U(m_L,m_C);
    C_matrix<dataType> LMU(m_L,m_C);
    L = 0.0;
    U = 0.0;
    LMU = 0.0;


    //Matrice L
    for(int i=0;i<m_L;i++)
    {
        L(i,i)=1.0;
        for (int j=0;j<i; j++)
        {
            L(i,j) = (*this)(i,j);
        }
    }

    //Matrice L
    for(int i=0;i<m_L;i++)
    {
        for (int j=i;j<m_C; j++)
        {
            U(i,j) = (*this)(i,j);
        }
    }

    L.show();
    U.show();

    LMU = L * U;
    return LMU;

}


//**********************************
//Résolution du système MX=B avec LU
//**********************************
template<class dataType> C_vector<dataType> C_matrix<dataType>::LineAlgEq_LU(C_vector<dataType> &B)
{

    C_matrix<dataType> MLU = LU();
    C_vector<dataType> X(m_C);
    C_vector<dataType> Y(m_C);


    //Solve for LY=B
    for(int i=0;i<m_L;i++)
    {
        Y.set(i,B[i]);
        for(unsigned short j=0;j<i;j++)
        {
            Y.set(i,Y[i] - MLU((unsigned short)i,j)*Y[j]);
        }
    }

    //Solve for UX=Y
    for(int i=m_L-1;i>=0;i--)
    {
        X.set(i, Y[i]);
        for(int j=m_L-1;j>i;j--)
        {
            X.set(i,X[i] - MLU(i,j)*X[j]);
        }
        X.set(i,X[i]/MLU(i,i));
    }
    return X;
}
/**/

//http://fr.wikipedia.org/wiki/%C3%89limination_de_Gauss-Jordan
template<class dataType> C_matrix<dataType> C_matrix<dataType>::inv(void)
{
    if(m_L!=m_C)
    {
        throw "matrix must be square";
    }
    std::cout << "use me carefully, I will return a matrix even though the matrix is not invertible" << std::endl;

    //copy the current matrix push
    C_matrix<dataType> tmp(*this);

    //
    long r = -1; //row index

    C_matrix<dataType> I(m_L,m_C);
    I= (dataType) 0.0;
    for(long i=0 ; i<m_L ; i++)
        I(i,i) = (dataType) 1.0;

    for(long j=0 ; j<m_C ; j++) //go through all column
    {
        //search max abs value in column j starting at index r+1
        dataType piv = -1.0;
        long k = 0;
        for(long i=r+1 ; i<m_L ; i++)
        {
            if(piv<ABS(tmp(i,j)))
            {
                piv = tmp(i,j);
                k = i;
            }
        }

        //if piv is not null
        if(piv>SMALL_NUM_F)
        {
            //divide row k by piv
            for(long i=0 ; i<m_C ; i++)
            {
                tmp(k,i) = tmp(k,i)/piv;
                I(k,i) = I(k,i)/piv;
            }

            //
            r++;

            //swapp rows k and r
            for(long i=0 ; i<m_C ; i++)
            {
                SWAP(tmp(k,i),tmp(r,i));
                SWAP(I(k,i),I(r,i));
            }

            for(long i=0 ; i<m_C ; i++)
            {
                if(i!=r)
                {
                    //compute Li <- Li - m_A[i][j]*Lr
                    dataType tmpD = tmp(i,j);
                    for(long m=0 ; m<m_C ; m++)
                    {
                        tmp(i,m) = tmp(i,m) - tmpD*tmp(r,m);
                        I(i,m) = I(i,m) - tmpD*I(r,m);
                    }
                }
            }
        }
    }
    return I;
}

#endif // C_MATRIX_H
