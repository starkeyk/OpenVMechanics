#ifndef NURB_H
#define NURB_H

#include <vector>
#include <iomanip>
#include<math.h>

typedef std::vector<double> d_vec;
typedef std::vector<int> i_vec;
typedef std::vector<unsigned int> ui_vec;
typedef std::vector<std::vector<double> > d_mat;

class CNurbSurface;
class CNurbSurfaceGP;
class CNurbGP;

class CNurb{
    friend CNurbSurface;
    friend CNurbSurfaceGP;
    friend CNurbGP;
    //test functions
    friend int test_nurb_basis_n1_k2();
    friend int test_nurb_basis_n2_k3();
    friend int test_nurb_basis_n3_k3();
    friend int test_nurb_basis_n4_k3();
private:
    int _n;                                                 // # of segments, n+1 -> control points
    int _k;                                                 // order of polynomial    k = 2 -> linear, k = 3 -> qudratic , ...
    int _ndim;                                         // dimmension of control points being used

    void SetKnot();
    double BsplineRecursive(int i, int k, double u)  const;
    double BsplineIthBasisFunctionDerivative(int i, int k, double u)  const;
public:
    d_vec _Weight;
    d_mat _X;
    d_vec  _T;
    ui_vec _global_indices;

    CNurb(int n, int k, d_mat X,int ndim);
    CNurb(int n, int k, d_mat X, int ndim, ui_vec index);
    CNurb();
    void SetKnot(const d_vec &T);
    void SetValues(int n, int k, d_mat X, int ndim, ui_vec index);
    void SetValues(int n, int k, d_mat X,int ndim);

    void SetWeights(d_vec h);
    void NurbEval_rogers(int c,double t,int npts, d_vec & n);
    void NurbEvalControlPoint(const double u, d_vec &N, d_vec &return_val)const ;
    void NurbEval(const double u, d_vec &N, const d_vec &F, double *return_val)const ;
    void  NurbCoordinateMappingDerivativeEval(const double u, d_vec &dN,  d_vec &return_val)  const;
    void  NurbDerivativeEvalXYZ(const double u, const d_vec &F, d_vec  &dN,  d_vec &return_val)  const;
    void BsplineDerivative(const double u, d_vec &dN, d_vec &return_val) const;
    void BsplineEval(const double u, d_vec &dN, d_vec &return_val)const;
    ~CNurb();
    void PrintClass();
};
#include"Nurb.C"
#endif
