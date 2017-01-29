#ifndef NURBSURFACE_H
#define NURBSURFACE_H
#include "Types.H"
#include <vector>
#include <iomanip>
#include<math.h>
#include"Nurb.H"

class CNurbSurface{
    friend CNurbSurfaceGP;
private:
    std::vector<CNurb> _nurb_surface;
    int _ndim_parameter;
    int _ndim;
    d_vec _weights;
    ui_vec  _global_indices;
    i_vec _n,_k;
    d_vec _dXdui,_dYdui,_dZdui;// hold the derivatives of [ x wrt u,v,w] [ y wrt u,v,w] [ z wrt u,v,w] after call to NurbCoordinateMappingDerivativeEval
    double GetControlPoint( int i, int j, int k, int dim)const;
    void ErrorCheckingOnConstructor(i_vec n, i_vec k, d_mat cp_u,d_mat cp_v,d_mat cp_w,int ndim,int ndim_parameter)const;
    void ErrorCheckingOnNurbEval(d_vec uvw, d_vec &N, d_vec F, double *return_val)const;
     int GetIndexForWeights(int i,int j,int k) const;
public:
    CNurbSurface(i_vec n, i_vec k, d_mat cp_u,d_mat cp_v,d_mat cp_w,int ndim,int ndim_parameter, ui_vec  global_indices);
    CNurbSurface(i_vec n, i_vec k, d_mat cp_u,d_mat cp_v,d_mat cp_w,int ndim,int ndim_parameter);
    void NurbEval(d_vec uvw, d_vec F, d_vec &N, double *return_val) const;
    void NurbEvalControlPoint(d_vec uvw,  d_vec &return_val) const;
    void NurbCoordinateMappingDerivativeEval(d_vec uvw, d_vec &dN,  d_vec &return_val);
    void NurbDerivativeEvalXYZ(d_vec uvw, d_vec F,d_mat &dN,  d_vec &return_val);
    void  PrintClass();
    ~CNurbSurface();
};
#include"NurbSurface.C"
#endif
