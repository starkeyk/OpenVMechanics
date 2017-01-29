#define pi  3.14159

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

CNurb::CNurb(int n, int k, d_mat X, int ndim){
    int i,j;
    _n = n;
    _ndim = ndim;
    _X.resize(_n+1);
    for(i=0;i<n+1;i++){
	_X[i].resize(_ndim);
    }
    _k = k;
    if(_n+1<k){
	std::cerr<<"not enough control points for spline order"<<std::endl;
	exit(0);
    }
    SetKnot();
    _Weight.resize(_n+1);
    for(i=0;i<_n+1;i++){                                    // default to B-spline
	_Weight[i] = 1;
	for(j=0;j<_ndim;j++){
	    _X[i][j] = X[i][j];
	}
    }
}

CNurb::CNurb(int n, int k, d_mat X, int ndim, ui_vec index){
    int i,j;
    _n = n;
    _ndim = ndim;
    _X.resize(_n+1);
    _global_indices.resize(n+1);
    for(i=0;i<n+1;i++){
	_X[i].resize(_ndim);
	_global_indices[i] = index[i];
    }
    _k = k;
    if(_n+1<k){
	std::cerr<<"not enough control points for spline order"<<std::endl;
	exit(0);
    }
    SetKnot();
    _Weight.resize(_n+1);
    for(i=0;i<_n+1;i++){                                    // default to B-spline
	_Weight[i] = 1;
	for(j=0;j<_ndim;j++){
	    _X[i][j] = X[i][j];
	}
    }
}

CNurb::CNurb(){


}

void CNurb::NurbEval_rogers(int c,double t,int npts, d_vec & n){
    int nplusc;
    int i,k;
    double  d,e;
    d_vec temp,x;

    nplusc = npts + c;
    temp.resize(nplusc+1);
    n.resize(npts);
    x.resize(nplusc+1);
    for (i = 1; i <= nplusc; i++){
	x[i] = _T[i-1];
    }
    /* calculate the first order basis functions n[i][1]	*/

    for (i = 1; i<= nplusc-1; i++){
	if (( t >= x[i]) && (t < x[i+1]))
	    temp[i] = 1;
	else
	    temp[i] = 0;
    }

    /* calculate the higher order basis functions */

    for (k = 2; k <= c; k++){
	for (i = 1; i <= nplusc-k; i++){
	    if (temp[i] != 0)    /* if the lower order basis function is zero skip the calculation */
		d = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
	    else
		d = 0;

	    if (temp[i+1] != 0)     /* if the lower order basis function is zero skip the calculation */
		e = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
	    else
		e = 0;

	    temp[i] = d + e;
	}
    }

    if (t == (float)x[nplusc]){		/*    pick up last point	*/
	temp[npts] = 1;
    }

    /* put in n array	*/

    for (i = 1; i <= npts; i++) {
	n[i-1] = temp[i];
    }
}

void CNurb::SetValues(int n, int k, d_mat X,int ndim){
    int i,j;
    _n = n;
    _ndim = ndim;
    _X.resize(_n+1);
    for(i=0;i<n+1;i++){
	_X[i].resize(_ndim);
    }
    _k = k;
    if(_n+1<k){
	std::cerr<<"not enough control points for spline order"<<std::endl;
	exit(0);
    }
    SetKnot();
    _Weight.resize(_n+1);
    for(i=0;i<_n+1;i++){                                    // default to B-spline
	_Weight[i] = 1;
	for(j=0;j<_ndim;j++){
	    _X[i][j] = X[i][j];
	}
    }
}

void CNurb::SetValues(int n, int k, d_mat X, int ndim, ui_vec index){
     int i,j;
    _n = n;
    _ndim = ndim;
    _X.resize(_n+1);
    _global_indices.resize(n+1);
    for(i=0;i<n+1;i++){
	_X[i].resize(_ndim);
	_global_indices[i] = index[i];
    }
    _k = k;
    if(_n+1<k){
	std::cerr<<"not enough control points for spline order"<<std::endl;
	exit(0);
    }
    SetKnot();
    _Weight.resize(_n+1);
    for(i=0;i<_n+1;i++){                                    // default to B-spline
	_Weight[i] = 1;
	for(j=0;j<_ndim;j++){
	    _X[i][j] = X[i][j];
	}
    }
}

void CNurb::SetWeights(d_vec h){
    int i=0;
    if(h.size() != _Weight.size()){
	std::cerr << "internal class weight size does not match input weights size"<<std::endl;
    }
    for(i=0;i<_n+1;i++){                                    // default to B-spline
	_Weight[i] = h[i];
    }
}

void CNurb::SetKnot(){
    int i=0;
    _T.resize(_n +_k+1 );
    for(i=0;i<_k;i++){
	_T[i] = 0;
    }
    for(i=_k;i<_n+1;i++){
	_T[i] = i - _k+1;
    }
    for(i=_n+1;i<_n+_k+1;i++){
	_T[i] = _n+2-_k;
    }

    for(i=0;i<_n+_k+1;i++){
	_T[i] = _T[i]/(_n+2-_k);//normalize
    }
}

void CNurb::SetKnot(const d_vec &T){
    int i=0;
    if (T.size()!=_n +_k+1 ){
	std::cerr<<"Input knot vector is the wrong size"<<std::endl;
	exit(0);
    }
    _T.resize(_n +_k+1 );
    for(i=0;i<_n +_k+1;i++){
	_T[i] = T[i];
    }
}

void CNurb::NurbEvalControlPoint(const double u,  d_vec &N, d_vec &return_val)const{
    int j,i;
    int num_points = _X.size();
    d_vec Xi(num_points,0);
    N.resize(num_points);
    for(j=0;j<_ndim;j++){
	for(i=0;i<num_points;i++){
	    Xi[i] = _X[i][j];
	}
	NurbEval( u,  N, Xi, &return_val[j]);
    }
    //     double sum=0;
    //     for(j=0;j<_ndim;j++){
    // 	return_val[j] = 0;
    //     }
    //
    //     for(i=0;i<num_points;i++){
    // 	double tmp=0;
    // 	tmp = BsplineRecursive(i, _k, u);
    // 	tmp*=_Weight[i];
    // 	sum+= tmp;
    // 	for(j=0;j<_ndim;j++){
    // 	    return_val[j]+=tmp*_X[i][j];
    // 	}
    // 	N[i] = tmp;
    //     }
    //     for(j=0;j<_ndim;j++){
    // 	return_val[j]/=sum;
    //     }
}

void CNurb::NurbEval(double u,  d_vec &N,const  d_vec &F, double *return_val)const{
    if(u>_T[_n+_k]){
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"u is too large, index out of range, u = "<< u<<", max index = "<<_T[_n+_k]<<std::endl;
	exit(0);
    }
    int i;
    int num_points = _X.size();
    double sum=0;
    return_val[0] = 0;
    N.resize(num_points);
    // compute extrapolation
    for(i=0;i<num_points;i++){
	double tmp=0;
	tmp = BsplineRecursive(i, _k, u);
	tmp*=_Weight[i];
	sum+= tmp;
	return_val[0] += F[i]*tmp;
	N[i] = tmp;
    }
    for(i=0;i<num_points;i++){
	N[i]/=sum;
    }
    return_val[0] /=sum;

}

void CNurb::NurbCoordinateMappingDerivativeEval(double u, d_vec &dR,  d_vec &return_val)const{
    int i, j;
    int num_points = _X.size();
    dR.resize(_n+1);
    return_val.resize(_n+1);
    d_vec N,N_sum;
    d_vec dN,dN_sum;
    d_vec tmp_return_val;
    BsplineDerivative(u,  dN,  tmp_return_val);
    BsplineEval(u,  N,  tmp_return_val);
    double W_prime = 0;
    double W = 0;
    N_sum.resize(_ndim);
    dN_sum.resize(_ndim);
    for(j=0;j<_ndim;j++){
	N_sum[j]=0; dN_sum[j]=0;return_val[j]=0;
    }
    for(i=0;i<num_points;i++){
	W_prime += _Weight[i] *dN[i];
	W += _Weight[i] *N[i];
	// 	std::cout<<std::setw(10)<<N[i]<<"  "<<std::setw(10)<<dN[i]<<std::endl;
    }
    // std::cout<<std::setw(10)<<W<<" "<<std::setw(10)<<W_prime<<std::endl;
    for(i=0;i<num_points;i++){
	for(j=0;j<_ndim;j++){
	    N_sum[j]+=N[i]*_Weight[i]*_X[i][j];
	    dN_sum[j]+=dN[i]*_Weight[i]*_X[i][j];
	}
	dR[i] = ( dN[i]*W - N[i]*W_prime ) *_Weight[i] / (W * W);
    }

    for(j=0;j<_ndim;j++){
	return_val[j]=(   dN_sum[j]*W  - N_sum[j]*W_prime   )  /   (W*W);
    }
}

void CNurb::NurbDerivativeEvalXYZ(double u, const d_vec &F, d_vec &dR,  d_vec &return_val)const{
    int i, j;
    int num_points = _X.size();
    dR.resize(_n+1);
    return_val.resize(_n+1);
    for(j=0;j<_ndim;j++){
	return_val[j]=0;
    }
    double dxidu = 0;
    d_vec dN;
    d_vec tmp_return_val;
    NurbCoordinateMappingDerivativeEval(u,  dN,  tmp_return_val);
    double tol = 1e-12;
    for(j = 0;j<_ndim;j++){
	dxidu = tmp_return_val[j];
// 		std::cout<<dxidu<<std::endl;
	for(i=0;i<num_points;i++){
	    if(fabs(dxidu)>tol){
		dR[i] = dN[i]/dxidu;
		return_val[j]+=F[i]*dN[i]/dxidu;
	    }
	}
    }
    //     std::cout<<return_val[_ndim-1]<<std::endl;
}

double CNurb::BsplineRecursive(int i, int k, double u) const{
    double Nik = 0;
    if(i>_n+_k+1){
	std::cerr<<"trying to evaluate at an index that is out of the parameters range"<<std::endl;
    }
    if(k>1){
	double b=0,num=0,denom=0;
	num = u - _T[i];
	denom = _T[i+k-1] - _T[i];
	if (denom != 0){   	                                       // indeterminate forms 0/0 are deemed to be zero
	    b = BsplineRecursive(i,k-1,u);
	    Nik = Nik + b*(num/denom);
	}
	num = _T[i+k] - u;
	denom = _T[i+k] - _T[i+1];
	if(denom != 0){
	    b = BsplineRecursive(i+1,k-1,u);
	    Nik = Nik + b*(num/denom);
	}

    }
    else if (_T[i+1] <_T[_n+_k]){
	if ((_T[i]<= u) && (u < _T[i+1])){
	    Nik = 1;
	}
    }
    else{
	if (_T[i]<= u){
	    Nik = 1;
	}
    }
    return Nik;
}

double CNurb::BsplineIthBasisFunctionDerivative(int i, int k, double u) const{
    double dNik = 0;
    double b=0,num=0,denom=0;
    b = BsplineRecursive(i,k-1,u);
    num = k-1;
    denom = _T[i+k-1] - _T[i];
    if (denom != 0){   	                                    // indeterminate forms 0/0 are deemed to be zero
	// 	std::cout<<std::endl;
	// 	std::cout<<std::setw(10)<<b<<" "<<std::setw(10)<<num<<" "<<std::setw(10)<<denom<<std::endl;
	dNik = dNik + b*(num/denom);
    }
    b = BsplineRecursive(i+1,k-1,u);
    num = k-1;
    denom = _T[i+k] - _T[i+1];
    if(denom != 0){
	// 	std::cout<<std::setw(10)<<b<<" "<<std::setw(10)<<num<<" "<<std::setw(10)<<denom<<std::endl;
	dNik = dNik - b*(num/denom);
    }
    return dNik;
}

void CNurb::BsplineEval(double u, d_vec &N, d_vec &return_val) const{
    int i,j,num_pts =_n+1;
    N.resize(_n+1);
    return_val.resize(_ndim);
    for(j=0;j<_ndim;j++){
	return_val[j]=0;
    }
    for(i=0;i<num_pts;i++){
	N[i] = BsplineRecursive(i,_k, u);
	for(j=0;j<_ndim;j++){
	    return_val[j]+=N[i]*_X[i][j];
	}
    }
}

void CNurb::BsplineDerivative(double u, d_vec &dN, d_vec &return_val) const{
    int i,j,num_pts =_n+1;
    dN.resize(_n+1);
    return_val.resize(_ndim);
    for(j=0;j<_ndim;j++){
	return_val[j]=0;
    }
    for(i=0;i<num_pts;i++){
	dN[i] = BsplineIthBasisFunctionDerivative(i,_k, u);
	for(j=0;j<_ndim;j++){
	    return_val[j]+=dN[i]*_X[i][j];
	}
    }
}

void CNurb::PrintClass(){
    int i,j, width=10;
    std::cout<<"Order is "<<_k<<", Number of segments is  "<<_n<<std::endl;
    std::cout.precision(3);
    std::cout<<"Knot Vector = "<<std::endl;
    for(i=0;i<_n + _k + 1;i++){
	std::cout<<std::setw(width)<<_T[i]<<std::endl;
    }
    std::cout<<"ndim is "<<_ndim<<std::endl<<"Control Points:  "<<_n+1<<"      Weights = "<<std::endl;
    for(i=0;i<_n + 1;i++){
	for (j=0;j<_ndim;j++){
	    std::cout<<std::setw(width)<<_X[i][j]<<" ";
	}
	std::cout<<std::setw(width)<<_Weight[i];
	std::cout<<std::endl;
    }
}

CNurb::~CNurb(){
    // 	PrintClass();

}
