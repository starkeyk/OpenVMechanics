void CNurbSurface::ErrorCheckingOnConstructor(i_vec n, i_vec k, d_mat cp_u,d_mat cp_v,d_mat cp_w,int ndim,int ndim_parameter) const{
    // data in control point is not set but either n or k is
    if(  (cp_u.size()==0) && ((n[0]!=0) || (k[0]!=0))   ){// data in control point is not set but either n or k is
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"data for control point are not set but either n or k for those control points are set "<<std::endl;
	exit(0);
    }
    if(  (cp_v.size()==0) && ((n[1]!=0) || (k[1]!=0))   ){// data in control point is not set but either n or k is
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"data for control point are not set but either n or k for those control points are set "<<std::endl;
	exit(0);
    }
    if(  (cp_w.size()==0) && ((n[2]!=0) || (k[2]!=0))   ){// data in control point is not set but either n or k is
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"data for control point are not set but either n or k for those control points are set "<<std::endl;
	exit(0);
    }

    //data for control points is set but either n or k is not given
    if(  (cp_u.size()!=0) && ((n[0]!=0) ^ (k[0]!=0))   ){// data in control point is not set but either n or k is
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"data for control point are set but either n or k for those control points are not set "<<std::endl;
	exit(0);
    }
    if(  (cp_v.size()!=0) && ((n[1]!=0) ^ (k[1]!=0))   ){// data in control point is not set but either n or k is
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"data for control point are  set but either n or k for those control points are not set "<<std::endl;
	exit(0);
    }
    if(  (cp_w.size()!=0) && ((n[2]!=0) ^ (k[2]!=0))   ){// data in control point is not set but either n or k is
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"data for control point are set but either n or k for those control points are not set "<<std::endl;
	exit(0);
    }
    if( ndim_parameter==3  &&  (n[2]==0  ||  k[2]==0)  ){
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"the values(# segments or k) for the 3rd parameter w are not set event though ndim_parameter=3 "<<std::endl;
	exit(0);
    }
    if( ndim_parameter==2  &&  (n[2]!=0  ||  k[2]!=0)  ){
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"ndim_parameter=2 but you have initialized values for n and k corresponding to the 3rd parameter "<<std::endl;
	exit(0);
    }
}

void CNurbSurface::ErrorCheckingOnNurbEval(d_vec uvw, d_vec &N, d_vec F, double *return_val) const{
   if(uvw[0]<0||uvw[1]<0){
       std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"u or v is negative "<<std::endl;
	exit(0);
   }

    if (N.size()!=( (_n[0]+1)*(_n[1]+1)*(_n[2]+1) ) ){
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"N is wrong size"<<std::endl;
	exit(0);
    }
    if (F.size()!=( (_n[0]+1)*(_n[1]+1)*(_n[2]+1) ) ){
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"F is wrong size"<<std::endl;
	exit(0);
    }
    if(uvw[0]>_nurb_surface[0]._T[_nurb_surface[0]._n+_nurb_surface[0]._k]){
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"u is too large, index out of range"<<std::endl;
	exit(0);
    }
    if(_ndim_parameter>1){
	if(uvw[1]>_nurb_surface[1]._T[_nurb_surface[1]._n+_nurb_surface[1]._k]){
	    std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	    std::cout<<"v is too large, index out of range"<<std::endl;
	    exit(0);
	}
    }
    if(_ndim_parameter>2){
	if(uvw[2]>_nurb_surface[2]._T[_nurb_surface[2]._n+_nurb_surface[2]._k]){
	    std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	    std::cout<<"w is too large, index out of range"<<std::endl;
	    exit(0);
	}
    }
    if(_ndim_parameter>_ndim){
	std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	std::cout<<"_ndim_parameter > _ndim"<<std::endl;
	exit(0);

    }
}

int CNurbSurface::GetIndexForWeights(int i,int j,int k) const{
    if(_ndim_parameter==2){
	k=0;
    }
    //     return i*(_n[0]+1) *(_n[1]+1) +  j*(_n[0]+1) + k ;
    return i +  j*(_n[0]+1) + k*(_n[0]+1)*(_n[1]+1) ;
}

void CNurbSurface::NurbEval(d_vec uvw, d_vec F,  d_vec &N,double *return_val) const{
    ErrorCheckingOnNurbEval( uvw, N,  F,  return_val);
    int i,j,k;
    double u=0,v=0,w=0;
    d_vec N_u(_n[0]+1,0); d_vec F_u(_n[0]+1,0);
    d_vec N_v(_n[1]+1,0); d_vec F_v(_n[1]+1,0);
    d_vec N_w(_n[2]+1,0);d_vec F_w(_n[2]+1,0);
    double  return_u=0;
    double  return_v=0;
    double  return_w=0;
    return_val[0] = 0;
    // performing next step to get N_u,N_v,N_w so i can compute the outer product
    for(i=0;i<_ndim_parameter;i++){
	switch ( i ) {
	    case 0:
		u=uvw[0];
		_nurb_surface[i].NurbEval(u,N_u,F_u,&return_u);
		break;
	    case 1:
		v=uvw[1];
		_nurb_surface[i].NurbEval(v,N_v,F_v,&return_v);
		break;
	    case 2:
		w=uvw[2];
		_nurb_surface[i].NurbEval(w,N_w,F_w,&return_w);
		break;
	    default:
		std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
		exit(0);
		break;
	}
    }


    // perform outer product to compute nurb  value at (u,v,w)
    double outer_product_sum = 0;
    if(_ndim_parameter==3){
	for(i=0;i<_n[0]+1;i++){ // get total sum first
	    for(j=0;j<_n[1]+1;j++){
		for(k = 0;k<_n[2]+1;k++){
		    int weight_index = GetIndexForWeights(i,j,k);
		    outer_product_sum+=_weights[weight_index]*N_u[i]*N_v[j]*N_w[k];
		}
	    }
	}
    }else{
	for(i=0;i<_n[0]+1;i++){ // get total sum first
	    for(j=0;j<_n[1]+1;j++){
		int weight_index = GetIndexForWeights(i,j,0);
		outer_product_sum+=_weights[weight_index]*N_u[i]*N_v[j];
	    }
	}

    }
    if(_ndim_parameter==3){
	for(i=0;i<_n[0]+1;i++){// compute value
	    for(j=0;j<_n[1]+1;j++){
		for(k = 0;k<_n[2]+1;k++){
		    int weight_index = GetIndexForWeights(i,j,k);
		    N[weight_index] = _weights[weight_index]*N_u[i]*N_v[j]*N_w[k]/outer_product_sum;
		    return_val[0]+=_weights[weight_index]*N_u[i]*N_v[j]*N_w[k]*F[weight_index]/outer_product_sum;
		}
	    }
	}
    }else{
	for(i=0;i<_n[0]+1;i++){// compute value
	    for(j=0;j<_n[1]+1;j++){
		int weight_index = GetIndexForWeights(i,j,0);
		N[weight_index] = _weights[weight_index]*N_u[i]*N_v[j]/outer_product_sum;
		return_val[0]+=_weights[weight_index]*N_u[i]*N_v[j]*F[weight_index]/outer_product_sum;
		// 				std::cout<<return_val[0]<<" "<<N[weight_index]<<std::endl;
	    }
	}
    }
}

void CNurbSurface::NurbEvalControlPoint(d_vec uvw,  d_vec &return_val) const{
    int i,j,k,l;
    double u=0,v=0,w=0;
    d_vec N_u(_n[0]+1,0); d_vec F_u(_n[0]+1,0);
    d_vec N_v(_n[1]+1,0); d_vec F_v(_n[1]+1,0);
    d_vec N_w(_n[2]+1,0);d_vec F_w(_n[2]+1,0);
    double  return_u=0;
    double  return_v=0;
    double  return_w=0;
    return_val.resize(_ndim);
    for(i=0;i<_ndim;i++){
	return_val[i]=0;
    }
    // performing next step to get N_u,N_v,N_w so i can compute the outer product
    for(i=0;i<_ndim_parameter;i++){
	switch ( i ) {
	    case 0:
		u=uvw[0];
		_nurb_surface[i].NurbEval(u,N_u,F_u,&return_u);
		break;
	    case 1:
		v=uvw[1];
		_nurb_surface[i].NurbEval(v,N_v,F_v,&return_v);
		break;
	    case 2:
		w=uvw[2];
		_nurb_surface[i].NurbEval(w,N_w,F_w,&return_w);
		break;
	    default:
		std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
		exit(0);
		break;
	}
    }

    // perform outer product to compute nurb  value at (u,v,w)
    double outer_product_sum = 0;
    switch ( _ndim_parameter) {
	case 2:
	    for(i=0;i<_n[0]+1;i++){ // get total sum first
		for(j=0;j<_n[1]+1;j++){
		    int weight_index = GetIndexForWeights(i,j,0);
		    outer_product_sum+=_weights[weight_index]*N_u[i]*N_v[j];
		}
	    }
	    for(i=0;i<_n[0]+1;i++){// compute value
		for(j=0;j<_n[1]+1;j++){
		    int weight_index = GetIndexForWeights(i,j,0);
		    for(l=0;l<_ndim;l++){
			return_val[l]+=_weights[weight_index]*N_u[i]*N_v[j]*GetControlPoint(  i,  j,  0,  l);
		    }
		    // 		    std::cout<<_weights[weight_index]<<" "<<N_u[i]<<" "<<_nurb_surface[0]._X[i][0]<<" "<<N_v[j]<<" "<<_nurb_surface[1]._X[j][0]<<" "<<std::endl;
		}
	    }
	    break;
	case 3:
	    for(i=0;i<_n[0]+1;i++){ // get total sum first
		for(j=0;j<_n[1]+1;j++){
		    for(k = 0;k<_n[2]+1;k++){
			int weight_index = GetIndexForWeights(i,j,k);
			outer_product_sum+=_weights[weight_index]*N_u[i]*N_v[j]*N_w[k];
		    }
		}
	    }
	    for(i=0;i<_n[0]+1;i++){// compute value
		for(j=0;j<_n[1]+1;j++){
		    for(k = 0;k<_n[2]+1;k++){
			int weight_index = GetIndexForWeights(i,j,k);
			for(l=0;l<_ndim;l++){
			    return_val[l]+=_weights[weight_index]*N_u[i]*N_v[j]*N_w[k]*GetControlPoint(  i,  j,  k,  l);
			}
		    }
		}
	    }
	    break;
	default:
	    std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
	    exit(0);
	    break;
    }

    for(i=0;i<_ndim;i++){
	return_val[i]/=outer_product_sum;
    }
}

double CNurbSurface::GetControlPoint( int i, int j, int k, int dim)const{
    switch ( _ndim_parameter ) {
	case 2:
	    return _nurb_surface[0]._X[i][dim]+_nurb_surface[1]._X[j][dim];
	    break;
	case 3:
	    return _nurb_surface[0]._X[i][dim]+_nurb_surface[1]._X[j][dim]+_nurb_surface[2]._X[k][dim];
	    break;

    }

}

void CNurbSurface::NurbCoordinateMappingDerivativeEval(d_vec uvw,d_vec &dN,  d_vec &return_val){
    int i,j,k,l;
    double u=0,v=0,w=0;

    d_vec N_u(_n[0]+1,0);
    d_vec N_v(_n[1]+1,0);
    d_vec N_w(_n[2]+1,0);

    d_vec dN_u(_n[0]+1,0); d_vec F_u(_n[0]+1,0);
    d_vec dN_v(_n[1]+1,0); d_vec F_v(_n[1]+1,0);
    d_vec dN_w(_n[2]+1,0); d_vec F_w(_n[2]+1,0);
    d_vec  return_dudxi(_ndim,0);
    d_vec  return_dvdxi(_ndim,0);
    d_vec  return_dwdxi(_ndim,0);
    d_vec  tmp_return_val(_ndim,0);

    double  W_prime_u = 0;
    double  W_prime_v = 0;
    double  W_prime_w = 0;

    for(i=0;i<_ndim_parameter;i++){
	switch ( i ) {
	    case 0:
		u=uvw[0];
		_nurb_surface[i].BsplineEval(u,  N_u,  tmp_return_val);  // just need to fill N_u
		_nurb_surface[i].BsplineDerivative(u,  dN_u,  tmp_return_val);
		break;
	    case 1:
		v=uvw[1];
		_nurb_surface[i].BsplineEval(v,  N_v,  tmp_return_val);  // just need to fill N_u
		_nurb_surface[i].BsplineDerivative(v,  dN_v,  tmp_return_val);
		break;
	    case 2:
		w=uvw[2];
		_nurb_surface[i].BsplineEval(w,  N_w,  tmp_return_val);  // just need to fill N_u
		_nurb_surface[i].BsplineDerivative(w,  dN_w,  tmp_return_val);
		break;
	    default:
		std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
		exit(0);
		break;
	}
    }
    double W = 0;
    for(i=0;i<_n[0]+1;i++){
	for(j=0;j<_n[1]+1;j++){
	    for(k=0;k<_n[2]+1;k++){
		int index =  GetIndexForWeights(i,j,k);
		switch ( _ndim_parameter ) {
		    case 2:
			W +=N_u[i]*N_v[j]*_weights[index];
			W_prime_u += dN_u[i]*N_v[j]*_weights[index];
			W_prime_v += dN_v[i]*N_u[j]*_weights[index];
			break;
		    case 3:
			W +=N_u[i]*N_v[j]*N_w[k]*_weights[index];
			W_prime_u += dN_u[i]*N_v[j]*N_w[k]*_weights[index];
			W_prime_v += dN_v[j]*N_u[i]*N_w[k]*_weights[index];
			W_prime_w += dN_w[k]*N_v[j]*N_u[i]*_weights[index];
			break;
		    default:
			std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
			exit(0);
			break;
		}
	    }
	}
    }
    _dXdui.resize(3);
    _dYdui.resize(3);
    _dZdui.resize(3);
    for(i=0;i<3;i++){
	_dXdui[i]=0;_dYdui[i]=0;_dZdui[i]=0;
    }

    for(i=0;i<_n[0]+1;i++){
	for(j=0;j<_n[1]+1;j++){
	    for(k=0;k<_n[2]+1;k++){
		int index =  GetIndexForWeights(i,j,k);
		switch ( _ndim_parameter ) {
		    case 2:
			//diff wrt u
			_dXdui[0] += ( dN_u[i]*N_v[j]*W - N_u[i]*N_v[j]*W_prime_u )*_weights[index]*GetControlPoint(  i,  j,  k,  0);
			_dYdui[0] += ( dN_u[i]*N_v[j]*W - N_u[i]*N_v[j]*W_prime_u )*_weights[index]*GetControlPoint(  i,  j,  k,  1);
			// 			    std::cout<<dN_u[i]<<" "<<N_v[j]<<" "<<W<<" "<< N_u[i]<<" "<<W_prime_u <<" "<<_dXdui[0]<<std::endl;
			//diff wrt v
			_dXdui[1] += ( dN_v[j]*N_u[i]*W - N_u[i]*N_v[j]*W_prime_v )*_weights[index]*GetControlPoint(  i,  j,  k,  0);
			_dYdui[1] += ( dN_v[j]*N_u[i]*W - N_u[i]*N_v[j]*W_prime_v )*_weights[index]*GetControlPoint(  i,  j,  k,  1);
			break;
		    case 3:
			//diff wrt u
			_dXdui[0] += ( dN_u[i]*N_v[j]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_u )*_weights[index]*GetControlPoint(  i,  j,  k,  0);
			_dYdui[0] += ( dN_u[i]*N_v[j]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_u )*_weights[index]*GetControlPoint(  i,  j,  k,  1);
			_dZdui[0] += ( dN_u[i]*N_v[j]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_u )*_weights[index]*GetControlPoint(  i,  j,  k,  2);
			//diff wrt v
			_dXdui[1] += ( dN_v[j]*N_u[i]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_v )*_weights[index]*GetControlPoint(  i,  j,  k,  0);
			_dYdui[1] += ( dN_v[j]*N_u[i]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_v )*_weights[index]*GetControlPoint(  i,  j,  k,  1);
			_dZdui[1] += ( dN_v[j]*N_u[i]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_v )*_weights[index]*GetControlPoint(  i,  j,  k,  2);
			//diff wrt w
			_dXdui[2] += ( dN_w[k]*N_u[i]*N_v[j]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_w )*_weights[index]*GetControlPoint(  i,  j,  k,  0);
			_dYdui[2] += ( dN_w[k]*N_u[i]*N_v[j]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_w )*_weights[index]*GetControlPoint(  i,  j,  k,  1);
			_dZdui[2] += ( dN_w[k]*N_u[i]*N_v[j]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_w )*_weights[index]*GetControlPoint(  i,  j,  k,  2);
			break;
		    default:
			std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
			exit(0);
			break;
		}
	    }
	}
    }
    switch ( _ndim_parameter ) {
	case 2:
	    return_val.resize(4);
	    return_val[0] = _dXdui[0];	    return_val[2] = _dXdui[1];
	    return_val[3] = _dYdui[0];	    return_val[1] = _dYdui[1];
	    break;
	case 3:
	    return_val.resize(9);
	    return_val[0] = _dXdui[0];	    return_val[3] = _dXdui[1];	    return_val[4] = _dXdui[2];
	    return_val[6] = _dYdui[0];	    return_val[1] = _dYdui[1];	    return_val[5] = _dYdui[2];
	    return_val[8] = _dZdui[0];	    return_val[7] = _dZdui[1];	    return_val[2] = _dZdui[2];
	    break;
    }
}

void CNurbSurface:: NurbDerivativeEvalXYZ(d_vec uvw, d_vec F,d_mat &dN,  d_vec &return_val){
    // set dxdu,  dxdv,  dxdw,  dy,u,  dydv,  dydw,  dzdu,  dzdv,  dzdw
    d_vec tmp;
    NurbCoordinateMappingDerivativeEval( uvw,  tmp,  return_val);// this is to set the mapping derivitives to be used later
    double dRdx=0;
    double dRdy=0;
    double dRdz=0;
    double u=0,v=0,w=0;
    double dRdu=0,dRdv=0,dRdw=0;
    int i,j,k;

    d_vec N_u(_n[0]+1,0);
    d_vec N_v(_n[1]+1,0);
    d_vec N_w(_n[2]+1,0);

    d_vec dN_u(_n[0]+1,0); d_vec F_u(_n[0]+1,0);
    d_vec dN_v(_n[1]+1,0); d_vec F_v(_n[1]+1,0);
    d_vec dN_w(_n[2]+1,0);d_vec F_w(_n[2]+1,0);
    d_vec  return_dudxi(_ndim,0);
    d_vec  return_dvdxi(_ndim,0);
    d_vec  return_dwdxi(_ndim,0);
    d_vec  tmp_return_val(_ndim,0);

    double  W_prime_u = 0;
    double  W_prime_v = 0;
    double  W_prime_w = 0;

    for(i=0;i<_ndim_parameter;i++){
	switch ( i ) {
	    case 0:
		u=uvw[0];
		return_val[0] = 0;
		_nurb_surface[i].BsplineEval(u,  N_u,  tmp_return_val);  // just need to fill N_u
		_nurb_surface[i].BsplineDerivative(u,  dN_u,  tmp_return_val);
		break;
	    case 1:
		v=uvw[1];
		return_val[1] = 0;
		_nurb_surface[i].BsplineEval(v,  N_v,  tmp_return_val);  // just need to fill N_v
		_nurb_surface[i].BsplineDerivative(v,  dN_v,  tmp_return_val);
		break;
	    case 2:
		w=uvw[2];
		return_val[2] = 0;
		_nurb_surface[i].BsplineEval(w,  N_w,  tmp_return_val);  // just need to fill N_w
		_nurb_surface[i].BsplineDerivative(w,  dN_w,  tmp_return_val);
		break;
	    default:
		std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
		exit(0);
		break;
	}
    }
    double W = 0;
    for(i=0;i<_n[0]+1;i++){        									 		//    W_prime's
	for(j=0;j<_n[1]+1;j++){
	    for(k=0;k<_n[2]+1;k++){
		int index =  GetIndexForWeights(i,j,k);
		switch ( _ndim_parameter ) {
		    case 2:
			W +=N_u[i]*N_v[j]*_weights[index];
			W_prime_u += dN_u[i]*N_v[j]*_weights[index];
			W_prime_v += dN_v[i]*N_u[j]*_weights[index];
			break;
		    case 3:
			W +=N_u[i]*N_v[j]*N_w[k]*_weights[index];
			W_prime_u += dN_u[i]*N_v[j]*N_w[k]*_weights[index];
			W_prime_v += dN_v[i]*N_u[j]*N_w[k]*_weights[index];
			W_prime_w += dN_w[k]*N_v[j]*N_u[i]*_weights[index];
			break;
		    default:
			std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
			exit(0);
			break;
		}
	    }
	}
    }

    d_vec dRdu_vec(_weights.size(),0); // used for creating dN
    d_vec dRdv_vec(_weights.size(),0);
    d_vec dRdw_vec(_weights.size(),0);
    //allocate and zero out dN
    dN.resize(_ndim_parameter);
    for(i=0;i<_ndim_parameter;i++){
	dN[i].resize(_weights.size());
	for(j = 0;j<_weights.size();j++){
	    dN[i][j] = 0;
	}
    }
    for(i=0;i<_n[0]+1;i++){        											// dRdu, dRdv, dRdw
	for(j=0;j<_n[1]+1;j++){
	    for(k=0;k<_n[2]+1;k++){
		int index =  GetIndexForWeights(i,j,k);
		switch ( _ndim_parameter ) {
		    case 2:
			//diff wrt u
			dRdu_vec[index] =  ( dN_u[i]*N_v[j]*W   -   N_u[i]*N_v[j]*W_prime_u  )*_weights[index];
			dRdu += ( dN_u[i]*N_v[j]*W   -   N_u[i]*N_v[j]*W_prime_u  )*_weights[index]*F[index];
			//diff wrt v
			dRdv_vec[index] = ( dN_v[j]*N_u[i]*W   -   N_u[i]*N_v[j]*W_prime_v )*_weights[index];
			dRdv += ( dN_v[j]*N_u[i]*W   -   N_u[i]*N_v[j]*W_prime_v )*_weights[index]*F[index];
			break;
		    case 3:
			//diff wrt u
			dRdu_vec[index] =  ( dN_u[i]*N_v[j]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_u )*_weights[index];
			dRdu += ( dN_u[i]*N_v[j]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_u )*_weights[index]*F[index];
			//diff wrt v
			dRdv_vec[index] = ( dN_v[j]*N_u[i]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_v )*_weights[index];
			dRdv += ( dN_v[j]*N_u[i]*N_w[k]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_v )*_weights[index]*F[index];
			//diff wrt w
			dRdw_vec[index] = ( dN_w[k]*N_u[i]*N_v[j]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_w )*_weights[index];
			dRdw += ( dN_w[k]*N_u[i]*N_v[j]*W - N_u[i]*N_v[j]*N_w[k]*W_prime_w )*_weights[index]*F[index];
			break;
		    default:
			std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
			exit(0);
			break;
		}
	    }
	}
    }

    double dxdu = _dXdui[0],  dxdv = _dXdui[1],  dxdw = _dXdui[2];
    double dydu = _dYdui[0],  dydv = _dYdui[1],  dydw = _dYdui[2];
    double dzdu = _dZdui[0],  dzdv = _dZdui[1],  dzdw = _dZdui[2];

    double tol = 1e-12;
    switch ( _ndim_parameter ) { // check for the divide by 0 cases, these terms are neglected by definition
	case 2:
	    // dR  =   dR du  + dR dv
	    // dx         du dx      dv dx
	    return_val[0] =0;
	    if(fabs(dxdu)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[0][i]+=dRdu_vec[i]/dxdu;
		}
		return_val[0] += dRdu/dxdu;
	    }
	    if(fabs(dxdv)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[0][i]+=dRdv_vec[i]/dxdv;
		}
		return_val[0] +=dRdv/dxdv;
	    }
	    // dR  =   dR du  + dR dv
	    // dy         du dy      dv dy
	    return_val[1] = 0;
	    if(fabs(dydu)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[1][i]+=dRdu_vec[i]/dydu;
		}
		return_val[1] += dRdu/dydu;
	    }
	    if(fabs(dydv)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[1][i]+=dRdv_vec[i]/dydv;
		}
		return_val[1] += dRdv/dydv;
	    }
	    break;
	case 3:
	    // dR  =   dR du  + dR dv  +  dR dw
	    // dx         du dx      dv dx       dw dx
	    return_val[0] =0;
	    if(fabs(dxdu)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[0][i]+=dRdu_vec[i]/dxdu;
		}
		return_val[0] += dRdu/dxdu;
	    }
	    if(fabs(dxdv)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[0][i]+=dRdv_vec[i]/dxdv;
		}
		return_val[0] +=dRdv/dxdv;
	    }
	    if(fabs(dxdw)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[0][i]+=dRdw_vec[i]/dxdw;
		}
		return_val[0] += dRdw/dxdw;
	    }
	    // dR  =   dR du  + dR dv  +  dR dw
	    // dy         du dy      dv dy       dw dy
	    return_val[1] =0;
	    if(fabs(dydu)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[1][i]+=dRdu_vec[i]/dydu;
		}
		return_val[1] += dRdu/dydu;
	    }
	    if(fabs(dydv)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[1][i]+=dRdv_vec[i]/dydv;
		}
		return_val[1] +=dRdv/dydv;
	    }
	    if(fabs(dydw)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[1][i]+=dRdw_vec[i]/dydw;
		}
		return_val[1] += dRdw/dydw;
	    }
	    // dR  =   dR du  + dR dv  +  dR dw
	    // dz         du dz      dv dz       dw dz
	    return_val[2] =0;
	    if(fabs(dzdu)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[2][i]+=dRdu_vec[i]/dzdu;
		}
		return_val[2] += dRdu/dzdu;
	    }
	    if(fabs(dzdv)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[2][i]+=dRdv_vec[i]/dzdv;
		}
		return_val[2] +=dRdv/dzdv;
	    }
	    if(fabs(dzdw)>=tol){
		for(i=0;i<_weights.size();i++){
		    dN[2][i]+=dRdw_vec[i]/dzdw;
		}
		return_val[2] += dRdw/dzdw;
	    }
	    break;
	    default:
		std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
		exit(0);
		break;
    }
//     std::cout<<"_________________________________________"<<std::endl;
//     std::cout<<"F = "<<std::endl;
//     for(i=0;i<F.size();i++){
// 	std::cout<<std::setw(10)<<F[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dN[0] = "<<dN.size()<<","<<dN[0].size()<<std::endl;
//     for(i=0;i<dN[0].size();i++){
// 	std::cout<<std::setw(10)<<dN[0][i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dN[1] = "<<std::endl;
//     for(i=0;i<dN[1].size();i++){
// 	std::cout<<std::setw(10)<<dN[1][i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dN[2] = "<<std::endl;
//     for(i=0;i<dN[2].size();i++){
// 	std::cout<<std::setw(10)<<dN[2][i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dRdu_vec = "<<std::endl;
//     for(i=0;i<dRdu_vec.size();i++){
// 	std::cout<<std::setw(10)<<dRdu_vec[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dRdv_vec= "<<std::endl;
//     for(i=0;i<dRdv_vec.size();i++){
// 	std::cout<<std::setw(10)<<dRdv_vec[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dRdw_vec= "<<std::endl;
//     for(i=0;i<dRdw_vec.size();i++){
// 	std::cout<<std::setw(10)<<dRdw_vec[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dN_u = "<<std::endl;
//     for(i=0;i<dN_u.size();i++){
// 	std::cout<<std::setw(10)<<dN_u[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dN_v = "<<std::endl;
//     for(i=0;i<dN_v.size();i++){
// 	std::cout<<std::setw(10)<<dN_v[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"dN_w = "<<std::endl;
//     for(i=0;i<dN_w.size();i++){
// 	std::cout<<std::setw(10)<<dN_w[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"N_u = "<<std::endl;
//     for(i=0;i<N_u.size();i++){
// 	std::cout<<std::setw(10)<<N_u[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"N_v = "<<std::endl;
//     for(i=0;i<N_v.size();i++){
// 	std::cout<<std::setw(10)<<N_v[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"N_w = "<<std::endl;
//     for(i=0;i<N_w.size();i++){
// 	std::cout<<std::setw(10)<<N_w[i];
//     }
//     std::cout<<std::endl;
//     std::cout<<"W_prime_u = "<<W_prime_u<<std::endl;
//     std::cout<<"W_prime_v = "<<W_prime_v<<std::endl;
//     std::cout<<"W_prime_w = "<<W_prime_w<<std::endl;
//     std::cout<<"dRdu = "<<std::setw(10)<<dRdu<<std::endl;
//     std::cout<<"dRdv = "<<std::setw(10)<<dRdv<<std::endl;
//     std::cout<<"dRdw = "<<std::setw(10)<<dRdw<<std::endl;
//     std::cout<<"dxdu = "<<std::setw(10)<<dxdu<<std::endl;
//     std::cout<<"dxdv = "<<std::setw(10)<<dxdv<<std::endl;
//     std::cout<<"dxdw = "<<std::setw(10)<<dxdw<<std::endl;
//     std::cout<<"dydu = "<<std::setw(10)<<dydu<<std::endl;
//     std::cout<<"dydv = "<<std::setw(10)<<dydv<<std::endl;
//     std::cout<<"dydw = "<<std::setw(10)<<dydw<<std::endl;
//     std::cout<<"dzdu = "<<std::setw(10)<<dzdu<<std::endl;
//     std::cout<<"dzdv = "<<std::setw(10)<<dzdv<<std::endl;
//     std::cout<<"dzdw = "<<std::setw(10)<<dzdw<<std::endl;
//     std::cout<<"_________________________________________"<<std::endl;
    //         int width = 10;
    //         std::cout<<"dx = " <<std::setw(width)<<dxdu<<std::setw(width)<<dxdv<<std::setw(width)<<dxdw<<std::endl;
    //         std::cout<<"dy = " <<std::setw(width)<<dydu<<std::setw(width)<<dydv<<std::setw(width)<<dydw<<std::endl;
    //         std::cout<<"dz = " <<std::setw(width)<<dzdu<<std::setw(width)<<dzdv<<std::setw(width)<<dzdw<<std::endl;
    //         std::cout<<"dRdu = "<<std::setw(width)<<"  "<<dRdu<<" dRdv"<<std::setw(width)<<dRdv<<" dRdw"<<std::setw(width)<<dRdw<<std::endl;
}

CNurbSurface::CNurbSurface(i_vec n, i_vec k, d_mat cp_u,d_mat cp_v,d_mat cp_w,int ndim,int ndim_parameter){
    ErrorCheckingOnConstructor( n,  k, cp_u, cp_v, cp_w, ndim, ndim_parameter);
    int i,j,m;
    int weight_size = (n[0]+1)*(n[1]+1)*(n[2]+1);
    _weights.resize(weight_size);
    for(i=0;i<weight_size;i++){
	_weights[i] = 1;
    }
    _ndim = ndim;
    _ndim_parameter = ndim_parameter;
    _nurb_surface.resize(ndim_parameter);
    _n.resize(n.size());
    _k.resize(k.size());
    for(i=0;i<n.size();i++){
	_n[i] = n[i];
    }
    for(i=0;i<k.size();i++){
	_k[i] = k[i];
    }
    for(i=0;i<ndim_parameter;i++){
	switch ( i ) {
	    case 0:
		_nurb_surface[0].SetValues(_n[0],_k[0],cp_u,ndim);
		break;
	    case 1:
		_nurb_surface[1].SetValues(_n[1],_k[1],cp_v,ndim);
		break;
	    case 2:
		_nurb_surface[2].SetValues(_n[2],_k[2],cp_w,ndim);
		break;
	    default:
		std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
		exit(0);
		break;
	}
    }
}

CNurbSurface::CNurbSurface(i_vec n, i_vec k, d_mat cp_u,d_mat cp_v,d_mat cp_w,int ndim,int ndim_parameter, ui_vec  global_indices){
    ErrorCheckingOnConstructor( n,  k, cp_u, cp_v, cp_w, ndim, ndim_parameter);
    int i,j,m;
    int weight_size = (n[0]+1)*(n[1]+1)*(n[2]+1);
    _weights.resize(weight_size);
    _global_indices.resize(weight_size);
    for(i=0;i<weight_size;i++){
	_weights[i] = 1;
	_global_indices[i] = global_indices[i];
    }
    _ndim = ndim;
    _ndim_parameter = ndim_parameter;
    _nurb_surface.resize(ndim_parameter);
    _n.resize(n.size());
    _k.resize(k.size());
    for(i=0;i<n.size();i++){
	_n[i] = n[i];
    }
    for(i=0;i<k.size();i++){
	_k[i] = k[i];
    }
    for(i=0;i<ndim_parameter;i++){
	switch ( i ) {
	    case 0:
		_nurb_surface[0].SetValues(_n[0],_k[0],cp_u,ndim);
		break;
	    case 1:
		_nurb_surface[1].SetValues(_n[1],_k[1],cp_v,ndim);
		break;
	    case 2:
		_nurb_surface[2].SetValues(_n[2],_k[2],cp_w,ndim);
		break;
	    default:
		std::cerr<<"error at "<<__LINE__<<" in "<< __func__<<std::endl;
		exit(0);
		break;
	}
    }

}

void  CNurbSurface::PrintClass(){
    int i;
    std::cout<<"_____________________________________________________________________________________"<<std::endl;
    for(i=0;i<_ndim_parameter;i++){
	std::cout<<"************************************  Nurb " <<i+1<<" ************************************"<<std::endl;
	_nurb_surface[i].PrintClass();
    }
    std::cout<<"*********************************************************************************"<<std::endl;
    int j,m,l;
    int width= 10;
    std::cout<<" |           Control Surface          | Index  |    Weights    | "<<std::endl;
    for(i=0;i<_n[0]+1;i++){
	for(j=0;j<_n[1]+1;j++){
	    for(m=0;m<_n[2]+1;m++){
		int index =  GetIndexForWeights(i,j,m);
		for(l=0;l<_ndim;l++){
		    std::cout<<std::setw(width)<<GetControlPoint(  i,  j,  m,  l)<<" ";
		}
		std::cout<<std::setw(width)<<index<<std::setw(width)<<_weights[index]<<std::endl;
	    }
	}
    }

    std::cout<<"_____________________________________________________________________________________"<<std::endl;
}

CNurbSurface::~CNurbSurface(){


}

