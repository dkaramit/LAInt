#ifndef LAT_HEAD
#define LAT_HEAD

#include<vector>
#include<cmath>
#include<functional>


template<class LD>
class LATrapz{
	//Local adaptive trapezoidal integration
	private:
	using vector=std::vector<LD>;
	using funcType=std::function<LD(LD)>;
	vector X;
	funcType integrand;
	
	public:
	~LATrapz()=default;
	LATrapz()=default;
	
	LATrapz(funcType integrand, const vector &X):integrand(integrand),X(X){}
	
	LD partialIntegral(size_t cell){
		LD x0,x1,dx,y0,y1;
		LD result1,result2;
		x0=X[cell-1];
		x1=X[cell];
		dx=x1-x0;
		y0=integrand(x0);
		y1=integrand(x1);
		result1 = 0.5*dx*(y1+y0);
		
		return result1;
	}
	
	LD integrate(){
		LD y0,y1,result=0;
		for(size_t i=1; i<X.size();++i){result+=partialIntegral(i);}
		
		return result;	
	}
	
	
	
	
};

#endif