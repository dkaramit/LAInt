#ifndef LAI_HEAD
#define LAI_HEAD

#include<vector>
#include<cmath>
#include<functional>


namespace lai{

template<class LD>
LD GaussLobattoRule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
	/*5-point Gauss-Lobatto quadrature rule*/
	LD dx=x1-x0;
	std::vector<LD> nodes={ 0.,-1.,1.,-std::sqrt(3/7.),std::sqrt(3/7.)};
	std::vector<LD> weights={32/45.,0.1,0.1,49/90.,49/90.};

	LD result=0;
	
	for(size_t i=0; i<5; ++i){
		result+=dx/2*weights[i]*integrand( dx/2*nodes[i] + x0+dx/2);
	}

	return result;

}

template<class LD>
LD GaussLegendreRule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
	/*5-point Gauss-Legendre quadrature rule*/
	LD dx=x1-x0;
	std::vector<LD> nodes={ 0, 1/3. *std::sqrt( 5-2*std::sqrt(10./7.) ),-1/3.*std::sqrt( 5-2*std::sqrt(10./7.) ),1/3. *std::sqrt( 5+2*std::sqrt(10./7.) ),-1/3.*std::sqrt( 5+2*std::sqrt(10./7.) )};
	std::vector<LD> weights={128/225.,(322+13*std::sqrt(70.))/900.,(322+13*std::sqrt(70.))/900.,(322-13*std::sqrt(70.))/900.,(322-13*std::sqrt(70.))/900.};
	
	LD result=0;
	
	for(size_t i=0; i<5; ++i){
		result+=dx/2*(weights[i]*integrand( dx/2*nodes[i] + x0+dx/2)  );
	}

	return result;

}



template<class LD>
LD SimpsonRule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
		LD dx=x1-x0;
		return dx/6. * ( integrand(x0)+4*integrand(x0+dx/2)+integrand(x0+dx)  );
}

template<class LD>
LD TrapezoidalRule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
		LD dx=x1-x0;
		return 0.5*dx*( integrand(x0)+integrand(x0+dx)  );
}

template<class LD>
LD MidpointRule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
		LD dx=x1-x0;
		return dx*integrand(x0+dx/2);
}

template<class LD>
LD LHSRule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
		LD dx=x1-x0;
		return dx*integrand(x0);
}

template<class LD>
LD RHSRule(std::function<LD(LD)> integrand, const LD & x0, const LD & x1){
		LD dx=x1-x0;
		return dx*integrand(x1);
}



template<class LD>
class LocalAdaptive{
	//Local adaptive trapezoidal integration
	private:
	using vector=std::vector<LD>;
	using funcType=std::function<LD(LD)>;
	using ruleType=std::function<LD(std::function<LD(LD)> integrand, const LD & x0, const LD & x1)>;



	static LD partialIntegral(funcType integrand, const LD & x0, const LD & x1, const LD & atol, const LD & rtol, const unsigned int & maxDivisions, ruleType rule){
		LD x=x0,dx=x1-x0;
		LD result1=0,result2=0,relDiff=0;

		result1 = rule(integrand,x0,x1);
		result2=result1;


		LD scale=0;
		unsigned int steps=0;
		unsigned int totalDiv=0;


		while(true){
			if(totalDiv>=maxDivisions){break;}

			x=x0;
			dx=dx*0.5;
			result2=0;

			
			steps=(x1-x0)/dx;			

			for(unsigned int count=0 ; count<steps ; ++count){
				result2 += rule(integrand,x,x+dx);
				
				x+=dx;
			}
			scale= atol + rtol*std::abs(result2);
			relDiff=std::abs( (result1-result2)/scale );

			// std::cout<<x0<<"\t"<<x1<<"\t"<<dx<<"\t"<<result2<<"\n";
			
			totalDiv++;
			result1=result2;
			if(relDiff<1){break;}

		}


		return result2;
	}
	
	public:
	~LocalAdaptive()=delete;
	LocalAdaptive()=delete;

	static LD integrate(funcType integrand, const vector &X , const LD & atol, const LD & rtol, const unsigned int & maxDivisions,
	ruleType rule){
		LD result=0;
		for(size_t i=1; i<X.size();++i){result+=partialIntegral(integrand,X[i-1],X[i],atol,rtol,maxDivisions,rule);}
		return result;	
	}
};


}




#endif