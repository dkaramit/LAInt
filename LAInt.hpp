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
	//Local adaptive integration.
	/*The integrate method is static, so  you don't need an instance. In fact, right now, you can't 
	declare an instance!

	I made this thing a class, because I might add some error estimation or other things that someone might need, and a class will make it easier to implement.

	The integration rules can simply be functions, or other collable objects, depending on the use case.
	*/
	private:
	using vector=std::vector<LD>;
	using funcType=std::function<LD(LD)>;
	using integrationRuleType=std::function<LD(std::function<LD(LD)> integrand, const LD & x0, const LD & x1)>;



	static LD partialIntegral(funcType integrand, const LD & x0, const LD & x1, const LD & atol, const LD & rtol, const unsigned int & maxDivisions, integrationRuleType integrationRule){
		LD x=x0,dx=x1-x0;
		LD result_prev=0,result=0,relDiff=0;

		result_prev = integrationRule(integrand,x0,x1);
		result=result_prev;


		LD scale=0;
		unsigned int NSubCells=0;
		unsigned int totalDiv=0;


		while(true){
			if(totalDiv>=maxDivisions){break;}

			x=x0;
			dx=dx*0.5;//divide the cell
			result=0;

			
			NSubCells=(x1-x0)/dx;//since the cell is divided, the integral in this cell will be the sum of all sub-cells (there are NSubCells of them)			

			for(unsigned int count=0 ; count<NSubCells ; ++count){
				result += integrationRule(integrand,x,x+dx);//sum take the integral of all subcells
				
				x+=dx;
			}

			/*When relDiff < 1, the requirements set by atol and rtol are satisfied*/
			scale= atol + rtol*std::abs(result); 
			relDiff=std::abs( (result_prev-result)/scale );

			
			totalDiv++;
			result_prev=result;
			if(relDiff<1){break;}

		}


		return result;
	}
	
	public:
	~LocalAdaptive()=delete;
	LocalAdaptive()=delete;


	/*
	Integrate integrand using the cells defined X. Each cell in X is subdivided until  atol and rtol are satisfied or until the maximum number of subdivisions (maxDivisions) 
	is reached. For the integration in each cell,  integrationRule is used (see the functions outside of this class).
	*/
	static LD integrate(funcType integrand, const vector &X , const LD & atol, const LD & rtol, const unsigned int & maxDivisions,
	integrationRuleType integrationRule){
		LD result=0;
		for(size_t i=1; i<X.size();++i){result+=partialIntegral(integrand,X[i-1],X[i],atol,rtol,maxDivisions,integrationRule);}
		return result;	
	}
};


}




#endif