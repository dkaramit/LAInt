#include<iostream>
#include<vector>
#include<cmath>
#include<functional>

#include"LAInt.hpp"
#include"misc.hpp"



using LD=double;

LD BreitWigner(LD s, LD m, LD Gamma){
	return 1/(std::pow(s-m,2) + std::pow(m*Gamma,2));
}

LD NarrowWidthApprox(std::function<LD(LD s)> csx, LD m, LD Gamma){
	return  M_PI/(m*Gamma) *csx(4*m*m);
}



int main(){

	LD m=1e3, Gamma=1e-5;
	std::function<LD(LD)> integrand=[&m,&Gamma](LD x){
		if(x==0){return static_cast<LD>(0);}

		LD s=x/(1-x);
		LD dsdx=1/(1-x)/(1-x);

		return BreitWigner(s,m,Gamma)*dsdx;

	};

	/*The results should agree with this*/
	std::cout << NarrowWidthApprox([&m,&Gamma](LD s){return 1;} , m, Gamma  )<< "\n";

	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,1,60),1e-8,1e-8,50,lai::LHSRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,1,60),1e-8,1e-8,50,lai::RHSRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,1,60),1e-8,1e-8,50,lai::MidpointRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,1,60),1e-8,1e-8,50,lai::TrapezoidalRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,1,60),1e-8,1e-8,50,lai::SimpsonRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,1,60),1e-8,1e-8,50,lai::GaussLegendreRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,1,60),1e-8,1e-8,50,lai::GaussLobattoRule<LD>)<< "\n";
	

	return 0;
}