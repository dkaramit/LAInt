#include<iostream>
#include<vector>
#include<cmath>
#include<functional>

#include"LAInt.hpp"
#include"misc.hpp"



using LD=double;



int main(){

	std::function<LD(LD)> integrand=[](LD x){return std::exp(-std::pow(x-10,2)*1e5);};

	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,50,500),1e-9,1e-9,100,lai::LHSRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,50,500),1e-9,1e-9,100,lai::RHSRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,50,500),1e-9,1e-9,100,lai::MidpointRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,50,500),1e-9,1e-9,100,lai::TrapezoidalRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,50,500),1e-9,1e-9,100,lai::SimpsonRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,50,500),1e-9,1e-9,100,lai::GaussLegendreRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,misc::linspace<LD>(0,50,500),1e-9,1e-9,100,lai::GaussLobattoRule<LD>)<< "\n";
	

	return 0;
}