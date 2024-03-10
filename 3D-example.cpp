#include<iostream>
#include<vector>
#include<cmath>
#include<functional>

#include"LAInt.hpp"
#include"misc.hpp"



using LD=double;


int main(){
    
    std::function integrand = [](LD x, LD y, LD z) { return std::sin(x/z)*y/z; };
	std::vector<std::vector<LD>> X={misc::linspace<LD>(0,1,10),misc::linspace<LD>(-2,4,5),misc::linspace<LD>(1,33,15)};

	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,X,1e-8,1e-8,5,lai::LHSRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,X,1e-8,1e-8,5,lai::RHSRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,X,1e-8,1e-8,5,lai::MidpointRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,X,1e-8,1e-8,5,lai::TrapezoidalRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,X,1e-8,1e-8,5,lai::SimpsonRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,X,1e-8,1e-8,5,lai::GaussLegendreRule<LD>)<< "\n";
	std::cout << lai::LocalAdaptive<LD>::integrate(integrand,X,1e-8,1e-8,5,lai::GaussLobattoRule<LD>)<< "\n";



return 0;
}